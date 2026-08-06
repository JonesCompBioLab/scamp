import pandas as pd
import os
import collections
from typing import List, Tuple
import numba as nb
import sklearn
import torch
from torch import nn
import numpy as np

class SCAMP(nn.Module):
    """Count-based MLP

    A simple multi-layer perceptron model that will predict ecDNA status
    from count distributions.

    Args:
        n_layers: Number of layers
        n_hidden: Hidden dimensions per layer
        activation: Activation function for each layer (default 'relu')
        dropout_rate: Dropout rate per layer (default 0.1)
    """

    def __init__(
        self,
        n_layers: int = 1,
        n_hidden: int = 10,
        activation: str = "relu",
        dropout_rate: float = 0.1,
    ):

        super().__init__()

        if activation == "relu":
            activation_fn = nn.ReLU
        elif activation == "tanh":
            activation_fn = nn.Tanh
        elif activation == "sigmoid":
            activation_fn = nn.Sigmoid
        elif activation == None:
            if n_layers > 0:
                raise Exception(
                    "If n_layers > 1 you must specify a valid activation function."
                )
        else:
            raise Exception(
                "Please specify one of activation functions: `relu`"
                " `sigmoid` or `tanh`"
            )

        self.trained = False
        self.statistics = {}

        n_in, n_out = 14, 2
        layers_dim = [n_in] + (n_layers) * [n_hidden] + [n_out]

        self.fc_layers = nn.Sequential(
            collections.OrderedDict(
                [
                    (
                        f"Layer {i}",
                        nn.Sequential(
                            nn.Linear(
                                n_in,
                                n_out,
                            ),
                            activation_fn() if activation != None else None,
                            (
                                nn.Dropout(p=dropout_rate)
                                if (dropout_rate > 0) and (i != len(layers_dim))
                                else None
                            ),
                        ),
                    )
                    for i, (n_in, n_out) in enumerate(
                        zip(layers_dim[:-1], layers_dim[1:], strict=True)
                    )
                ]
            )
        )

        # save module kwargs
        self._module_kwargs = {
            "n_layers": n_layers,
            "n_hidden": n_hidden,
            "dropout_rate": dropout_rate,
            "activation": activation,
        }

    def forward(self, x):
        """Forward computation on `x`

        Args:
            x: Tensor of size n_in (here 14)

        Returns:
            Softmax probability.
        """

        for i, layers in enumerate(self.fc_layers):
            for layer in layers:
                if layer is None:
                    continue
                x = layer(x)
        return x

    def fit(
        self,
        X,
        y,
        n_epochs: int = 100,
        learning_rate: float = 1e-3,
        batch_size: int = 128,
        optimizer: str = "Adam",
        num_workers: int = 1,
        verbose=True,
        reporting_freq=1000,
    ):
        """Train model.

        Train model based on paired count distributions and ecDNA/HSR
        annotations.

        Args:
            X_train: Features for training set.
            y_train: Labels for training set.
            n_epochs: Number of epochs.
            learning_rate: Learning rate for optimizer.
            batch_size: Batch size.
            optimizer: Algorithm for optimization.
            num_workers: Number of workers to use.
        """

        # initialize training objects
        if optimizer == "Adam":
            optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)

        # loss_fn = nn.BCELoss()
        loss_fn = nn.CrossEntropyLoss()

        X_train, y_train = torch.Tensor(X), torch.Tensor(y).view(
            y.shape[0], 2
        ).to(torch.float)

        # create batchified dataset and loader
        scamp_dataset = torch.utils.data.TensorDataset(X_train, y_train)
        dataloader = torch.utils.data.DataLoader(
            scamp_dataset,
            batch_size=batch_size,
            shuffle=True,
            num_workers=num_workers,
        )

        # log statistics
        self.statistics["Loss"] = []

        # put into train mode
        self.train()

        for epoch in range(n_epochs):

            running_loss = 0.0
            for _, data in enumerate(dataloader, 0):

                X_batch, y_batch = data
                _y_batch = self.forward(X_batch)
                loss = loss_fn(_y_batch, y_batch)

                # backward pass
                loss.backward()

                # updates
                optimizer.step()

                # zero gradients
                optimizer.zero_grad()

                # track loss
                running_loss += loss.item()

            self.statistics["Loss"].append(running_loss)

            if verbose and (epoch % reporting_freq == 0):

                preds = self.proba(X_train).round()

                # put back into train mode
                self.train()

                accuracy = sklearn.metrics.accuracy_score(
                    preds.detach().numpy(), y_train.detach().numpy()
                )
                print(
                    f"[epoch:{epoch}]: The loss value for training part={running_loss}, accuracy={accuracy}"
                )

        self.trained = True

    @torch.no_grad()
    def proba(self, x):

        self.eval()
        _x = self.forward(torch.Tensor(x))
        return torch.softmax(_x, dim=1)

    def save(self, model_save_dir: str, overwrite: bool = False):
        """Save model.

        Args:
            path_dir: File path to save directory.
            overwrite: Overwrite existing saved model.
        """

        if self.trained:

            if (os.path.exists(model_save_dir)) and (not overwrite):
                raise Exception(
                    "Model save directory already exists and overwrite is False."
                )

            os.makedirs(model_save_dir, exist_ok=overwrite)

            model_state_dict = self.state_dict()
            model_kwargs = self._module_kwargs

            torch.save(
                {
                    "model_state_dict": model_state_dict,
                    "attr_dict": model_kwargs,
                },
                f"{model_save_dir}/model.pt",
            )
            print(f"Saved model at {model_save_dir}/model.pt")
        else:
            raise Exception("Please train model before saving.")

    @classmethod
    def load(cls, model_save_dir):

        saved_model_state = torch.load(f"{model_save_dir}/model.pt")

        model_params = saved_model_state["model_state_dict"]
        model_kwargs = saved_model_state["attr_dict"]

        loaded_model = cls(**model_kwargs)
        loaded_model.load_state_dict(model_params)
        
        loaded_model.trained = True
        loaded_model.eval()

        return loaded_model

    def prepare_copy_numbers(
        self,
        copy_numbers,
        genes,
        min_copy_number,
        max_percentile,
        filter_copy_number=3,
    ) -> Tuple[np.array, np.array]:
        """Prepare copy-number data.

        Prepares copy-number data for models. Currently summarizes copy-number
        distribution according to mean, variance, dispersion, inter-quartile
        range and quantiles.

        Args:
            copy_numbers: Copy-number file path.
            min_copy_number: Ignore copy-numbers below this value for computing
                amplification statistics.
            max_percentile: Max percentile to cap distribution.
            filter_copy_number: Throw out genes whose mean copy-number is less
                than a pre-defined value (default = 3).

        Returns:
            Numpy array storing counts per cell of genes.
        """

        # filter by mean copy number
        means = np.mean(copy_numbers, axis=0)
        kii = np.where(means >= filter_copy_number)[0]
        genes = genes[kii]
        copy_numbers = copy_numbers[:, kii]

        X = np.zeros((len(genes), 14))

        for g_i in range(len(genes)):

            vals = copy_numbers[:, g_i]

            cap = np.percentile(vals, max_percentile)
            vals[vals > cap] = np.nan
            mu, sig, dispersion, percentiles, iqr = (
                self.summarize_count_distribution(
                    np.array(vals).astype(float), min_count=min_copy_number
                )
            )
            X[g_i, :] = [mu, sig, dispersion, iqr] + percentiles

        mask = ~np.isnan(X).any(axis=1)
        genes = genes[mask]
        X = X[mask, :]

        return X, genes
    
    def summarize_count_distribution(self, counts, min_count=2):

        sub_counts = counts[counts > min_count]
        if len(sub_counts) < 2:
            return np.nan, np.nan, np.nan, [np.nan]*10, np.nan
        
        mu = np.nanmean(sub_counts)
        sig = np.nanvar(sub_counts)
        dispersion = sig / mu
        quantiles = [np.nanpercentile(sub_counts, p) for p in np.arange(0, 100, step=10)]

        iqr = np.nanpercentile(sub_counts, 75) - np.nanpercentile(sub_counts, 25)
        
        return mu, sig, dispersion, quantiles, iqr

