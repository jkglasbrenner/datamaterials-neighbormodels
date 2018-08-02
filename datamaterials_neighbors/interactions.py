# -*- coding: utf-8 -*-

from typing import Dict, Union

import numpy as np
from pandas import DataFrame

MagneticPatterns = Dict[str, Union[int, float]]


def build_magnetic_patterns_data_frame(
    magnetic_patterns: MagneticPatterns,
) -> DataFrame:
    """Converts dictionary of magnetic patterns into the pandas ``DataFrame``
    format.

    :param magnetic_configurations:
    :return: A pandas ``DataFrame`` labeling and specifying magnetic configurations.
    """
    magnetic_patterns_df: DataFrame = (
        DataFrame(data=magnetic_patterns)
            .reset_index(level=0)
            .rename(columns={"index": "site"})
            .melt(id_vars=["site"],
                  value_vars=list(magnetic_patterns.keys()),
                  var_name="pattern",
                  value_name="spin")
            .loc[:, ["pattern", "site", "spin"]]
    )  # yapf: disable

    return magnetic_patterns_df


def get_pairwise_interactions(
    magnetic_patterns_df: DataFrame,
) -> DataFrame:
    """Builds data frame of pairwise interaction terms.

    :param magnetic_patterns_df: A pandas ``DataFrame`` labeling and specifying
        magnetic configurations.
    :return: A pandas ``DataFrame`` of pairwise interaction terms.
    """
    pairwise_interactions_df: DataFrame = (
        magnetic_patterns_df
            .merge(right=magnetic_patterns_df, how="inner", on="pattern",
                   suffixes=("_i", "_j"))
            .rename(columns={"site_i": "i", "site_j": "j"})
            .loc[:, ["pattern", "i", "j", "spin_i", "spin_j"]]
            .assign(interaction=lambda x: x["spin_i"] * x["spin_j"])
            .query("i <= j")
    )  # yapf: disable

    return pairwise_interactions_df


def get_interaction_coefficients(neighbors, pairwise_interactions_df):
    coefficients_df: DataFrame = (
        neighbors
            .merge(pairwise_interactions_df, on=["i", "j"])
            .sort_values("pattern")
            .assign(interaction=lambda x: x["interaction"] * x["n"])
            .groupby(["pattern", "distance_bin", "subspecies_i", "subspecies_j"])
            .apply(lambda x: np.sum(x[["interaction"]]))
    )  # yapf: disable

    return coefficients_df
