# -*- coding: utf-8 -*-

from typing import Dict, Union

import numpy as np
from pandas import DataFrame

MagneticPatterns = Dict[str, Union[int, float]]


def get_interaction_coefficients(
    neighbors_df: DataFrame,
    pairwise_interactions_df: DataFrame,
    bin_to_distance_map: DataFrame,
    num_sites: int,
) -> DataFrame:
    """Computes interaction coefficients using the neighbors and interactions data
    frames.

    :param neighbors_df: A pandas ``DataFrame`` of neighbor counts aggregated over
        site-index pairs and separation distances.
    :param pairwise_interactions_df: A pandas ``DataFrame`` of pairwise interaction terms.
    :param bin_to_distance_map:
    :param num_sites:
    :return: A pandas ``DataFrame`` of interaction coefficients.
    """
    coefficients_df: DataFrame = (
        neighbors_df
            .merge(pairwise_interactions_df, on=["i", "j"])
            .sort_values("pattern")
            .assign(coefficient=lambda x: x["interaction"] * x["n"])
            .assign(rank=lambda x: x.groupby(["pattern", "subspecies_i",
                                              "subspecies_j"])
                                    .apply(lambda x: x.loc[:, ["distance_bin"]]
                                                      .rank(method="dense")))
            .assign(rank=lambda x: pd.to_numeric(x["rank"], downcast="integer"))
            .groupby(["pattern", "distance_bin", "rank", "subspecies_i", "subspecies_j"])
            .apply(lambda x: np.sum(x[["coefficient"]]) / num_sites)
            .reset_index()
            .merge(bin_to_distance_map, on=["distance_bin"])
            .sort_values(by=["pattern", "subspecies_i", "subspecies_j",
                             "distance_bin"])
            .loc[:, ["pattern", "subspecies_i", "subspecies_j", "rank",
                     "distance_ij", "coefficient"]]
    )  # yapf: disable

    return coefficients_df


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


def get_pairwise_interactions(magnetic_patterns_df: DataFrame, ) -> DataFrame:
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

