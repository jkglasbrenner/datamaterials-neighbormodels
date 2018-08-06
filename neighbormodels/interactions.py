# -*- coding: utf-8 -*-

from typing import Dict, Union

import numpy as np
import pandas as pd
from pandas import DataFrame, Series

from neighbormodels.neighbors import NeighborData

MagneticPatterns = Dict[str, Union[int, float]]


def build_model(
    neighbor_data: NeighborData,
    magnetic_patterns: MagneticPatterns,
) -> DataFrame:
    """Builds and returns a data frame describing a pairwise interaction model. The
    intended use-case for the model is fitting magnetic energies taken from density
    functional theory calculations and extracting exchange parameters.

    :param neighbor_data: A named tuple with three field names:

        ``data_frame``
            A pandas ``DataFrame`` of neighbor counts aggregated over site-index pairs
            and separation distances.

        ``bins_data_frame``
            A pandas ``DataFrame`` of neighbor distances mapped to unique bin
            intervals.

        ``structure``
            A copy of the ``Structure`` object defining the crystal structure.
    :param magnetic_patterns: A dictionary of magnetic patterns to be mapped onto
        the crystal structure and used to compute the interaction coefficients of the
        model.
    :return: A pandas ``DataFrame`` of the interaction parameter names and coefficients
        for the pairwise interaction model.
    """
    magnetic_patterns_df: DataFrame = build_magnetic_patterns_data_frame(
        magnetic_patterns=magnetic_patterns
    )

    return magnetic_patterns_df \
        .pipe(compute_interaction_signs) \
        .pipe(compute_model_coefficients, neighbor_data=neighbor_data) \
        .pipe(spread_parameter_name_column)


def compute_interaction_signs(magnetic_patterns_df: DataFrame) -> DataFrame:
    """Computes the signs of the pairwise interactions for the magnetic model.

    :param magnetic_patterns_df: A pandas ``DataFrame`` labeling and specifying
        magnetic patterns.
    :return: A pandas ``DataFrame`` of the signs of the pairwise interactions.
    """
    return magnetic_patterns_df \
        .merge(right=magnetic_patterns_df,
               how="inner",
               on="pattern",
               suffixes=("_i", "_j")) \
        .rename(columns={"site_i": "i", "site_j": "j"}) \
        .loc[:, ["pattern", "i", "j", "spin_i", "spin_j"]] \
        .assign(sign=lambda x: x["spin_i"] * x["spin_j"]) \
        .query("i <= j")


def compute_model_coefficients(
    interaction_signs_df: DataFrame,
    neighbor_data: NeighborData,
) -> DataFrame:
    """Computes the model coefficients by aggregating over the dot product of the
    neighbor counts and interaction signs.

    :param interaction_signs_df: A pandas ``DataFrame`` of the signs of the pairwise
        interactions.
    :param neighbor_data: A named tuple with three field names:

        ``data_frame``
            A pandas ``DataFrame`` of neighbor counts aggregated over site-index pairs
            and separation distances.

        ``bins_data_frame``
            A pandas ``DataFrame`` of neighbor distances mapped to unique bin
            intervals.

        ``structure``
            A copy of the ``Structure`` object defining the crystal structure.
    :return: A pandas ``DataFrame`` of the model coefficients.
    """
    return neighbor_data.data_frame \
        .pipe(multiply_interaction_signs_and_neighbor_count,
              interaction_signs_df=interaction_signs_df) \
        .pipe(group_subspecie_pairs_and_rank_by_distance) \
        .pipe(label_interaction_parameters) \
        .pipe(aggregate_interaction_coefficients,
              num_sites=neighbor_data.structure.num_sites)


def multiply_interaction_signs_and_neighbor_count(
    data_frame: DataFrame,
    interaction_signs_df: DataFrame,
) -> DataFrame:
    """Adds coefficient column to data frame. Compatible with the pandas ``pipe()``
    method.

    :param data_frame: A pandas ``DataFrame`` of neighbor counts aggregated over
        site-index pairs and separation distances.
    :param interaction_signs_df: A pandas ``DataFrame`` of the signs of the pairwise
        interactions.
    :return: A copy of input ``data_frame`` with the coefficient column added.
    """
    return data_frame \
        .merge(interaction_signs_df, on=["i", "j"]) \
        .sort_values("pattern") \
        .assign(coefficient=lambda x: x["sign"] * x["n"])


def group_subspecie_pairs_and_rank_by_distance(data_frame: DataFrame) -> DataFrame:
    """Adds rank column to data frame that sorts subspecie pairs by neighbor distance.
    Compatible with the pandas ``pipe()`` method.

    :param data_frame: A pandas ``DataFrame`` of neighbor counts aggregated over
        site-index pairs and separation distances.
    :return: A copy of input ``data_frame`` with the rank column added.
    """
    df: DataFrame = data_frame.copy()
    df["rank"] = df \
        .groupby(["subspecies_i", "subspecies_j"]) \
        .apply(lambda x: x[["distance_bin"]].rank(method="dense")) \
        .apply(lambda x: pd.to_numeric(arg=x, downcast="integer"))

    return df


def label_interaction_parameters(data_frame: DataFrame) -> DataFrame:
    """Adds parameter_names column to data frame that labels the unique parameters
    of the interaction model. Compatible with the pandas ``pipe()`` method.

    :param data_frame: A pandas ``DataFrame`` of neighbor counts aggregated over
        site-index pairs and separation distances ranked by subspecies pairs.
    :return: A copy of input ``data_frame`` with the parameter_name column added.
    """
    df = data_frame.copy()

    parameter_names: Series = "J" + df["rank"].astype(str)

    single_specie_check: Series = df["subspecies_i"] == df["subspecies_j"]
    single_specie: bool = single_specie_check.all()

    if not single_specie:
        parameter_names += "_" + df["subspecies_i"] + df["subspecies_j"]

    df["parameter_name"] = parameter_names

    return df


def aggregate_interaction_coefficients(
    data_frame: DataFrame,
    num_sites: int,
) -> DataFrame:
    """Aggregates interaction coefficients by summing them within groups defined by
    the magnetic patterns and interaction parametres. Compatible with the pandas
    ``pipe()`` method.

    :param data_frame: A pandas ``DataFrame`` of interaction coefficients and parameters
        grouped by magnetic patterns, subspecies pairs, and neighbor distances.
    :param num_sites: The total number of magnetic sites in the unit cell.
    :return: A data frame where the interaction coefficients were summed within
        groups defined by the magnetic patterns and interaction parameters.
    """
    return data_frame \
        .groupby(["pattern", "parameter_name"]) \
        .apply(lambda x: np.sum(x[["coefficient"]]) / num_sites) \
        .reset_index() \
        .sort_values(by=["pattern", "parameter_name"]) \
        .loc[:, ["pattern", "parameter_name", "coefficient"]]


def build_magnetic_patterns_data_frame(
    magnetic_patterns: MagneticPatterns,
) -> DataFrame:
    """Converts dictionary of magnetic patterns into the pandas ``DataFrame``
    format.

    :param magnetic_patterns: A dictionary of magnetic patterns to be mapped onto
        the crystal structure and used to compute the interaction coefficients of the
        model.
    :return: A pandas ``DataFrame`` labeling and specifying magnetic patterns.
    """
    return DataFrame(data=magnetic_patterns) \
        .reset_index(level=0) \
        .rename(columns={"index": "site"}) \
        .melt(id_vars=["site"],
              value_vars=list(magnetic_patterns.keys()),
              var_name="pattern",
              value_name="spin") \
        .loc[:, ["pattern", "site", "spin"]]


def spread_parameter_name_column(data_frame: DataFrame) -> DataFrame:
    """Spreads interaction parameter names into their own columns with the interaction
    coefficient as rows. Compatible with the pandas ``pipe()`` method.

    :param data_frame: A data frame where the interaction coefficients were summed
        within groups defined by the magnetic patterns and interaction parameters.
    :return: A data frame with the parameter names pivoted into their own columns.
    """
    df: DataFrame = data_frame \
        .pivot(index="pattern", columns="parameter_name", values="coefficient") \
        .reset_index()

    df.columns.name = ""

    return df
