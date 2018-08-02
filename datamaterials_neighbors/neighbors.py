# -*- coding: utf-8 -*-

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pandas import DataFrame, IntervalIndex, Series
from pandas.core.groupby import DataFrameGroupBy
from pymatgen import PeriodicSite, Structure

from datamaterials_neighbors.structure import label_subspecies

Neighbor = Tuple[PeriodicSite, float, int]
SiteNeighbors = List[Optional[Neighbor]]
AllNeighborDistances = List[SiteNeighbors]
NeighborDistances = Dict[str, Union[List[str], List[float], List[int]]]


def count_neighbors_grouped_by_site_index_pairs_and_distance(
    cell_structure: Structure,
    r: float,
    unordered_pairs: bool = True,
) -> DataFrame:
    """Builds a data frame containing neighbor counts grouped over site-index pairs
    and separation distances.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param r: Radius of sphere.
    :param unordered_pairs: Treat site-index pairs as unordered if True (default True).
    :return: A pandas ``DataFrame`` of neighbor counts aggregated over site-index pairs
        and separation distances.
    """
    cell_structure = cell_structure.copy()
    add_subspecie_labels_if_missing(cell_structure=cell_structure)

    neighbor_distances_df: DataFrame = get_neighbor_distances_data_frame(
        cell_structure=cell_structure,
        r=r,
        unordered_pairs=unordered_pairs,
    )

    grouped_distances: DataFrameGroupBy = \
        group_neighbors_within_site_index_pairs_by_distance(
            neighbor_distances_df=neighbor_distances_df,
        )

    distances_summary: DataFrame = (
        grouped_distances
            .describe()
            .loc[:, "distance_ij"]
            .loc[:, ["max", "count"]]
            .reset_index()
            .loc[:, ["i", "j", "subspecies_i", "subspecies_j", "max", "count"]]
            .rename(columns={"max": "distance_ij", "count": "n"})
            .assign(n=lambda x: pd.to_numeric(x["n"], downcast='integer'))
    )  # yapf: disable

    return distances_summary


def group_neighbors_within_site_index_pairs_by_distance(
    neighbor_distances_df: DataFrame,
) -> DataFrameGroupBy:
    """Iterate over all sites, grouping by site-index pairs and binning by distance.

    :param neighbor_distances_df: A pandas ``DataFrame`` containing all pairwise
        neighbor distances.
    :return: A pandas data frame grouped by site-index pairs and neighbor distances.
    """
    unique_distances: np.ndarray = find_unique_distances(
        distance_ij=neighbor_distances_df["distance_ij"]
    )

    bin_intervals: IntervalIndex = define_bin_intervals(
        unique_distances=unique_distances
    )

    binned_distances: Series = pd.cut(
        x=neighbor_distances_df["distance_ij"],
        bins=bin_intervals,
    )

    return neighbor_distances_df.groupby([
        "i",
        "j",
        "subspecies_i",
        "subspecies_j",
        binned_distances,
    ])


def find_unique_distances(distance_ij: Series) -> np.ndarray:
    """Find the unique distances that define the neighbor groups.

    :param distance_ij: A pandas ``Series`` of pairwise neighbor distances.
    :return: An array of unique neighbor distances.
    """
    unique_floats: np.ndarray = np.sort(distance_ij.unique())
    next_distance_not_close: np.ndarray = np.logical_not(
        np.isclose(unique_floats[1:], unique_floats[:-1])
    )
    unique_distances: np.ndarray = np.concatenate((
        unique_floats[:1],
        unique_floats[1:][next_distance_not_close],
    ))

    return unique_distances


def define_bin_intervals(unique_distances: np.ndarray) -> IntervalIndex:
    """Constructs intervals with bins centered around the neighbor distances.

    This binning procedure provides a robust method for grouping data based on a
    variable with a float data type.

    :param unique_distances: An array of neighbor distances returned by asking
        pandas to return the unique distances.
    :return: A pandas ``IntervalIndex`` that defines bins centered around the
        neighbor distances.
    """
    bin_centers: np.ndarray = np.concatenate(([0], unique_distances))
    bin_edges: np.ndarray = np.concatenate([
        bin_centers[:-1] + (bin_centers[1:] - bin_centers[:-1]) / 2,
        bin_centers[-1:] + (bin_centers[-1:] - bin_centers[-2:-1]) / 2,
    ])
    bin_intervals: IntervalIndex = IntervalIndex.from_breaks(breaks=bin_edges)

    return bin_intervals


def get_neighbor_distances_data_frame(
    cell_structure: Structure,
    r: float,
    unordered_pairs: bool,
) -> DataFrame:
    """Get data frame of pairwise neighbor distances for each atom in the unit cell,
    out to a distance ``r``.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param r: Radius of sphere.
    :param unordered_pairs: Treat site-index pairs as unordered if True.
    :return: A pandas ``DataFrame`` of pairwise neighbor distances.
    """
    all_neighbors: AllNeighborDistances = cell_structure.get_all_neighbors(
        r=r,
        include_index=True,
    )

    neighbor_distances: NeighborDistances = extract_neighbor_distance_data(
        cell_structure=cell_structure,
        all_neighbors=all_neighbors,
        unordered_pairs=unordered_pairs,
    )

    return DataFrame(data=neighbor_distances)


def extract_neighbor_distance_data(
    cell_structure: Structure,
    all_neighbors: AllNeighborDistances,
    unordered_pairs: bool,
) -> NeighborDistances:
    """Extracts the site indices, site species, and neighbor distances for each pair
    and stores it in a dictionary.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param all_neighbors: A list of lists containing the neighbors for each site in
        the structure.
    :param unordered_pairs: Treat site-index pairs as unordered if True.
    :return: A dictionary of site indices, site species, and neighbor distances for
        each pair.
    """
    neighbor_distances: NeighborDistances = {
        "i": [],
        "j": [],
        "subspecies_i": [],
        "subspecies_j": [],
        "distance_ij": [],
    }

    site_i_index: int
    site_i_neighbors: SiteNeighbors
    for site_i_index, site_i_neighbors in enumerate(all_neighbors):
        append_site_i_neighbor_distance_data(
            site_i_index=site_i_index,
            site_i_neighbors=site_i_neighbors,
            cell_structure=cell_structure,
            neighbor_distances=neighbor_distances,
            unordered_pairs=unordered_pairs,
        )

    return neighbor_distances


def append_site_i_neighbor_distance_data(
    site_i_index: int,
    site_i_neighbors: SiteNeighbors,
    cell_structure: Structure,
    neighbor_distances: NeighborDistances,
    unordered_pairs: bool,
) -> None:
    """Helper function to append indices, species, and distances in the
    ``neighbor_distances`` dictionary.

    :param site_i_index:
    :param site_i_neighbors:
    :param cell_structure: A pymatgen ``Structure`` object.
    :param neighbor_distances: A dictionary of site indices, site species, and neighbor
        distances for each pair.
    :param unordered_pairs: Treat site-index pairs as unordered if True.
    """
    site_j: Neighbor
    for site_j in site_i_neighbors:
        subspecies_pair: List[str] = sorted([
            cell_structure[site_i_index].properties["subspecie"],
            cell_structure[site_j[2]].properties["subspecie"],
        ])

        index_pair: List[str] = [site_i_index, site_j[2]]
        if unordered_pairs:
            index_pair.sort()

        neighbor_distances["i"].append(index_pair[0])
        neighbor_distances["j"].append(index_pair[1])
        neighbor_distances["subspecies_i"].append(subspecies_pair[0])
        neighbor_distances["subspecies_j"].append(subspecies_pair[1])
        neighbor_distances["distance_ij"].append(site_j[1])


def add_subspecie_labels_if_missing(cell_structure: Structure) -> None:
    """Checks if ``cell_structure`` has the subspecie site property. If not, then
    label each site using the site's atomic specie name.

    :param cell_structure: A pymatgen ``Structure`` object.
    """
    if "subspecie" not in cell_structure.site_properties:
        label_subspecies(
            cell_structure=cell_structure,
            site_indices=[],
        )
