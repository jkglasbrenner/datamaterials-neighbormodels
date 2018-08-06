# -*- coding: utf-8 -*-

from typing import Dict, List, NamedTuple, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pandas import Categorical, DataFrame, IntervalIndex, Series
from pandas.core.groupby import DataFrameGroupBy
from pymatgen import PeriodicSite, Structure

from neighbormodels.structure import label_subspecies

Neighbor = Tuple[PeriodicSite, float, int]
SiteNeighbors = List[Optional[Neighbor]]
AllNeighborDistances = List[SiteNeighbors]
NeighborDistances = Dict[str, Union[List[str], List[float], List[int]]]


class NeighborData(NamedTuple):
    data_frame: DataFrame
    bins_data_frame: DataFrame
    structure: Structure


def count_neighbors(
    cell_structure: Structure,
    r: float,
    unordered_pairs: bool = True,
) -> NeighborData:
    """Builds a data frame containing neighbor counts grouped over site-index pairs
    and separation distances.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param r: Radius of sphere.
    :param unordered_pairs: Treat site-index pairs as unordered if True (default True).
    :return: A named tuple with three field names:

        ``data_frame``
            A pandas ``DataFrame`` of neighbor counts aggregated over site-index pairs
            and separation distances.

        ``bins_data_frame``
            A pandas ``DataFrame`` of neighbor distances mapped to unique bin
            intervals.

        ``structure``
            A copy of the ``Structure`` object defining the crystal structure.
    """
    cell_structure = add_subspecie_labels_if_missing(cell_structure=cell_structure)

    neighbor_distances_df: DataFrame = get_neighbor_distances_data_frame(
        cell_structure=cell_structure,
        r=r,
        unordered_pairs=unordered_pairs,
    )

    distance_bins_df: DataFrame = neighbor_distances_df \
        .pipe(define_bins_to_group_and_sort_by_distance)

    neighbor_count_df: DataFrame = neighbor_distances_df \
        .pipe(group_site_index_pairs_by_distance, distance_bins_df=distance_bins_df) \
        .pipe(count_neighbors_within_distance_groups)

    return NeighborData(
        data_frame=neighbor_count_df,
        bins_data_frame=distance_bins_df,
        structure=cell_structure,
    )


def count_neighbors_within_distance_groups(
    grouped_distances: DataFrameGroupBy,
) -> DataFrame:
    """Count number of neighbors within each group of same-distance site-index pairs.

    :param grouped_distances: A data frame grouped over site-index pairs, subspecies
        pairs, and bin intervals.
    :return: A pandas ``DataFrame`` of neighbor counts aggregated over site-index pairs
        and separation distances.
    """
    return grouped_distances \
        .apply(lambda x: pd.to_numeric(arg=x["distance_ij"].count(),
                                       downcast="integer")) \
        .rename("n") \
        .reset_index()


def group_site_index_pairs_by_distance(
    neighbor_distances_df: DataFrame,
    distance_bins_df: DataFrame,
) -> DataFrameGroupBy:
    """Iterate over all sites, grouping by site-index pairs, subspecies pairs, and
    bin intervals.

    :param neighbor_distances_df: A pandas ``DataFrame`` containing all pairwise
        neighbor distances.
    :param distance_bins_df: A pandas ``DataFrame`` of neighbor distances mapped to
        unique bin intervals.
    :return: A data frame grouped over site-index pairs, subspecies pairs, and
        bin intervals.
    """
    binned_distances: Series = \
        pd.cut(x=neighbor_distances_df["distance_ij"], bins=distance_bins_df.index) \
          .rename("distance_bin")

    return neighbor_distances_df \
        .groupby(["i", "j", "subspecies_i", "subspecies_j", binned_distances])


def define_bins_to_group_and_sort_by_distance(
    neighbor_distances_df: DataFrame,
) -> DataFrame:
    """Defines bin intervals to group and sort neighbor pairs by distance.

    :param neighbor_distances_df: A pandas ``DataFrame`` of pairwise neighbor
        distances.
    :return: A pandas ``DataFrame`` of neighbor distances mapped to unique bin
        intervals.
    """
    unique_distances: np.ndarray = find_unique_distances(
        distance_ij=neighbor_distances_df["distance_ij"]
    )

    bin_intervals: IntervalIndex = define_bin_intervals(
        unique_distances=unique_distances
    )

    return DataFrame(
        data={
            "distance_bin": Categorical(values=bin_intervals, ordered=True),
            "distance_ij": Categorical(values=unique_distances, ordered=True),
        },
        index=bin_intervals,
    )


def find_unique_distances(distance_ij: Series) -> np.ndarray:
    """Finds the unique distances that define the neighbor groups.

    :param distance_ij: A pandas ``Series`` of pairwise neighbor distances.
    :return: An array of unique neighbor distances.
    """
    unique_floats: np.ndarray = np.sort(distance_ij.unique())

    next_distance_not_close: np.ndarray = np.logical_not(
        np.isclose(unique_floats[1:], unique_floats[:-1])
    )

    return np.concatenate((
        unique_floats[:1],
        unique_floats[1:][next_distance_not_close],
    ))


def define_bin_intervals(unique_distances: np.ndarray) -> IntervalIndex:
    """Constructs bin intervals used to group over neighbor distances.

    This binning procedure provides a robust method for grouping data based on a
    variable with a float data type.

    :param unique_distances: An array of neighbor distances returned by asking
        pandas to return the unique distances.
    :return: A pandas ``IntervalIndex`` defining bin intervals can be used to sort
        and group neighbor distances.
    """
    bin_centers: np.ndarray = np.concatenate(([0], unique_distances))

    bin_edges: np.ndarray = np.concatenate([
        bin_centers[:-1] + (bin_centers[1:] - bin_centers[:-1]) / 2,
        bin_centers[-1:] + (bin_centers[-1:] - bin_centers[-2:-1]) / 2,
    ])

    return IntervalIndex.from_breaks(breaks=bin_edges)


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

    :param site_i_index: Site index of first site in neighbor pair.
    :param site_i_neighbors: A list of site i's neighbors.
    :param cell_structure: The pymatgen ``Structure`` object that defines the crystal
        structure.
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


def add_subspecie_labels_if_missing(cell_structure: Structure) -> Structure:
    """Makes a copy of ``cell_structure`` and then checks if ``cell_structure`` has
    the subspecie site property. If it does, then return the copy as-is, otherwise
    label each site of the copy using the site's atomic specie name and then return
    it.

    :param cell_structure: A pymatgen ``Structure`` object.
    :return: An exact copy of the input ``cell_structure`` object with subspecie
        labels added, if missing.
    """
    cell_structure = cell_structure.copy()

    if "subspecie" not in cell_structure.site_properties:
        label_subspecies(cell_structure=cell_structure, site_indices=[])

    return cell_structure
