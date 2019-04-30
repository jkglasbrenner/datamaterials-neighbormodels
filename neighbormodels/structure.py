# -*- coding: utf-8 -*-

from collections import Counter
from typing import List, NamedTuple, Tuple, Union

from pymatgen import Lattice, Structure


class StructureParameters(NamedTuple):
    abc: Tuple[float, float, float]
    ang: Tuple[float, float, float]
    spacegroup: int
    species: List[str]
    coordinates: List[List[float]]


def from_parameters(structure_parameters: StructureParameters) -> Structure:
    """Generates a pymatgen ``Structure`` object using a material's structural
    parameters.

    :param structure_parameters: A ``StructureParameters`` tuple that specifies a
        material's crystal structure.
    :return: A pymatgen ``Structure`` object.
    """
    cell_lattice: Lattice = Lattice.from_lengths_and_angles(
        abc=structure_parameters.abc, ang=structure_parameters.ang
    )

    cell_structure: Structure = Structure.from_spacegroup(
        sg=structure_parameters.spacegroup,
        lattice=cell_lattice,
        species=structure_parameters.species,
        coords=structure_parameters.coordinates,
    )

    return cell_structure


def from_file(structure_file: str) -> Structure:
    """Generates a pymatgen ``Structure`` object from a supported file format.

    :param structure_file: Path to the structure file. Supported formats include CIF,
        POSCAR/CONTCAR, CHGCAR, LOCPOT, vasprun.xml, CSSR, Netcdf, and pymatgen's
        serialized structures.
    :return: A pymatgen ``Structure`` object.
    """
    cell_structure: Structure = Structure.from_file(
        filename=structure_file, primitive=False, sort=False, merge_tol=0.01
    )

    return cell_structure


def label_subspecies(
    cell_structure: Structure, site_indices: Union[List[int], int] = []
) -> None:
    """Toggles subspecies grouping on the specified site indices. Sites not found in
    the list are labeled with the atomic species name.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param site_indices: A site index or list of site indices all present in
        ``cell_structure`` (default []).
    """
    if isinstance(site_indices, int):
        site_indices = [site_indices]

    site_properties_subspecies: List[str] = get_subspecies_labels(
        cell_structure=cell_structure.copy(), site_indices=site_indices
    )

    cell_structure.add_site_property(
        property_name="subspecie", values=site_properties_subspecies
    )


def get_subspecies_labels(
    cell_structure: Structure, site_indices: List[int]
) -> List[Union[str, None]]:
    """Generates subspecies labels using the provided site indices. Sites not found in
    the list are labeled with the atomic species name.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param site_indices: A list of site indices.
    :return: A list of subspecies labels.
    """
    species_counter: Counter = Counter()
    site_properties_subspecies: List[str] = []

    for site_index, site in enumerate(cell_structure):
        specie_name: str = site.specie.name

        if site_index in site_indices:
            species_counter[specie_name] += 1
            site_properties_subspecies.append(
                f"{specie_name}{species_counter[specie_name]}"
            )

        else:
            site_properties_subspecies.append(f"{specie_name}")

    return site_properties_subspecies
