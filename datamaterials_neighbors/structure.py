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
        abc=structure_parameters.abc,
        ang=structure_parameters.ang,
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
        filename=structure_file,
        primitive=False,
        sort=False,
        merge_tol=0.01,
    )

    return cell_structure


def set_sublattice_parameters(
    cell_structure: Structure,
    species: Union[List[str], str],
) -> None:
    """Applies unique sublattice labels to sites matching the specified atomic species.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param species: An atomic specie or list of atomic species all present in
        ``cell_structure``.
    """
    if isinstance(species, str):
        species = [species]

    site_properties_sublattices: List[Union[str, None]] = enumerate_sublattices(
        cell_structure=cell_structure.copy(),
        species=species,
    )

    cell_structure.add_site_property(
        property_name="sublattice",
        values=site_properties_sublattices,
    )


def enumerate_sublattices(
    cell_structure: Structure,
    species: List[str],
) -> List[Union[str, None]]:
    """Enumerates sublattice labels to sites matching the specified atomic species into
    a list.

    :param cell_structure: A pymatgen ``Structure`` object.
    :param species: A list of atomic species.
    :return: A list of unique sublattice labels.
    """
    species_counter: Counter = Counter()
    site_properties_sublattices: List[Union[str, None]] = []

    for site in cell_structure:
        specie: str = site.species_string

        if specie in species:
            species_counter[specie] += 1
            site_properties_sublattices.append(f"{specie}{species_counter[specie]}")

        else:
            site_properties_sublattices.append(None)

    return site_properties_sublattices
