# Python 2/3 compatibility
# from __future__ import absolute_import, division, print_function
# from builtins import (bytes, str, open, super, range, zip, round, input, int, pow, object)

from enum import Enum

available_descriptors = ['bispectrum_so4',
                         'bispectrum_so3',
                         'behler',
                         'distance_2b',
                         'coordination',
                         'angle_3b',
                         'co_angle_3b',
                         'co_distance_2b',
                         'cosnx',
                         'trihis',
                         'water_monomer',
                         'water_dimer',
                         'A2_dimer',
                         'AB_dimer',
                         'bond_real_space',
                         'atom_real_space',
                         'power_so3',
                         'power_so4',
                         'soap',
                         'AN_monomer',
                         'general_monomer',
                         'general_dimer',
                         'general_trimer',
                         'rdf',
                         'as_distance_2b',
                         'molecule_lo_d',
                         'alex',
                         'com_dimer',
                         'distance_Nb']


class MissingParameter(Exception):
    def __init__(self, parameter=None):
        if parameter:
            message = f"The {parameter} parameter is missing."
        else:

            message = f"A parameter is missing."

        # Call the base class constructor with the parameters it needs
        super().__init__(message)


class sparse_method(Enum):
    """sparse_method options are:

    * **RANDOM**: default, chooses n_sparse random datapoints
    * **PIVOT**: based on the full covariance matrix finds the n_sparse "pivoting" points
    * **CLUSTER**: based on the full covariance matrix performs a k-medoid clustering into n_sparse clusters, returning the medoids
    * **UNIFORM**: makes a histogram of the data based on n_sparse and returns a data point from each bin
    * **KMEANS**: k-means clustering based on the data points
    * **COVARIANCE**: greedy data point selection based on the sparse covariance matrix, to minimise the GP variance of all datapoints
    * **UNIQ**: selects unique datapoints from the dataset
    * **FUZZY**: fuzzy k-means clustering
    * **FILE**: reads sparse points from a file
    * **INDEX_FILE**: reads indices of sparse points from a file
    * **CUR_COVARIANCE**: CUR, based on the full covariance matrix
    * **CUR_POINTS**: CUR, based on the datapoints
    """
    RANDOM = 'RANDOM'
    PIVOT = 'PIVOT'
    CLUSTER = 'CLUSTER'
    UNIFORM = 'UNIFORM'
    KMEANS = 'KMEANS'
    COVARIANCE = 'COVARIANCE'
    UNIQ = 'UNIQ'
    FUZZY = 'FUZZY'
    FILE = 'FILE'
    INDEX_FILE = 'INDEX_FILE'
    CUR_COVARIANCE = 'CUR_COVARIANCE'
    CUR_POINTS = 'CUR_POINTS'


class covariance_type(Enum):
    """Type of covariance function to use. Available:

    * **ARD_SE**,
    * **DOT_PRODUCT**,
    * **BOND_REAL_SPACE**,
    * **PP**: piecewise polynomial
    """
    ARD_SE = 'ARD_SE'
    DOT_PRODUCT = 'DOT_PRODUCT'
    BOND_REAL_SPACE = 'BOND_REAL_SPACE'
    PP = 'PP'


class Descriptor(object):
    pass


class DescriptorNew(object):
    """ GAP options
    This subroutine parses the options given in the gap string, for each GAP.

    Args:
        covariance_type (covariance_type): Type of covariance function to use. Available: ARD_SE, DOT_PRODUCT, BOND_REAL_SPACE, PP (piecewise polynomial)
        delta (float): Set the standard deviation of the Gaussian process. Typically this would be set to the standard deviation (i.e. root mean square) of the function that is approximated with the Gaussian process.
        add_species (bool): Create species-specific descriptor, using the descriptor string as a template.
        config_type_n_sparse (str): Number of sparse points in each config type. Format: {type1:50:type2:100}
        f0 (float): Set the mean of the Gaussian process. Defaults to 0.
        mark_sparse_atoms (bool): Reprints the original xyz file after sparsification process. sparse propery added, true for atoms associated with a sparse point.
        n_sparse (int): Number of sparse points to use in the sparsification of the Gaussian process
        print_sparse_index (str): If given, after determinining the sparse points, their 1-based indices are appended to this file
        sparse_file (str): Sparse points from a file. Integers, in single line.
        sparse_method (sparse_method): Sparsification method. RANDOM(default), PIVOT, CLUSTER, UNIFORM, KMEANS, COVARIANCE, NONE, FUZZY, FILE, INDEX_FILE, CUR_COVARIANCE, CUR_POINTS
        theta_fac (str): Set the width of Gaussians for the ARD_SE and PP kernel by multiplying the range of each descriptor by theta_fac. Can be a single number or different for each dimension. For multiple theta_fac separate each value by whitespaces.
        theta_file (str): Set the width of Gaussians for the ARD_SE kernel from a file. There should be as many real numbers as the number of dimensions, in a single line
        theta_uniform (float): Set the width of Gaussians for the ARD_SE and PP kernel, same in each dimension.
        unique_descriptor_tolerance (float): Descriptor tolerance when filtering out duplicate data points
        unique_hash_tolerance (float): Hash tolerance when filtering out duplicate data points
        zeta (float): Exponent of soap type dot product covariance kernel

    """

    def __init__(self,
                 covariance_type=None,
                 delta=None,
                 Name='Default',
                 add_species=False,
                 config_type_n_sparse=None,
                 f0=0.0,
                 mark_sparse_atoms=False,
                 n_sparse=0,
                 print_sparse_index=None,
                 sparse_file=None,
                 sparse_method=sparse_method.RANDOM,
                 theta_fac=1.0,
                 theta_file=None,
                 theta_uniform=0.0,
                 unique_descriptor_tolerance=1.0e-10,
                 unique_hash_tolerance=1.0e-10,
                 zeta=1.0,
                 ):

        print(delta)
        if covariance_type is None:
            raise MissingParameter()
        if delta is None:
            raise MissingParameter('delta')

        self.covariance_type = covariance_type
        self.delta = delta
        self.Name = Name
        self.add_species = add_species
        self.config_type_n_sparse = config_type_n_sparse
        self.f0 = f0
        self.mark_sparse_atoms = mark_sparse_atoms
        self.n_sparse = n_sparse
        self.print_sparse_index = print_sparse_index


        self.sparse_file = sparse_file
        self.sparse_method = sparse_method
        self.theta_fac = theta_fac
        self.theta_file = theta_file
        self.theta_uniform = theta_uniform
        self.unique_descriptor_tolerance = unique_descriptor_tolerance
        self.unique_hash_tolerance = unique_hash_tolerance
        self.zeta = zeta


class fourier_so4(Descriptor):
    """fourier_so4

    Args:
        cutoff (float): Cutoff for SO4 bispectrum
        z0_ratio (float): Ratio of radius of 4D projection sphere times PI and the cutoff.
        j_max (int): Max of expansion of bispectrum, i.e. resulution
        Z (int): Atomic number of central atom
        n_species (int): Number of species for the descriptor
        species_Z (int __OR__ ?): Atomic number of species
        w (float __OR__ ?): Weight associated to each atomic type
    """

    def __init__(self, cutoff=2.75, z0_ratio=0.0, j_max=4, Z=0, n_species=1, species_Z=None, w=None):
        pass


class bispectrum_so4(Descriptor):
    """bispectrum_so4

    Args:
    """

    def __init__(self):
        pass


class bispectrum_so3(Descriptor):
    """bispectrum_so3

    Args:
        cutoff (float): Cutoff for bispectrum_so3-type descriptors
        min_cutoff (float): Cutoff for minimal distances in bispectrum_so3-type descriptors
        l_max (int): L_max for bispectrum_so3-type descriptors
        n_max (int): N_max for bispectrum_so3-type descriptors
        Z (int): Atomic number of central atom
        n_species (int): Number of species for the descriptor
        species_Z (int __OR__ ?)): Atomic number of species
        w (float __OR__ ?)): Weight associated to each atomic type
    """

    def __init__(self, cutoff=0.00, min_cutoff=0.00, l_max=4, n_max=4, Z=0, n_species=1, species_Z=None,
                 w=None):
        pass


class behler(Descriptor):
    """behler

    Args:
        behler_cutoff (float): Cutoff for Behler-type descriptors
    """

    def __init__(self, behler_cutoff=2.75):
        pass


class distance_2b(DescriptorNew):
    """distance_2b

    Args:
        cutoff (float): Cutoff for distance_2b-type descriptors
        cutoff_transition_width (float): Transition width of cutoff for distance_2b-type descriptors
        Z1 (int): Atom type #1 in bond
        Z2 (int): Atom type #2 in bond
        resid_name (str): Name of an integer property in the atoms object giving the residue id of the molecule to which the atom belongs.
        only_intra (bool): Only calculate INTRAmolecular pairs with equal residue ids (bonds)
        only_inter (bool): Only apply to INTERmolecular pairs with different residue ids (non-bonded)
    """

    def __init__(self, cutoff=0.00, cutoff_transition_width=0.5, Z1=0, Z2=0, resid_name='', only_intra=False,
                 only_inter=False, **kwargs):
        print(kwargs)
        super().__init__(**kwargs)


class coordination(Descriptor):
    """coordination

    Args:
        cutoff (float): Cutoff for coordination-type descriptors
        transition_width (float): Width of transition region from 1 to 0
        Z (int): Atomic number of central atom
    """

    def __init__(self, cutoff=0.00, transition_width=0.20, Z=0):
        pass


class angle_3b(Descriptor):
    """angle_3b

    Args:
        cutoff (float): Cutoff for angle_3b-type descriptors
        Z (int): Atomic number of central atom
        Z1 (int): Atomic number of neighbour #1
        Z2 (int): Atomic number of neighbour #2
    """

    def __init__(self, cutoff=0.00, Z=0, Z1=0, Z2=0):
        pass


class co_angle_3b(Descriptor):
    """co_angle_3b

    Args:
        cutoff (float): Cutoff for co_angle_3b-type descriptors
        coordination_cutoff (float): Cutoff for coordination function in co_angle_3b-type descriptors
        coordination_transition_width (float): Transition width for co_angle_3b-type descriptors
        Z (int): Atomic number of central atom
        Z1 (int): Atomic number of neighbour #1
        Z2 (int): Atomic number of neighbour #2
    """

    def __init__(self, cutoff=0.00, coordination_cutoff=0.00, coordination_transition_width=0.00, Z=0, Z1=0, Z2=0):
        pass


class co_distance_2b(Descriptor):
    """co_distance_2b

    Args:
        cutoff (float): Cutoff for co_distance_2b-type descriptors
        transition_width (float): Transition width of cutoff for co_distance_2b-type descriptors
        coordination_cutoff (float): Cutoff for coordination function in co_distance_2b-type descriptors
        coordination_transition_width (float): Transition width for co_distance_2b-type descriptors
        Z1 (int): Atom type #1 in bond
        Z2 (int): Atom type #2 in bond
    """

    def __init__(self, cutoff=0.00, transition_width=0.50, coordination_cutoff=0.00, coordination_transition_width=0.00,
                 Z1=0, Z2=0):
        pass


class cosnx(Descriptor):
    """cosnx

    Args:
        cutoff (float): Cutoff for cosnx-type descriptors
        min_cutoff (float): Cutoff for minimal distances in cosnx-type descriptors
        l_max (int): L_max for cosnx-type descriptors
        n_max (int): N_max for cosnx-type descriptors
        Z (int): Atomic number of central atom
        n_species (int): Number of species for the descriptor
        species_Z (int __OR__ ?)): Atomic number of species
        w (float __OR__ ?)): Weight associated to each atomic type
    """

    def __init__(self, cutoff=0.00, min_cutoff=0.00, l_max=4, n_max=4, Z=0, n_species=1, species_Z=None,
                 w=None):
        pass


class trihis(Descriptor):
    """trihis

    Args:
        cutoff (float): Cutoff for trihis-type descriptors
        n_gauss (int): Number of Gaussians for trihis-type descriptors
        trihis_gauss_centre: Number of Gaussians for trihis-type descriptors
        trihis_gauss_width: Number of Gaussians for trihis-type descriptors
    """

    def __init__(self, cutoff=0.00, n_gauss=0, trihis_gauss_centre=None, trihis_gauss_width=None):
        pass


class water_monomer(Descriptor):
    """water_monomer

    Args:
        cutoff (float): Cutoff for water_monomer-type descriptors
    """

    def __init__(self, cutoff=0.00):
        pass


class water_dimer(Descriptor):
    """water_dimer

    Args:
        cutoff (float): Cutoff for water_dimer-type descriptors
        cutoff_transition_width (float): Width of smooth cutoff region for water_dimer-type descriptors
        monomer_cutoff (float): Monomer cutoff for water_dimer-type descriptors
        OHH_ordercheck (bool): T: find water molecules. F: use default order OHH
        power (float): Power of distances to be used in the kernel
    """

    def __init__(self, cutoff=0.00, cutoff_transition_width=0.50, monomer_cutoff=1.50, OHH_ordercheck=True,
                 power=1.0):
        pass


class A2_dimer(Descriptor):
    """A2_dimer

    Args:
        cutoff (float): Cutoff for A2_dimer-type descriptors
        monomer_cutoff (float): Monomer cutoff for A2_dimer-type descriptors
        atomic_number (int): Atomic number in A2_dimer-type descriptors
    """

    def __init__(self, cutoff=0.00, monomer_cutoff=1.50, atomic_number=1):
        pass


class AB_dimer(Descriptor):
    """AB_dimer

    Args:
        cutoff (float): Cutoff for AB_dimer-type descriptors
        monomer_cutoff (float): Monomer cutoff for AB_dimer-type descriptors
        atomic_number1 (int): Atomic number of atom 1 in AB_dimer-type descriptors
        atomic_number2 (int): Atomic number of atom 2 in AB_dimer-type descriptors
    """

    def __init__(self, cutoff=0.00, monomer_cutoff=1.50, atomic_number1=1, atomic_number2=9):
        pass


class bond_real_space(Descriptor):
    """bond_real_space

    Args:
        bond_cutoff (float): Bond cutoff for bond_real_space-type descriptors
        bond_transition_width (float): Bond transition width for bond_real_space-type descriptors
        cutoff (float): Space cutoff for bond_real_space-type descriptors
        transition_width (float): Space transition width for bond_real_space-type descriptors
        atom_sigma (float): Atom sigma for bond_real_space-type descriptors
        max_neighbours (int): Maximum number of neighbours
    """

    def __init__(self, bond_cutoff=0.00, bond_transition_width=0.00, cutoff=0.00, transition_width=0.00,
                 atom_sigma=0.00, max_neighbours=0):
        pass


class atom_real_space(Descriptor):
    """atom_real_space

    Args:
        cutoff (float): Space cutoff for atom_real_space-type descriptors
        cutoff_transition_width (float): Space transition width for atom_real_space-type descriptors
        l_max (int): Cutoff for spherical harmonics expansion
        alpha (float): Width of atomic Gaussians
        zeta (float): Exponent of covariance function
    """

    def __init__(self, cutoff=0.00, cutoff_transition_width=0.00, l_max=0, alpha=1.0, zeta=1.0):
        pass


class power_so3(Descriptor):
    """power_so3

    Args:
        cutoff (float): Cutoff for power_so3-type descriptors
        min_cutoff (float): Cutoff for minimal distances in power_so3-type descriptors
        l_max (int): L_max for power_so3-type descriptors
        n_max (int): N_max for power_so3-type descriptors
        Z (int): Atomic number of central atom
        n_species (int): Number of species for the descriptor
        species_Z (int __OR__ ?)): Atomic number of species
        w (float) __OR__ ?): Weight associated to each atomic type

    """

    def __init__(self, cutoff=0.00, min_cutoff=0.00, l_max=4, n_max=4, Z=0, n_species=1, species_Z=None, w=None):
        pass


class power_so4(Descriptor):
    """power_so4

    Args:
    """

    def __init__(self): pass

    class soap(Descriptor):
        """soap

        Args:
            species_Z: Atomic number of species
            Z (int): Atomic number of central atom, 0 is the wild-card __OR__ Atomic numbers to be considered for central atom, must be a list
            cutoff: Cutoff for soap-type descriptors
            atom_sigma: Width of atomic Gaussians for soap-type descriptors
            n_max: N_max (number of radial basis functions) for soap-type descriptors
            l_max: L_max (spherical harmonics basis band limit) for soap-type descriptors
            n_Z (int): How many different types of central atoms to consider
            n_species (int): Number of species for the descriptor
            cutoff_transition_width (float): Cutoff transition width for soap-type descriptors
            cutoff_dexp (int): Cutoff decay exponent
            cutoff_scale (float): Cutoff decay scale
            cutoff_rate (float): Inverse cutoff decay rate
            central_weight (float): Weight of central atom in environment
            central_reference_all_species (bool): Place a Gaussian reference for all atom species densities. By default (F) only consider when neighbour is the same species as centre
            covariance_sigma0 (float): sigma_0 parameter in polynomial covariance function
            xml_version (int): Version of GAP the XML potential file was created
            basis_error_exponent (float): 10^(-basis_error_exponent) is the max difference between the target and the expanded function
            average (bool): Whether to calculate averaged SOAP - one descriptor per atoms object. If false (default) atomic SOAP is returned.
            diagonal_radial (bool): Only return the n1=n2 elements of the power spectrum.
            normalise (bool): Normalise descriptor so magnitude is 1. In this case the kernel of two equivalent environments is 1.

        """

        def __init__(self, cutoff, n_max, l_max, atom_sigma, cutoff_transition_width=0.50, cutoff_dexp=0,
                     cutoff_scale=1.0, cutoff_rate=1.0, central_weight=1.0, central_reference_all_species=False,
                     average=False, diagonal_radial=False, covariance_sigma0=0.0, normalise=True,
                     basis_error_exponent=10.0, n_Z=1, n_species=1, xml_version=1426512068, species_Z=None, Z=0):
            pass


class AN_monomer(Descriptor):
    """AN_monomer

    Args:
        cutoff (float): Cutoff for AN_monomer-type descriptors
        atomic_number (int): Atomic number in AN_monomer-type descriptors
        N (int): Number of atoms in cluster
        do_atomic (bool): Descriptors are cluster based or atom-based
    """

    def __init__(self, cutoff=0.00, atomic_number=1, N=4, do_atomic=True):
        pass


class general_monomer(Descriptor):
    """general_monomer

    Args:
        cutoff (float): Cutoff for general_monomer-type descriptors
        signature: Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}
        atom_ordercheck (bool): T: find molecules. F: go by order of atoms
        strict (bool): Raise error if not all atoms assigned to monomer
        power (float): Power of distances to be used in the kernel
    """

    def __init__(self, cutoff=0.00, signature=None, atom_ordercheck=True, strict=True, power=1.0):
        pass


class com_dimer(Descriptor):
    """com_dimer

    Args:
        cutoff (float): Cutoff(intermolecular) for com_dimer-type descriptors
        monomer_one_cutoff (float): Cutoff(mono1) for com_dimer-type descriptors
        monomer_two_cutoff (float): Cutoff(mono2) for com_dimer-type descriptors
        cutoff_transition_width (float): Width of smooth cutoff region for com_dimer-type descriptors
        atom_ordercheck (bool): T: find molecules. F: go by order of atoms
        strict (bool): Raise error if not all atoms assigned to monomer or if no monomer pairs found
        mpifind (bool): Use find_monomer_pairs_MPI
        signature_one: Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}
        signature_two: Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}
    """

    def __init__(self, cutoff=0.00, monomer_one_cutoff=0.00, monomer_two_cutoff=0.00, cutoff_transition_width=0.50,
                 atom_ordercheck=True, strict=True, mpifind=False, signature_one=None, signature_two=None):
        pass


class general_dimer(Descriptor):
    """general_dimer

    Args:
        cutoff (float): Cutoff(intermolecular) for general_dimer-type descriptors
        monomer_one_cutoff (float): Cutoff(mono1) for general_dimer-type descriptors
        monomer_two_cutoff (float): Cutoff(mono2) for general_dimer-type descriptors
        cutoff_transition_width (float): Width of smooth cutoff region for general_dimer-type descriptors
        internal_swaps_only (bool): F: energies will be symmetrised over swaps of nuclei between monomers
        atom_ordercheck (bool): T: find molecules. F: go by order of atoms
        double_count (bool): T: double count when constructing the dimers, for compatibility with water dimer descriptor, default False
        strict (bool): Raise error if not all atoms assigned to monomer or if no monomer pairs found
        use_com (bool): Use COM instead of COG
        mpifind (bool): Use find_monomer_pairs_MPI
        signature_one: Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}
        signature_two: Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}
    """

    def __init__(self, cutoff=0.00, monomer_one_cutoff=0.00, monomer_two_cutoff=0.00, cutoff_transition_width=0.50,
                 internal_swaps_only=True, atom_ordercheck=True, double_count=False, strict=True, use_com=False,
                 mpifind=False, signature_one=None, signature_two=None):
        pass


class general_trimer(Descriptor):
    """general_trimer

    Args:
        cutoff (float): Cutoff(intermolecular) for general_trimer-type descriptors
        monomer_one_cutoff (float): Cutoff(mono1) for general_trimer-type descriptors
        monomer_two_cutoff (float): Cutoff(mono2) for general_trimer-type descriptors
        monomer_three_cutoff (float): Cutoff(mono3) for general_trimer-type descriptors
        cutoff_transition_width (float): Width of smooth cutoff region for general_trimer-type descriptors
        internal_swaps_only (bool): F: energies will be symmetrised over swaps of nuclei between monomers
        atom_ordercheck (bool): T: find molecules. F: go by order of atoms
        strict (bool): Raise error if not all atoms assigned to monomer or if no monomer pairs found
        use_com (bool): Use COM instead of COG
        mpifind (bool): Use find_monomer_triplets_MPI
        signature_one: Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}
        signature_two: Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}
        signature_three: Atomic numbers of monomer three, format {Z1 Z2 Z3 ...}
        power (float): Power of distances to be used in the kernel
    """

    def __init__(self, cutoff=0.00, monomer_one_cutoff=0.00, monomer_two_cutoff=0.00, monomer_three_cutoff=0.00,
                 cutoff_transition_width=0.50, internal_swaps_only=True, atom_ordercheck=True, strict=True,
                 use_com=False, mpifind=False, signature_one=None, signature_two=None, signature_three=None, power=1.0):
        pass


class rdf(Descriptor):
    """rdf

    Args:
        cutoff (float): Cutoff for rdf-type descriptors
        transition_width (float): Width of transition region from 1 to 0
        Z (int): Atomic number of central atom
        r_min (float): Atomic number of central atom
        r_max (float): Atomic number of central atom
        n_gauss (int): Atomic number of central atom
        w_gauss (float): Atomic number of central atom
    """

    def __init__(self, cutoff=0.00, transition_width=0.20, Z=0, r_min=0.0, r_max=0.0, n_gauss=10, w_gauss=0.0):
        pass


class as_distance_2b(Descriptor):
    """as_distance_2b

    Args:
        min_cutoff (float): Lower cutoff for as_distance_2b-type descriptors
        max_cutoff: Higher cutoff for as_distance_2b-type descriptors
        as_cutoff: Cutoff of asymmetricity
        overlap_alpha (float): Cutoff of asymmetricity
        min_transition_width (float): Transition width of lower cutoff for as_distance_2b-type descriptors
        max_transition_width (float): Transition width of higher cutoff for as_distance_2b-type descriptors
        as_transition_width (float): Transition width of asymmetricity cutoff for as_distance_2b-type descriptors
        coordination_cutoff: Cutoff for coordination function in as_distance_2b-type descriptors
        coordination_transition_width (float): Transition width for as_distance_2b-type descriptors
        Z1 (int): Atom type #1 in bond
        Z2 (int): Atom type #2 in bond
    """

    def __init__(self, min_cutoff=0.00, max_cutoff=None, as_cutoff=None, overlap_alpha=0.50, min_transition_width=0.50,
                 max_transition_width=0.50, as_transition_width=0.10, coordination_cutoff=None,
                 coordination_transition_width=0.50, Z1=0, Z2=0):
        pass


class molecule_lo_d(Descriptor):
    """molecule_lo_d

    Args:
        cutoff (float): Cutoff for molecule_lo_d-type descriptors
        atoms_template_string: Atoms object which serves as a template - written to a string
        neighbour_graph_depth (int): Ignore distances between atoms separated by more than this number of bonds
        signature (str): Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}
        symmetry_property_name (str): Integer arrays specifying symmetries - see header of make_permutations_v2.f95 for format
        append_file (str): Pairs of atoms for which we want the distance to be additionally included in the descriptor
        atom_ordercheck (bool): T: basic check that atoms in same order as in template F: assume all xyz frames have atoms in same order
        desctype (int): 0: distance matrix, 1: inverse distance matrix, 2: Coulomb matrix, 3: exponential
    """

    def __init__(self, cutoff=0.00, atoms_template_string='', neighbour_graph_depth=2, signature='',
                 symmetry_property_name='symm', append_file='', atom_ordercheck=True, desctype=0):
        pass


class alex(Descriptor):
    """alex

    Args:
        cutoff (float): Cutoff for alex-type descriptors
        Z (int): Atomic number of central atom
        power_min (int): Minimum power of radial basis for the descriptor
        power_max (int): Maximum power of the radial basis for the descriptor
        n_species (int): Number of species for the descriptor
        species_Z: Atomic number of species
    """

    def __init__(self, cutoff=0.00, Z=0, power_min=5, power_max=10, n_species=1, species_Z=None):
        pass


class distance_Nb(Descriptor):
    """distance_Nb

    Args:
        cutoff: Cutoff for distance_Nb-type descriptors
        cutoff_transition_width (float): Transition width of cutoff for distance_Nb-type descriptors
        order: Many-body order, in terms of number of neighbours
        compact_clusters (bool): If true, generate clusters where the atoms have at least one connection to the central atom. If false, only clusters where all atoms are connected are generated.
        Z: Atomic type of neighbours
        atom_mask_name=None: Name of a logical property in the atoms object. For atoms where this property is true descriptors are calculated.
        xml_version (int): Version of GAP the XML potential file was created
        do_transfer (bool): Enable transfer function
        transfer_factor (float): Transfer function: stretch factor
        transfer_width (float): Transfer function: transition width
        transfer_r0 (float): Transfer function: transition distance
    """

    def __init__(self, order, cutoff, cutoff_transition_width=0.5, compact_clusters=False, Z=None,
                 atom_mask_name=None, xml_version=1423143769, do_transfer=False, transfer_factor=5.0,
                 transfer_width=1.0, transfer_r0=3.0):
        pass


# class Descriptors(object):
#     pass
#
#
# class distance_2b(Descriptors):
#     """2-body descriptor
#
#     Args:
#         cutoff: Cutoff for distance_2b-type descriptors
#         cutoff_transition_width: Transition width of cutoff for
#         Z1: Atom type #1 in bond
#         Z2: Atom type #2 in bond
#         resid_name:Name of an integer property in
#         only_intra: Only calculate INTRAmolecular pairs with equal residue ids (
#         only_inter: Only apply to INTERmolecular pairs with different residue ids (
#     """
#
#     def __init__(self,
#                  cutoff=0.00,
#                  cutoff_transition_width=0.5,
#                  Z1=0,
#                  Z2=0,
#                  resid_name='',
#                  only_intra=False,
#                  only_inter=False,
#                  ):
#         pass
#
#
# # Type of descriptor is soap
# class Soap(Descriptors):
#     """
#     SOAP descriptor
#
#     Args:
#         cutoff: Cutoff for soap-type descriptors
#         cutoff_transition_width: Cutoff transition width for
#         cutoff_dexp: Cutoff decay exponent
#         cutoff_scale: Cutoff decay scale
#         cutoff_rate: Inverse cutoff decay rate
#         l_max: L_max (spherical harmonics basis band limit) for soap-type
#         n_max: N_max (number of radial basis functions) for soap-type
#         atom_sigma: Width of atomic Gaussians for soap-type
#         central_weight: Weight of central atom in environment
#         central_reference_all_species: Place a Gaussian reference for all atom species densities. By default (F) only consider when neighbour is the same species as centre
#         average: Whether to calculate averaged SOAP - one descriptor per atoms object.
#         diagonal_radial: Only return the n1=n2 elements of the power
#         covariance_sigma0: sigma_0 parameter in polynomial covariance
#         normalise: Normalise descriptor so magnitude is 1. In this case the kernel
#         basis_error_exponent: 10^(-basis_error_exponent) is the max
#         n_Z: How many different types of central atoms to consider
#         n_species: Number of species for the
#         species_Z: Atomic number of species
#         xml_version: Version of GAP the XML potential file was created
#         species_Z: Atomic number of species
#         Z: Atomic number of central atom, 0 is the wild-card __OR__ Atomic numbers to be considered for central atom, must be a list
#     """
#
#     def __init__(self,
#                  cutoff,
#                  n_max,
#                  l_max,
#                  atom_sigma,
#                  cutoff_transition_width=0.50,
#                  cutoff_dexp=0,
#                  cutoff_scale=1.0,
#                  cutoff_rate=1.0,
#                  central_weight=1.0,
#                  central_reference_all_species=False,
#                  average=False,
#                  diagonal_radial=False,
#                  covariance_sigma0=0.0,
#                  normalise=True,
#                  basis_error_exponent=10.0,
#                  n_Z=1,
#                  n_species=1,
#                  xml_version=1426512068,
#                  species_Z=None,
#                  Z=None
#                  ):
#         pass


if __name__ == '__main__':
    d = distance_2b(delta=0.1)
    print(d)
