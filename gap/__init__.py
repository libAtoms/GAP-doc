# -*- coding: utf-8 -*-

# Main options
# quippy.teach_sparse_parse_command_line(*args, **kwargs)
# This subroutine parses the main command line options.
import os
from enum import Enum
from gap.descriptors import Descriptor


class Verbose(Enum):
    """ Verbosity control. Options:

    * **NORMAL**,
    * **VERBOSE**,
    * **NERD**,
    * **ANAL**.
    """
    Normal = 'NORMAL'
    Verbose = 'VERBOSE'
    Nerd = 'NERD'
    Anal = 'ANAL'


class GAP(object):
    """GAP Wrapper
    This subroutine parses the main command line options.

    Args:
        gap (list of Descriptors): Initialisation string for GAPs
        default_sigma (float): Error in [energies forces virials hessians]
        config_type_sigma (str): What sigma values to choose for each type of data. Format:
            {type:energy:force:virial:hessian}
        core_ip_args (str): QUIP init string for a potential to subtract from data (and added back after prediction)
        core_param_file (str): QUIP XML file for a potential to subtract from data (and added back after prediction)
        do_e0_avg (bool): Method of calculating e0 if not explicitly specified. If true, computes the average atomic
            energy in input data. If false, sets e0 to the lowest atomic energy in the input data.
        do_ip_timing (bool): To enable or not timing of the interatomic potential.
        e0 (str): Atomic energy value to be subtracted from energies before fitting (and added back on after
            prediction). Specifiy a single number (used for all species) or by species: {Ti:-150.0:O:-320}.
            energy = core + GAP + e0
        e0_offset (float): Offset of baseline. If zero, the offset is the average atomic energy of the input data or
            the e0 specified manually.
        hessian_delta (float): Delta to use in numerical differentiation when obtaining second derivative for the
            Hessian covariance
        sigma_parameter_name (str): Sigma parameters (error hyper) for a given configuration in the database. Overrides
            the command line sigmas. In the XYZ, it must be prepended by energy\_, force\_, virial\_ or hessian\_
        sigma_per_atom (bool): Interpretation of the energy and virial sigmas specified in >>default_sigma<< and
            >>config_type_sigma<<. If >>T<<, they are interpreted as per-atom errors, and the variance will be scaled
            according to the number of atoms in the configuration. If >>F<< they are treated as absolute errors and no
            scaling is performed. NOTE: sigmas specified on a per-configuration basis (see >>sigma_parameter_name<<)
            are always absolute.
        sparse_jitter (float): Intrisic error of atomic/bond energy, used to regularise the sparse covariance matrix
        sparse_use_actual_gpcov (bool): Use actual GP covariance for sparsification methods
        template_file (str): Template XYZ file for initialising object
        verbosity (Verbose): Verbosity control.
    """

    def __init__(self, gap,
                 default_sigma,
                 config_type_sigma=None,
                 core_ip_args=None,
                 core_param_file='quip_params.xml',
                 do_e0_avg=True,
                 do_ip_timing=False,
                 e0='0',
                 e0_offset=0.0,
                 hessian_delta=1.0,
                 sigma_parameter_name='sigma',
                 sigma_per_atom=True,
                 sparse_jitter=1.0,
                 sparse_use_actual_gpcov=False,
                 template_file='template.xml',
                 verbosity=Verbose.Normal,
                 ):
        self.gap = gap
        self.default_sigma = default_sigma

        self.config_type_sigma = config_type_sigma
        self.core_ip_args = core_ip_args
        self.core_param_file = core_param_file
        self.do_e0_avg = do_e0_avg
        self.do_ip_timing = do_ip_timing
        self.e0 = e0
        self.e0_offset = e0_offset
        self.gap = gap
        self.hessian_delta = hessian_delta
        self.sigma_parameter_name = sigma_parameter_name
        self.sigma_per_atom = sigma_per_atom
        self.sparse_jitter = sparse_jitter
        self.sparse_use_actual_gpcov = sparse_use_actual_gpcov
        self.template_file = template_file
        self.verbosity = verbosity

        self.at_file = None
        self.gp_file = None
        self.config_type_parameter_name = None
        self.energy_parameter_name = None
        self.force_parameter_name = None
        self.hessian_parameter_name = None
        self.virial_parameter_name = None
        self.do_copy_at_file = None
        self.sparse_separate_file = None
        self.rnd_seed = None

    def teach(self, at_file,
              gp_file='gp_new.xml',
              config_type_parameter_name='config_type',
              energy_parameter_name='energy',
              force_parameter_name='force',
              hessian_parameter_name='hessian',
              virial_parameter_name='virial',
              do_copy_at_file=True,
              sparse_separate_file=True,
              rnd_seed=-1,
              ):
        """

        Args:
            at_file (str): XYZ file with teaching configurations
            gp_file (str): Output XML file
            config_type_parameter_name (str): Identifier of property determining the type of input data in the at_file
            do_copy_at_file (bool): Do copy the at_file into the GAP XML file (should be set to False for NetCDF input).
            energy_parameter_name (str): Name of energy property in the at_file that describes the data
            force_parameter_name (str): Name of force property in the at_file that describes the data
            hessian_parameter_name (str): Name of hessian property in the at_file that describes the data
            rnd_seed (int): Random seed.
            sparse_separate_file (bool): Save sparse coordinates data in separate file
            virial_parameter_name (str): Name of virial property in the at_file that describes the data
        """

        self.at_file = at_file
        self.gp_file = gp_file
        self.config_type_parameter_name = config_type_parameter_name
        self.energy_parameter_name = energy_parameter_name
        self.force_parameter_name = force_parameter_name
        self.hessian_parameter_name = hessian_parameter_name
        self.virial_parameter_name = virial_parameter_name
        self.do_copy_at_file = do_copy_at_file
        self.sparse_separate_file = sparse_separate_file
        self.rnd_seed = rnd_seed

        print(self.cmd())

    def cmd(self):
        cmd = ' \{}'.format(os.linesep).join([
            'teach_sparse',
            # 'at_file={}'.format(at_file),
            # 'gp_file={}'.format(gp_file),
            'gap={}'.format(self.gap),
            '{}'.format(self)
        ])

        return cmd


if __name__ == '__main__':
    # g = GAP(gap=None, default_sigma=None, asdjkhdsg=234, asdgfj=234)
    # g.teach(at_file='input.xyz')

    # !teach_sparse at_file=tmp.xyz \
    # gap={distance_Nb order=2 \
    #                  cutoff=5.0 \
    #                  covariance_type=ARD_SE \
    #                  theta_uniform=1.0 \
    #                  n_sparse=15 \
    #                  delta=1.0:\
    #      distance_Nb order=3 \
    #                  cutoff=4.0 \
    #                  covariance_type=ARD_SE \
    #                  theta_uniform=1.0 \
    #                  n_sparse=50 \
    #                  delta=0.004} \
    # e0=-29.716948405885105 \
    # default_sigma={0.005 0.5 0.0 0.0} \
    # do_copy_at_file=F sparse_separate_file=F \
    # gp_file=gap_3b.xml

    from gap.descriptors import *

    distance_Nb(order=2, cutoff=5.0)

    g = GAP(gap=None, default_sigma=None)

    g.teach(at_file='tmp.xyz', gp_file='gap_3b.xml', )
