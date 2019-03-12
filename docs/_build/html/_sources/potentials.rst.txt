
Using Potentials
================

QUIP
----

QUantum mechanics and Interatomic Potentials (QUIP)

Potentials implemented in ``QUIP``
++++++++++++++++++++++++++++++++++

Classical interatomic potentials:

- BKS (silica)
- Brenner (carbon)
- EAM (fcc)
- Fanourgakis-Xantheas
- Finnis-Sinclair (bcc)
- Flikkema-Bromley
- GAP (general many-body)
- Guggenheim-McGlashan
- Lennard-Jones
- Morse
- Partridge-Schwenke (water monomer)
- Si-MEAM (silicon)
- Stillinger-Weber (carbon, silicon, germanium)
- Stillinger-Weber + Vashishta (silicon/silica interfaces)
- Sutton-Chen
- Tangney-Scandolo (silica, titania etc)
- Tersoff (silicon, carbon)

Plus several tight binding parameterisations (Bowler, DFTB, GSP,
NRL-TB, ...)

External packages:

- ``CASTEP``--- DFT, planewaves, ultrasoft pseudopotentials
- ``CP2K`` --- DFT, mixed Gaussian/planewave basis set, various pseudopotentials.
  Our ``CP2K`` Driver supports QM, MM and QM/MM.
- ``MOLPRO`` --- All electron quantum chemistry code. DFT, CCSD(T), MP2
- ``VASP`` --- DFT, planewaves, PAW or ultrasoft pseudopotentials
- Interface to `OpenKIM <http://www.openkim.org>`_ project
- Relatively easy to add new codes

``QUIP`` also has a full interface to the Atomic Simulation
Environment, `ASE <https://wiki.fysik.dtu.dk/ase>`_

- ``ASE`` adds support for several more codes e.g. ``ABINIT``, ``Elk``,
  ``Exciting``, ``GPAW``, ``SIESTA``, ...

- ``ASE`` also works with the `CMR <https://wiki.fysik.dtu.dk/cmr>`_ database system

Performing calculations
+++++++++++++++++++++++

* As well as preparing structures and post-processing results,
  ``quippy`` allows calculations to be run

* In ``QUIP`` and ``quippy``, all calculations are performed with a
  Potential object (very similar to the
  :class:`~ase.calculators.interface.Calculator` concept in ``ASE``)

* Types of potential

  - *Internal*: interatomic potential or tight binding
  - *External*: file-based communication with external code or callback-based communication with a Python function
  - Plus flexible combinations of other potentials

* *Internal* potentials use XML parameter strings
* *External* potentials use template parameter files

Creating a Potential
++++++++++++++++++++

Internal potential::

  >>> sw_pot = Potential('IP SW')

External potential::

  >>> castep = Potential('FilePot',
  ...                    command='./castep-driver.sh')

Driver script can be a shell script, an executable program using
``QUIP`` or a ``quippy`` script. It can even invoke code on a remote
machine.



.. automodule:: quippy.potential
    :synopsis: Evaluate interatomic potentials

.. autoclass:: Potential
    :members:
    :inherited-members:
    :show-inheritance:

.. autoclass:: ForceMixingPotential
    :members:
    :show-inheritance:

.. autoclass:: Minim
    :members:

.. autofunction:: force_test

GAP
---

Example:
  sw_pot = Potential('IP SW')

XML functions
+++++++++++++

.. automodule:: quippy.qpxml
   :synopsis: Functions for manipulating GAP xml
   :members: