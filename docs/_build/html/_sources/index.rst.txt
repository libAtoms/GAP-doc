.. GAP Wrapper documentation master file, created by
   sphinx-quickstart on Thu Jan 10 17:45:59 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: quippy

QUIP and quippy documentation
=============================

The ``QUIP`` package (`GitHub <https://github.com/libAtoms/QUIP>`_) is a
collection of software tools to carry out molecular dynamics
simulations. It implements a variety of interatomic potentials and
tight binding quantum mechanics, and is also able to call external
packages, and serve as plugins to other software such as `LAMMPS
<http://lammps.sandia.gov>`_, `CP2K <http://www.cp2k.org>`_ and also
the python framework `ASE <https://wiki.fysik.dtu.dk/ase>`_.  Various
hybrid combinations are also supported in the style of QM/MM, with a
particular focus on materials systems such as metals and
semiconductors.

``quippy`` is a Python interface to ``QUIP`` that provides deep access to
most of the Fortran types and routines. The quippy interface is principally
maintained by `James Kermode <http://www.warwick.ac.uk/jrkermode>`_.

Long term support of the package is ensured by:
 - Noam Bernstein (Naval Research Laboratory)
 - Gabor Csanyi (University of Cambridge)
 - James Kermode (University of Warwick)

Portions of this code were written by: Albert Bartok-Partay, Livia
Bartok-Partay, Federico Bianchini, Anke Butenuth, Marco Caccin,
Silvia Cereda, Gabor Csanyi, Alessio Comisso, Tom Daff, ST John,
Chiara Gattinoni, Gianpietro Moras, James Kermode, Letif Mones,
Alan Nichol, David Packwood, Lars Pastewka, Giovanni Peralta, Ivan
Solt, Oliver Strickson, Wojciech Szlachta, Csilla Varnai, Steven
Winfield.

Copyright 2006-2016.

Most of the publicly available version is released under the GNU
General Public license, version 2, with some portions in the public
domain.

Overview of ``libAtoms`` and ``QUIP``
-------------------------------------

- The `libAtoms <http://www.libatoms.org>`_ package is a software
  library written in Fortran 95 for the purposes of carrying out
  molecular dynamics simulations.

- The ``QUIP`` (**QU**\ antum mechanics and **I**\ nteratomic
  **P**\ otentials) package, built on top of ``libAtoms``, implements a
  wide variety of interatomic potentials and tight binding quantum
  mechanics, and is also able to call external packages.

- Various hybrid combinations are also supported in the style of
  QM/MM, including `Learn on the Fly` scheme [LOTF]_

- `quippy <http://www.jrkermode.co.uk/quippy>`_ is a Python interface
  to libAtoms and QUIP.

Features
--------

The following interatomic potentials are presently coded or linked in QUIP:

* EAM (fcc metals)
* Fanourgakis-Xantheas (water)
* Finnis-Sinclair (bcc metals)
* Flikkema-Bromley
* GAP (Gaussian Approximation Potentials: general many-body)
* Guggenheim-!McGlashan
* Brenner (carbon)
* OpenKIM (general interface)
* Lennard-Jones
* Morse
* Partridge-Schwenke (water monomer)
* Stillinger-Weber (carbon, silicon, germanium)
* SiMEAM (silicon)
* Sutton-Chen
* Tangney-Scandolo (silica, titania etc)
* Tersoff (silicon, carbon)

The following tight-binding functional forms and parametrisations are implemented:

* Bowler
* DFTB
* GSP
* NRL-TB

The following external packages can be called:

* CASTEP
* VASP
* CP2K
* ASAP
* ASE (recent version, 3.11+, recommended)
* Molpro


.. toctree::
   :maxdepth: 4
   :caption: Contents:

   introduction.rst
   gap.rst
   descriptors.rst
   potentials.rst
   tutorials.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
----------

.. [LOTF] Csányi, G., Albaret, T., Payne, M., & De Vita,
   A. 'Learn on the Fly': A Hybrid Classical and Quantum-Mechanical
   Molecular Dynamics Simulation. Physical Review Letters,
   93(17), 175503. (2004) http://prl.aps.org/abstract/PRL/v93/i17/e175503>