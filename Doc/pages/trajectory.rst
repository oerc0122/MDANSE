
.. _trajectory-converters:

Trajectory Converters
=====================

Below is the list of converters present in MDANSE. These converters allow for the
outputs of a variety of MD simulation packages to be used in
MDANSE by converting the various file formats to the MDT format. These trajectory
converters convert positions and other vital
information such as the unit cell parameters.

The MDANSE MD engine specific converters were developed to convert
specific trajectory formats outputted by those MD engines. If you are
unable to find the MD package you used or if you had dumped the
trajectory to a format not support by that particular
converter, we recommend trying the ASE, MDAnalysis or MDTraj converters
which may support the file format you need to convert. These general
trajectory converters utilise their respective libraries to read MD
trajectory files and convert them to MDT.

- ASE
- CASTEP
- DCD
- CP2K
- Forcite
- Discover
- DFTB
- DMol
- DL_POLY
- Gromacs
- LAMMPS
- MDAnalysis
- MDTraj
- VASP
- CHARMM
- NAMD
- XPLOR
