# LAMMPS-Membrane

## membrane_polymerization.cpp
C++ script which generates a LAMMPS script `membrane.in` and data file `membrane.dat`. Parameters are read in from `params_membrane.txt`.

## membrane.dat
[Data file](http://lammps.sandia.gov/doc/read_data.html) which stores the initial positions of the lipids in a bilayer membrane. Includes positions, bonds, and angles. Also includes the initial configuration of a specified number of rigid filaments which will polymerize.

## membrane.in
LAMMPS script modeling a bilipid membrane in a periodic simulation box filled with DPD fluid. `membrane.dat` must be in the same directory to run. This script runs with the version of LAMMPS compiled Aug. 10, 2015.

_Outputs_:
* `xpos.out` - Average membrane X position over time
* `ypos.out` - Average membrane unwrapped Y position over time
* `zpos.out` - Average membrane unwrapped Z position over time
* `membrane.lammpstrj` - LAMMPS trajectory file for visualization. Loadable in, e.g., [VMD](http://www.ks.uiuc.edu/Research/vmd/). Does not include fluid
