# Random Velocity Generator for Rigid Water Models
## Purpose
Rigid water models (SPCE, TIP4P, etc.) are widely used in molecular simulation. Many molecular dynamics simulation softwares has their own functions that generates random initial velocities based on a normal distribution at simulation temperature. However, since rigid models only have six degrees of freedom (nine for flexible bodies), the initial velocities created by the above method will produce an initial temperature 1/3 higher than the simulation temperature. This will push the system away from equilibration. 

This script generates initial velocities according to six degrees of freedom. The translational and rotational velocities for oxygen and hydrogen atoms are created separately while keeping the total momentum of water molecule fixed. The strenching force between atoms in a molecule is also eliminated. 

## Dependencies
Numpy 

Pandas

## Example
The main functions for creating initial velocities are included in velocity.py. The lammps_water_example.py file presents an example of generating a SPCE water input file for lammps. The initial velocities created keep the system at 500 K. 
