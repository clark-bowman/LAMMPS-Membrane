/*
LAMMPS Script Maker - Bilipid Membrane with Polymerizing Filaments, DPD Fluid (clark_bowman@brown.edu)
For LAMMPS version Aug. 10, 2015
*/

#include <fstream>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265359

double randDouble()
{
    return double(rand()) / double(RAND_MAX);
}

int main()
{
    srand(5000);
    using namespace std;

    double box_width, membrane_center, intra_spacing, filament_spacing, box_length, membrane_force, temperature, v_lipid_bond, v_filament_bond, v_lipid_angle, v_lattice_density, v_monomer_density, v_timestep, v_polymer_rmin, v_polymer_prob;
    int target_lipids, tail_length, num_tails, num_heads, n_filaments_oned, filament_length, v_run_length, v_dump_interval;

    // Read parameters from text file
    {
        ifstream file_in ("params_membrane.txt");
        file_in >> box_width; // Y, Z widths of simulation box (LJ units)
        file_in >> box_length; // X length of simulation box
        file_in >> membrane_center; // Starting X position of membrane
        file_in >> intra_spacing; // Distance in between atoms of lipid tails
        file_in >> target_lipids; // Target number of lipids in the membrane
        file_in >> tail_length; // Length of each lipid tail
        file_in >> num_tails; // Number of hydrophobic tails per lipid
        file_in >> num_heads; // Number of hydrophilic heads per lipid
        file_in >> n_filaments_oned; // Number of polymerizing filaments (1D) - e.g., 3 means a 3x3 grid of filaments
        file_in >> filament_length; // Length of filaments
        file_in >> filament_spacing; // Distance in between atoms of filaments
        file_in >> membrane_force; // Magnitude of force opposing forward motion of membrane (applied to each atom)
        file_in >> temperature; // LJ temperature of system
        file_in >> v_lipid_bond; // Harmonic bond strength in lipid
        file_in >> v_filament_bond; // Harmonic bond strength in filaments
        file_in >> v_lipid_angle; // Harmonic angle strength in lipid tails
        file_in >> v_lattice_density; // Numerical density of fluid initialization
        file_in >> v_monomer_density; // Numerical density of free monomers
        file_in >> v_timestep; // LJ time per simulation step
        file_in >> v_run_length; // Number of steps to rain data-gathering simulation
        file_in >> v_dump_interval; // Dump visualization data every # steps
        file_in >> v_polymer_rmin; // Distance within which polymerization can occur
        file_in >> v_polymer_prob; // Probability to create bonds with free monomers at each step
    }

    double inter_spacing = box_width / sqrt(double(target_lipids / 2));
    int n_iter = int(ceil(box_width / inter_spacing));
    int n_filaments = n_filaments_oned * n_filaments_oned;

    // Generate membrane initial condition
    {
        ofstream fileOut ("membrane.dat");
        fileOut.precision(5);
        fileOut << "# Membrane data file" << endl << endl;

        // Header with number of atoms, bonds, angles, dihedrals, and simulation information (dihedrals are not used here)
        fileOut << target_lipids * (num_heads + num_tails * tail_length) + n_filaments * filament_length << " atoms" << endl;
        fileOut << target_lipids * (num_tails * tail_length + num_heads - 1) + n_filaments * (filament_length - 1) << " bonds" << endl;
        fileOut << target_lipids * num_tails * (tail_length - 1) + n_filaments * (filament_length - 2) << " angles" << endl;
        fileOut << 0 << " dihedrals" << endl << endl;
        fileOut << 8 << " atom types" << endl << "2 bond types" << endl;
        fileOut << "2 angle types" << endl << "0 dihedral types" << endl << endl;
        fileOut << 0. << " " << box_length << " xlo xhi" << endl;
        fileOut << 0. << " " << box_width << " ylo yhi" << endl;
        fileOut << 0. << " " << box_width << " zlo zhi" << endl << endl;

        // Extra bonds and angles are needed for polymerization
        fileOut << "2 extra bond per atom" << endl << "3 extra angle per atom" << endl << endl;

        // Masses for each particle type
        fileOut << "Masses" << endl << endl;
        fileOut << "1 1.0" << endl << "2 1.0" << endl << "3 1.0" << endl << "4 1.0" << endl << "5 1.0" << endl << "6 1.0" << endl << "7 1.0" << endl << "8 1.0" << endl;

        // Loop to create atoms
        fileOut << endl << "Atoms" << endl << endl;

        int atom_no = 1;
        int mol_no = 1;
        int atom_no2 = 1;
        double ang, xpos, myx, myy, myz, yoff, zoff, ang_perp, ang_perp_inc;
        double ypos = 0.;

        // Looping over all Y positions for lipids...
        for (int yctr = 0; yctr < n_iter && mol_no <= target_lipids + n_filaments; ypos += inter_spacing)
        {
            double zpos = 0.;

            // Looping over all Z positions for lipids...
            for (int zctr = 0; zctr < n_iter && mol_no <= target_lipids + n_filaments; zpos += inter_spacing)
            {
                // Start with X position at the membrane center, then subtract off the length of a tail
                xpos = membrane_center - (double(tail_length) + 0.5) * intra_spacing;
                yoff = 0.;
                zoff = 0.;

                // For each head particle per lipid...
                for (int hctr = 0; hctr < num_heads; hctr++)
                {
                    // Create the head atom
                    fileOut << atom_no + hctr << " " << mol_no << " 1 " << xpos << " " << ypos + yoff << " " << zpos + zoff << endl;

                    // If enough tails haven't been made yet...
                    if (hctr < num_tails)
                    {
                        // Initialize the tail location
                        myx = xpos;
                        myy = ypos + yoff;
                        myz = zpos + zoff;

                        // Create each tail atom, moving along the length of the tail in the X direction
                        for (int j = 0; j < tail_length; j++)
                        {
                            myx += intra_spacing;
                            fileOut << atom_no + num_heads + hctr * tail_length + j << " " << mol_no << " 2 " << myx << " " << myy << " " << myz << endl;
                        }
                    }

                    // Choose a random angle and offset the next head (if present) in that direction
                    ang = randDouble() * 2 * PI;
                    yoff += cos(ang) * intra_spacing;
                    zoff += sin(ang) * intra_spacing;
                }

                // Update atom count for number of atoms in one lipid, increment molecule count by 1
                atom_no += num_heads + num_tails * tail_length;
                mol_no++;

                // Now do the same thing for the opposite layer of the membrane:
                // Start with X position at the membrane center, then add the length of a tail
                xpos = membrane_center + (double(tail_length) + 0.5) * intra_spacing;
                yoff = 0.;
                zoff = 0.;
                for (int hctr = 0; hctr < num_heads; hctr++)
                {
                    fileOut << atom_no + hctr << " " << mol_no << " 1 " << xpos << " " << ypos + yoff << " " << zpos + zoff << endl;
                    if (hctr < num_tails)
                    {
                        myx = xpos;
                        myy = ypos + yoff;
                        myz = zpos + zoff;
                        for (int j = 0; j < tail_length; j++)
                        {
                            myx -= intra_spacing;
                            fileOut << atom_no + num_heads + hctr * tail_length + j << " " << mol_no << " 2 " << myx << " " << myy << " " << myz << endl;
                        }
                    }
                    ang = randDouble() * 2 * PI;
                    yoff += cos(ang) * intra_spacing;
                    zoff += sin(ang) * intra_spacing;
                }
                atom_no += num_heads + num_tails * tail_length;
                mol_no++;
                zctr++;
            }
            yctr++;
        }

        // Test for shallowest angle at which filaments will fit
        for (ang_perp_inc = 0.; true; ang_perp_inc += 0.001)
        {
            ang_perp = 0.;
            double tempx = 0.1;
            for (int j = 1; j < filament_length; j++)
            {
                tempx += filament_spacing * cos(ang_perp);
                ang_perp = min(ang_perp + ang_perp_inc, PI * 0.5);
            }
            if (tempx < membrane_center - double(tail_length + 1) * intra_spacing)
                break;
        }

        // Looping over the desired number of filaments...
        for (int i = 0; i < n_filaments; i++)
        {
            // Start on the left side, at evenly distributed X and Y positions
            xpos = 0.1;
            ypos = box_width * (0.5 + (i % n_filaments_oned)) / double(n_filaments_oned);
            double zpos = box_width * (0.5 + (i / n_filaments_oned)) / double(n_filaments_oned);

            // At which angle the filament will bend if it is too long
            ang = randDouble() * 2. * PI;

            // Current perpendicular angle
            ang_perp = 0.;

            // Create the pinned atom at the base
            fileOut << atom_no << " " << mol_no << " 5 " << xpos << " " << ypos << " " << zpos << endl;
            atom_no++;

            // Add the rest of the atoms along the filament according to the angles
            for (int j = 1; j < filament_length; j++)
            {
                xpos += filament_spacing * cos(ang_perp);
                ypos += filament_spacing * sin(ang_perp) * cos(ang);
                zpos += filament_spacing * sin(ang_perp) * sin(ang);
                fileOut << atom_no << " " << mol_no << " " << 4 + int(j == 1) + 3 * int(j == filament_length - 1) << " " << xpos << " " << ypos << " " << zpos << endl;
                atom_no++;
                ang_perp = min(ang_perp + ang_perp_inc, PI * 0.5);
            }
            mol_no++;
        }

        // Loop to create bonds
        fileOut << endl << "Bonds" << endl << endl;
        int bond_no = 1;
        atom_no = 1;

        // For each lipid...
        for (int ctr = 0; ctr < target_lipids; ctr++)
        {
            // The heads are bonded
            for (int i = 0; i < num_heads - 1; i++)
            {
                fileOut << bond_no << " 1 " << atom_no + i << " " << atom_no + i + 1 << endl;
                bond_no++;
            }
            for (int i = 0; i < num_tails; i++)
            {
                // Each tail is bonded to a head
                fileOut << bond_no << " 1 " << atom_no + i << " " << atom_no + num_heads + i * tail_length << endl;
                bond_no++;

                // Atoms within a tail are bonded
                for (int j = 0; j < tail_length - 1; j++)
                {
                    fileOut << bond_no << " 1 " << atom_no + num_heads + i * tail_length + j << " " << atom_no + num_heads + i * tail_length + j + 1 << endl;
                    bond_no++;
                }
            }
            atom_no += (num_heads + num_tails * tail_length);
        }

        // For each filament...
        for (int i = 0; i < n_filaments; i++)
        {
            // Atoms within a filament are bonded
            for (int j = 0; j < filament_length - 1; j++)
            {
                fileOut << bond_no << " 2 " << atom_no + j << " " << atom_no + j + 1 << endl;
                bond_no++;
            }
            atom_no += filament_length;
        }

        // Loop to create angles
        fileOut << endl << "Angles" << endl << endl;
        int ang_no = 1;
        atom_no = 1;

        // For each lipid...
        for (int ctr = 0; ctr < target_lipids; ctr++)
        {
            // Add angles along the length of each tail (no angle potential is imposed among heads, even if there are 3 or more)
            for (int i = 0; i < num_tails; i++)
            {
                fileOut << ang_no << " 1 " << atom_no + i << " " << atom_no + num_heads + i * tail_length << " " << atom_no + num_heads + i * tail_length + 1 << endl;
                ang_no++;
                for (int j = 0; j < tail_length - 2; j++)
                {
                    fileOut << ang_no << " 1 " << atom_no + num_heads + i * tail_length + j << " " << atom_no + num_heads + i * tail_length + j + 1 << " " << atom_no + num_heads + i * tail_length + j + 2 << endl;
                    ang_no++;
                }
            }
            atom_no += (num_heads + num_tails * tail_length);
        }

        // For each filament...
        for (int i = 0; i < n_filaments; i++)
        {
            // Add angles along the length of each filament
            for (int j = 0; j < filament_length - 2; j++)
            {
                fileOut << ang_no << " 2 " << atom_no + j << " " << atom_no + j + 1 << " " << atom_no + j + 2 << endl;
                ang_no++;
            }

            // Deprecated: long-distance angles to impose more rigid filaments
            /*for (int j = 0; j < filament_length - 10; j++)
            {
                fileOut << ang_no << " 3 " << atom_no + j << " " << atom_no + j + 5 << " " << atom_no + j + 10 << endl;
                ang_no++;
            }*/
            atom_no += filament_length;
        }
    }

    // Create LAMMPS script
    {
        ofstream fileOut ("membrane.in");
        fileOut.precision(5);

        fileOut << "# LAMMPS Script - Bilipid Membrane with Polymerizing Filaments, DPD Fluid (clark_bowman@brown.edu)" << endl;
        fileOut << "# LAMMPS version Aug. 10, 2015" << endl << endl << endl << endl << endl << endl;

        fileOut << "# SEC: This section initializes the simulation." << endl << endl;
        fileOut << "# 	Define units, atom style, log path, and neighbor settings; and read configuration data for the polymer." << endl;
        fileOut << "# 	Configuration data for the polymer in polymer.dat must be generated separately." << endl;
        fileOut << "# 	Final line communicates ghost data and is necessary for DPD parallelizing. On some versions of LAMMPS, use instead `communicate single vel yes'" << endl << endl;
        fileOut << "units lj" << endl;
        fileOut << "atom_style molecular" << endl;
        fileOut << "log membrane.log" << endl;
        fileOut << "read_data membrane.dat" << endl;
        fileOut << "neighbor 0.3 bin" << endl;
        fileOut << "neigh_modify delay 3" << endl;
        fileOut << "comm_modify vel yes" << endl << endl;

        fileOut << "# 	Define bond, angle, pairwise interactions." << endl;
        fileOut << "# 	For pairwise interactions: type 1 is lipid head, 2 is lipid tail, 3 is water," << endl;
        fileOut << "#   4 is filament, 5 is pinned filament, 6 is free monomers, 7 is terminal monomers, 8 is added filament." << endl;
        fileOut << "#   43872 is a temperature seed; change for different randomness." << endl << endl;
        fileOut << "bond_style harmonic" << endl;
        fileOut << "bond_coeff 1 " << v_lipid_bond << " " << intra_spacing << endl;
        fileOut << "bond_coeff 2 " << v_filament_bond << " " << filament_spacing << endl;
        fileOut << "angle_style cosine/delta" << endl;
        fileOut << "angle_coeff 1 " << v_lipid_angle << " 180.0" << endl;
        fileOut << "angle_coeff 2 100.0 180.0" << endl;
        fileOut << "dihedral_style none" << endl << endl;
        fileOut << "pair_style dpd " << temperature << " 1.0 43872" << endl;
        fileOut << "pair_coeff * * 25.0 4.5 1.0" << endl;
        fileOut << "pair_coeff 1 3*8 35.0 4.5 1.0" << endl;
        fileOut << "pair_coeff 2 3*8 75.0 20.0 1.0" << endl;
        fileOut << "pair_coeff 1 2 50.0 9.0 1.0" << endl;
        fileOut << "pair_coeff 2 2 15.0 4.5 1.0" << endl << endl;

        fileOut << "# SEC: This section initializes the geometry." << endl << endl;
        fileOut << "#   Define safe regions where fluid may be placed. This is the union of two regions outside the membrane." << endl;
        fileOut << "# 	The simulation box is " << box_length << " x " << box_width << " x " << box_width << "." << endl;
        fileOut << "#   Also creates diffuse monomers on the filament side." << endl << endl;
        fileOut << "region sim_box block 0 " << box_length << " -2 " << box_width + 2. << " -2 " << box_width + 2. << " units box side in" << endl << endl;
        fileOut << "region safe1 block 0.1 " << membrane_center - double(tail_length + 1) * intra_spacing << " 0.1 " << box_width << " 0.1 " << box_width << " units box" << endl;
        fileOut << "region safe2 block " << membrane_center + double(tail_length + 1) * intra_spacing << " " << box_length << " 0.1 " << box_width << " 0.1 " << box_width << " units box" << endl;
        fileOut << "region safe union 2 safe1 safe2" << endl;
        fileOut << "lattice fcc " << v_lattice_density << endl;
        fileOut << "create_atoms 3 region safe" << endl;
        fileOut << "lattice fcc " << v_monomer_density << endl;
        fileOut << "create_atoms 6 region safe1" << endl << endl;

        fileOut << "# SEC: This section defines LAMMPS groups and computes that will be used in the simulation." << endl << endl;
        fileOut << "# 	Groups are named indicatively of their membership." << endl;
        fileOut << "# 	Computes include the x, y, and z positions of the center of mass of the membrane." << endl << endl;
        fileOut << "group fluid type 3" << endl;
        fileOut << "group polymer type 4" << endl;
        fileOut << "group pinned type 5" << endl;
        fileOut << "group membrane subtract all fluid polymer pinned" << endl;
        fileOut << "group bounded subtract all fluid pinned" << endl;
        fileOut << "group free subtract all pinned" << endl;
        fileOut << "group bondable type 6 7" << endl;
        fileOut << "compute xpos membrane reduce ave x" << endl;
        fileOut << "compute ypos_a membrane property/atom yu" << endl;
        fileOut << "compute zpos_a membrane property/atom zu" << endl;
        fileOut << "compute ypos membrane reduce ave c_ypos_a" << endl;
        fileOut << "compute zpos membrane reduce ave c_zpos_a" << endl << endl;

        fileOut << "# SEC: This section initializes the particle velocities." << endl << endl;
        fileOut << "velocity free create " << temperature << " 21456" << endl << endl;

        fileOut << "# SEC: This section defines fixes to impose forces in the simulation." << endl << endl;
        fileOut << "# 	NVE integration, with limit for initialization." << endl;
        fileOut << "# 	Fix 2 prevents bounded particles (e.g., monomers) from wrapping around in the X direction." << endl;
        fileOut << "# 	Fix 3 pins the bases of the filaments." << endl << endl;
        fileOut << "fix 1 all nve/limit 0.05" << endl;
        fileOut << "fix 2 bounded wall/region sim_box lj126 0.001 0.3 1.0" << endl;
        fileOut << "fix 3 pinned setforce 0.0 0.0 0.0" << endl << endl;

        fileOut << "# SEC: This section runs the simulation." << endl << endl;
        fileOut << "# 	Simulation timestep (LJ time units)." << endl;
        fileOut << "timestep " << v_timestep << endl << endl;
        fileOut << "# 	How often to output thermo data and first run phase." << endl;
        fileOut << "thermo 1000" << endl;
        fileOut << "run 500" << endl << endl;

        fileOut << "# 	Release limit on integrator." << endl;
        fileOut << "unfix 1" << endl;
        fileOut << "fix 1 all nve" << endl;
        fileOut << "run 10000" << endl << endl;

        fileOut << "#   At this point, begin imposing the force on the membrane." << endl;
        fileOut << "# 	Fixes 4, 4a, 4b dump the computes at specified short intervals." << endl;
        fileOut << "# 	Dump 1 is a fluid-less simulation dump at the specified interval." << endl << endl;
        fileOut << "fix 5 membrane addforce -" << membrane_force << " 0 0" << endl;
        fileOut << "fix 4 membrane ave/time 10 25 250 c_xpos file xpos.out" << endl;
        fileOut << "fix 4a membrane ave/time 10 25 250 c_ypos file ypos.out" << endl;
        fileOut << "fix 4b membrane ave/time 10 25 250 c_zpos file zpos.out" << endl;
        fileOut << "fix 6 bondable bond/create 1 6 7 " << v_polymer_rmin << " 2 iparam 1 7 jparam 2 8 prob " << v_polymer_prob << " 23478 atype 2" << endl;
        fileOut << "dump 1 bounded atom " << v_dump_interval << " membrane.lammpstrj" << endl << endl;

        fileOut << "# 	Run data-collecting simulation." << endl;
        fileOut << "run " << v_run_length << endl;
    }
    return 0;
}
