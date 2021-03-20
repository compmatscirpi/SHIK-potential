# SHIK-potential
SHIK interaction potential files for LAMMPS

Introduction:

The SHIK (Sundarararaman, Huang, Ispas, Kob) potential files for oxide glasses are a system of self consistent interaction potential files for molecular dynamics simulations. These potentials were parameterized to accurately reproduce various properties over a wide range of compositions and currently include the network formers silicon, boron and aluminium and the network modifiers sodium, potassium, lithium, magnesium, and calcium.

Using the SHIK potential in LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator):

Here, you will find the potential files (pair_SHIK_wolf.cpp and pair_SHIK_wolf.h) along with a sample input script for lammps (input_SHIK_equil.NVE) to use the potential. To use the SHIK potential, you will need to put the potential files (pair_SHIK_wolf.cpp and pair_SHIK_wolf.h) in the src directory of your lammps installation and recompile. The pair style that you would need to use for the SHIK potential is shik/wolf which should be followed by the short range cutoff and the long range cutoff respectively. pair_style shik/wolf 8.0 10.0 The pair coefficient format is as follows pair_coeff atom_type atom_type A B C D Gamma_shortrange Gamma_longrange (Please see reference 1 for details of functional form and details of parameters)

A sample input file (input_SHIK_equil.NVE) for running a simple MD simulation of silica glass has been included here for reference.

These potential files have been currently tested upto the LAMMPS Mar 2018 package but should work in the current version of LAMMPS too (2020). Potential parameters for the various systems can be found in the references below. If you have any questions, please contact me at ssundararaman@lbl.gov or siddharth.s410@gmail.com.

References:

For pure silica glass: Sundararaman Siddharth, Liping Huang, Simona Ispas, and Walter Kob. "New optimization scheme to obtain interaction potentials for oxide glasses." The Journal of Chemical Physics 148, no. 19 (2018): 194504.
For Alkali and Alkaline-earth silicate glasses: Sundararaman Siddharth, Liping Huang, Simona Ispas, and Walter Kob. "New interaction potentials for alkali and alkaline-earth aluminosilicate glasses." The Journal of chemical physics 150, no. 15 (2019): 154505.
For Alkali borate glasses: Sundararaman Siddharth, Liping Huang, Simona Ispas, and Walter Kob. "New interaction potentials for borate glasses with mixed network formers." The Journal of chemical physics 152, no. 10 (2020): 104501.
For Alkaline-earth silicate and borate glasses: Yueh-Ting Shih, Sundararaman Siddharth, Simona Ispas, and Liping Huang. "New interaction potentials for alkaline earth silicate and borate glasses." Journal of non-crystalline solids (Under review).
