# SHIK-potential
SHIK interaction potential files for LAMMPS

Introduction:

The SHIK (Sundarararaman, Huang, Ispas, Kob) potential files for oxide glasses are a system of self consistent interaction potential files for molecular dynamics simulations. These potentials were parameterized to accurately reproduce various properties over a wide range of compositions and currently include the network formers silicon, boron and aluminium and the network modifiers sodium, potassium, lithium, magnesium, and calcium.

Using the SHIK potential in LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator):

Here, you will find the potential files (pair_SHIK_wolf.cpp and pair_SHIK_wolf.h) along with several sample input scripts for lammps (input_SHIK_equil.NVT) to use the potential. To use the SHIK potential, you will need to put the potential files (pair_SHIK_wolf.cpp and pair_SHIK_wolf.h) in the src directory of your lammps installation and recompile. The pair style that you would need to use for the SHIK potential is shik/wolf which should be followed by the short range cutoff and the long range cutoff respectively. pair_style shik/wolf 8.0 10.0 The pair coefficient format is as follows pair_coeff atom_type atom_type A B C D Gamma_shortrange Gamma_longrange (Please see reference 1 for details of functional form and details of parameters)

Sample input files for running simple MD simulation of different glass system have been included here for reference.
File "input_SHIK_equil_sio2.NVT" is for silica glass; "input_SHIK_equil_na2osio2.NVT" and "input_SHIK_equil_caoosio2.NVT" are for sodium silicate and calcium silicate respectively; "input_SHIK_equil_na2ob2o3.NVT" and "input_SHIK_equil_caob2o3.NVT" are for sodium borate and calcium borate respectively.

These potential files have been currently tested upto the LAMMPS Mar 2018 package but should work in the current version of LAMMPS too (2020). Potential parameters for the various systems can be found in the references below. If you have any questions, please contact Siddharth Sundararaman at ssundararaman@lbl.gov or siddharth.s410@gmail.com or Yueh-Ting (Tim) Shih at ytshih@mail.ntut.edu.tw.

Update notes: 

09192022: The file pair_SHIK_wolf_2022 is compatible with the 2022 LAMMPS version 

References:

1. For pure silica glass: Siddharth Sundararaman, Liping Huang, Simona Ispas, and Walter Kob. "New optimization scheme to obtain interaction potentials for oxide glasses." The Journal of Chemical Physics 148, no. 19 (2018): 194504.
2. For Alkali and Alkaline-earth silicate glasses: Siddharth Sundararaman, Liping Huang, Simona Ispas, and Walter Kob. "New interaction potentials for alkali and alkaline-earth aluminosilicate glasses." The Journal of chemical physics 150, no. 15 (2019): 154505.
3. For Alkali borate glasses: Siddharth Sundararaman, Liping Huang, Simona Ispas, and Walter Kob. "New interaction potentials for borate glasses with mixed network formers." The Journal of chemical physics 152, no. 10 (2020): 104501.
4. For Alkaline-earth silicate and borate glasses: Yueh-Ting Shih, Siddharth Sundararaman, Simona Ispas, and Liping Huang. "New interaction potentials for alkaline earth silicate and borate glasses." Journal of non-crystalline solids 565 (2021): 120853.
