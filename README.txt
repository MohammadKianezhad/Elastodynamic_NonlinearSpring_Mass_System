# Elasto-dynamic Simulation of Nonlinear Spring-Mass System with Explicit Time Integration

This code has been developed for a second-round interview of a PhD position at ETH Zurich.

Notes:
   - The formulation used in this code is taken from book "DYNAMICS OF
     STRUCTURES" written by Anil K. Chopra.
   - in this program the units are in meter (m), Second (s), Kilogram(Kg), and
     Newtons (N).
   - in this program "no*" stands for "number of *".
   - in this program "elmn" stands for "element".
   - in this program "dof" stands for "degree of freedom".
   - in this program "u" stands for "displacement".
   - in this program "ud" stands for "velocity" as the first derivative of "u".
   - in this program "udd" stands for "accelaration" as the Second derivative of "u".

This program consists of various functions, which are briefly introduced below:
   - Main: the main body of the program
   - Mesh_Square: this function produces a mesh with "rank" iterations of unitcellin X- and Y-direction.
   - Mesh_Triangle: this function produces a mesh with "rank" iterations of triangle unitcell.
   - dt_max_CFLcondition: this function calculates the largest time step with which the simulation will be stable.
   - MakeStiffness: this function calculates the stiffness matrix of the whole system.
   - Solver: function for 0-1 method Solver.
   - Stress: this function calculates the stress of each element.
   - VideoMaker: This function makes a video of the structure's deformation throughout the analysis.
   - VideoMaker_Stress: this function makes a video of the structure's deformation along with its stress contour throughout the analysis.
   - Draw: this function draws the structure in its initial configuration.
   - Drawdef: this function draws the structure in its deformed configuration.
   - Drawdef_Stress: his function draws the structure in its deformed configuration along with the stress contour on elements.