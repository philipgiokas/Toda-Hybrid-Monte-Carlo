# Toda-Hybrid-Monte-Carlo

This code is for the implementation of a two dimensional affine Toda quantum field theory hybrid Monte Carlo Simulation and the fitting of resulting wall correlation functions to their respective parameters (the principal one of interest one being mass), also the statistical errors of the fitted parameters are calculated.

The files toda_hybrid.hpp/cpp build the necessary class for the simulation, parameter fits and parameter statistical error estimations.

The files mass_curve_fit.hpp/cpp provide an interface to the GSL, GNU Scientific Library, parameter fitting library.

The main file runs a Monte Carlo simualtion whose parameters are set from the main arguments.

The code has been tested with the G2 affine algebra and good agreement has been found with the previous work of Prof G.Watts of King's College London.  Further testing with the D4 affine algebra is in process.
