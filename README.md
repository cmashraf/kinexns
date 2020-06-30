## kinexns (An Open Source Python Package for Chemical Kinetic Modeling, Sensitivity Analysis and Optimization)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/cmashraf/kinexns.svg?branch=master)](https://travis-ci.org/cmashraf/kinexns)
[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/cmashraf/kinexns/badges/quality-score.png?b=master)](https://scrutinizer-ci.com/g/cmashraf/kinexns/?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/cmashraf/kinexns/badge.svg?branch=master)](https://coveralls.io/github/cmashraf/kinexns?branch=master)

Chemical kinetic models are usually composed of a set of stiff ordinary differential equations (ODEs), where a number of rate parameters have to be estimated to fit with experimental observations. This presents a twofold challenge—firstly, to solve the stiff set of ODEs accurately and efficiently; and secondly, to estimate the model parameters precisely by optimizing the objective function. In recent years, stochastic optimization methods for parameter estimation have gained popularity over the classical optimization methods as the former do not require a reasonable initial guess and have the capability to escape local minima. Here, we developed an open source python package kinexns to efficiently solve the kinetic model using CVode solver, perform sensitivity analysis to determine important model parameters, and optimize the model parameters by using the different stochastic methods. The different algorithms we considered are: Monte Carlo (MC), Latin-Hypercube Sampling (LHS), Maximum Likelihood Estimation (MLE), Markov-Chain Monte-Carlo (MCMC), Shuffled Complex Evolution Algorithm (SCE-UA), Simulated Annealing (SA), RObust Parameter Estimation (ROPE), Artificial Bee Colony (ABC), Fitness Scaled Chaotic Artificial Bee Colony (FSCABC), and Dynamically Dimensioned Search algorithm (DDS). 

This package can generate the set of ODEs for a kinetic model either through parsing Chemkin input files or by parsing though some user defined files that contains the reactions in SMILES format. Also, this package can parse Chemkin thermo files or get the thermodynamic information through a user provided file to calculate the rate constants for the reactions.

Complete example of handling **Chemkin input files**, building model, solving the set of ODES, performing sensitivity analysis and model optimization procedures can be found at 
**[Notebooks/example_chemkin](https://github.com/cmashraf/kinexns/tree/master/Notebooks/example_chemkin)**

Complete example of handling **user deined input files**, building model, solving the set of ODES, performing sensitivity analysis and model optimization procedures can be found at 
**[Notebooks/example_xylose](https://github.com/cmashraf/kinexns/tree/master/Notebooks/example_xylose)**

-------
### Software dependencies and license information

**Programming language:**
Python version 3.7.7 (https://www.python.org)

**Python packages needed:**
**kinexns** uses basic python packages like NumPy, Pandas etc. Also, to perform tasks like solving odes, performing sensitivity analysis and model optimization, kinexns uses other open source python packages like assimulo, SALib and spotpy.

To efficiently solve the system of stiff ODEs, kinexns uses **[assimulo](https://jmodelica.org/assimulo/)** through which it can access the **CVOde class** which provides a connection to the **[Sundials](https://computation.llnl.gov/casc/sundials/main.html)** solver CVode. CVOde is a highly efficient stiff ODE solver developed in Lawrence Livermore National Lab.

To perform sensitivity analysis, kinexns uses **[SALib](https://salib.readthedocs.io/en/latest/)** package that uses Sobol Sensitivity Analysis.

Finaly for parameter optimization, kinexns uses **[spotpy](https://pypi.org/project/spotpy/)**, which is another open source python package that makes various stochastic algorithms available.

The ODE solver that we used for our research is a modified version of DDASAC that is unfortunately not open source.  We chose this solver because it performed the best on the stiff set of ODEs in this model, but future users can modify the code (by replacing our `ddasac_utils.py` module) to use other solvers, such as those in the python package scipy.integrate.

**License information:**
The MIT License (MIT)
