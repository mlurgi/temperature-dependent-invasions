# temperature-dependent-invasions

Start here!

Author: Dr Miguel Lurgi
Lecturer in Biosciences (Computational Ecology)
Computational Ecology Lab - Department of Biosciences
Swansea University, UK
 
and

Centre for Biodiversity Theory and Modelling
Theoretical and Experimental Ecology Station, CNRS, France

Date Created: 19-12-2019

Copyright (c) Miguel Lurgi, 2019
Email: miguel.lurgi@swansea.ac.uk

This repository contains the source code used to produce the results presented in our paper entitled:

'Warming indirectly incrases invasion success in food webs'

Uploaded to BioRXiv. https://doi.org/10.1101/2020.07.20.211516'

In this file I describe the contents of the archive (i.e. files contained here) and provide a brief description of what source file does when executed.

To execute the scripts successfully, place all the scripts contained in this repository within the same folder.

The main script file is:

1.- 'invasion_experiments_for_temp_dep.r' : Run this script to perform the temperature dependent invasion experiments as described in the paper.

Complementary scripts comprise:

2.- 'temp-dep-food-web.r' : This script implements the network dynamical model of food webs based on the bioenergetic model and using the niche model for constructing realistic food webs. It incorporates temperature scaling of different parameters as proposed in Binzer et al. (2016) Interactive effects of warming, eutrophication and size structure: Impacts on biodiversity and food-web structure. Global Change Biology 22, 220–227, doi: 10.1111/gcb.13086

3.- 'SEMs.r' : This script contains all the statistical analyses performed as part of the structural equation models, including their fits.

4.- 'calculate_size_effects.r' : This script contains all the statistical analyses performed to obtain the sizes of the effects of the invasions on the food webs across a temperature gradient.

In addition to the scripts mentioned above, another source code file containing helper functions is provided to support calculations of food web properties: utils.r. Also two data files are provided to facilitate running of the scripts concerning the analysis of the output from the experiments (to save the time of running all the simulations): 'data-for-sems.rda' and 'output-processing.RData'.

I hope you enjoy the code!

If you have any questions / run into any issues, please contact me at miguel.lurgi@swansea.ac.uk






