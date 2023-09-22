# Statistical approaches for the integration of external controls in a cystic fibrosis clinical trial- a simulation and an application code

This repository contains the code that corresponds with the manuscript Statistical approaches for the integration of external controls in a cystic fibrosis clinical trial- a simulation and an application by Mark N. Warden, Sonya L. Heltshe, Noah Simon, Stephen J. Mooney, Nicole Mayer-Hamblett, and Amalia S. Magaret. It was submitted to the American Journal of Epidemiology.

Details regarding the structure of the simulation to assess the three selected statistical methods as well as information about the data application can be found in the manuscript.

## Folders
- **Application Code**: Folder that contains the code that implements the methods in the real dataset
-- **Epic_Optimize_Implementation.R**- R script that implements the data application
-- **JAGS_model_A3.bug** - Jags/Winbugs model file used to specific the propsensity score based power prior model for the application
-- **JAGS_model_A3_opt_only.bug** - Jags/Winbugs model file used to specific the propsensity score based power prior model that uses  ONLY the Optimize trial data (i.e. does not use historical controls)
-- **JAGS_model_Epic_Optimize_CP.bug** - Jags/Winbugs model file used to specific the commensurate prior model for the application
-- **JAGS_model_Epic_Optimize_CP_opt_only.bug** - Jags/Winbugs model file used to specific the commensurate prior model that uses ONLY the Optimize trial data (i.e. does not use historical controls)
- **Output**: Folder that stores simulation results
-- **Parameter set results** - Folder to store simulation results (empty in this repository)
-- **Merging_node_results_simulation.R** - R script to merge node results from the simulation

## Functions (R scripts) and Files
- **D1.R** - Trial data (historical and prospective) generating function based on inputed parameters
- **A0.R** - Fits a propensity score for being in the prospective trial and saves it as a column variable
- **A1.R** - Implements inverse probability weighting and estimates a treatment effect
- **A2.R** - Implements Bayesian commensurate prior heirarchical modeling and estimates a treatment effect
- **A3.R** - Implements Bayesian propensity score based power priors and estimates a treatment effect
- **L1.R** - Loops over one parameter set, implementing A0-A3 and saving the treatment effect estimates
- **L2.R** - Prepares for parallel programming and then Loops L1 over ALL parameter sets in parallel
- **S1.R** - Summarizes the simulation results over 1 parameter set
- **S2.R** - Loops S1 over all parameter sets and saves the outcome
- **autoexec_simulation.R** - Written similar to SAS programs, prepares the workspace for the simulation
- **Cluster taskmaster linear programming.R** - Used for linear programming to find the optimal way to break up the computing job into tasks and distribute them to the various cores in the AWS cloud computing servers.
- **JAGS_model_A2.bug** - the Jags/Winbugs model specification for the commensurate prior method
- **JAGS_model_A3.bug** - the Jags/Winbugs model specification for the propensity score based power prior method


## Documentation
Scripts are commented with a header to help understand their function.

## Contact information
Mark Warden can be reached at mark.n.warden@gmail.com

## License

MIT
Copyright (c) [2023] [Mark Warden]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [dill]: <https://github.com/joemccann/dillinger>

