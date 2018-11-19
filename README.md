# Steak-Cooking
Code to simulate a steak using Flory-Rehner Theory and implementing local shrinkage scheme.

1. Set a parameter specifying struct, defined by 'setdefaultparams_steak.m' to defined physical and numerical parameters such as heat capacity and time step.

e.g. P=setdefaultparams_steak();

2. Run the main file associate with the boundary conditions of interest

e.g. steak_sim_dirichlet(P) to impose dirichlet Temperature boundary conditions as if the meat were pan seared.

Note 

This code has seen several revisions and is not actively maintained, but a snapshot of the code written during the summer of 2018 as part of the Research Experiences for Undergraduate (REU) program at James Madison University, funded by National Science Foundation grant NSF-DMS 1560151. Moreover, the code does not represent its final revision, but a snapshot after much progress and the majority of coding was complete.
