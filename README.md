# BESMo - Bedload Scenario Model
This is a bedload transport model that was created to:
- Simulate many different sediment transport scenarios in parallel (Monte-Carlo style)
- Simulate extreme sediment supply conditions (large, episodic supply events)

The two main applications were
1. Simulating episodic sediment in the flume lab at UBC. Paper:
Mueller & Hassan (2018)
2. Simulate different options in management of a potential San Clemente Dam removal
Paper: in preperation
AGU 2020 talk:

Matlab versions used were 2015b and 2018b

## Learning about sediment transport estimations
Garry Parker's ebook is most helpful. Many of the solutions are based on that.


## Current issues
* Only Wilcock and Crow works. Ashida Michue is not tested thoroughly, but files are there..
* Reaches are not implemented in this version, but sediment feed can be implemented at any node.
* There are a bunch of unused/old scripts in the repo. Best follow the execution, halt the program if you want to implement something. Set breakpoints etc.
* I refactored some of the variable names a while back which might be a bit confusing (some unused old variable names might still be active). 

## Cleaned up
I deleted a bunch of configurations for the 508 projects of Alex and Conor.
The current configuration executes a version of the permuted pulse sequence from Maria's experiments.

## Run the model
### Generate a runconfig file
Generate a .mat file (binary workspace snapshot) using a configuration script: `genconf_508Flume.m`
Result is a `runconfig_*.mat` file that has all variables stored needed for execution.
If there are multiple scenarios of for example grain size distribution or feedrate, 
one `runconfig_*.mat` file is generated per scenario. Each scenario can additionally be
executed for a number of runs each. For example I iterate through grain size distributions
per scenario and then iterate through event frequencies per run.
Using scenario-iterations is recommended over using run-iterations.

### Run the model
In all `run_*.m` scripts the configuration is referred to by name (Set at the top of `genconf_*.m`).
If there are multiple scenarios they have to be iterated through by their file endings.
The iteration should be copied from the `runconfig_*.mat` file generation at the end of the  `genconf_*.m` script

### Run analysis

## Analyzing model output
Starting an analysis is done similarly to executing a simulation. Iterate through 
`runconfig_*.mat`, generate folders and do the plotting.

This is a bit rough as of writing (2020-11-25). The analysis scripts are horrible, sorry for that.
Problems: Handling heterogeneous model data, normalizing output
Things that work well: loading the data


