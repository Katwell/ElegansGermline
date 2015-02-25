# ElegansGermline
ElegansGermline is an add-on project for the Chaste C++ biological modelling library (https://chaste.cs.ox.ac.uk/trac). Chaste stands for Cancer Heart and Soft Tissue Environment and is a library developed at the University of Oxford. Amongst other things, it has been used to simulate the heart, colonic crypt and solid tumor growth. This particular add-on project allows the user to run 3D mechanics simulations of the _C. elegans_ germ line, an experimental system used by biologists to research the cell proliferation and fate decisions. The ElegansGermline project includes code for a tube shaped boundary condition that forms along the path of a leader cell, and support statechart cell cycle models. Both these features could be useful in other developmental biology modelling projects.

## Installation
Before using the ElegansGermline code, it is necessary to install Chaste on your machine. A guide to Chaste installation can be found at https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/InstallGuide, covering Ubuntu, other Linux systems and Mac OS X. Installing Chaste on Windows 8.1 is also possible, but the process is experimental and is still being improved. The current Windows install guide is available at https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/Windows, and troubleshooting notes based on our own experiences are included in the docs folder. 

Before proceeding, it is useful to run the Chaste test suite to check that the installation was successful. From your main Chaste directory, type:
```
    scons
```
in the command line to run all tests. This can take several hours, so for a quicker check try:
```
    scons cell_based/test/tutorial/TestCellBasedDemoTutorial.hpp
```    
to just run a short test of the main cell mechanics functionality. 

Next, place the contents of this repository in the projects folder, which is located under the main Chaste directory. You should now have _Chaste/projects/ElegansGermline/_. 

Finally, a small change to the Chaste library is required to be able to run statechart cell cycle models. Go to the folder _Chaste/cell_based/src/population/cell_cycle/_ and open the file _AbstractCellCycleModel.hpp_. Add the keyword "virtual" in front of the function "SetCell" (at around line 186). The code should now read: 
```
    virtual void SetCell(CellPtr pCell); 
```
At this point the add-on material in ElegansGermline should be available to use in Chaste (see instructions below).

## Usage
The ElegansGermline project can be used in two ways. Either you can run our existing model of the _C. elegans_ germ line (provided in the test folder), or you can use the individual source files as a library to develop your own model. For more information on how to write models using Chaste, take a look at this tutorial: https://chaste.cs.ox.ac.uk/trac/wiki/CellBasedChastePractical. You may also find helpful the more detailed description of each source file provided in the docs folder. 

To simply run our existing model, type the following in the command line from the main chaste directory:
```
    scons co=1 b=GccOpt projects/ElegansGermline/test/TestElegansGermline.hpp
```
This compiles the test file called _TestElegansGermline.hpp_. _co=1_ specifies "compile but do not run", while _b=GccOpt_ requests a high degree of optimisation for speed (recommended). An executable called _TestElegansGermlineRunner_ will be created in the directory _Chaste/projects/ElegansGermline/build/optimised_. From that directory type:
```
    ./TestElegansGermlineRunner "Baseline.txt"
```
to run a simulation. Baseline.txt here is a parameter input file, one of several example files provided in the data directory.

## About parameter files
As our model has a large number of parameters and options, we have endeavored to allow them all to be set in one place - a parameter input file that is read in by the model as a command line argument. Parameter files must be placed in the directory  _Chaste/projects/ElegansGermline/data_; several examples are already there for you. The format of the file is:
- Line 1: your choice of output directory name, followed by a newline character.
- Subsequent lines: a parameter of type double, a tab character, a comment (optional, must not contain tab or newline), then a newline character.

The comments in the parameter files provided explain which parameters our current model expects and in what order. Providing parameters in a different order will generate garbage. Note that the comment field is only for human-readability and does not tell the computer how to interpret each double provided; that simply depends on the order of the input.

## Visualising the data
Simulation output is placed in the testoutput directory. By default, this will be: _tmp/(YOUR USERNAME)/testoutput/_. A different output directory can be specified by setting the environment variable CHASTE_TEST_OUTPUT. In the testoutput directory should be a folder with the name specified in the parameter file; for the example above that would be _Baseline_. Inside that should be two text files: _GonadData.txt_ and _TrackingData.txt_, as well as a folder _results_from_time_0_ containing a number of .vtu files.

_GonadData.txt_ contains tab delimited measurements of various germ line properties, recorded each simulated hour. The columns of this file should be interpreted as:

1. Hours since the start of the simulation
2. Organ length
3. Effective cell cycle length, adding on time in arrest
4. Sperm count
5. Proliferative cell count
6. Cell death rate
7. Total cell count
8. Micron distance from the DTC to the last proliferative cell
9. Micron distance from the DTC to the first meiotic cell
10. G1 cell count
11. S phase cell count
12. G2 cell count
13. M phase cell count
14. Meiotic S phase cell count
15. Estimate of the first cell row in which two meiotic cells appear
16. Estimate of the last cell row in which a proliferative cell appears

A summary of some of this data can be plotted in R using the script _plotGonadData.R_, included in the RScripts directory. Before trying to run this script, open it in a text editor and set the path to your output, following the directions provided in the comments.

_TrackingData.txt_ outputs the position of cells over time. Its columns should be interpreted as:

1. Hours since the start of the simulation
2. Cell ID
3. x coordinate
4. y coordinate
5. z coordinate

We have not provided an R script to plot this data directly. However it is used in the script _GetTreePath.R_, which generates a tree of tracks for one cell and all its daughters, that can then be visualized in Paraview.

The .vtu files in _results_from_time_0_ are 3D snapshots of the simulation, captured each hour. These files can be opened using the visualisation software Paraview, available from http://www.paraview.org/. On opening, click apply in the left hand panel, then select "3D Glyphs" as the data representation (initial value is Surface). In the new Glyph's properties menu, choose the following options:

- Glyph Type = Sphere
- Radius = 1
- Active Attributes, Scalars = Radius
- Scaling, Scale mode = scalar
- Scaling, Scale factor = 1
- Masking, Glyph Mode = All Points

After clicking Apply, the screen should show each cell as a sphere of the appropriate radius (you may need to adjust the zoom). Cells can be colored acording to their properties using the Glyph's coloring dropdown menu, and subsets of cells sharing a common property can be selected by using a Threshold Filter (icon of a cube containing a T, toward the top left of the screen).

## Source file descriptions
This project contains the following source code files:

- _test/TestElegansGermline.hpp_
- _src/boundary_condition/DTCMovementModel.hpp(cpp)_
- _src/boundary_condition/LeaderCellBoundaryCondition.hpp(cpp)_
- _src/cell_removal/Fertilisation.hpp(cpp)_
- _src/cell_removal/RandomCellKillerByType.hpp(cpp)_
- _src/data_input/GlobalParameterStruct.hpp(cpp)_
- _src/data_output/CellTrackingOutput.hpp(cpp)_
- _src/data_output/GonadArmDataOutput.hpp(cpp)_
- _src/force_law/RepulsionForceSizeCorrected.hpp(cpp)_
- _src/statechart/AbstractStatechartCellCycleModel.hpp_
- _src/statechart/StatechartCellCycleModel.hpp_
- _src/statechart/StatechartInterface.hpp_
- _src/statechart/FateDecisionCoupledToCycle.hpp(cpp)_
- _src/statechart/FateUpcoupledFromCycle.hpp(cpp)_

A full description of each is given in the docs, in the file _SourceCodeDetails_. ElegansGermline also contains the following R scripts, with descriptions in comments at the top of each script:

- _RScripts/plotGonadData.R_
- _RScripts/plotGrowthRates.R_
- _RScripts/deathRateComparison.R_
- _RScripts/GetTreePath.R_

## Reproducing our results
If all goes to plan, a paper concerning this germ line model will soon be published (watch this space for the link). Instructions for reproducing the specific figures in that paper are provided in the docs, in the file _FigureInstructions_. The same file also explains a bit about the various parameter files included in the data directory.
