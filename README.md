# ElegansGermline
This is an add-on project for the Chaste C++ biological modelling library https://chaste.cs.ox.ac.uk/trac. Chaste stands for Cancer Heart and Soft Tissue Environment. This particular project allows the user to run 3D mechanics simulations of the _C. elegans_ germ line, an experimental system used by biologists to research cell proliferation and differentiation. The ElegansGermline project includes code for enforcing a tube shaped boundary condition formed along the path of a leader cell, and support for cell cycle models controlled by a statechart, both of which could be useful in other modelling contexts.

## Installation
Before using the ElegansGermline code, it is necessary to install Chaste on your machine. A guide to Chaste installation can be found at https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/InstallGuide, covering Ubuntu, other Linux systems and Mac OS X. Installing Chaste on Windows 8.1 is also possible, but the process is experimental and is still being improved. The current Windows install guide is available at https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/Windows, and troubleshooting notes based on Katwell's experiences are included in the docs folder. 

Before proceeding, it is useful to run the Chaste test suite to check that the installation was successful. From the main Chaste directory, type:
    scons
in the command line to run all tests. This may take several hours, so for a quicker check try:
    scons cell_based/test/tutorial/TestCellBasedDemoTutorial.hpp

Next, place the contents of this repository in the projects folder, which is located under the main Chaste directory. You should now have _Chaste/projects/ElegansGermline/_. The add-on material is now available to use (see instructions below).

## Usage
The ElegansGermline project can be used in two ways. Either you can run our existing model of the _C. elegans_ germ line (provided in the test folder), or you can use the individual source files like a library to write your own model. For more information on how to write your own models, see this tutorial including helpful instructions for writing Chaste tests https://chaste.cs.ox.ac.uk/trac/wiki/CellBasedChastePractical; see also the more detailed description of the source files below. 

To simply run our existing model, type the following in the command line from the main chaste directory:
    scons co=1 b=GccOpt projects/kathryna/test/TestElegansGermline.hpp
This compiles the model file called _TestElegansGermline.hpp_. _co=1_ specifies "compile but do not run", while _b=GccOpt_ requests a high level of optimisation for speed (recommended). An executable called _TestElegansGermlineRunner_ will be created in the directory _Chaste/projects/ElegansGermline/build/optimised_. In that directory type:
    ./TestElegansGermlineRunner "Baseline.txt"
Baseline being a parameter input file, found in the data directory. The model should then execute.

## About parameter files
Since our model has a large number of parameters and options, we have endeavoured to allow them all to be set in one place - a parameter input file taken in by the model as a command line argument. Parameter files must be placed in _Chaste/projects/ElegansGermline/data_; several examples are already there for you. The format of the file is:
1. On line one - the name of the output directory.
2. Subsequent lines - a parameter of type double, a tab character, a comment (optional, must not contain tab or newline), a newline.
The comments in the existing parameter files tell you the parameters that are expected and in what order. Providing parameters in a different order will generate garbage. The comment is just for user readability and does not tell the computer how to interpret each double provided.

## Visualising the data
Output is placed in your testoutput directory. By default, this is: _tmp/<your-username>_/testoutput/_. A different output directory can be specified by setting the environment variable CHASTE_TEST_OUTPUT. In the testoutput directory should be a folder with the name given in the parameter file, in the example "Baseline". Inside is a text file "GonadData.txt" and a folder results_from_time_0 containing a number of .vtu files and one .pvd file.

GonadData.txt contains tab delimited measurements of various germ line properties, recorded every simulated hour. The columns are:
1. Hours since start of simulation
2. Organ length
3. Effective cell cycle length, accounting for arrests
4. Sperm count
5. Proliferative cell count
6. Cell death rate
7. Total cell count
8. Micron distance from DTC to last proliferative cell
9. Micron distance from DTC to first cell in meiosis
10. G1 cell count
11. S phase cell count
12. G2 cell count
13. M phase cell count
14. Meiotic S phase cell count
15. First cell row in which two meiotic cells appear
16. Last cell row in which a proliferative cell appears
A summary of some of this data can be plotted in R using the script "plotGonadData.R", included in the RScripts directory. Before trying to run this script, open it in a text editor and set the path to your output, following the directions provided in the comments.

The .vtu and .pvd files are 3D snapshots of the simulation, captured each hour. The .pvd file provides a summary of the .vtu files, each of which represents a single frame. The .pvd file can be opened using the visualisation software paraview, available from http://www.paraview.org/. On opening, click apply, then select "3D Glyphs" as the data representation (initial value is Surface). In the new Glyph's properties menu, select the following options:
- Glyph Type = Sphere
- Radius = 1
- Active Attributes, Scalars = Radius
- Scaling, Scale mode = scalar
- Scaling, Scale factor = 1
- Masking, Glyph Mode = All Points
After clicking Apply, the screen should show each cell as a sphere of the appropriate radius (you may need to adjust the zoom). Cells can be colored acording to their properties using the Glyph's Coloring dropdown menu, and subsets of cells sharing a certain property can be selected by using a Threshold Filter (icon of a cube containing a T shape, toward the top left).


## Source file descriptions
This project contains the following files:


## Reproducing our figures
If all goes to plan, our paper concerning this germ line model will be published in the very near future (watch this space for the link). To reproduce the specific figures in that paper, the following steps should be taken.

