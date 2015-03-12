/*
Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TESTELEGANSGERMLINE_HPP_
#define TESTELEGANSGERMLINE_HPP_

//Chaste and system headers
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellAncestorWriter.hpp"
#include "RandomNumberGenerator.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <ctime>

//Elegans specific headers
#include "GlobalParameterStruct.hpp"                // parameter storage and read-in from file
#include "DTCMovementModel.hpp"                     // leader cell movement
#include "LeaderCellBoundaryCondition.hpp"          // boundary condition
#include "RepulsionForceSizeCorrected.hpp"          // force law
#include "GonadArmDataOutput.hpp"                   // data recording
#include "CellTrackingOutput.hpp"                   // cell tracking
#include "OocyteFatedCellApoptosis.hpp"             // apoptosis
#include "Fertilisation.hpp"                        // fertilisation
#include "StatechartCellCycleModel.hpp"             // statechart wrapper class
#include "ElegansDevStatechartCellCycleModel.hpp"   // elegans specific changes in cell cycle length
#include "FateUncoupledFromCycle.hpp"               // statechart model of cell behaviour 


/*
* Simulates the larval development and adult maintenance of the C. elegans gonad, which forms behind
* the moving Distal Tip Cell. The gonad is filled with dividing germ cells that behave according to a
* given statechart model.
* Run as ./TestElegansGermlineRunner "Baseline.txt"
*/

class TestElegansWithLeaderCell : public AbstractCellBasedTestSuite
{

public:

    void TestLarvalDevelopment() throw(Exception){

        //0) Seed random number generator with the system time----------------------

        int seed = (int)std::time(NULL);
        std::cout << seed << std::endl;
        RandomNumberGenerator::Instance()->Reseed(seed);      
    

        //1) Read in parameter set from the config file specified in the first command 
        //line argument---------------------------------------------------------------
        
        std::cout << std::endl << "Selected parameters file: " << 
                     (*(CommandLineArguments::Instance()->p_argv))[1] << std::endl;      
        std::string myParameterFilesDirectory = "./projects/ElegansGermline/data/";
        
        GlobalParameterStruct* parameters = GlobalParameterStruct::Instance();
        parameters -> ConfigureFromFile( (*(CommandLineArguments::Instance()->p_argv))[1], myParameterFilesDirectory);
        
        // HELPER CODE FOR PARAMETER SWEEPING! When uncommented, the program takes a second 
        // command line argument: the name of the results output directory. It also accepts
        // arbitrarily many subsequent pairs of command line arguments that can be used to make
        // small alterations to the parameter set. So for instance:
        //
        // ./TestElegansGermlineRunner "Baseline.txt" "MyOutput1"  5  3.4 
        //
        // would run with most parameters as specified in Baseline.txt; the output directory 
        // would be set to MyOutput1, and parameter number 5 would be reset to equal 3.4.
        //
        parameters->ResetDirectoryName( (*(CommandLineArguments::Instance()->p_argv))[2] );   
        std::cout << "New output directory name: " << (*(CommandLineArguments::Instance()->p_argv))[2] << std::endl;
        
        int nArgs = (*(CommandLineArguments::Instance()->p_argc));
        
        for (int i=3; i<nArgs-1; i+=2){
            parameters->ResetParameter(atof( (*(CommandLineArguments::Instance()->p_argv))[i] ) ,
                                       atof( (*(CommandLineArguments::Instance()->p_argv))[i+1] ));
        
            std::cout << "Param number: " << 
            atof( (*(CommandLineArguments::Instance()->p_argv))[i] ) <<
            " New value: " <<
            atof( (*(CommandLineArguments::Instance()->p_argv))[i+1] ) << std::endl;
        }

        //---------------------------------------------------------------------------
    


        // 2) Create a node corresponding to each starting cell----------------------

        std::vector< Node<3>* > nodes;
        int cellIndex = 0;
        for (double i = parameters->GetParameter(0); i >= 0; i--){ //parameters[0] is number of starting cells
            Node<3>* newNode;
            newNode = new Node<3>(cellIndex, false, 0, -parameters->GetParameter(7), 1.8*i); //cellID, _, x, y, z
            nodes.push_back(newNode);
            cellIndex++;
        }
        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>; //make a 3D mesh from the nodes
        p_mesh->ConstructNodesWithoutMesh(nodes, 25);    //specify a max interaction distance of 25, > 2x max cell radius
    
        //----------------------------------------------------------------------------



        // 3) Initialise some cells from the nodes, and select a cell cycle model-----

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        //choice of model is "ElegansDevStatechartCellCycleModel" with the statechart "FateUncoupledFromCycle":
        CellsGenerator< ElegansDevStatechartCellCycleModel < FateUncoupledFromCycle >, 3> cells_generator;
        
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<3> cell_population(*p_mesh, cells); // <- finally, we have a cell population!
    
        //----------------------------------------------------------------------------



        // 4) Setup the properties of the cells---------------------------------------       

        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End(); ++cell_iter)
        {
            Node<3>* node = cell_population.GetNode(cell_population.GetLocationIndexUsingCell(*cell_iter));
            if (node->GetIndex() == 0){
                cell_iter->GetCellData()->SetItem("IsDTC", 1.0);  //Make the cell with node index 0 the DTC
                cell_iter->SetCellProliferativeType(p_diff_type); //The DTC is terminally differentiated
            }
            else{
                cell_iter->GetCellData()->SetItem("IsDTC", 0.0);
            }
            cell_iter->GetCellData()->SetItem("DistanceAwayFromDTC", 0.0); //Set other cell data
            cell_iter->GetCellData()->SetItem("RowNumber", 0.0);
            cell_iter->GetCellData()->SetItem("Radius", parameters->GetParameter(9));  // parameters[9] = initial radius
            cell_iter->GetCellData()->SetItem("MaxRadius", parameters->GetParameter(9));   
            cell_iter->GetCellData()->SetItem("SpermFated", 0.0);
            cell_iter->GetCellData()->SetItem("OocyteFated", 0.0);  
            cell_iter->GetCellData()->SetItem("Differentiation_Sperm", 0.0);
            cell_iter->GetCellData()->SetItem("Differentiation_Oocyte", 0.0);
            cell_iter->GetCellData()->SetItem("Fertilised", 0.0);
            cell_iter->GetCellData()->SetItem("PreviousClosestPointIndex", -1);
            cell_iter->GetCellData()->SetItem("ArrestedFor", 0.0);
            cell_iter->GetCellData()->SetItem("InProximalArm", 1.0);
        }
        cell_population.SetAbsoluteMovementThreshold(2.5);                      // Max cell movement in one timestep
        cell_population.SetDampingConstantNormal(parameters->GetParameter(12)); // Baseline cell drag coefficient
        cell_population.SetUseVariableRadii(true);                              // Different sized cells will be used
        cell_population.SetCellAncestorsToLocationIndices();                    // Request cell lineage tracking
        cell_population.AddCellWriter<CellAncestorWriter>();

        //----------------------------------------------------------------------------



        // 5) Setup the simulation----------------------------------------------------

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(parameters->GetDirectory().c_str());    // Set output directory name
        simulator.SetSamplingTimestepMultiple(parameters->GetParameter(36)); // How frequently to output a snapshot
        simulator.SetDt(1.0/parameters->GetParameter(36));                   // Length of a timestep (parameters[36])
        simulator.SetEndTime(parameters->GetParameter(35));                  // End time (parameters[35])
        
        //----------------------------------------------------------------------------


    
        // 6) Add a force between cells-----------------------------------------------

        MAKE_PTR(RepulsionForceSizeCorrected<3>, p_force);
        p_force->SetMeinekeSpringStiffness(parameters->GetParameter(13));   //Set force strength (parameters[13])
        simulator.AddForce(p_force);
    
        //----------------------------------------------------------------------------


    
        // 7) Boundary Condition------------------------------------------------------

        double MidlinePointSpacing = 2;
        double initialGonadLength = 32;
        std::vector< c_vector<double, 3> > MidlinePointCollection; // Points on the DTC path
        std::vector< int > MidlinePointTypes;                      // Records whether a point is part of a straight or the turn
        
        //initialise the points on the gonad midline at start of the simulation: 
        double newPointZ = 0;
        c_vector<double, 3> aPoint;
        while (newPointZ < initialGonadLength){
            aPoint[0] = 0;
            aPoint[1] = -parameters->GetParameter(7);   //y coord of the middle of the proximal straight
            aPoint[2] = newPointZ;
            MidlinePointCollection.push_back(aPoint);
            MidlinePointTypes.push_back(0);             //0 codes for "part of the proximal straight"
            newPointZ += MidlinePointSpacing;
        }
        //add code that moves the DTC
        MAKE_PTR_ARGS(DTCMovementModel<3>, dtcMovement, (false, false, 0.0, MidlinePointCollection, MidlinePointTypes, aPoint, MidlinePointSpacing));
        simulator.AddSimulationModifier(dtcMovement);
        //add a leader cell based boundary condition
        MAKE_PTR_ARGS(LeaderCellBoundaryCondition<3>, boundaryCondition, (&cell_population, dtcMovement, parameters->GetParameter(28)) ); //Parameters[28]: initial gonad radius
        simulator.AddCellPopulationBoundaryCondition(boundaryCondition);
    
        //---------------------------------------------------------------------------


    
        // 8) Cell volume tracking, required to apply contact inhibition-------------

        MAKE_PTR(VolumeTrackingModifier<3>, volumeTrackingForContactInhibition);
        simulator.AddSimulationModifier(volumeTrackingForContactInhibition);

        //---------------------------------------------------------------------------

    
    
        // 9) Cell removal, by fertilization and apoptosis---------------------------

        double lengthOfOvulationRegion = 20.0; //How close to the gonad's proximal end must a cell be before it can be removed
        MAKE_PTR_ARGS(Fertilisation<3>, removalByFertilisation, (&cell_population, lengthOfOvulationRegion));
        simulator.AddCellKiller(removalByFertilisation);
        MAKE_PTR_ARGS(OocyteFatedCellApoptosis<3>, removalByApoptosis, (&cell_population, parameters->GetParameter(21))); // parameters[21] = cell death rate
        simulator.AddCellKiller(removalByApoptosis);

        //---------------------------------------------------------------------------

    
    
        // 10) Add some data output--------------------------------------------------

        MAKE_PTR_ARGS(GonadArmDataOutput<3>, dataRecording, (parameters->GetParameter(36))); // parameters[36] = timesteps per hour
        simulator.AddSimulationModifier(dataRecording);
        MAKE_PTR_ARGS(CellTrackingOutput<3>, positionRecording, (parameters->GetParameter(36), 1));
        simulator.AddSimulationModifier(positionRecording);
    
        //----------------------------------------------------------------------------

    

        // 11) Run simulation and save the final state--------------------------------
        
        simulator.Solve();
        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);
    
        //----------------------------------------------------------------------------



        // 12) Delete nodes-----------------------------------------------------------
        
        for (unsigned i = 0; i<nodes.size(); i++)
        {
            delete nodes[i];
        } 

        //----------------------------------------------------------------------------
    }
};

#endif /* TESTELEGANSGERMLINE_HPP_ */
