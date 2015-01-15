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

/*Tests larval development of the C. elegans gonad, which forms behind an explicitly modelled Distal Tip Cell*/
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
#include "CellBasedSimulationArchiver.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>

//Elegans specific headers
#include "DTCMovementModel.hpp"
#include "LeaderCellBoundaryCondition.hpp"
#include "RepulsionForceSizeCorrected.hpp"
#include "RepulsionForce.hpp"
#include "StatechartCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GonadArmDataOutput.hpp"
#include "CellTrackingOutput.hpp"
#include "GlobalParameterStruct.hpp"
#include "RandomCellKillerByType.hpp"
#include "Fertilisation.hpp"
#include "CellBasedEventHandler.hpp"


class TestElegansWithLeaderCell : public AbstractCellBasedTestSuite
{

public:

    void TestLarvalDevelopment() throw(Exception){
    
        //---------------------------Set parameters----------------------------------
        int count;
        // Display each command-line argument.
        std::cout << "\nCommand-line arguments:" << *(CommandLineArguments::Instance()->p_argc) << "\n";
        for (count = 0; count < *(CommandLineArguments::Instance()->p_argc); count++)
            std::cout << "  argv[" << count << "]   "
            << (*(CommandLineArguments::Instance()->p_argv))[count] << "\n";
        
        std::string myParametersFilesDirectory = "/Users/kathryn/Documents/PhDYear2Term1/Chaste/projects/kathryna/data/";
        GlobalParameterStruct::Instance((*(CommandLineArguments::Instance()->p_argv))[1], myParametersFilesDirectory);
        double startingLength = 32;
    
        //------------Initialise simulation by setting up nodes----------------------     
        std::vector<Node<3>*> nodes;
        int index = 0;
        for (double i = GlobalParameterStruct::Instance()->GetParameter(0); i >= 0; i--){
            Node<3>* myNode;
            myNode = new Node<3>(index, false, 0, -GlobalParameterStruct::Instance()->GetParameter(7), 1.8*i);
            nodes.push_back(myNode);
            index++;
        }
        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 25);
    
        //---------------------Initialise cells from the nodes------------------------ 
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StatechartCellCycleModelSerializable, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<3> cell_population(*p_mesh, cells);
    
        //Setup the properties of the cells        
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End(); ++cell_iter)
        {
            Node<3>* node = cell_population.GetNode(cell_population.GetLocationIndexUsingCell(*cell_iter));
            if (node->GetIndex() == 0){
                cell_iter->GetCellData()->SetItem("IsDTC", 1.0);
                cell_iter->SetCellProliferativeType(p_diff_type);
            }
            else{
                cell_iter->GetCellData()->SetItem("IsDTC", 0.0);
            }
            cell_iter->GetCellData()->SetItem("DistanceAwayFromDTC", 0.0);
            cell_iter->GetCellData()->SetItem("RowNumber", 0.0);
            cell_iter->GetCellData()->SetItem("Radius", GlobalParameterStruct::Instance()->GetParameter(9));
            cell_iter->GetCellData()->SetItem("MaxRadius", GlobalParameterStruct::Instance()->GetParameter(9));     
            cell_iter->GetCellData()->SetItem("Differentiation_Sperm", 0.0);
            cell_iter->GetCellData()->SetItem("Differentiation_Oocyte", 0.0);
            cell_iter->GetCellData()->SetItem("SpermFated", 0.0);
            cell_iter->GetCellData()->SetItem("OocyteFated", 0.0);
            cell_iter->GetCellData()->SetItem("Fertilised", 0.0);
            cell_iter->GetCellData()->SetItem("PreviousClosestPointIndex", -1);
            cell_iter->GetCellData()->SetItem("ArrestedFor", 0.0);
        }
        cell_population.SetAbsoluteMovementThreshold(2.5);
        cell_population.SetDampingConstantNormal(GlobalParameterStruct::Instance()->GetParameter(12));
        cell_population.SetUseVariableRadii(true);

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.AddCellWriter<CellAncestorWriter>();


        //-------------------------Setup simulation----------------------------
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(GlobalParameterStruct::Instance()->GetDirectory().c_str());
        simulator.SetSamplingTimestepMultiple(GlobalParameterStruct::Instance()->GetParameter(36));
        simulator.SetDt(1.0/GlobalParameterStruct::Instance()->GetParameter(36));
        simulator.SetEndTime(GlobalParameterStruct::Instance()->GetParameter(35));
        
    
        //-------------------------------Forces-------------------------------- 
        MAKE_PTR(RepulsionForceSizeCorrected<3>, p_force);
        p_force->SetMeinekeSpringStiffness(GlobalParameterStruct::Instance()->GetParameter(13));
        simulator.AddForce(p_force);
    
    
        //-----------------------Boundary Condition Modifiers--------------------  
        double MidlinePointSpacing = 2;
        std::vector< c_vector<double, 3> > MidlinePointCollection;
        std::vector< int > MidlinePointTypes;
        
        //initialise midline points for the gonad at start of simulation 
        double newPointZ = 0;
        c_vector<double, 3> aPoint;
        while (newPointZ < startingLength){
            aPoint[0] = 0;
            aPoint[1] = -GlobalParameterStruct::Instance()->GetParameter(7);
            aPoint[2] = newPointZ;
            MidlinePointCollection.push_back(aPoint);
            MidlinePointTypes.push_back(0);
            newPointZ += MidlinePointSpacing;
        }
        
        MAKE_PTR_ARGS(DTCMovementModel<3>, dtcMovement, (MidlinePointCollection, MidlinePointTypes, aPoint, MidlinePointSpacing, false, false, 0.0));
        simulator.AddSimulationModifier(dtcMovement);
        MAKE_PTR_ARGS(LeaderCellBoundaryCondition<3>, boundaryCondition, (&cell_population, dtcMovement, GlobalParameterStruct::Instance()->GetParameter(28)));
        simulator.AddCellPopulationBoundaryCondition(boundaryCondition);
    
    
        //----------------------------Other Modifiers------------------------------ 
        MAKE_PTR(VolumeTrackingModifier<3>, volumeTrackingForContactInhibition);
        simulator.AddSimulationModifier(volumeTrackingForContactInhibition);
    
    
        //------------------------------Data output--------------------------------  
        MAKE_PTR_ARGS(GonadArmDataOutput<3>, dataRecording, (GlobalParameterStruct::Instance()->GetParameter(36)) );
        simulator.AddSimulationModifier(dataRecording);
        //MAKE_PTR_ARGS(CellTrackingOutput<3>, positionRecording, (GlobalParameterStruct::Instance()->GetParameter(36),5) );
        //simulator.AddSimulationModifier(positionRecording);
    
    
        //------------------------------Cell killers-------------------------------         
        double lengthOfSpermatheca = 20.0;
        MAKE_PTR_ARGS(Fertilisation<3>, removalByFertilisation, (&cell_population, lengthOfSpermatheca));
        simulator.AddCellKiller(removalByFertilisation);
        MAKE_PTR_ARGS(RandomCellKillerByType<3>, removalByApoptosis, (&cell_population, GlobalParameterStruct::Instance()->GetParameter(21)));
        simulator.AddCellKiller(removalByApoptosis);
    
    
        //----------------------------------Run-------------------------------------  
        simulator.Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);
    
        //----------------------------Garbage collection----------------------------  
        for (unsigned i = 0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    
    }

};

#endif /* TESTELEGANSGERMLINE_HPP_ */
