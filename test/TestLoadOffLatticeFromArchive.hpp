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

#ifndef TESTLOADOFFLATTICEFROMARCHIVE_HPP_
#define TESTLOADOFFLATTICEFROMARCHIVE_HPP_

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
#include "VolumeTrackingModifier.hpp"
#include "CellAncestorWriter.hpp"
#include "CellBasedSimulationArchiver.hpp"

//Elegans specific headers
#include "GlobalParameterStruct.hpp"
#include "DTCMovementModel.hpp"
#include "LeaderCellBoundaryCondition.hpp"
#include "RepulsionForceSizeCorrected.hpp"
#include "RepulsionForce.hpp"
#include "GonadArmDataOutput.hpp"
#include "CellTrackingOutput.hpp"
#include "OocyteFatedCellApoptosis.hpp"
#include "Fertilisation.hpp"
#include "StatechartCellCycleModel.hpp"
#include "ElegansDevStatechartCellCycleModel.hpp"
#include "FateUncoupledFromCycle.hpp"

/* 
* Tests loading a C. elegans simulation from a saved file (unarchiving). 
*
* startTime = the time at which the simulation being loaded finished
* endTime = the desired end time for this simulation
* directory = the directory name for the loaded simulation (under testoutput) -
* currently set to Baseline
*/

class TestLoadOffLatticeFromArchive : public AbstractCellBasedTestSuite
{

public:

  void TestUnarchive() throw(Exception){

    double startTime = 1.0;
    double endTime = 1.5;
    std::string directory = std::string("Baseline");
	
  	OffLatticeSimulation<3>* p_simulator = 
            CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load(directory, startTime);
  
  	p_simulator->SetEndTime(endTime);
  	p_simulator->Solve();

  };

};

#endif /* TESTLOADOFFLATTICEFROMARCHIVE_HPP_ */
