/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "CellTrackingOutput.hpp"
#include "GlobalParameterStruct.hpp"

template<unsigned DIM>
CellTrackingOutput<DIM>::CellTrackingOutput(int samplingInterval, int cellIdInterval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputFile(NULL)
{
  interval = samplingInterval;
  idInterval = cellIdInterval;
  OutputFileHandler rOutputFileHandler(GlobalParameterStruct::Instance()->GetDirectory(), false);
  OutputFile = rOutputFileHandler.OpenOutputFile("TrackingData.txt");
}

template<unsigned DIM>
CellTrackingOutput<DIM>::~CellTrackingOutput()
{
}

//At each timestep, loop through all cells and output IDs and positions along gonad arm.
template<unsigned DIM>
void CellTrackingOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
  //std::cout << SimulationTime::Instance()->GetTime() << "\t" << "Tracking" << "\t" << time(NULL) << std::endl;

  if (SimulationTime::Instance()->GetTimeStepsElapsed() % interval == 0){

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
      cell_iter != rCellPopulation.End(); ++cell_iter)
    {
      Node<DIM>* node = rCellPopulation.GetNode(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
      unsigned id = cell_iter->GetCellId();

      if (id%idInterval == 0 && id!=0){
        c_vector<double, DIM> location = node->rGetLocation();
        *OutputFile << SimulationTime::Instance()->GetTime() << "\t" << id << "\t" <<  location[0] << "\t" <<  location[1] << "\t" <<  location[2] << "\n";
      }
    }
   
  OutputFile->flush();

  }

  if(SimulationTime::Instance()->IsFinished()){
    //Close output file
    OutputFile->close();
  }
}

template<unsigned DIM>
void CellTrackingOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}


template<unsigned DIM>
void CellTrackingOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  // No parameters to output
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellTrackingOutput<1>;
template class CellTrackingOutput<2>;
template class CellTrackingOutput<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellTrackingOutput)
