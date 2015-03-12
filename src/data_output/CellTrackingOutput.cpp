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

#include "CellTrackingOutput.hpp"
#include "GlobalParameterStruct.hpp"


//Constructor 
template<unsigned DIM>
CellTrackingOutput<DIM>::CellTrackingOutput(int samplingInterval, int cellIdInterval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputFile(NULL),
      mSamplingInterval(samplingInterval),
      mCellIdInterval(cellIdInterval){}


//Destructor
template<unsigned DIM>
CellTrackingOutput<DIM>::~CellTrackingOutput(){};


//Prepare for solve by opening an output file in the appropriate directory.
template<unsigned DIM>
void CellTrackingOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
std::string outputDirectory){
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  OutputFile = rOutputFileHandler.OpenOutputFile("TrackingData.txt");
};



/*
* Getters for the private member variables
*/
template<unsigned DIM>
int CellTrackingOutput<DIM>::GetSamplingInterval() const
{
  return mSamplingInterval;
};
template<unsigned DIM>
int CellTrackingOutput<DIM>::GetCellIdInterval() const
{
  return mCellIdInterval;
};



/*
* Actual data recording function. 
* At each timestep, if it's time for a new datapoint to be recorded, loop through all the cells and then output IDs
* and positions for those with IDs divisible by mCellIdInterval.
*/
template<unsigned DIM>
void CellTrackingOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
  
  //If it's an output timestep
  if (SimulationTime::Instance()->GetTimeStepsElapsed() % GetSamplingInterval() == 0){

    //Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End(); ++cell_iter)
    {

      //Get cell ID from the cell's associated node
      Node<DIM>* node = rCellPopulation.GetNode(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
      unsigned id = cell_iter->GetCellId();

      //If the cell meets the cellIDInterval requirement and is not ID = 0, output the position.
      //Prohibition on the ID=0 cell is because this cell is the DTC, not a germ cell in C. elegans germ line
      // simulations. That condition can be removed.  
      if (id % GetCellIdInterval() == 0 && id != 0){
        c_vector<double, DIM> location = node->rGetLocation();
        *OutputFile << SimulationTime::Instance()->GetTime() << "\t" << id 
        << "\t" <<  location[0] << "\t" <<  location[1] << "\t" <<  location[2] << "\n";
      }
    }
  
  //Output the data to file immediately, so if the simulation crashes you will have some data saved
  OutputFile->flush();
  }

  //If simulation is finished, close the output file.
  if(SimulationTime::Instance()->IsFinished()){
    OutputFile->close();
  }
}


//Output properties of this cell data recorder to file.
template<unsigned DIM>
void CellTrackingOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<SamplePositionsEveryXTimesteps>" << mSamplingInterval << "</SamplePositionsEveryXTimesteps>\n";
  *rParamsFile << "\t\t\t<SampleEveryNthCellPosition>" << mCellIdInterval << "</SampleEveryNthCellPosition>\n";
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
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
