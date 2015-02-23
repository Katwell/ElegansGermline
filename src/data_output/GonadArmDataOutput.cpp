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

#include "GonadArmDataOutput.hpp"
#include "GlobalParameterStruct.hpp"
#include "MeshBasedCellPopulation.hpp"

template<unsigned DIM>
GonadArmDataOutput<DIM>::GonadArmDataOutput(int samplingInterval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputPositionFile(NULL)
{
  interval = samplingInterval;
  OutputFileHandler rOutputFileHandler(GlobalParameterStruct::Instance()->GetDirectory(), false);
  OutputPositionFile = rOutputFileHandler.OpenOutputFile("GonadData.txt");
}

template<unsigned DIM>
GonadArmDataOutput<DIM>::~GonadArmDataOutput()
{
  //if (OutputPositionFile!=NULL){
    OutputPositionFile->close();
  //}
}

//At each timestep, loop through all cells and output IDs and positions along gonad arm.
template<unsigned DIM>
void GonadArmDataOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  if(SimulationTime::Instance()->GetTimeStepsElapsed()%interval==0){

    double cellCycleDuration = (rCellPopulation.Begin())->GetCellCycleModel()->GetSDuration()
                      +(rCellPopulation.Begin())->GetCellCycleModel()->GetG2Duration()
                      +(rCellPopulation.Begin())->GetCellCycleModel()->GetTransitCellG1Duration()
                      +(rCellPopulation.Begin())->GetCellCycleModel()->GetMDuration();

    double gonadLength=0;

    double lastProliferativeCell=0;
    double firstMeioticCell=DBL_MAX;
    
    int totalCells = 0;
    int spermCount = 0;
    int prolifCount = 0;
    int G1count=0;
    int Scount=0;
    int G2count = 0;
    int Mcount=0;
    int MeioticS=0;

    double TimeArrested=0;

    //Work out width of one cell row, based on separations:
    int numberOfMeasurements = 0;
    double meanSeparation = 0;
    double meioticCountArray[128] = {0};
    double mitoticCountArray[128] = {0};
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();++cell_iter)
    { //consider distal arm only
      if(cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC")<75.0 && cell_iter->GetCellData()->GetItem("IsDTC")==0.0){
        double vol = cell_iter->GetCellData()->GetItem("volume");
        double diam = 2*pow(vol*0.239,0.3333);
        meanSeparation += diam;
        numberOfMeasurements++;
      }
    }
    meanSeparation = meanSeparation/numberOfMeasurements;


    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();++cell_iter)
    {

      double dist = cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC");
      if(dist<150.0){
       double row = round(dist/meanSeparation);
       cell_iter->GetCellData()->SetItem("RowNumber",row);
       if(cell_iter->GetCellData()->GetItem("CellCyclePhase")<0){
        meioticCountArray[(int)row] = meioticCountArray[(int)row] + 1;
       }else{
        mitoticCountArray[(int)row] = mitoticCountArray[(int)row] + 1;
       }
      }

      //Establish gonad length
      if(cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC")>gonadLength){
        gonadLength = cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC");
      }

      //Count sperm 
      if(cell_iter->GetCellData()->GetItem("Differentiation_Sperm")==1.0){
        spermCount++;
      }

      //count proliferative cells and length of zone
      if(cell_iter->GetCellData()->GetItem("CellCyclePhase")!=-1.0){
        prolifCount++;
      if (cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC")>lastProliferativeCell){
        lastProliferativeCell = cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC");
        }
      }else{
        if (cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC")<firstMeioticCell){
        firstMeioticCell = cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC");
        }
      }

      //Count number in each phase
      if(cell_iter->GetCellData()->GetItem("CellCyclePhase")==1.0){
          G1count++;
          TimeArrested+=cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(cell_iter->GetCellData()->GetItem("CellCyclePhase")==2.0){
        Scount++;
        TimeArrested += cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(cell_iter->GetCellData()->GetItem("CellCyclePhase")==3.0){
        G2count++;
          TimeArrested+=cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(cell_iter->GetCellData()->GetItem("CellCyclePhase")==4.0){
        Mcount++;
        TimeArrested += cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(cell_iter->GetCellData()->GetItem("CellCyclePhase")==2.5){
          MeioticS++;
      }

      //Count total cells
      totalCells++;
    }

    int firstMeioticRow = 128;
    int lastMitoticRow = 0;
    for(int i=0; i<128; i++){
      if(meioticCountArray[i]>1 && i<firstMeioticRow){
        firstMeioticRow = i;
      }
      if(mitoticCountArray[i]>0 && i>lastMitoticRow){
        lastMitoticRow = i;
      }
    }

    double deathRate;
    if (SimulationTime::Instance()->GetTime() > 17.0){
      deathRate = GlobalParameterStruct::Instance()->GetParameter(21);
    }
    else{
      deathRate = 0.0;
    }
    
    cellCycleDuration += (TimeArrested) / (Mcount + Scount + G2count + G1count);

    *OutputPositionFile << SimulationTime::Instance()->GetTime() << "\t" << gonadLength << "\t" << cellCycleDuration << "\t" << spermCount << "\t" <<
        prolifCount << "\t" << deathRate << "\t" << totalCells << "\t" << lastProliferativeCell << "\t" << firstMeioticCell << "\t" <<
            G1count << "\t" << Scount << "\t" << G2count << "\t" << Mcount << "\t" << MeioticS << "\t" << firstMeioticRow << "\t" << lastMitoticRow  << "\n";
    OutputPositionFile->flush();
  }
}

template<unsigned DIM>
void GonadArmDataOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}


template<unsigned DIM>
void GonadArmDataOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  // No parameters to output
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class GonadArmDataOutput<1>;
template class GonadArmDataOutput<2>;
template class GonadArmDataOutput<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GonadArmDataOutput)
