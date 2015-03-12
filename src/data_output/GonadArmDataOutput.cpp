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

#include "GonadArmDataOutput.hpp"
#include "NodeBasedCellPopulation.hpp"


//Constructor, initialises sampling interval and sets output file to null
template<unsigned DIM>
GonadArmDataOutput<DIM>::GonadArmDataOutput(int interval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputFile(NULL),
      mInterval(interval)
{}


//Empty destructor 
template<unsigned DIM>
GonadArmDataOutput<DIM>::~GonadArmDataOutput(){}


//Getter methods for private members
template<unsigned DIM>
int GonadArmDataOutput<DIM>::GetInterval() const
{
  return mInterval;
};


//Open an output file GonadData.txt in the simulation directory
template<unsigned DIM>
void GonadArmDataOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  OutputFile = rOutputFileHandler.OpenOutputFile("GonadData.txt");
}


//At each timestep, if the time is a sampling time, loop through all cells and compile some general gonad data.
//Output that data to file.
template<unsigned DIM>
void GonadArmDataOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  //If it's a sampling time, start gathering some useful data
  if(SimulationTime::Instance()->GetTimeStepsElapsed() % GetInterval() ==0){


    //Estimates the length in microns of one cell row, based on the cell separations in the first
    //75 microns of the distal zone:
    int    numberOfMeasurements = 0;
    double meanSeparation = 0;
    //loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End();++cell_iter)
    { //consider distal arm germ cells only
      if(cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC")<75.0 &&
         cell_iter->GetCellData()->GetItem("IsDTC")==0.0)
      {
        //Get cell's compressed volume
        double vol = cell_iter->GetCellData()->GetItem("volume");
        //Work out effective diameter and add to meanSeparation running total
        double diam = 2*pow(vol*0.239,0.3333);
        meanSeparation += diam;
        numberOfMeasurements++;
      }
    }
    meanSeparation = meanSeparation/numberOfMeasurements;


    //Initialize variables to store the data we'll generate while looping over all cells
    double cellCycleDuration = (rCellPopulation.Begin())->GetCellCycleModel()->GetSDuration()
                              +(rCellPopulation.Begin())->GetCellCycleModel()->GetG2Duration()
                              +(rCellPopulation.Begin())->GetCellCycleModel()->GetTransitCellG1Duration()
                              +(rCellPopulation.Begin())->GetCellCycleModel()->GetMDuration();
    double gonadLength = 0;
    double lastProliferativeCell = 0;
    double firstMeioticCell = DBL_MAX;
    int totalCells = 0;
    int spermCount = 0;
    int prolifCount = 0;
    int G1count = 0;
    int Scount = 0;
    int G2count = 0;
    int Mcount = 0;
    int MeioticS = 0;
    double TimeArrested = 0;
    double meioticCountArray[128] = {0}; // <- These arrays store info on whether there is
    double mitoticCountArray[128] = {0}; // a prolif cell / 2 meiotic cells in each row for 
                                         // the first 128 rows.

    //Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End();++cell_iter)
    {

      //Work out a cell's row (counting from the DTC), based on the mean compressed cell diameter
      double dist = cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC");
      double row = round(dist / meanSeparation);
      cell_iter->GetCellData()->SetItem("RowNumber",row);
      double phase = cell_iter->GetCellData()->GetItem("CellCyclePhase");
      //For the first 128 rows:
      if(row < 128.0){
       if(phase < 0){
        //If cell is in meiosis, record the fact in meioticCountArray
        meioticCountArray[(int)row] = meioticCountArray[(int)row] + 1;
       }else{
        //otherwise record a proliferative cell in mitoticCountArray
        mitoticCountArray[(int)row] = mitoticCountArray[(int)row] + 1;
       }
      }

      //Establish the gonad length
      if(dist > gonadLength){
        gonadLength = dist;
      }

      //Count sperm 
      if(cell_iter->GetCellData()->GetItem("Differentiation_Sperm") == 1.0){
        spermCount++;
      }

      //count proliferative cells and calculate the positions of the closes meiotic cell to the
      //DTC and furthest mitotic cell from DTC
      if(phase > 0.0){
        prolifCount++;
        if (dist >lastProliferativeCell){
          lastProliferativeCell = dist;
        }
      }else{
        if (dist < firstMeioticCell){
          firstMeioticCell = dist;
        }
      }

      //Count number of cells in each phase, and track how long cells are spending in arrest
      if(phase == 1.0){
        G1count++;
         TimeArrested += cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(phase == 2.0){
        Scount++;
        TimeArrested += cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(phase == 3.0){
        G2count++;
          TimeArrested += cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(phase == 4.0){
        Mcount++;
        TimeArrested += cell_iter->GetCellData()->GetItem("ArrestedFor");
      }else if(phase == 2.5){
        MeioticS++;
      }

      //Count total cells
      totalCells++;
    }

    //Loop through our stored info on the cells each row contains. Work out which is
    //the first meiotic row (with 2 meiotic cells) and which is the last mitotic row
    //(1 proliferative cell)
    int firstMeioticRow = 128;
    int lastMitoticRow = 0;
    for(int i=0; i<128; i++){
      if(meioticCountArray[i] > 1 && i < firstMeioticRow){
        firstMeioticRow = i;
      }
      if(mitoticCountArray[i] > 0 && i > lastMitoticRow){
        lastMitoticRow = i;
      }
    }

    //Get the cell death rate parameter, just for reference really. TODO: replace with
    //a count of the actual number of apoptotic cells? Would be more informative.
    double deathRate = GlobalParameterStruct::Instance()->GetParameter(21);
    
    //Add mean time spent in arrest to the programmed in cell cycle duration
    cellCycleDuration += (TimeArrested) / (Mcount + Scount + G2count + G1count);

    //Write data
    *OutputFile << SimulationTime::Instance()->GetTime() << "\t" 
                << gonadLength << "\t" 
                << cellCycleDuration << "\t" 
                << spermCount << "\t" 
                << prolifCount << "\t" 
                << deathRate << "\t" 
                << totalCells << "\t" 
                << lastProliferativeCell << "\t" 
                << firstMeioticCell << "\t" 
                << G1count << "\t" 
                << Scount << "\t" 
                << G2count << "\t" 
                << Mcount << "\t" 
                << MeioticS << "\t" 
                << firstMeioticRow << "\t" 
                << lastMitoticRow  << "\n";
    
    //Flush the output file to record data as soon as possible
    OutputFile->flush();

    //Output time for reference
    std::cout << "Time = " << SimulationTime::Instance()->GetTime() << std::endl;
  }

  //If the simulation is finished, close the output file.
  if(SimulationTime::Instance()->IsFinished()){
    OutputFile->close();
  }
}


//Output this class's parameters to a log file
template<unsigned DIM>
void GonadArmDataOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<SampleDataEveryXTimesteps>" << mInterval << "</SampleDataEveryXTimesteps>\n";
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile); 
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
