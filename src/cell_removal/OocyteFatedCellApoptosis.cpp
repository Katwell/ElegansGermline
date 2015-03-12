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
   this list of conditiçons and the following disclaimer in the documentation
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


#include "OocyteFatedCellApoptosis.hpp"


//Constructor, initialises mHourlyProbabilityOfDeath
template<unsigned DIM>
OocyteFatedCellApoptosis<DIM>::OocyteFatedCellApoptosis(AbstractCellPopulation<DIM>* pCellPopulation, double HourlyProbabilityOfDeath)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mHourlyProbabilityOfDeath(HourlyProbabilityOfDeath){
}


//Empty destructor
template<unsigned DIM>
OocyteFatedCellApoptosis<DIM>::~OocyteFatedCellApoptosis(){
}


//Getter for mHourlyProbabilityOfDeath
template<unsigned DIM>
double OocyteFatedCellApoptosis<DIM>::GetHourlyProbabilityOfDeath() const
{
    return mHourlyProbabilityOfDeath;
}


/*
* Actual cell killer method. Also sets the property Apoptosis to 1.0, incase another class wants
* to count the number of dying cells at any given moment.
*/
template<unsigned DIM>
void OocyteFatedCellApoptosis<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
  if (!pCell->HasApoptosisBegun()){
    pCell->StartApoptosis();
    pCell->GetCellData()->SetItem("Apoptosis", 1.0);
  }
}


//Method that determines which cells should diie
template<unsigned DIM>
void OocyteFatedCellApoptosis<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{

  //Loop over the cell population
  for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
  cell_iter != this->mpCellPopulation->End(); ++cell_iter)
  {

    //Identify oocyte fated cells that are still less than 250 microns from the DTC
    if (cell_iter->GetCellData()->GetItem("OocyteFated") == 1.0 &&
      cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC") < 250.0){

      // Random number generator used to determine if this cell should undergo apoptosis. The probability
      // of death used here is based on the Chaste class "RandomCellKiller"
      if (RandomNumberGenerator::Instance()->ranf() < 
         (1.0 - pow((1.0 - mHourlyProbabilityOfDeath), SimulationTime::Instance()->GetTimeStep()) ) ){
         CheckAndLabelSingleCellForApoptosis(*cell_iter);
      }
    }
  }
}


//Output mHourlyProbabilityOfDeath parameter to a log file
template<unsigned DIM>
void OocyteFatedCellApoptosis<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<HourlyProbabilityOfDeath>" << mHourlyProbabilityOfDeath << "</HourlyProbabilityOfDeath>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OocyteFatedCellApoptosis<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OocyteFatedCellApoptosis)
