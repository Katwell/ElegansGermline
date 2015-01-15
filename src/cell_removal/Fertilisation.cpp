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

#include "Fertilisation.hpp"

template<unsigned DIM>
Fertilisation<DIM>::Fertilisation(AbstractCellPopulation<DIM>* pCellPopulation, double spermathecaLength)
    : AbstractCellKiller<DIM>(pCellPopulation),
    mSpermathecaLength(spermathecaLength){
}


template<unsigned DIM>
Fertilisation<DIM>::~Fertilisation(){
}


template<unsigned DIM>
double Fertilisation<DIM>::GetSpermathecaLength() const
{
    return mSpermathecaLength;
}


template<unsigned DIM>
void Fertilisation<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    if (!pCell->HasApoptosisBegun()){
        pCell->StartApoptosis();
    }
}


template<unsigned DIM>
void Fertilisation<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{

    bool stop = false;

    //If the worm has reached adulthood...
    if (SimulationTime::Instance()->GetTime() > 17.0){

        //Determine the final length of the gonad, so we can work out how far a cell must be from the DTC
        //to be in the spermatheca
        double gonadLength = 0;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End();
            ++cell_iter)
        {
            if (cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC") > gonadLength){
                gonadLength = cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC");
            }
        }


        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End();
            ++cell_iter)
        {
            
            //Loop over all cells and only consider mature oocytes in the spermatheca 
            if (cell_iter->GetCellData()->GetItem("Differentiation_Oocyte") == 1.0 
                && gonadLength - cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC") <= mSpermathecaLength 
                && cell_iter->GetCellData()->GetItem("Radius") > 11.0 
                && cell_iter->HasApoptosisBegun() == false
                && stop==false){

                //Now loop through in search of fertilising sperm
                for (typename AbstractCellPopulation<DIM>::Iterator cell_iter2 = this->mpCellPopulation->Begin();
                    cell_iter2 != this->mpCellPopulation->End(); ++cell_iter2)
                {
                    if (cell_iter2->GetCellData()->GetItem("Differentiation_Sperm") == 1.0
                        && cell_iter2->HasApoptosisBegun() == false
                        && gonadLength - cell_iter2->GetCellData()->GetItem("DistanceAwayFromDTC") <= mSpermathecaLength
                        && stop==false){

                        CheckAndLabelSingleCellForApoptosis(*cell_iter2);
                        CheckAndLabelSingleCellForApoptosis(*cell_iter);
                        stop = true;
                    }
                }

            }
        }

    }
}


template<unsigned DIM>
void Fertilisation<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<LengthOfSpermatheca>" << mSpermathecaLength << "</LengthOfSpermatheca>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class Fertilisation<1>;
template class Fertilisation<2>;
template class Fertilisation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Fertilisation)
