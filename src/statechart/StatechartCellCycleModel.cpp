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

#include "StatechartCellCycleModel.hpp"
#include "GlobalParameterStruct.hpp"

StatechartCellCycleModel::StatechartCellCycleModel(bool LoadingFromArchive): AbstractCellCycleModel(){
	
	mLoadingFromArchive=LoadingFromArchive;
	TempVariableStorage=std::vector<double>();
	TempStateStorage=0;
	mDimension=3;
	mCurrentCellCyclePhase=G_ONE_PHASE;
    MAKE_PTR(CellStatechart,newStatechart);
    pStatechart=newStatechart;

	//Set initial cell cycle phase durations
	mG1Duration			   = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
	mTransitCellG1Duration = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
	mStemCellG1Duration    = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
	mSDuration			   = GlobalParameterStruct::Instance()->GetParameter(17)*GlobalParameterStruct::Instance()->GetParameter(14);
	mG2Duration		   	   = GlobalParameterStruct::Instance()->GetParameter(18)*GlobalParameterStruct::Instance()->GetParameter(14);
	mMDuration			   = GlobalParameterStruct::Instance()->GetParameter(19)*GlobalParameterStruct::Instance()->GetParameter(14);
};


AbstractCellCycleModel* StatechartCellCycleModel::CreateCellCycleModel(){
	//Create a new cell cycle model
	StatechartCellCycleModel* newStatechartCellCycleModel = new StatechartCellCycleModel();
	//Ensure values are inhereted from parent as appropriate
	newStatechartCellCycleModel->SetBirthTime(mBirthTime);
	newStatechartCellCycleModel->mG1Duration = mG1Duration;
    newStatechartCellCycleModel->SetMinimumGapDuration(mMinimumGapDuration);
    newStatechartCellCycleModel->SetStemCellG1Duration(mStemCellG1Duration);
    newStatechartCellCycleModel->SetTransitCellG1Duration(mTransitCellG1Duration);
    newStatechartCellCycleModel->SetSDuration(mSDuration);
    newStatechartCellCycleModel->SetG2Duration(mG2Duration);
    newStatechartCellCycleModel->SetMDuration(mMDuration);
	newStatechartCellCycleModel->SetDimension(mDimension);
	newStatechartCellCycleModel->mLoadingFromArchive=mLoadingFromArchive;
	//Create a new statechart.
	MAKE_PTR(CellStatechart, newStatechart);
	//Set its cell pointer to the parent cell to avoid it being null when state constructors are called.
	newStatechart->SetCell(mpCell);
	//Copy the state of the parent chart. Give result to the daughter cell cycle model.
	newStatechartCellCycleModel->pStatechart=pStatechart->Copy(newStatechart);
	//Return the new cell cycle model. The cell pointer will be changed to point to the daughter when SetCell is called.
	return newStatechartCellCycleModel;
 };


void StatechartCellCycleModel::Initialise(){
	pStatechart->initiate();
 	//Handles advancing the cell to some point in the cell cycle so that the population DOES NOT start out synchronised.

	bool synchronisation = false;

 	//1) Get the cell's randomly generated age to work out what phase it's in and how long it's been there.
	double randStartingAge;
	if(synchronisation){
		randStartingAge = 2.5;
	}else{
		randStartingAge = mpCell->GetAge();
	}
	double intPart;
	modf(randStartingAge/(mG1Duration+mG2Duration+mMDuration+mSDuration),&intPart);
	double remainder=randStartingAge-intPart*(mG1Duration+mG2Duration+mMDuration+mSDuration);

	//Fire an event to force the statechart into the appropriate phase
	CellCyclePhase_ phase;
	if(remainder<mG1Duration){
		pStatechart->process_event(EvGoToCellCycle_Mitosis_G1());
		pStatechart->TimeInPhase=remainder;
	}else if(remainder>mG1Duration && remainder<mG1Duration+mSDuration){
		pStatechart->process_event(EvGoToCellCycle_Mitosis_S());
		pStatechart->TimeInPhase=remainder-mG1Duration;
	}
	else if(remainder>mG1Duration+mSDuration && remainder<mG1Duration+mSDuration+mG2Duration){
		pStatechart->process_event(EvGoToCellCycle_Mitosis_G2());		
		pStatechart->TimeInPhase=remainder-(mG1Duration+mSDuration);
	}else{
		pStatechart->process_event(EvGoToCellCycle_Mitosis_M());		
		pStatechart->TimeInPhase=remainder-(mG1Duration+mSDuration+mG2Duration);
	}
};


void StatechartCellCycleModel::UpdateCellCyclePhase(){
	//To update the phase, just update the statechart
	pStatechart->process_event(EvCheckCellData());
};


//Update cell cycle phase durations over time. Updates only occur when a phase duration is actually needed.
double StatechartCellCycleModel::GetG1Duration(){
	double adultG1  = GlobalParameterStruct::Instance()->GetParameter(31)*GlobalParameterStruct::Instance()->GetParameter(15);
	double larvalG1 = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
	double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5; 
	double duration = 4.5;
	double currentG1=larvalG1;
	if (SimulationTime::Instance()->GetTime()>delay && SimulationTime::Instance()->GetTime() <= delay + duration){
		currentG1=(larvalG1 + (SimulationTime::Instance()->GetTime() - delay)*((adultG1 - larvalG1) / duration));
	}
	else if (SimulationTime::Instance()->GetTime() >= delay + duration){
		currentG1 = (adultG1);
	}
	return currentG1;
};
double StatechartCellCycleModel::GetStemCellG1Duration(){
	double adultG1 = GlobalParameterStruct::Instance()->GetParameter(31)*GlobalParameterStruct::Instance()->GetParameter(15);
	double larvalG1 = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
	double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5;
	double duration = 4.5;
	double currentG1 = larvalG1;
	if (SimulationTime::Instance()->GetTime()>delay && SimulationTime::Instance()->GetTime() <= delay + duration){
		currentG1 = (larvalG1 + (SimulationTime::Instance()->GetTime() - delay)*((adultG1 - larvalG1) / duration));
	}
	else if (SimulationTime::Instance()->GetTime() >= delay + duration){
		currentG1 = (adultG1);
	}
	this->SetStemCellG1Duration(currentG1);
	return currentG1;
};
double StatechartCellCycleModel::GetTransitCellG1Duration(){
	double adultG1 = GlobalParameterStruct::Instance()->GetParameter(31)*GlobalParameterStruct::Instance()->GetParameter(15);
	double larvalG1 = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
	double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5;
	double duration = 4.5;
	double currentG1 = larvalG1;
	if (SimulationTime::Instance()->GetTime()>delay && SimulationTime::Instance()->GetTime() <= delay + duration){
		currentG1 = (larvalG1 + (SimulationTime::Instance()->GetTime() - delay)*((adultG1 - larvalG1) / duration));
	}
	else if (SimulationTime::Instance()->GetTime() >= delay + duration){
		currentG1 = (adultG1);
	}
	this->SetTransitCellG1Duration(currentG1);
	return currentG1;
};
double StatechartCellCycleModel::GetSDuration(){
	double adultS = GlobalParameterStruct::Instance()->GetParameter(32)*GlobalParameterStruct::Instance()->GetParameter(15);
	double larvalS = GlobalParameterStruct::Instance()->GetParameter(17)*GlobalParameterStruct::Instance()->GetParameter(14);
	double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5;
	double duration = 4.5;
	double currentS = larvalS;
	if (SimulationTime::Instance()->GetTime()>delay && SimulationTime::Instance()->GetTime() <= delay + duration){
		currentS = (larvalS + (SimulationTime::Instance()->GetTime() - delay)*((adultS - larvalS) / duration));
	}
	else if (SimulationTime::Instance()->GetTime() >= delay + duration){
		currentS = (adultS);
	}
	this->SetSDuration(currentS);
	return currentS;
};
double StatechartCellCycleModel::GetG2Duration(){
	double adultG2 = GlobalParameterStruct::Instance()->GetParameter(33)*GlobalParameterStruct::Instance()->GetParameter(15);
	double larvalG2 = GlobalParameterStruct::Instance()->GetParameter(18)*GlobalParameterStruct::Instance()->GetParameter(14);
	double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5;
	double duration = 4.5;
	double currentG2 = larvalG2;
	if (SimulationTime::Instance()->GetTime()>delay && SimulationTime::Instance()->GetTime() <= delay + duration){
		currentG2 = (larvalG2 + (SimulationTime::Instance()->GetTime() - delay)*((adultG2 - larvalG2) / duration));
	}
	else if (SimulationTime::Instance()->GetTime() >= delay + duration){
		currentG2 = (adultG2);
	}
	this->SetG2Duration(currentG2);
	return currentG2;
};
double StatechartCellCycleModel::GetMDuration(){
	double adultM = GlobalParameterStruct::Instance()->GetParameter(34)*GlobalParameterStruct::Instance()->GetParameter(15);
	double larvalM = GlobalParameterStruct::Instance()->GetParameter(19)*GlobalParameterStruct::Instance()->GetParameter(14);
	double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5;
	double duration = 4.5;
	double currentM = larvalM;
	if (SimulationTime::Instance()->GetTime()>delay && SimulationTime::Instance()->GetTime() <= delay + duration){
		currentM = (larvalM + (SimulationTime::Instance()->GetTime() - delay)*((adultM - larvalM) / duration));
	}
	else if (SimulationTime::Instance()->GetTime() >= delay + duration){
		currentM = (adultM);
	}
	this->SetMDuration(currentM);
	return currentM;
};
double StatechartCellCycleModel::GetSG2MDuration(){
	return (this->GetSDuration() + this->GetG2Duration() + this->GetMDuration());
};


void StatechartCellCycleModel::SetCell(CellPtr pCell){
	//Hold a pointer to the cell in this class
	mpCell = pCell;
	//Set the statechart's cell pointer to point to the same cell.
	pStatechart->SetCell(mpCell);
	//If we're loading from an archive, now is an appropriate time to initiate the statechart and set its state
	//and variables from stored values. 
	if (mLoadingFromArchive == true){
		pStatechart->initiate();
		pStatechart->SetState(TempStateStorage);
		pStatechart->SetVariables(TempVariableStorage);
		mLoadingFromArchive = false;
	}
};


bool StatechartCellCycleModel::ReadyToDivide(){
	if (!mReadyToDivide)	//If not dividing, update the statechart
	{
		UpdateCellCyclePhase();
	}
	return mReadyToDivide;  //Return whether ready to divide now
};


void StatechartCellCycleModel::ResetForDivision(){
	//To reset, change the mReadyToDivide flag to false: the message has been received
	mReadyToDivide = false;
};


//Setter methods that allow the statechart to set members of the cell cycle model
void StatechartCellCycleModel::SetCellCyclePhase(CellCyclePhase_ Phase){
	mCurrentCellCyclePhase = Phase;
};
void StatechartCellCycleModel::SetReadyToDivide(bool Ready){
	mReadyToDivide = Ready;
};


//Standard outputting of parameters associated with the cell cycle model to a params file for storage
//Haven't added any new parameters in particular.
void StatechartCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
	AbstractCellCycleModel::OutputCellCycleModelParameters( rParamsFile);
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StatechartCellCycleModel)
