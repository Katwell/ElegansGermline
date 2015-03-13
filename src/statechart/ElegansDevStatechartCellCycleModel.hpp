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

#ifndef ELEGANSDEVSTATECHARTCELLCYCLEMODEL_HPP_
#define ELEGANSDEVSTATECHARTCELLCYCLEMODEL_HPP_

#include "StatechartCellCycleModel.hpp"
#include "GlobalParameterStruct.hpp" 
#include "ChasteSerialization.hpp"

/*
* A child class of StatechartCellCycleModel for running C. elegans germ line development simulations. 
* The main change is that the getter methods for cell cycle phase durations have been overriden. 
* When, say, the length of G2 is requested, this class calculates an appropriate answer based on the
* simulation parameters and the current time. This allows the cell cycle length to change over the 
* course of the worm's larval development. 
*
* As in the parent class, CELLSTATECHART is a statechart model of cell behaviour, see 
* FateUncoupledFromCycle for an example.
*/

template< typename CELLSTATECHART >
class ElegansDevStatechartCellCycleModel: public StatechartCellCycleModel< CELLSTATECHART >{

private:

    //Standard serialization block:
    friend class boost::serialization::access;   
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {   
        SerializableSingleton<GlobalParameterStruct>* s_wrapper = GlobalParameterStruct::Instance()->GetSerializationWrapper();
        archive & s_wrapper;    
        archive & boost::serialization::base_object<StatechartCellCycleModel<CELLSTATECHART> >(*this); 
    }

public:

	ElegansDevStatechartCellCycleModel(bool LoadingFromArchive = false, bool SynchronisedCells = false):
	   StatechartCellCycleModel< CELLSTATECHART >(LoadingFromArchive, SynchronisedCells){

        this->SetG1Duration( GlobalParameterStruct::Instance()->GetParameter(16)*
                             GlobalParameterStruct::Instance()->GetParameter(14) );
        this->SetTransitCellG1Duration( GlobalParameterStruct::Instance()->GetParameter(16)*
                                        GlobalParameterStruct::Instance()->GetParameter(14) );
        this->SetStemCellG1Duration( GlobalParameterStruct::Instance()->GetParameter(16)*
                                     GlobalParameterStruct::Instance()->GetParameter(14) );
        this->SetSDuration( GlobalParameterStruct::Instance()->GetParameter(17)*
                            GlobalParameterStruct::Instance()->GetParameter(14) );
        this->SetG2Duration( GlobalParameterStruct::Instance()->GetParameter(18)*
                             GlobalParameterStruct::Instance()->GetParameter(14) );
        this->SetMDuration( GlobalParameterStruct::Instance()->GetParameter(19)*
                            GlobalParameterStruct::Instance()->GetParameter(14) );
	}

    ~ElegansDevStatechartCellCycleModel(){};



    /**
    * Builder method that creates new instances of this cell-cycle model for daughter cells.
    */
    virtual AbstractCellCycleModel* CreateCellCycleModel(){        

        //THIS is the key point. The rest is identical to StatechartCellCycleModel
        ElegansDevStatechartCellCycleModel* newStatechartCellCycleModel = new ElegansDevStatechartCellCycleModel();

        newStatechartCellCycleModel->SetG1Duration(GetG1Duration());
        newStatechartCellCycleModel->SetStemCellG1Duration(GetG1Duration());
        newStatechartCellCycleModel->SetTransitCellG1Duration(GetG1Duration());
        newStatechartCellCycleModel->SetSDuration(GetSDuration());
        newStatechartCellCycleModel->SetG2Duration(GetG2Duration());
        newStatechartCellCycleModel->SetMDuration(GetMDuration());
        newStatechartCellCycleModel->SetDimension(AbstractCellCycleModel::mDimension);
        newStatechartCellCycleModel->SetBirthTime(AbstractCellCycleModel::mBirthTime);
        newStatechartCellCycleModel->SetMinimumGapDuration(AbstractCellCycleModel::mMinimumGapDuration);
        newStatechartCellCycleModel->mLoadingFromArchive = false;
        
        MAKE_PTR(CELLSTATECHART, newStatechart);
        newStatechart->SetCell(AbstractCellCycleModel::mpCell);
        newStatechartCellCycleModel->pStatechart = StatechartCellCycleModel<CELLSTATECHART>::pStatechart->CopyInto(newStatechart);

        return newStatechartCellCycleModel;
    }; 




    /*
    * Method that ensures cells are not synchronised at the start of a simulation, and initialises
    * the statechart
    */
    virtual void Initialise(){

        //Initiate can only be called AFTER the cell pointer has been set.
        StatechartCellCycleModel<CELLSTATECHART>::pStatechart->initiate();
    
        double randStartingAge;
        if(StatechartCellCycleModel<CELLSTATECHART>::mSynchronisedCells){
            randStartingAge = 2.5;
        }else{
            randStartingAge = AbstractCellCycleModel::mpCell->GetAge();
        }
        double nCyclesCompleted;
        modf(randStartingAge/(GetG1Duration() + GetG2Duration() + GetMDuration() + GetSDuration()), &nCyclesCompleted);
        double remainder = randStartingAge - nCyclesCompleted * (GetG1Duration() + GetG2Duration() + GetMDuration() + GetSDuration());
    
        //Fire an event to force the statechart into the appropriate starting phase
        CellCyclePhase_ phase;
        if(remainder < GetG1Duration()){
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->process_event(EvGoToCellCycle_Mitosis_G1());
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->TimeInPhase = remainder;
        }else if(remainder > GetG1Duration() && remainder < GetG1Duration() + GetSDuration()){
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->process_event(EvGoToCellCycle_Mitosis_S());
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->TimeInPhase = remainder - GetG1Duration();
        }
        else if(remainder > GetG1Duration()+ GetSDuration() && remainder < GetG1Duration() +
            GetSDuration() + GetG2Duration()){
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->process_event(EvGoToCellCycle_Mitosis_G2());       
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->TimeInPhase = remainder - (GetG1Duration() + GetSDuration());
        }else{
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->process_event(EvGoToCellCycle_Mitosis_M());        
            StatechartCellCycleModel<CELLSTATECHART>::pStatechart->TimeInPhase = remainder - (GetG1Duration() + GetSDuration() + GetG2Duration());
        }
    };





	/*
    * Override all phase duration getter methods, so that the appropriate phase duration 
    * for a worm of the current age is calculated when the chart requests it.
    */
    virtual double GetG1Duration(){
        double adultG1  = GlobalParameterStruct::Instance()->GetParameter(31)*GlobalParameterStruct::Instance()->GetParameter(15);
        double larvalG1 = GlobalParameterStruct::Instance()->GetParameter(16)*GlobalParameterStruct::Instance()->GetParameter(14);
        double delay = GlobalParameterStruct::Instance()->GetParameter(20) - 18.5; 
        double duration = 4.5;
        double currentG1 = larvalG1;
        if (SimulationTime::Instance()->GetTime() > delay && SimulationTime::Instance()->GetTime() <= delay + duration){
            currentG1=(larvalG1 + (SimulationTime::Instance()->GetTime() - delay)*((adultG1 - larvalG1) / duration));
        }
        else if (SimulationTime::Instance()->GetTime() >= delay + duration){
            currentG1 = (adultG1);
        }
        return currentG1;
    };

    virtual double GetStemCellG1Duration(){
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

    virtual double GetTransitCellG1Duration(){
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

    virtual double GetSDuration(){
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

    virtual double GetG2Duration(){
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

    virtual double GetMDuration(){
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

    virtual double GetSG2MDuration(){
        return (this->GetSDuration() + this->GetG2Duration() + this->GetMDuration());
    };


};

#endif /*ELEGANSDEVSTATECHARTCELLCYCLEMODEL_HPP_*/