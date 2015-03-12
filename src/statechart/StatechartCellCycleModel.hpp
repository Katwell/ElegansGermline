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

#ifndef STATECHARTCELLCYCLEMODEL_HPP_
#define STATECHARTCELLCYCLEMODEL_HPP_

#define MAX_STATE_COUNT 32

#include "AbstractCellCycleModel.hpp"
#include "AbstractStatechartCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include <boost/statechart/event.hpp>
namespace sc = boost::statechart;

#include "SmartPointers.hpp"
#include "Identifiable.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/bitset.hpp>

/*This class is a wrapper around a statechart model of cell behaviour <CELLSTATECHART>. 
* It ensures that to the rest of Chaste, the statechart appears as a normal cell cycle model.
*
* In fact, the whole class is basically a typical Chaste cell cycle model, plus an extra member
* that is a pointer to a statechart. The statechart updates whenever UpdateCellCyclePhase
* is called, and takes over responsibility for setting the ReadyToDivide and CellCyclePhase flags.
*
* The StatechartCellCycleModel itself remains responsible for setting the G1,S,G2 and M durations.
* These values are queried by the statechart using getter methods.
* Because Boost Statecharts doesn't support archiving, this wrapper also deals with saving
* the state of a chart and any variables associated with it as required. 
* The state is encoded as a bitset for archiving purposes, while any statechart associated 
* variables are stored as a std::vector. 
*/


/*
* 5 predeclared events that any Statechart model must be able to respond to, in order to be compatible
* with the wider Chaste code. EvCheckCellData prompts the chart to update with reference to the current 
* properties of its cell. The other 4 events allow the statechart to be forced into a particular cell 
* cycle phase at the start of a simulation.
*/ 
struct EvCheckCellData :  sc::event< EvCheckCellData > {};
struct EvGoToCellCycle_Mitosis_M : sc::event< EvGoToCellCycle_Mitosis_M > {};
struct EvGoToCellCycle_Mitosis_S : sc::event< EvGoToCellCycle_Mitosis_S > {};
struct EvGoToCellCycle_Mitosis_G2 : sc::event< EvGoToCellCycle_Mitosis_G2 > {};
struct EvGoToCellCycle_Mitosis_G1 : sc::event< EvGoToCellCycle_Mitosis_G1 > {};



template<typename CELLSTATECHART>
class StatechartCellCycleModel : public AbstractCellCycleModel, public AbstractStatechartCellCycleModel
{

private:

    //Standard serialization block, this doesn't deal with archiving the statechart:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;  
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    }


public:

    /*
    * Constructor:
    * 1) Makes a new AbstractCellCycleModel
    * 2) Leaves the statechart's pointer to its cell NULL
    */
    StatechartCellCycleModel(bool LoadingFromArchive = false, bool SynchronisedCells = false):
    AbstractCellCycleModel(){

        mLoadingFromArchive = LoadingFromArchive;
        mSynchronisedCells = SynchronisedCells;
        TempVariableStorage = std::vector<double>();
        TempStateStorage = 0;
        
        mDimension = 3;
        mCurrentCellCyclePhase = G_ONE_PHASE;
        
        MAKE_PTR(CELLSTATECHART, myStatechart);
        pStatechart = myStatechart;
    };  

    ~StatechartCellCycleModel(){};


    /* Pointer to a statechart model, of type CELLSTATECHART */
    boost::shared_ptr<CELLSTATECHART> pStatechart;   
    
    /* Variables associated with archiving a statechart */
    bool mLoadingFromArchive;
    bool mSynchronisedCells;               
    std::vector<double> TempVariableStorage;
    std::bitset< MAX_STATE_COUNT > TempStateStorage;




    /*
    * Because a cell cycle model doesn't get a pointer to its cell until AFTER construction, this method is
    * where we set the cell pointer and then pass it on to the statechart. The method SetCell must be made 
    * virtual in AbstractCellCycleModel so we can override it here. 
    */
    void SetCell(CellPtr pCell){ 
        mpCell = pCell;
        pStatechart->SetCell(mpCell);
        //If loading from an archive, now is an appropriate time to initiate the chart and set its state
        //and variables from stored values. 
        if (mLoadingFromArchive == true){
            pStatechart->initiate();
            pStatechart->SetState(TempStateStorage);
            pStatechart->SetVariables(TempVariableStorage);
            mLoadingFromArchive = false;
        }
    };
    


    /**
    * Builder method that creates new instances of the cell-cycle model for daughter cells.
    */
    virtual AbstractCellCycleModel* CreateCellCycleModel(){     

        StatechartCellCycleModel* newStatechartCellCycleModel = new StatechartCellCycleModel();
        //Ensure certain values are inhereted from the parent model:
        newStatechartCellCycleModel->SetBirthTime(mBirthTime);
        newStatechartCellCycleModel->SetG1Duration(GetG1Duration());
        newStatechartCellCycleModel->SetMinimumGapDuration(mMinimumGapDuration);
        newStatechartCellCycleModel->SetStemCellG1Duration(GetStemCellG1Duration());
        newStatechartCellCycleModel->SetTransitCellG1Duration(GetTransitCellG1Duration());
        newStatechartCellCycleModel->SetSDuration(GetSDuration());
        newStatechartCellCycleModel->SetG2Duration(GetG2Duration());
        newStatechartCellCycleModel->SetMDuration(GetMDuration());
        newStatechartCellCycleModel->SetDimension(mDimension);
        newStatechartCellCycleModel->mLoadingFromArchive = false;
        //Create a new statechart.
        MAKE_PTR(CELLSTATECHART, newStatechart);
        //Set its cell pointer to point at the parent cell. We'll change this later, but the cell pointer
        //can't be null when the state constructors are called for the first time.
        newStatechart->SetCell(mpCell);
        //Copy the state of the parent's chart into the daughter chart.
        newStatechartCellCycleModel->pStatechart = pStatechart->CopyInto(newStatechart);
        //Return the new cell cycle model. The cell pointer will be changed to point at the daughter cell 
        //when SetCell is called on the new model.
        return newStatechartCellCycleModel;
    };    
    


    /*
    * Method that ensures cells are not synchronised at the start of a simulation. Also initialises the
    * statechart, which can only be done AFTER the cell pointer has been set
    */
    virtual void Initialise(){

        //Initialise chart
        pStatechart->initiate();
    
        //Set cell starting state
        double startingAge;
        if(mSynchronisedCells){
            startingAge = 0;
        }else{
            startingAge = mpCell->GetAge();
        }
        double nCyclesCompleted;
        modf(startingAge/(GetG1Duration() + GetG2Duration() + GetMDuration() + GetSDuration()), &nCyclesCompleted);
        double remainder = startingAge - nCyclesCompleted * (GetG1Duration() + GetG2Duration() + GetMDuration() + GetSDuration());
    
        //Fire an event to force the statechart into the appropriate starting phase
        if(remainder < GetG1Duration()){
            pStatechart->process_event(EvGoToCellCycle_Mitosis_G1());
            pStatechart->TimeInPhase = remainder;
        }else if(remainder > GetG1Duration() && remainder < GetG1Duration() + GetSDuration()){
            pStatechart->process_event(EvGoToCellCycle_Mitosis_S());
            pStatechart->TimeInPhase = remainder - GetG1Duration();
        }
        else if(remainder > GetG1Duration()+GetSDuration() && remainder < GetG1Duration() + GetSDuration() + GetG2Duration()){
            pStatechart->process_event(EvGoToCellCycle_Mitosis_G2());       
            pStatechart->TimeInPhase = remainder - (GetG1Duration() + GetSDuration());
        }else{
            pStatechart->process_event(EvGoToCellCycle_Mitosis_M());        
            pStatechart->TimeInPhase = remainder - (GetG1Duration() + GetSDuration() + GetG2Duration());
        }
    };



    /**
    * @return whether the cell is ready to divide. Set by the statechart.
    */
    bool ReadyToDivide(){
        if (!mReadyToDivide){
            UpdateCellCyclePhase();
        }
        return mReadyToDivide;
    };    
    
    /**
    * This method updates the statechart, and the statechart in turn updates the current phase and
    * sets ReadyToDivide as appropriate
    */
    void UpdateCellCyclePhase(){
        pStatechart->process_event(EvCheckCellData());
    };  

    /*
    * Only a slight change to this method from the usual one - we don't need to reset the cell cycle phase
    * after division because the statechart will handle that.
    */
    void ResetForDivision(){
        mReadyToDivide = false;
    };    



    //Setter methods for use by the statechart
    void SetCellCyclePhase(CellCyclePhase_ Phase){
        mCurrentCellCyclePhase = Phase;
    };
    void SetReadyToDivide(bool Ready){
        mReadyToDivide = Ready;
    };

    //A G1 setter method for use by child classes. Doesn't exist in AbstractCellCycleModel for some reason.
    void SetG1Duration(double duration){
        mG1Duration = duration;
    }



    //Standard parameter logging. TODO: Is there anything I could usefully be adding here?
    void OutputCellCycleModelParameters(out_stream& rParamsFile){
        AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
    };

};



//ARCHIVING METHODS FOR THE STATECHART:
namespace boost
{
namespace serialization
{

    template<class Archive, typename CELLSTATECHART>
    inline void save_construct_data(
    Archive & ar, const StatechartCellCycleModel<CELLSTATECHART>* t, const BOOST_PFTO unsigned int file_version)
    {
        std::bitset<MAX_STATE_COUNT> state = t->pStatechart->GetState();
        ar << state;
        std::vector<double> v = t->pStatechart->GetVariables();
        int numberOfVars=v.size();
        ar << numberOfVars;
        for(int i=0; i<(int)v.size(); i++){
            ar<<v.at(i);
        }
    }
    
    template<class Archive, typename CELLSTATECHART>
    inline void load_construct_data(
        Archive & ar, StatechartCellCycleModel<CELLSTATECHART>* t, const unsigned int file_version)
    {
        std::bitset<MAX_STATE_COUNT> state;
        ar >> state;
        int numberOfVars;
        ar >> numberOfVars;
        std::vector<double> v;
        for(int i=0; i<numberOfVars; i++){
            double value;
            ar >> value;
            v.push_back(value);
        }

        // Construct a new cell cycle model and set the statechart's state and variable values.
        ::new(t)StatechartCellCycleModel<CELLSTATECHART>(true, false);
        t->TempStateStorage = state;
        t->TempVariableStorage = v;
    }
}
} // namespace ...

#endif  /*STATECHARTCELLCYCLEMODEL_HPP_*/

/* 
*  NOTE: We do not export the StatechartCellCycleModel class and make the serializer aware of it here, 
*  because the type is not complete. Instead the EXPORT_TEMPLATE_CLASS1 statement comes at the end of 
*  each statechart model file. There we export StatechartCellCycleModel<MyModel>.
*/ 
