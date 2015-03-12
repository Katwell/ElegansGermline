#include <StatechartInterface.hpp>
#include <FateUncoupledFromCycle.hpp>

/*
* This file actually implements the statechart's functions. 
*
* It begins by implementing a series of methods required by 
* StatechartCellCycleModel, which must be present in every model. 
*
* It then moves on to functions that are called by individual states 
* in response to certain events. These reaction functions really define the
* model's behaviour.
*/



//FIRST, SOME HARDCODED AND RARELY ALTERED DETAILS. SPECIFIC TO THIS MODEL.
double stochasticity         = 0.1;
bool   contactInhibitionInG1 = false;
bool   contactInhibitionInG2 = true;



// 1) IMPLEMENT THE STATECHART'S FUNCTIONS:

//Constructor. Sets the cell pointer to null and any associated variables to 0
FateUncoupledFromCycle::FateUncoupledFromCycle(){
        pCell=boost::shared_ptr<Cell>(); /*!REQUIRED!*/
        TimeInPhase=0;                   /*!REQUIRED!*/
        SpermatocyteDivisions=0;
        SpermDevelopmentDelay=0;
};

//Setter method for the pointer to this chart's cell 
void FateUncoupledFromCycle::SetCell(CellPtr newCell){ /*!REQUIRED!*/
     assert(newCell!=NULL);
     pCell=newCell;
};

//Gets a vector containing all the chart's associated variables
std::vector<double> FateUncoupledFromCycle::GetVariables(){ /*!REQUIRED!*/
    std::vector<double> variables;
    variables.push_back(TimeInPhase);
    variables.push_back(SpermatocyteDivisions);
    variables.push_back(SpermDevelopmentDelay);
    return variables;
}

//Sets the values of all chart associated variables from an input vector
void FateUncoupledFromCycle::SetVariables(std::vector<double> variables){ /*!REQUIRED!*/
    TimeInPhase = variables.at(0);
    SpermatocyteDivisions = variables.at(1);
    SpermDevelopmentDelay = variables.at(2);
}

//Get an encoding of the current state in bitset form
std::bitset<MAX_STATE_COUNT> FateUncoupledFromCycle::GetState(){ /*!REQUIRED!*/
    std::bitset<MAX_STATE_COUNT> state;
    if(state_cast<const GLP1_Unbound*>()!=0){ 
        state.set(1,1);
    }
    if(state_cast<const GLP1_Bound*>()!=0){ 
        state.set(2,1);
    }
    if(state_cast<const GLP1_Absent*>()!=0){ 
        state.set(3,1);
    }
    if(state_cast<const LAG1_Inactive*>()!=0){ 
        state.set(4,1);
    }
    if(state_cast<const LAG1_Active*>()!=0){ 
        state.set(5,1);
    }
    if(state_cast<const GLD1_Inactive*>()!=0){ 
        state.set(6,1);
    }
    if(state_cast<const GLD1_Active*>()!=0){ 
        state.set(7,1);
    }
    if(state_cast<const GLD2_Inactive*>()!=0){ 
        state.set(8,1);
    }
    if(state_cast<const GLD2_Active*>()!=0){ 
        state.set(9,1);
    }
    if(state_cast<const CellCycle_Mitosis_G1*>()!=0){ 
        state.set(10,1);
    }
    if(state_cast<const CellCycle_Mitosis_S*>()!=0){ 
        state.set(11,1);
    }
    if(state_cast<const CellCycle_Mitosis_G2*>()!=0){ 
        state.set(12,1);
    }
    if(state_cast<const CellCycle_Mitosis_M*>()!=0){ 
        state.set(13,1);
    }
    if(state_cast<const CellCycle_ExitedProlif_G1*>()!=0){ 
        state.set(14,1);
    }
    if(state_cast<const CellCycle_ExitedProlif_MeioticS*>()!=0){ 
        state.set(15,1);
    }
    if(state_cast<const CellCycle_ExitedProlif_Meiosis*>()!=0){ 
        state.set(16,1);
    }
    if(state_cast<const Differentiation_Precursor*>()!=0){ 
        state.set(17,1);
    }
    if(state_cast<const Differentiation_SpermFated*>()!=0){ 
        state.set(18,1);
    }
    if(state_cast<const Differentiation_OocyteFated*>()!=0){ 
        state.set(19,1);
    }
    if(state_cast<const Differentiation_Sperm*>()!=0){ 
        state.set(20,1);
    }
    if(state_cast<const Differentiation_Oocyte*>()!=0){ 
        state.set(21,1);
    }
    return state;
}

//Set the current state from a bitset
void FateUncoupledFromCycle::SetState(std::bitset<MAX_STATE_COUNT> state){ /*!REQUIRED!*/
    if(state[1]==1){ 
        process_event(EvGoToGLP1_Unbound());
    }
    if(state[2]==1){ 
        process_event(EvGoToGLP1_Bound());
    }
    if(state[3]==1){ 
        process_event(EvGoToGLP1_Absent());
    }
    if(state[4]==1){ 
        process_event(EvGoToLAG1_Inactive());
    }
    if(state[5]==1){ 
       process_event(EvGoToLAG1_Active());
    }
    if(state[6]==1){ 
        process_event(EvGoToGLD1_Inactive());
    }
    if(state[7]==1){ 
        process_event(EvGoToGLD1_Active());
    }
    if(state[8]==1){ 
        process_event(EvGoToGLD2_Inactive());
    }
    if(state[9]==1){ 
        process_event(EvGoToGLD2_Active());
    }
    if(state[10]==1){ 
        process_event(EvGoToCellCycle_Mitosis_G1());
    }
    if(state[11]==1){ 
        process_event(EvGoToCellCycle_Mitosis_S());
    }
    if(state[12]==1){ 
        process_event(EvGoToCellCycle_Mitosis_G2());
    }
    if(state[13]==1){ 
        process_event(EvGoToCellCycle_Mitosis_M());
    }
    if(state[14]==1){ 
        process_event(EvGoToCellCycle_ExitedProlif_G1());
    }
    if(state[15]==1){ 
        process_event(EvGoToCellCycle_ExitedProlif_MeioticS());
    }
    if(state[16]==1){ 
        process_event(EvGoToCellCycle_ExitedProlif_Meiosis());
    }
    if(state[17]==1){ 
        process_event(EvGoToDifferentiation_Precursor());
    }
    if(state[18]==1){ 
        process_event(EvGoToDifferentiation_SpermFated());
    }
    if(state[19]==1){ 
        process_event(EvGoToDifferentiation_OocyteFated());
    }
    if(state[20]==1){ 
        process_event(EvGoToDifferentiation_Sperm());
    }
    if(state[21]==1){ 
        process_event(EvGoToDifferentiation_Oocyte());
    }
}


//Copy this statechart's state into a new chart (used during division) /*!REQUIRED!*/
boost::shared_ptr<FateUncoupledFromCycle> FateUncoupledFromCycle::CopyInto(boost::shared_ptr<FateUncoupledFromCycle> myNewStatechart){
    myNewStatechart->initiate();
    if(state_cast<const GLP1_Unbound*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Unbound());
    }
    if(state_cast<const GLP1_Bound*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Bound());
    }
    if(state_cast<const GLP1_Absent*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Absent());
    }
    if(state_cast<const LAG1_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToLAG1_Inactive());
    }
    if(state_cast<const LAG1_Active*>()!=0){
        myNewStatechart->process_event(EvGoToLAG1_Active());
    }
    if(state_cast<const GLD1_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToGLD1_Inactive());
    }
    if(state_cast<const GLD1_Active*>()!=0){
        myNewStatechart->process_event(EvGoToGLD1_Active());
    }
    if(state_cast<const GLD2_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToGLD2_Inactive());
    }
    if(state_cast<const GLD2_Active*>()!=0){
        myNewStatechart->process_event(EvGoToGLD2_Active());
    }
    if(state_cast<const CellCycle_Mitosis_G1*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_G1());
    }
    if(state_cast<const CellCycle_Mitosis_S*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_S());
    }
    if(state_cast<const CellCycle_Mitosis_G2*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_G2());
    }
    if(state_cast<const CellCycle_Mitosis_M*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_M());
    }
    if(state_cast<const CellCycle_ExitedProlif_G1*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_ExitedProlif_G1());
    }
    if(state_cast<const CellCycle_ExitedProlif_MeioticS*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_ExitedProlif_MeioticS());
    }
    if(state_cast<const CellCycle_ExitedProlif_Meiosis*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_ExitedProlif_Meiosis());
    }
    if(state_cast<const Differentiation_Precursor*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_Precursor());
    }
    if(state_cast<const Differentiation_SpermFated*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_SpermFated());
    }
    if(state_cast<const Differentiation_OocyteFated*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_OocyteFated());
    }    
    if(state_cast<const Differentiation_Sperm*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_Sperm());
    }
    if(state_cast<const Differentiation_Oocyte*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_Oocyte());
    }    
    myNewStatechart->SpermatocyteDivisions = this->SpermatocyteDivisions;
    myNewStatechart->SpermDevelopmentDelay = this->SpermDevelopmentDelay;
    return (myNewStatechart);
};





// 2) DEFINE THE CASCADE OF EVENTS THAT OCCURS WHEN THE CHART IS PROMPTED TO UPDATE

sc::result Running::react( const EvCheckCellData & ){
    /*!REQUIRED! - Chaste will call EvCheckCellData to prompt the chart to update. In
    response, the state Running posts an update event to each orthogonal region, then
    discards the EvCheckCellData event */
    post_event(EvCellCycleUpdate());
    post_event(EvDifferentiationUpdate());
    post_event(EvGLD2Update());
    post_event(EvGLD1Update());
    post_event(EvLAG1Update());
    post_event(EvGLP1Update());
    return discard_event();
};





// 3) DECLARE CONSTRUCTORS FOR STATES WITH CHILDREN.
//
// In this model, we have no actions on entry for any of these states, so all the 
// constructors are empty
GLP1::GLP1( my_context ctx ):my_base( ctx ){};
LAG1::LAG1( my_context ctx ):my_base( ctx ){};
GLD1::GLD1( my_context ctx ):my_base( ctx ){};
GLD2::GLD2( my_context ctx ):my_base( ctx ){};
CellCycle::CellCycle( my_context ctx ):my_base( ctx ){};
CellCycle_Mitosis::CellCycle_Mitosis( my_context ctx ):my_base( ctx ){};
CellCycle_ExitedProlif::CellCycle_ExitedProlif( my_context ctx ):my_base( ctx ){};
Differentiation::Differentiation( my_context ctx ):my_base( ctx ){};




// 4) DEFINE THE RESPONSE TO EVENTS GIVEN BY EACH LEAF STATE (i.e. each
// state with no children). This really determines the model's behaviour.

//--------------------------------------------------------------------------
//---------------------------GLP1_Unbound-----------------------------------

GLP1_Unbound::GLP1_Unbound( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result GLP1_Unbound::react( const EvGLP1Update & ){ //On update...

    CellPtr myCell=context<FateUncoupledFromCycle>().pCell; //Get a pointer to the cell.

    if(GetDistanceFromDTC(myCell) < 35){ //If cell is more than 35 microns from DTC... 
        return transit<GLP1_Bound>();    //Transition to the bound state
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//----------------------------GLP1_Bound------------------------------------

GLP1_Bound::GLP1_Bound( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result GLP1_Bound::react( const EvGLP1Update & ){  //On update...

    CellPtr myCell=context<FateUncoupledFromCycle>().pCell; //Get a pointer to the cell.
    double DTCSignallingLength = GlobalParameterStruct::Instance()->GetParameter(24); //Check this parameter

   if(GetDistanceFromDTC(myCell) > DTCSignallingLength){ //If cell out of range of DTC signal
       return transit<GLP1_Absent>();                    //Transition GLP1 signal to absent
   }
   return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//---------------------------GLP1_Absent------------------------------------
GLP1_Absent::GLP1_Absent( my_context ctx ):my_base( ctx ){ }; //Empty constructor

sc::result GLP1_Absent::react( const EvGLP1Update & ){ // This state is unresponsive
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//-------------------------LAG1_Inactive------------------------------------
LAG1_Inactive::LAG1_Inactive( my_context ctx ):my_base( ctx ){ }; //Empty constructor

sc::result LAG1_Inactive::react( const EvLAG1Update & ){ //On update...
    
    if(state_cast<const GLP1_Bound*>()!=0){ //If GLP1 is Bound...
        return transit<LAG1_Active>();      //Transition to LAG1 Active
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//-------------------------LAG1_Active--------------------------------------
LAG1_Active::LAG1_Active( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result LAG1_Active::react( const EvLAG1Update & ){ //On update...

    if(state_cast<const GLP1_Bound*>()==0){ //If GLP1 is in any other state BUT Bound...
        return transit<LAG1_Inactive>();    //Transition to LAG1 Inactive
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------GLD1_Inactive-----------------------------------
GLD1_Inactive::GLD1_Inactive( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result GLD1_Inactive::react( const EvGLD1Update & ){ //On update...

    if(state_cast<const LAG1_Inactive*>()!=0){ //If LAG1 is Inactive...
        return transit<GLD1_Active>();         //Transition to GLD1 Active
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------GLD1_Active-------------------------------------
GLD1_Active::GLD1_Active( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result GLD1_Active::react( const EvGLD1Update & ){ //On update...

    if(state_cast<const LAG1_Active*>()!=0){ //If LAG1 is Active...
        return transit<GLD1_Inactive>();     //Transition to GLD1 Inactive
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//----------------------------------------------------------------------
//-----------------------------GLD2-------------------------------------
GLD2_Inactive::GLD2_Inactive( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result GLD2_Inactive::react( const EvGLD2Update & ){ //On update...

    if(state_cast<const LAG1_Inactive*>()!=0){ //If LAG1 is Inactive...
        return transit<GLD2_Active>();         //Transition to GLD2 Active
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//---------------------------GLD2_Active------------------------------------
GLD2_Active::GLD2_Active( my_context ctx ):my_base( ctx ){}; //Empty constructor

sc::result GLD2_Active::react( const EvGLD2Update & ){ //On update...

    if(state_cast<const LAG1_Active*>()!=0){ //If LAG1 is Active...
        return transit<GLD2_Inactive>();     //Transition to GLD2 Inactive
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------CellCycle_Mitosis_G1----------------------------
CellCycle_Mitosis_G1::CellCycle_Mitosis_G1( my_context ctx ):my_base( ctx ){ //Constructor with entry events

    CellPtr myCell = context<FateUncoupledFromCycle>().pCell;                   //Get cell pointer
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();               
    double CurrentG1 = GetG1Duration(myCell);                                   //Get G1 duration at current time
    Duration = p_gen->NormalRandomDeviate(CurrentG1, stochasticity*CurrentG1);  //Add some random variation
    context<FateUncoupledFromCycle>().TimeInPhase = 0.0;                        //Time in phase initially 0
    
    SetCellCyclePhase(myCell, G_ONE_PHASE);                                     //Set cell cycle phase to 1.0/G1
    myCell->GetCellData()->SetItem("CellCyclePhase",1.0);                       //in cell cycle model and cellData    
    myCell->GetCellData()->SetItem("DNAContent",1.0);
    
    if(contactInhibitionInG1 && GetTime()>17){                                  //If this is an adult worm, get the
        CompressionThresh=GlobalParameterStruct::Instance()->GetParameter(29);  //threshold compression volume for
        myCell->GetCellData()->SetItem("ArrestedFor",0.0);                      //contact inhibition. 
    }
}

sc::result CellCycle_Mitosis_G1::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    
    if (contactInhibitionInG1 && GetTime()>17){                                 //If the worm is an adult and contactInhibition is applied in G1... 
        double rad = myCell->GetCellData()->GetItem("Radius");                  //Compare compressed volume and relaxed volume. 
        if(myCell->GetCellData()->GetItem("volume") < CompressionThresh * 4.18879 * rad*rad*rad && CompressionThresh < 1.0){   //4.18879 = 4/3 pi
            myCell->GetCellData()->SetItem("ArrestedFor", myCell->GetCellData()->GetItem("ArrestedFor") + GetTimestep());   //If heavily compressed, increment time in arrest 
                                                                                                                            //and make no progress through G1
        }else{
            context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();
        }
    }else{
        context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();   //If not contact inhibited, increment time in phase
    }

    if(context<FateUncoupledFromCycle>().TimeInPhase >= Duration){        //If time elapsed greater than phase duration...
        return transit<CellCycle_Mitosis_S>();                            //Transit into S phase
    }
    if(GetTime()>1 && (state_cast<const GLD2_Active*>()!=0 || state_cast<const GLD1_Active*>()!=0)){ //If GLD1 / GLD2 on...
        return transit<CellCycle_ExitedProlif_G1>();                                                 //Transit into CellCycle_ExitedProlif_G1
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//------------------------CellCycle_Mitosis_S-------------------------------
CellCycle_Mitosis_S::CellCycle_Mitosis_S( my_context ctx ):my_base( ctx ){ //Constructor with entry events
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    Duration = GetSDuration(myCell);                                       //Get the current length of S phase
    context<FateUncoupledFromCycle>().TimeInPhase = 0.0;                   //Time in S phase so far = 0

    SetCellCyclePhase(myCell,S_PHASE);                                     //Set S phase labels.
    myCell->GetCellData()->SetItem("CellCyclePhase",2.0);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
}

sc::result CellCycle_Mitosis_S::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;

    myCell->GetCellData()->SetItem("DNAContent",1.0 + context<FateUncoupledFromCycle>().TimeInPhase/Duration); //Increment DNA content

    context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();         //Increment time in S phase
    if(context<FateUncoupledFromCycle>().TimeInPhase >= Duration){          //If time elapsed > S phase duration 
        return transit<CellCycle_Mitosis_G2>();                             //Transit into G2
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//----------------------CellCycle_Mitosis_G2--------------------------------
CellCycle_Mitosis_G2::CellCycle_Mitosis_G2( my_context ctx ):my_base( ctx ){   //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();       
    double CurrentG2 = GetG2Duration(myCell);                                  //Get current G2 duration
    Duration = p_gen->NormalRandomDeviate(CurrentG2, stochasticity*CurrentG2); //Add random noise to cell cycle phase length
    context<FateUncoupledFromCycle>().TimeInPhase = 0.0;                       //Current time in phase = 0.0

    SetCellCyclePhase(myCell,G_TWO_PHASE);                                     //Set cell cycle labels to G2 (3.0)
    myCell->GetCellData()->SetItem("CellCyclePhase",3.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);

    if(contactInhibitionInG2 && GetTime()>17){                                 //If worm is adult and contact inhibition enabled in G2:
        CompressionThresh=GlobalParameterStruct::Instance()->GetParameter(29); //Get compression threshold and set time arrested = 0.0
        myCell->GetCellData()->SetItem("ArrestedFor",0.0);
    }

}

sc::result CellCycle_Mitosis_G2::react( const EvCellCycleUpdate & ){    //On update...
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    
    if (contactInhibitionInG2 && GetTime()>17){                         //If worm is adult, and CI in G2 enabled
        double rad=myCell->GetCellData()->GetItem("Radius");            //Compare compressed and relaxed volumes and for heavily compressed
                                                                        //cells arrest instead on making progress through G2.
        if(myCell->GetCellData()->GetItem("volume") < CompressionThresh * 4.18879 * rad*rad*rad && CompressionThresh<1.0){    
           myCell->GetCellData()->SetItem("ArrestedFor", myCell->GetCellData()->GetItem("ArrestedFor")+GetTimestep());
        }else{
            context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();  //If cell not too compressed, progress through G2
        }
    }else{
        context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep(); //If worm is larval or no CI enabled, progress through G2 
    }

    if(context<FateUncoupledFromCycle>().TimeInPhase >= Duration){  //If time in G2 is over...
        if(myCell->GetCellData()->GetItem("IsDTC")==0.0){           //If not the DTC:
            SetReadyToDivide(myCell,true);                          //Call for a cell division from Chaste
        }
        return transit<CellCycle_Mitosis_M>();                      //Then transit into M phase
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------CellCycle_Mitosis_M-----------------------------
CellCycle_Mitosis_M::CellCycle_Mitosis_M( my_context ctx ): my_base( ctx ){    //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    Duration = GetMDuration(myCell);                                           //Get duration appropriate for M phase
    context<FateUncoupledFromCycle>().TimeInPhase = 0.0;                       //Current time in phase = 0.0
    SetCellCyclePhase(myCell, M_PHASE);                                        //Set M phase labels
    myCell->GetCellData()->SetItem("CellCyclePhase",4.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);
}

sc::result CellCycle_Mitosis_M::react( const EvCellCycleUpdate & ){             //On update...
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();             //Increment time in phase...
    myCell->GetCellData()->SetItem("DNAContent",2.0-context<FateUncoupledFromCycle>().TimeInPhase/Duration); //Decrease DNA content
    if(context<FateUncoupledFromCycle>().TimeInPhase >= Duration){              //If time in phase > M duration
        return transit<CellCycle_Mitosis_G1>();                                 //Transit into G1
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------CellCycle_ExitedProlif_G1-----------------------
CellCycle_ExitedProlif_G1::CellCycle_ExitedProlif_G1( my_context ctx ):my_base( ctx ){  //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    context<FateUncoupledFromCycle>().TimeInPhase = 0.0;                 //Time in phase is initially 0
    Duration = GetG1Duration(myCell);                                    //Get length appropriate for G1
    myCell->GetCellData()->SetItem("CellCyclePhase",1.0);                //Set cell cycle label to 1.0 G1
    myCell->GetCellData()->SetItem("DNAContent",1.0);
};
sc::result CellCycle_ExitedProlif_G1::react( const EvCellCycleUpdate & ){ //On update ...
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();       //Increment time in phase
    if(context<FateUncoupledFromCycle>().TimeInPhase > Duration ){        //If time elapsed > G1 Duration  
        return transit<CellCycle_ExitedProlif_MeioticS>();                //Transit into Meiotic S
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------CellCycle_ExitedProlif_MeioticS-----------------
CellCycle_ExitedProlif_MeioticS::CellCycle_ExitedProlif_MeioticS( my_context ctx ):my_base( ctx ){  //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    context<FateUncoupledFromCycle>().TimeInPhase = 0.0;        //Time in phase is initially 0
    Duration = GetSDuration(myCell);                            //Get length appropriate for S phase
    myCell->GetCellData()->SetItem("CellCyclePhase",2.5);       //Set cell cycle phase label to 2.5 (Meiotic S)
    myCell->GetCellData()->SetItem("DNAContent",1.0);
};
sc::result CellCycle_ExitedProlif_MeioticS::react( const EvCellCycleUpdate & ){ //On update ...
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    context<FateUncoupledFromCycle>().TimeInPhase += GetTimestep();              //Increment time in phase
    myCell->GetCellData()->SetItem("DNAContent",1.0 + context<FateUncoupledFromCycle>().TimeInPhase/Duration); //Increase DNA content
    if(context<FateUncoupledFromCycle>().TimeInPhase > Duration ){               //If time elapsed > S Duration  
        return transit<CellCycle_ExitedProlif_Meiosis>();                        //Transit into Meiosis
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//--------------------------CellCycle_ExitedProlif_Meiosis------------------
CellCycle_ExitedProlif_Meiosis::CellCycle_ExitedProlif_Meiosis( my_context ctx ):my_base( ctx ){ //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    myCell->GetCellData()->SetItem("CellCyclePhase",-1.0);                  //Set cell cycle phase label to -1 (meiosis)
    myCell->GetCellData()->SetItem("DNAContent",2.0);
};
sc::result CellCycle_ExitedProlif_Meiosis::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    if(state_cast<const Differentiation_Sperm*>()==0 && state_cast<const Differentiation_Oocyte*>()==0){    //If the cell has not undergone sex determination
        UpdateRadiusMeiotic(myCell);                                                                        //Update radius according to meiotic cell rules.
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//---------------------Differentiation_Precursor----------------------------
Differentiation_Precursor::Differentiation_Precursor(my_context ctx) :my_base(ctx){}; //Empty constructor

sc::result Differentiation_Precursor::react(const EvDifferentiationUpdate &){         //On update...
     CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    if (GetTime()>1 &&  GetDistanceFromDTC(myCell) > 200){                            //If distance from DTC > 200 microns
        if( GetTime() < GlobalParameterStruct::Instance()->GetParameter(22)){         //And time < threshold time
            return transit<Differentiation_SpermFated>();                             //Become sperm fated
        }else{
            return transit<Differentiation_OocyteFated>();                            //Otherwise become oocyte fated
        }
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//---------------------Differentiation_SpermFated----------------------------
Differentiation_SpermFated::Differentiation_SpermFated(my_context ctx) :my_base(ctx){ //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    myCell->GetCellData()->SetItem("SpermFated", 1.0);                                //Set sperm fated label to 1.0
    context<FateUncoupledFromCycle>().SpermDevelopmentDelay = 0.0;                    //Time spent becoming a sperm set to 0
};

sc::result Differentiation_SpermFated::react(const EvDifferentiationUpdate &){        //On update...
    CellPtr myCell = context<FateUncoupledFromCycle>().pCell;
    context<FateUncoupledFromCycle>().SpermDevelopmentDelay += GetTimestep();         //Increment time spent becoming a mature sperm
    if (context<FateUncoupledFromCycle>().SpermDevelopmentDelay > GlobalParameterStruct::Instance()->GetParameter(23)){ //If enough time has elapsed
        return transit<Differentiation_Sperm>();                                      //Transit into sperm state
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//---------------------Differentiation_OocyteFated----------------------------
Differentiation_OocyteFated::Differentiation_OocyteFated(my_context ctx) :my_base(ctx){ //Constructor
    CellPtr myCell=context<FateUncoupledFromCycle>().pCell;
    myCell->GetCellData()->SetItem("OocyteFated", 1.0);                                 //Set oocyte fated label to 1.0
};

sc::result Differentiation_OocyteFated::react(const EvDifferentiationUpdate &){         //On update...
    CellPtr myCell = context<FateUncoupledFromCycle>().pCell;
    if(GetDistanceFromDTC(myCell) > 250){                                               //If cell > 250 microns from DTC
        UpdateRadiusOocyte(myCell);                                                     //Grow oocyte
    }
    if(myCell->GetCellData()->GetItem("Radius")>10.0){                                  //If cell radius > 10 microns
       return transit<Differentiation_Oocyte>();                                        //Transit to oocyte state    
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//------------------------Differentiation_Sperm-----------------------------
Differentiation_Sperm::Differentiation_Sperm(my_context ctx):my_base(ctx){    //Constructor
    CellPtr myCell = context<FateUncoupledFromCycle>().pCell;
    myCell->GetCellData()->SetItem("Differentiation_Sperm", 1.0);             //Set sperm label to 1.0
};

sc::result Differentiation_Sperm::react(const EvDifferentiationUpdate &){     //On update...
    CellPtr myCell = context<FateUncoupledFromCycle>().pCell;
    if (context<FateUncoupledFromCycle>().SpermatocyteDivisions == 1){        //If 1 sperm division has occurred so far...
        SetRadius(myCell, 1.5);                                               //Set the cell radius to sperm value
        SetReadyToDivide(myCell, true);                                       //Divide one more time
        context<FateUncoupledFromCycle>().SpermatocyteDivisions++;            //Increment sperm divisions count
    }
    if (context<FateUncoupledFromCycle>().SpermatocyteDivisions == 0){        //If 0 sperm divisions have occurred so far...
        double radius = myCell->GetCellData()->GetItem("Radius");             //Halve the cell volume
        SetRadius(myCell, radius/1.26); 
        SetReadyToDivide(myCell, true);                                       //Undergo a sperm division
        context<FateUncoupledFromCycle>().SpermatocyteDivisions++;            //Increment sperm divisions count
    }
    return discard_event();
};
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//-----------------------Differentiation_Oocyte-----------------------------
Differentiation_Oocyte::Differentiation_Oocyte(my_context ctx):my_base(ctx){ //Constructor
    CellPtr myCell = context<FateUncoupledFromCycle>().pCell;
    myCell->GetCellData()->SetItem("Differentiation_Oocyte", 1.0);          //Set Oocyte label to 1.0 in CellData
};

sc::result Differentiation_Oocyte::react(const EvDifferentiationUpdate &){  //Unresponsive state
    return discard_event();
};
//--------------------------------------------------------------------------



// 5) FINALLY, DECLARE THAT StatechartCellCycleModel AND ElegansDevStatechartCellCycleModel
// CAN TAKE THIS AS A TEMPLATE PARAMETER

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(StatechartCellCycleModel, FateUncoupledFromCycle)
EXPORT_TEMPLATE_CLASS1(ElegansDevStatechartCellCycleModel, FateUncoupledFromCycle)

