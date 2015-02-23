#include <StatechartInterface.hpp>
#include <FateDecisionCoupledToCycle.hpp>

//--------------------STATECHART FUNCTIONS------------------------------

//HARDCODED + RARELY ALTERED DETAILS
double stochasticity         = 0.1;
bool   contactInhibitionInG1 = false;
bool   contactInhibitionInG2 = true;

FateDecisionCoupledToCycle::FateDecisionCoupledToCycle(){
        pCell=boost::shared_ptr<Cell>(); 
        //Set default value of all state associated variables
        TimeInPhase=0;
        SpermatocyteDivisions=0;
        SpermDevelopmentDelay=0;
};

//Set associated cell
void FateDecisionCoupledToCycle::SetCell(CellPtr newCell){
     assert(newCell!=NULL);
     pCell=newCell;
};

//Get a vector containing all state associated variables
std::vector<double> FateDecisionCoupledToCycle::GetVariables(){
    std::vector<double> variables;
    variables.push_back(TimeInPhase);
    variables.push_back(SpermatocyteDivisions);
    variables.push_back(SpermDevelopmentDelay);
    return variables;
}

//Set values of all state associated variables from an input vector
void FateDecisionCoupledToCycle::SetVariables(std::vector<double> variables){
    TimeInPhase=variables.at(0);
    SpermatocyteDivisions=variables.at(1);
    SpermDevelopmentDelay=variables.at(2);
}

//Get an encoding of the current state in bitset form
std::bitset<32> FateDecisionCoupledToCycle::GetState(){
  std::bitset<32> state;
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

//Set current state from a bitset vector
void FateDecisionCoupledToCycle::SetState(std::bitset<32> state){
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

//Copy an existing statechart -----> HERE!!!!!
boost::shared_ptr<FateDecisionCoupledToCycle> FateDecisionCoupledToCycle::Copy(boost::shared_ptr<FateDecisionCoupledToCycle> myNewStatechart){
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
    myNewStatechart->SpermatocyteDivisions=this->SpermatocyteDivisions;
    myNewStatechart->SpermDevelopmentDelay=this->SpermDevelopmentDelay;
    return (myNewStatechart);
};

//--------------------FIRST RESPONDER------------------------------
//Describes the events that are fired when an update of the chart is called for

sc::result Running::react( const EvCheckCellData & ){
    post_event(EvCellCycleUpdate());
    post_event(EvDifferentiationUpdate());
    post_event(EvGLD2Update());
    post_event(EvGLD1Update());
    post_event(EvLAG1Update());
    post_event(EvGLP1Update());
    return discard_event();
};


//Constructors for empty container states
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
GLP1::GLP1( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
LAG1::LAG1( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
GLD1::GLD1( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
GLD2::GLD2( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
CellCycle::CellCycle( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
CellCycle_Mitosis::CellCycle_Mitosis( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
CellCycle_ExitedProlif::CellCycle_ExitedProlif( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
Differentiation::Differentiation( my_context ctx ):my_base( ctx ){
};


//Leaf states that actually do something...
//--------------------------------------------------------------------------
//---------------------------GLP1_Unbound-----------------------------------
GLP1_Unbound::GLP1_Unbound( my_context ctx ):my_base( ctx ){};

sc::result GLP1_Unbound::react( const EvGLP1Update & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;

    if(GetDistanceFromDTC(myCell) < 35){
        return transit<GLP1_Bound>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//----------------------------GLP1_Bound------------------------------------
GLP1_Bound::GLP1_Bound( my_context ctx ):my_base( ctx ){};

sc::result GLP1_Bound::react( const EvGLP1Update & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    double t = GetTime();
    double prolifZoneLength = GlobalParameterStruct::Instance()->GetParameter(24);

   if(GetDistanceFromDTC(myCell) > prolifZoneLength){
       return transit<GLP1_Absent>();
   }
   return discard_event();
};
//--------------------------------------------------------------------------
//---------------------------GLP1_Absent------------------------------------
GLP1_Absent::GLP1_Absent( my_context ctx ):my_base( ctx ){ };

sc::result GLP1_Absent::react( const EvGLP1Update & ){
    return discard_event();
};
//--------------------------------------------------------------------------
//-------------------------LAG1_Inactive------------------------------------
LAG1_Inactive::LAG1_Inactive( my_context ctx ):my_base( ctx ){};

sc::result LAG1_Inactive::react( const EvLAG1Update & ){
    if(state_cast<const GLP1_Bound*>()!=0){
        return transit<LAG1_Active>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//-------------------------LAG1_Active--------------------------------------
LAG1_Active::LAG1_Active( my_context ctx ):my_base( ctx ){};

sc::result LAG1_Active::react( const EvLAG1Update & ){
    if(state_cast<const GLP1_Bound*>()==0){
        return transit<LAG1_Inactive>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------GLD1_Inactive-----------------------------------
GLD1_Inactive::GLD1_Inactive( my_context ctx ):my_base( ctx ){};

sc::result GLD1_Inactive::react( const EvGLD1Update & ){
    if(state_cast<const LAG1_Inactive*>()!=0){
        return transit<GLD1_Active>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------GLD1_Active-------------------------------------
GLD1_Active::GLD1_Active( my_context ctx ):my_base( ctx ){};

sc::result GLD1_Active::react( const EvGLD1Update & ){
    if(state_cast<const LAG1_Active*>()!=0){
        return transit<GLD1_Inactive>();
    }
    return discard_event();
};
//----------------------------------------------------------------------
//-----------------------------GLD2-------------------------------------
GLD2_Inactive::GLD2_Inactive( my_context ctx ):my_base( ctx ){};

sc::result GLD2_Inactive::react( const EvGLD2Update & ){
    if(state_cast<const LAG1_Inactive*>()!=0){
        return transit<GLD2_Active>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//---------------------------GLD2_Active------------------------------------
GLD2_Active::GLD2_Active( my_context ctx ):my_base( ctx ){};

sc::result GLD2_Active::react( const EvGLD2Update & ){
    if(state_cast<const LAG1_Active*>()!=0){
        return transit<GLD2_Inactive>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------CellCycle_Mitosis_G1----------------------------
CellCycle_Mitosis_G1::CellCycle_Mitosis_G1( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
    double CurrentG1 = GetG1Duration(myCell);
    Duration = p_gen->NormalRandomDeviate(CurrentG1, stochasticity*CurrentG1);
    context<FateDecisionCoupledToCycle>().TimeInPhase = 0.0;
    
    SetCellCyclePhase(myCell, G_ONE_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",1.0);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
    
    if(contactInhibitionInG1 && GetTime()>17){
        CompressionThresh=GlobalParameterStruct::Instance()->GetParameter(29);
        myCell->GetCellData()->SetItem("ArrestedFor",0.0);
    }
}

sc::result CellCycle_Mitosis_G1::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    if (contactInhibitionInG1 && GetTime()>17){
        double rad = myCell->GetCellData()->GetItem("Radius");
        if(myCell->GetCellData()->GetItem("volume") < CompressionThresh * 4.18879 * rad*rad*rad && CompressionThresh<1.0){   //4.18879 = 4/3 pi
            myCell->GetCellData()->SetItem("ArrestedFor", myCell->GetCellData()->GetItem("ArrestedFor") + GetTimestep());
        }else{
            context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
        }
    }else{
        context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
    }

    if(context<FateDecisionCoupledToCycle>().TimeInPhase >= Duration){
        return transit<CellCycle_Mitosis_S>();
    }
    if(GetTime()>1 && (state_cast<const GLD2_Active*>()!=0 || state_cast<const GLD1_Active*>()!=0)){
        return transit<CellCycle_ExitedProlif_G1>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//------------------------CellCycle_Mitosis_S-------------------------------
CellCycle_Mitosis_S::CellCycle_Mitosis_S( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    Duration = GetSDuration(myCell);
    context<FateDecisionCoupledToCycle>().TimeInPhase = 0.0;

    SetCellCyclePhase(myCell,S_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",2.0);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
}

sc::result CellCycle_Mitosis_S::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    myCell->GetCellData()->SetItem("DNAContent",1.0 + context<FateDecisionCoupledToCycle>().TimeInPhase/Duration);

    context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
    if(context<FateDecisionCoupledToCycle>().TimeInPhase >= Duration){
        return transit<CellCycle_Mitosis_G2>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//----------------------CellCycle_Mitosis_G2--------------------------------
CellCycle_Mitosis_G2::CellCycle_Mitosis_G2( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double CurrentG2 = GetG2Duration(myCell);
    Duration = p_gen->NormalRandomDeviate(CurrentG2, stochasticity*CurrentG2);
    context<FateDecisionCoupledToCycle>().TimeInPhase = 0.0;

    SetCellCyclePhase(myCell,G_TWO_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",3.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);

    if(contactInhibitionInG2 && GetTime()>17){
        CompressionThresh=GlobalParameterStruct::Instance()->GetParameter(29);
        myCell->GetCellData()->SetItem("ArrestedFor",0.0);
    }
}

sc::result CellCycle_Mitosis_G2::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    if (contactInhibitionInG2 && GetTime()>17){
        double rad=myCell->GetCellData()->GetItem("Radius");
        if(myCell->GetCellData()->GetItem("volume") < CompressionThresh * 4.18879 * rad*rad*rad && CompressionThresh<1.0){    
           myCell->GetCellData()->SetItem("ArrestedFor", myCell->GetCellData()->GetItem("ArrestedFor")+GetTimestep());
        }else{
            context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
        }
    }else{
        context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
    }

    if(context<FateDecisionCoupledToCycle>().TimeInPhase >= Duration){
        if(myCell->GetCellData()->GetItem("IsDTC")==0.0){
            SetReadyToDivide(myCell,true);
        }
        return transit<CellCycle_Mitosis_M>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------CellCycle_Mitosis_M-----------------------------
CellCycle_Mitosis_M::CellCycle_Mitosis_M( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    Duration = GetMDuration(myCell);
    context<FateDecisionCoupledToCycle>().TimeInPhase = 0.0;
    SetCellCyclePhase(myCell, M_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",4.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);
}

sc::result CellCycle_Mitosis_M::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
    myCell->GetCellData()->SetItem("DNAContent",2.0-context<FateDecisionCoupledToCycle>().TimeInPhase/Duration);
    if(context<FateDecisionCoupledToCycle>().TimeInPhase >= Duration){
        return transit<CellCycle_Mitosis_G1>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------CellCycle_ExitedProlif_G1-----------------------
CellCycle_ExitedProlif_G1::CellCycle_ExitedProlif_G1( my_context ctx ):my_base( ctx ){ 
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    context<FateDecisionCoupledToCycle>().TimeInPhase = 0.0;
    Duration = GetG1Duration(myCell);
    myCell->GetCellData()->SetItem("CellCyclePhase",1.0);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
};
sc::result CellCycle_ExitedProlif_G1::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
    if(context<FateDecisionCoupledToCycle>().TimeInPhase > Duration ){
        return transit<CellCycle_ExitedProlif_MeioticS>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------CellCycle_ExitedProlif_MeioticS-----------------
CellCycle_ExitedProlif_MeioticS::CellCycle_ExitedProlif_MeioticS( my_context ctx ):my_base( ctx ){ 
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    context<FateDecisionCoupledToCycle>().TimeInPhase = 0.0;
    Duration = GetSDuration(myCell);
    myCell->GetCellData()->SetItem("CellCyclePhase",2.5);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
};
sc::result CellCycle_ExitedProlif_MeioticS::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    context<FateDecisionCoupledToCycle>().TimeInPhase += GetTimestep();
    myCell->GetCellData()->SetItem("DNAContent",1.0 + context<FateDecisionCoupledToCycle>().TimeInPhase/Duration);
    if(context<FateDecisionCoupledToCycle>().TimeInPhase > Duration ){
        return transit<CellCycle_ExitedProlif_Meiosis>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//--------------------------CellCycle_ExitedProlif_Meiosis------------------
CellCycle_ExitedProlif_Meiosis::CellCycle_ExitedProlif_Meiosis( my_context ctx ):my_base( ctx ){ 
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    myCell->GetCellData()->SetItem("CellCyclePhase",-1.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);
};
sc::result CellCycle_ExitedProlif_Meiosis::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    if(state_cast<const Differentiation_Sperm*>()==0 && state_cast<const Differentiation_Oocyte*>()==0){
        UpdateRadiusMeiotic(myCell);
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//---------------------Differentiation_Precursor----------------------------
Differentiation_Precursor::Differentiation_Precursor(my_context ctx) :my_base(ctx){};

sc::result Differentiation_Precursor::react(const EvDifferentiationUpdate &){
     CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    if(state_cast<const CellCycle_ExitedProlif_Meiosis*>()!=0){
        if( GetTime() < GlobalParameterStruct::Instance()->GetParameter(22)){
            return transit<Differentiation_SpermFated>();
        }else{
            return transit<Differentiation_OocyteFated>();
        }
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//---------------------Differentiation_SpermFated----------------------------
Differentiation_SpermFated::Differentiation_SpermFated(my_context ctx) :my_base(ctx){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    myCell->GetCellData()->SetItem("SpermFated", 1.0);
    context<FateDecisionCoupledToCycle>().SpermDevelopmentDelay = 0.0;
};

sc::result Differentiation_SpermFated::react(const EvDifferentiationUpdate &){
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    context<FateDecisionCoupledToCycle>().SpermDevelopmentDelay += GetTimestep();
    if (context<FateDecisionCoupledToCycle>().SpermDevelopmentDelay > GlobalParameterStruct::Instance()->GetParameter(23)){
        return transit<Differentiation_Sperm>();
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//---------------------Differentiation_OocyteFated----------------------------
Differentiation_OocyteFated::Differentiation_OocyteFated(my_context ctx) :my_base(ctx){
    CellPtr myCell=context<FateDecisionCoupledToCycle>().pCell;
    myCell->GetCellData()->SetItem("OocyteFated", 1.0);
};

sc::result Differentiation_OocyteFated::react(const EvDifferentiationUpdate &){
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    if(//myCell->GetCellData()->GetItem("InProximalArm")==1.0){
        GetDistanceFromDTC(myCell) > 250){
        UpdateRadiusOocyte(myCell);
    }
    if(myCell->GetCellData()->GetItem("Radius")>10.0){
       return transit<Differentiation_Oocyte>(); 
    }
    return discard_event();
};
//--------------------------------------------------------------------------
//------------------------Differentiation_Sperm-----------------------------
Differentiation_Sperm::Differentiation_Sperm(my_context ctx) :
my_base(ctx){
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    myCell->GetCellData()->SetItem("Differentiation_Sperm", 1.0);
};

sc::result Differentiation_Sperm::react(const EvDifferentiationUpdate &){
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    if (context<FateDecisionCoupledToCycle>().SpermatocyteDivisions == 1){
        SetRadius(myCell, 1.5);
        SetReadyToDivide(myCell, true);
        context<FateDecisionCoupledToCycle>().SpermatocyteDivisions++;
    }
    if (context<FateDecisionCoupledToCycle>().SpermatocyteDivisions == 0){
        double radius = myCell->GetCellData()->GetItem("Radius");
        SetRadius(myCell, radius/1.26); //Half the volume
        SetReadyToDivide(myCell, true);
        context<FateDecisionCoupledToCycle>().SpermatocyteDivisions++;
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//-----------------------Differentiation_Oocyte-----------------------------
Differentiation_Oocyte::Differentiation_Oocyte(my_context ctx) :
my_base(ctx){
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    myCell->GetCellData()->SetItem("Differentiation_Oocyte", 1.0);
};

sc::result Differentiation_Oocyte::react(const EvDifferentiationUpdate &){
    CellPtr myCell = context<FateDecisionCoupledToCycle>().pCell;
    return discard_event();
};


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS1(StatechartCellCycleModel, FateDecisionCoupledToCycle)
