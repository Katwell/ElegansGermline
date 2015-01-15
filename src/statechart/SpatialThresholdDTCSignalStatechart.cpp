#include <StatechartInterface.hpp>
#include <SpatialThresholdDTCSignalStatechart.hpp>

//--------------------STATECHART FUNCTIONS------------------------------

//HARDCODED + RARELY ALTERED DETAILS
double stochasticity         = 0.1;
bool   contactInhibitionInG1 = false;
bool   contactInhibitionInG2 = false;

CellStatechart::CellStatechart(){
        pCell=boost::shared_ptr<Cell>(); 
        //Set default value of all state associated variables
        TimeInPhase=0;
        Spermatocyte_divisions=0;
        LateMeiosisDelay=0;
        Fate=0;
};

//Set associated cell
void CellStatechart::SetCell(CellPtr newCell){
     assert(newCell!=NULL);
     pCell=newCell;
};

//Get a vector containing all state associated variables
std::vector<double> CellStatechart::GetVariables(){
    std::vector<double> variables;
    variables.push_back(TimeInPhase);
    variables.push_back(Spermatocyte_divisions);
    variables.push_back(LateMeiosisDelay);
    variables.push_back(Fate);
    return variables;
}

//Set values of all state associated variables from an input vector
void CellStatechart::SetVariables(std::vector<double> variables){
    TimeInPhase=variables.at(0);
    Spermatocyte_divisions=variables.at(1);
    LateMeiosisDelay=variables.at(2);
    Fate=variables.at(3);
}

//Get an encoding of the current state in bitset form
std::bitset<100> CellStatechart::GetState(){
  std::bitset<100> state;
 if(state_cast<const GLD2_Inactive*>()!=0){ 
     state.set(1,1);
 }
 if(state_cast<const GLD2_Active*>()!=0){ 
     state.set(2,1);
 }
 if(state_cast<const GLD1_Inactive*>()!=0){ 
     state.set(3,1);
 }
 if(state_cast<const GLD1_Active*>()!=0){ 
     state.set(4,1);
 }
 if(state_cast<const LAG1_Inactive*>()!=0){ 
     state.set(5,1);
 }
 if(state_cast<const LAG1_Active*>()!=0){ 
     state.set(6,1);
 }
 if(state_cast<const GLP1_Absent*>()!=0){ 
     state.set(7,1);
 }
 if(state_cast<const GLP1_Active*>()!=0){ 
     state.set(8,1);
 }
 if(state_cast<const GLP1_Inactive*>()!=0){ 
     state.set(9,1);
 }
 if(state_cast<const GLP1_Bound*>()!=0){ 
     state.set(10,1);
 }
 if(state_cast<const GLP1_Unbound*>()!=0){ 
     state.set(11,1);
 }
 if(state_cast<const CellCycle_Meiosis*>()!=0){ 
     state.set(12,1);
 }
 if(state_cast<const CellCycle_Mitosis_M*>()!=0){ 
     state.set(13,1);
 }
 if(state_cast<const CellCycle_Mitosis_S*>()!=0){ 
     state.set(14,1);
 }
 if(state_cast<const CellCycle_Mitosis_G2*>()!=0){ 
     state.set(15,1);
 }
 if(state_cast<const CellCycle_Mitosis_G1*>()!=0){ 
     state.set(16,1);
 }
 if(state_cast<const Differentiation_Sperm*>()!=0){ 
     state.set(17,1);
 }
 if(state_cast<const Differentiation_Oocyte*>()!=0){ 
     state.set(18,1);
 }
 if(state_cast<const Differentiation_EarlyMeiosis*>()!=0){ 
     state.set(19,1);
 }
 if(state_cast<const Differentiation_Precursor*>()!=0){ 
     state.set(20,1);
 }
 if(state_cast<const OocyteEffector_Inactive*>()!=0){ 
     state.set(21,1);
 }
 if(state_cast<const OocyteEffector_Active*>()!=0){ 
     state.set(22,1);
 }
 if(state_cast<const SpermEffector_Inactive*>()!=0){ 
     state.set(23,1);
 }
 if(state_cast<const SpermEffector_Active*>()!=0){ 
     state.set(24,1);
 }
 if(state_cast<const Differentiation_LateMeiosis*>()!=0){ 
     state.set(25,1);
 }
 return state;
}

//Set current state from a bitset vector
void CellStatechart::SetState(std::bitset<100> state){
 if(state[1]==1){ 
     process_event(EvGoToGLD2_Inactive());
 }
 if(state[2]==1){ 
     process_event(EvGoToGLD2_Active());
 }
 if(state[3]==1){ 
     process_event(EvGoToGLD1_Inactive());
 }
 if(state[4]==1){ 
     process_event(EvGoToGLD1_Active());
 }
 if(state[5]==1){ 
     process_event(EvGoToLAG1_Inactive());
 }
 if(state[6]==1){ 
     process_event(EvGoToLAG1_Active());
 }
 if(state[7]==1){ 
     process_event(EvGoToGLP1_Absent());
 }
 if(state[8]==1){ 
     process_event(EvGoToGLP1_Active());
 }
 if(state[9]==1){ 
     process_event(EvGoToGLP1_Inactive());
 }
 if(state[10]==1){ 
     process_event(EvGoToGLP1_Bound());
 }
 if(state[11]==1){ 
     process_event(EvGoToGLP1_Unbound());
 }
 if(state[12]==1){ 
     process_event(EvGoToCellCycle_Meiosis());
 }
 if(state[13]==1){ 
     process_event(EvGoToCellCycle_Mitosis_M());
 }
 if(state[14]==1){ 
     process_event(EvGoToCellCycle_Mitosis_S());
 }
 if(state[15]==1){ 
     process_event(EvGoToCellCycle_Mitosis_G2());
 }
 if(state[16]==1){ 
     process_event(EvGoToCellCycle_Mitosis_G1());
 }
 if(state[17]==1){ 
     process_event(EvGoToDifferentiation_Sperm());
 }
 if(state[18]==1){ 
     process_event(EvGoToDifferentiation_Oocyte());
 }
 if(state[19]==1){ 
     process_event(EvGoToDifferentiation_EarlyMeiosis());
 }
 if(state[20]==1){ 
     process_event(EvGoToDifferentiation_Precursor());
 }
 if(state[21]==1){ 
     process_event(EvGoToOocyteEffector_Inactive());
 }
 if(state[22]==1){ 
     process_event(EvGoToOocyteEffector_Active());
 }
 if(state[23]==1){ 
     process_event(EvGoToSpermEffector_Inactive());
 }
 if(state[24]==1){ 
     process_event(EvGoToSpermEffector_Active());
 }
 if(state[25]==1){ 
     process_event(EvGoToDifferentiation_LateMeiosis());
 }
}

//Copy an existing statechart
boost::shared_ptr<CellStatechart> CellStatechart::Copy(boost::shared_ptr<CellStatechart> myNewStatechart){
    myNewStatechart->initiate();
    if(state_cast<const GLD2_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToGLD2_Inactive());
    }
    if(state_cast<const GLD2_Active*>()!=0){
        myNewStatechart->process_event(EvGoToGLD2_Active());
    }
    if(state_cast<const GLD1_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToGLD1_Inactive());
    }
    if(state_cast<const GLD1_Active*>()!=0){
        myNewStatechart->process_event(EvGoToGLD1_Active());
    }
    if(state_cast<const LAG1_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToLAG1_Inactive());
    }
    if(state_cast<const LAG1_Active*>()!=0){
        myNewStatechart->process_event(EvGoToLAG1_Active());
    }
    if(state_cast<const GLP1_Absent*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Absent());
    }
    if(state_cast<const GLP1_Active*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Active());
    }
    if(state_cast<const GLP1_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Inactive());
    }
    if(state_cast<const GLP1_Bound*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Bound());
    }
    if(state_cast<const GLP1_Unbound*>()!=0){
        myNewStatechart->process_event(EvGoToGLP1_Unbound());
    }
    if(state_cast<const CellCycle_Meiosis*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Meiosis());
    }
    if(state_cast<const CellCycle_Mitosis_M*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_M());
    }
    if(state_cast<const CellCycle_Mitosis_S*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_S());
    }
    if(state_cast<const CellCycle_Mitosis_G2*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_G2());
    }
    if(state_cast<const CellCycle_Mitosis_G1*>()!=0){
        myNewStatechart->process_event(EvGoToCellCycle_Mitosis_G1());
    }
    if(state_cast<const Differentiation_Sperm*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_Sperm());
    }
    if(state_cast<const Differentiation_Oocyte*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_Oocyte());
    }
    if(state_cast<const Differentiation_EarlyMeiosis*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_EarlyMeiosis());
    }
    if(state_cast<const Differentiation_Precursor*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_Precursor());
    }
    if(state_cast<const OocyteEffector_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToOocyteEffector_Inactive());
    }
    if(state_cast<const OocyteEffector_Active*>()!=0){
        myNewStatechart->process_event(EvGoToOocyteEffector_Active());
    }
    if(state_cast<const SpermEffector_Inactive*>()!=0){
        myNewStatechart->process_event(EvGoToSpermEffector_Inactive());
    }
    if(state_cast<const SpermEffector_Active*>()!=0){
        myNewStatechart->process_event(EvGoToSpermEffector_Active());
    }
    if(state_cast<const Differentiation_LateMeiosis*>()!=0){
        myNewStatechart->process_event(EvGoToDifferentiation_LateMeiosis());
    }
    myNewStatechart->Spermatocyte_divisions=this->Spermatocyte_divisions;
    myNewStatechart->LateMeiosisDelay=this->LateMeiosisDelay;
    myNewStatechart->Fate=this->Fate;
    return (myNewStatechart);
};

//--------------------FIRST RESPONDER------------------------------
//Describes the events that are fired when an update of the chart is called for

sc::result Running::react( const EvCheckCellData & ){
    post_event(EvCellCycleUpdate());
    post_event(EvDifferentiationUpdate());
    post_event(EvOocyteEffectorUpdate());
    post_event(EvSpermEffectorUpdate());
    post_event(EvGLD2Update());
    post_event(EvGLD1Update());
    post_event(EvLAG1Update());
    post_event(EvGLP1Update());
    return discard_event();
};



//Constructors for empty container states
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
GLD2::GLD2( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
GLD1::GLD1( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
LAG1::LAG1( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
GLP1::GLP1( my_context ctx ):my_base( ctx ){
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
Differentiation::Differentiation( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
OocyteEffector::OocyteEffector( my_context ctx ):my_base( ctx ){
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
SpermEffector::SpermEffector( my_context ctx ):my_base( ctx ){
};



//Leaf states that actually do something...

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

//--------------------------------------------------------------------------
//-------------------------LAG1_Inactive------------------------------------
LAG1_Inactive::LAG1_Inactive( my_context ctx ):my_base( ctx ){};

sc::result LAG1_Inactive::react( const EvLAG1Update & ){
    if(state_cast<const GLP1_Active*>()!=0){
        return transit<LAG1_Active>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//-------------------------LAG1_Active--------------------------------------
LAG1_Active::LAG1_Active( my_context ctx ):my_base( ctx ){};

sc::result LAG1_Active::react( const EvLAG1Update & ){
    if(state_cast<const GLP1_Active*>()==0){
        return transit<LAG1_Inactive>();
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
//----------------------------GLP1_Active-----------------------------------
GLP1_Active::GLP1_Active( my_context ctx ):my_base( ctx ){};

sc::result GLP1_Active::react( const EvGLP1Update & ){
   CellPtr myCell=context<CellStatechart>().pCell;
   double t = GetTime();
   double prolifZoneLength=GlobalParameterStruct::Instance()->GetParameter(24);

   if(GetDistanceFromDTC(myCell) > prolifZoneLength){
       return transit<GLP1_Absent>();
   }
   return discard_event();
};

//--------------------------------------------------------------------------
//--------------------------GLP1_Inactive-----------------------------------
GLP1_Inactive::GLP1_Inactive( my_context ctx ):my_base( ctx ){};

sc::result GLP1_Inactive::react( const EvGLP1Update & ){
    return transit<GLP1_Active>();
};

//--------------------------------------------------------------------------
//----------------------------GLP1_Bound------------------------------------
GLP1_Bound::GLP1_Bound( my_context ctx ):my_base( ctx ){};

sc::result GLP1_Bound::react( const EvGLP1Update & ){
    return transit<GLP1_Inactive>();
};

//--------------------------------------------------------------------------
//---------------------------GLP1_Unbound-----------------------------------
GLP1_Unbound::GLP1_Unbound( my_context ctx ):my_base( ctx ){};

sc::result GLP1_Unbound::react( const EvGLP1Update & ){
    CellPtr myCell=context<CellStatechart>().pCell;
    if(GetDistanceFromDTC(myCell)<35){
        return transit<GLP1_Bound>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//--------------------------CellCycle_Meiosis-------------------------------
CellCycle_Meiosis::CellCycle_Meiosis( my_context ctx ):my_base( ctx ){ 
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase = 0.0;
    CurrentG1 = GetG1Duration(myCell);
    CurrentS = GetSDuration(myCell);
};

sc::result CellCycle_Meiosis::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase+=GetTimestep();

    if(state_cast<const Differentiation_EarlyMeiosis*>()!=0){
        if(context<CellStatechart>().TimeInPhase < CurrentG1 ){
            myCell->GetCellData()->SetItem("CellCyclePhase",1.0);
            myCell->GetCellData()->SetItem("DNAContent",1.0);
        }
        if(context<CellStatechart>().TimeInPhase > CurrentG1 &&
            context<CellStatechart>().TimeInPhase < (CurrentG1+CurrentS) ){
            myCell->GetCellData()->SetItem("CellCyclePhase",2.5);
            myCell->GetCellData()->SetItem("DNAContent",1.0+context<CellStatechart>().TimeInPhase/CurrentS);
        }
        if(context<CellStatechart>().TimeInPhase > (CurrentG1+CurrentS) ){
            myCell->GetCellData()->SetItem("CellCyclePhase",-1.0);
            myCell->GetCellData()->SetItem("DNAContent",2.0);
        }
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//--------------------------CellCycle_Mitosis_M-----------------------------
CellCycle_Mitosis_M::CellCycle_Mitosis_M( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase = 0.0;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double CurrentM = GetMDuration(myCell);
    //HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Duration = p_gen->NormalRandomDeviate(CurrentM,stochasticity*CurrentM);
    Duration = CurrentM;
    //HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SetCellCyclePhase(myCell, M_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",4.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);
}

sc::result CellCycle_Mitosis_M::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase+=GetTimestep();
    myCell->GetCellData()->SetItem("DNAContent",2.0-context<CellStatechart>().TimeInPhase/Duration);

    if(context<CellStatechart>().TimeInPhase>=Duration){
        return transit<CellCycle_Mitosis_G1>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//------------------------CellCycle_Mitosis_S-------------------------------
CellCycle_Mitosis_S::CellCycle_Mitosis_S( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase = 0.0;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double CurrentS = GetSDuration(myCell);
    //HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Duration=p_gen->NormalRandomDeviate(CurrentS,stochasticity*CurrentS);
    Duration = CurrentS;
    //HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SetCellCyclePhase(myCell,S_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",2.0);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
}

sc::result CellCycle_Mitosis_S::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;
    myCell->GetCellData()->SetItem("DNAContent",1.0+context<CellStatechart>().TimeInPhase/Duration);

    context<CellStatechart>().TimeInPhase+=GetTimestep();
    if(context<CellStatechart>().TimeInPhase>=Duration){
        return transit<CellCycle_Mitosis_G2>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//----------------------CellCycle_Mitosis_G2--------------------------------
CellCycle_Mitosis_G2::CellCycle_Mitosis_G2( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase = 0.0;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double CurrentG2 = GetG2Duration(myCell);
    Duration=p_gen->NormalRandomDeviate(CurrentG2,stochasticity*CurrentG2);
    SetCellCyclePhase(myCell,G_TWO_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",3.0);
    myCell->GetCellData()->SetItem("DNAContent",2.0);
    CompressionThresh=GlobalParameterStruct::Instance()->GetParameter(29);
    myCell->GetCellData()->SetItem("ArrestedFor",0.0);
}

sc::result CellCycle_Mitosis_G2::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;

    if (contactInhibitionInG2 && GetTime()>17){
        double rad=myCell->GetCellData()->GetItem("Radius");
        //4.18879 = 4/3 * pi
        if(myCell->GetCellData()->GetItem("volume")<CompressionThresh*4.18879*rad*rad*rad && CompressionThresh<1.0){    
            myCell->GetCellData()->SetItem("ArrestedFor", myCell->GetCellData()->GetItem("ArrestedFor")+GetTimestep());
        }else{
            context<CellStatechart>().TimeInPhase+=GetTimestep();
        }
    }else{
        context<CellStatechart>().TimeInPhase += GetTimestep();
    }

    if(context<CellStatechart>().TimeInPhase>=Duration){
        if(myCell->GetCellData()->GetItem("IsDTC")==0.0){
            SetReadyToDivide(myCell,true);
        }
        return transit<CellCycle_Mitosis_M>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//--------------------------CellCycle_Mitosis_G1----------------------------
CellCycle_Mitosis_G1::CellCycle_Mitosis_G1( my_context ctx ):
my_base( ctx ){ 
    CellPtr myCell=context<CellStatechart>().pCell;
    context<CellStatechart>().TimeInPhase = 0.0;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
    double CurrentG1 = GetG1Duration(myCell);
    Duration=p_gen->NormalRandomDeviate(CurrentG1,stochasticity*CurrentG1);
    SetCellCyclePhase(myCell,G_ONE_PHASE);
    myCell->GetCellData()->SetItem("CellCyclePhase",1.0);
    myCell->GetCellData()->SetItem("DNAContent",1.0);
    CompressionThresh=GlobalParameterStruct::Instance()->GetParameter(29);
    myCell->GetCellData()->SetItem("ArrestedFor",0.0);
}

sc::result CellCycle_Mitosis_G1::react( const EvCellCycleUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;

    if (contactInhibitionInG1 && GetTime()>17){
        double rad=myCell->GetCellData()->GetItem("Radius");
        if(myCell->GetCellData()->GetItem("volume")<CompressionThresh*4.18879*rad*rad*rad && CompressionThresh<1.0){
            myCell->GetCellData()->SetItem("ArrestedFor", myCell->GetCellData()->GetItem("ArrestedFor")+GetTimestep());
        }else{
            context<CellStatechart>().TimeInPhase+=GetTimestep();
        }
    }else{
        context<CellStatechart>().TimeInPhase+=GetTimestep();
    }

    if(context<CellStatechart>().TimeInPhase>=Duration){
        return transit<CellCycle_Mitosis_S>();
    }
    if(GetTime()>1 && (state_cast<const GLD2_Active*>()!=0 || state_cast<const GLD1_Active*>()!=0)){
        return transit<CellCycle_Meiosis>();
    }
    return discard_event();
};


//--------------------------------------------------------------------------
//---------------------Differentiation_Precursor----------------------------
Differentiation_Precursor::Differentiation_Precursor(my_context ctx) :my_base(ctx){};

sc::result Differentiation_Precursor::react(const EvDifferentiationUpdate &){
    if (state_cast<const CellCycle_Meiosis*>() != 0){
        return transit<Differentiation_EarlyMeiosis>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//----------------------Differentiation_EarlyMeiosis------------------------
Differentiation_EarlyMeiosis::Differentiation_EarlyMeiosis( my_context ctx ):my_base( ctx ){};

sc::result Differentiation_EarlyMeiosis::react( const EvDifferentiationUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;
    if(myCell->GetCellData()->GetItem("CellCyclePhase")==-1.0){
        return transit<Differentiation_LateMeiosis>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//-------------------------Differentiation_LateMeiosis----------------------
Differentiation_LateMeiosis::Differentiation_LateMeiosis( my_context ctx ):my_base( ctx ){ 
    
    CellPtr myCell=context<CellStatechart>().pCell;
    myCell->GetCellData()->SetItem("CellCyclePhase",-1.0);
    context<CellStatechart>().LateMeiosisDelay = 0.0;
    
    context<CellStatechart>().Fate = 0;
    if(state_cast<const SpermEffector_Active*>()!=0){
        context<CellStatechart>().Fate = 1;
        myCell->GetCellData()->SetItem("SpermFated", 1.0);
    }else if(state_cast<const OocyteEffector_Active*>()!=0){
        context<CellStatechart>().Fate = 2;
        myCell->GetCellData()->SetItem("OocyteFated", 1.0);
    }
};

sc::result Differentiation_LateMeiosis::react( const EvDifferentiationUpdate & ){
    CellPtr myCell=context<CellStatechart>().pCell;
    
    context<CellStatechart>().LateMeiosisDelay += GetTimestep();

    //if(GetDistanceFromDTC(myCell)>250 && context<CellStatechart>().Fate==2){
    //    UpdateRadiusOocyte(myCell);
    //}else{
        UpdateRadiusLateMeiotic(myCell);
    //}

    double delay = GlobalParameterStruct::Instance()->GetParameter(23);
    if(context<CellStatechart>().Fate==1 && context<CellStatechart>().LateMeiosisDelay > delay){
        return transit<Differentiation_Sperm>();
    }
    if(context<CellStatechart>().Fate==2 && GetDistanceFromDTC(myCell)>250){
        return transit<Differentiation_Oocyte>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//-----------------------OocyteEffector_Inactive----------------------------
OocyteEffector_Inactive::OocyteEffector_Inactive(my_context ctx) :my_base(ctx){};

sc::result OocyteEffector_Inactive::react(const EvOocyteEffectorUpdate &){
    if (GetTime()>GlobalParameterStruct::Instance()->GetParameter(22)){
        return transit<OocyteEffector_Active>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//--------------------------OocyteEffector_Active---------------------------
OocyteEffector_Active::OocyteEffector_Active(my_context ctx) :my_base(ctx){};

sc::result OocyteEffector_Active::react(const EvOocyteEffectorUpdate &){
    if (GetTime()<GlobalParameterStruct::Instance()->GetParameter(22)){
        return transit<OocyteEffector_Inactive>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//-------------------------SpermEffector_Inactive---------------------------
SpermEffector_Inactive::SpermEffector_Inactive(my_context ctx) :my_base(ctx){};

sc::result SpermEffector_Inactive::react(const EvSpermEffectorUpdate &){
    if (GetTime()<GlobalParameterStruct::Instance()->GetParameter(22)){
        return transit<SpermEffector_Active>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//----------------------------SpermEffector_Active--------------------------
SpermEffector_Active::SpermEffector_Active(my_context ctx) :my_base(ctx){};

sc::result SpermEffector_Active::react(const EvSpermEffectorUpdate &){
    if (GetTime()>GlobalParameterStruct::Instance()->GetParameter(22)){
        return transit<SpermEffector_Inactive>();
    }
    return discard_event();
};

//--------------------------------------------------------------------------
//------------------------Differentiation_Sperm-----------------------------
Differentiation_Sperm::Differentiation_Sperm(my_context ctx) :
my_base(ctx){
    CellPtr myCell = context<CellStatechart>().pCell;
    myCell->GetCellData()->SetItem("Differentiation_Sperm", 1.0);
    myCell->GetCellData()->SetItem("CellCyclePhase", -1.0);
};

sc::result Differentiation_Sperm::react(const EvDifferentiationUpdate &){
    CellPtr myCell = context<CellStatechart>().pCell;
    double radius = myCell->GetCellData()->GetItem("Radius");

    if (context<CellStatechart>().Spermatocyte_divisions == 1){
        SetReadyToDivide(myCell, true);
        //Half the volume
        SetRadius(myCell, 1.5);
        context<CellStatechart>().Spermatocyte_divisions++;
    }
    if (context<CellStatechart>().Spermatocyte_divisions == 0){
        SetReadyToDivide(myCell, true);
        //Half the volume
        SetRadius(myCell, radius / 1.26);
        context<CellStatechart>().Spermatocyte_divisions++;
    }
    //if (context<CellStatechart>().Spermatocyte_divisions > 1 && radius>1.5){
    //  SetRadius(myCell, radius-0.5*GetTimestep());
    //}
    return discard_event();
};

//--------------------------------------------------------------------------
//-----------------------Differentiation_Oocyte-----------------------------
Differentiation_Oocyte::Differentiation_Oocyte(my_context ctx) :
my_base(ctx){
    CellPtr myCell = context<CellStatechart>().pCell;
    myCell->GetCellData()->SetItem("Differentiation_Oocyte", 1.0);
    myCell->GetCellData()->SetItem("CellCyclePhase", -1.0);
};

sc::result Differentiation_Oocyte::react(const EvDifferentiationUpdate &){
    CellPtr myCell = context<CellStatechart>().pCell;
    if (GetDistanceFromDTC(myCell)>250){
        UpdateRadiusOocyte(myCell);
    }
    return discard_event();
};
