#ifndef StatechartINTERFACE_HPP_
#define StatechartINTERFACE_HPP_

//Include any Chaste headers you need.
#include <Cell.hpp>
#include <AbstractCellPopulation.hpp>
#include <CellCyclePhases.hpp>
#include <StatechartCellCycleModel.hpp>

//Fill out all the functions that will get and set cell properties. Some common functions provided.

//Get cell cycle phase durations
double GetMDuration(CellPtr pCell){
    return pCell->GetCellCycleModel()->GetMDuration();
};

double GetSDuration(CellPtr pCell){
    return pCell->GetCellCycleModel()->GetSDuration();
};

double GetG1Duration(CellPtr pCell){
    return pCell->GetCellCycleModel()->GetG1Duration();
};

double GetG2Duration(CellPtr pCell){
    return pCell->GetCellCycleModel()->GetG2Duration();
};


//Setters for cell cycle model
void SetCellCyclePhase(CellPtr pCell, CellCyclePhase_ phase){
    AbstractCellCycleModel* model = pCell->GetCellCycleModel();
    dynamic_cast<StatechartCellCycleModel*>(model)->SetCellCyclePhase(phase);
}
void SetReadyToDivide(CellPtr pCell, bool Ready){
    AbstractCellCycleModel* model = pCell->GetCellCycleModel();
    dynamic_cast<StatechartCellCycleModel*>(model)->SetReadyToDivide(Ready);
};

//Misc
bool IsDead(CellPtr pCell){
     return pCell->IsDead();
};
double GetTimestep(){
     return SimulationTime::Instance()->GetTimeStep();
};
double GetTime(){
     return SimulationTime::Instance()->GetTime();
};


//Elegans Specific
void SetProliferationFlag(CellPtr pCell, double Flag){
    pCell->GetCellData()->SetItem("Proliferating",Flag);
};
void SetRadius(CellPtr pCell, double radius){
      pCell->GetCellData()->SetItem("Radius",radius);
};
double GetDistanceFromDTC(CellPtr pCell){
    return pCell->GetCellData()->GetItem("DistanceAwayFromDTC");
};
double GetMaxRadius(CellPtr pCell){
    return pCell->GetCellData()->GetItem("MaxRadius");
};

void UpdateRadiusOocyte(CellPtr pCell){
  double MaxRad = GetMaxRadius(pCell);
  double Rad = pCell->GetCellData()->GetItem("Radius");
  if(Rad<(MaxRad-0.05)){
    SetRadius(pCell,Rad+=GetTimestep()*GlobalParameterStruct::Instance()->GetParameter(11)); //1 micron per hour
  }
};
void UpdateRadiusLateMeiotic(CellPtr pCell){
  double MaxRad = GetMaxRadius(pCell);
  double Rad = pCell->GetCellData()->GetItem("Radius");
  if(Rad<fmin(MaxRad-0.05,4)){
    SetRadius(pCell,Rad+=GetTimestep()*GlobalParameterStruct::Instance()->GetParameter(12));  //0.5 micron per hour
  }
};
#endif
