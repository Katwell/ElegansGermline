#ifndef StatechartINTERFACE_HPP_
#define StatechartINTERFACE_HPP_

#include <Cell.hpp>
#include <AbstractCellPopulation.hpp>
#include <AbstractStatechartCellCycleModel.hpp>
#include <GlobalParameterStruct.hpp>
#include <CellCyclePhases.hpp>

/*
* Implements some common functions that may be needed by many statechart models of cell
* behaviour.
*/

//CONTROLLING THE CELL CYCLE:

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
void SetCellCyclePhase(CellPtr pCell, CellCyclePhase_ phase){
    AbstractCellCycleModel* model = pCell->GetCellCycleModel();
    dynamic_cast<AbstractStatechartCellCycleModel*>(model)->SetCellCyclePhase(phase);
}
void SetReadyToDivide(CellPtr pCell, bool Ready){
    AbstractCellCycleModel* model = pCell->GetCellCycleModel();
    dynamic_cast<AbstractStatechartCellCycleModel*>(model)->SetReadyToDivide(Ready);
};


//MISC
bool IsDead(CellPtr pCell){
     return pCell->IsDead();
};
double GetTimestep(){
     return SimulationTime::Instance()->GetTimeStep();
};
double GetTime(){
     return SimulationTime::Instance()->GetTime();
};


//C ELEGANS SPECIFIC

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
    return pCell->GetCellData()->GetItem("MaxRadius"); // Max radius that will fit in the gonad.
};

//grows cell, provided it is not going to end up too big to fit in the gonad.
void UpdateRadiusOocyte(CellPtr pCell){
  double MaxRad = GetMaxRadius(pCell);
  double Rad = pCell->GetCellData()->GetItem("Radius");
  if(Rad<(MaxRad-0.05)){
    SetRadius(pCell,Rad+=GetTimestep()*GlobalParameterStruct::Instance()->GetParameter(11)); //1 micron per hour
  }
};

//grows cell, provided it is not going to end up too big to fit in the gonad, or larger than the max meiotic
//cell radius (Parameter 38).
void UpdateRadiusMeiotic(CellPtr pCell){
  double MaxRad = GetMaxRadius(pCell);
  double Rad = pCell->GetCellData()->GetItem("Radius");
  if(Rad<fmin(MaxRad-0.05,GlobalParameterStruct::Instance()->GetParameter(38))){
    SetRadius(pCell,Rad+=GetTimestep()*GlobalParameterStruct::Instance()->GetParameter(10));  //1.0 micron per hour
  }
};
#endif
