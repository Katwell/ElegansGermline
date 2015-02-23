#ifndef ABSTRACTSTATECHARTCELLCYCLEMODEL_HPP_
#define ABSTRACTSTATECHARTCELLCYCLEMODEL_HPP_

class AbstractStatechartCellCycleModel{

public:

	//Two new setter methods - these expose two protected variables of AbstractCellCycleModel to the statechart:
    //the Phase- and ReadyToDivide flags; otherwise it can't set them.
    virtual void SetCellCyclePhase(CellCyclePhase_ Phase) = 0;
    virtual void SetReadyToDivide(bool Ready) = 0;  

};


#endif