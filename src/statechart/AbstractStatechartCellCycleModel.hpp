#ifndef ABSTRACTSTATECHARTCELLCYCLEMODEL_HPP_
#define ABSTRACTSTATECHARTCELLCYCLEMODEL_HPP_

/*
* This class provides a guarantee that StatechartCellCycleModel will implement two setter methods: 
* SetCellCyclePhase and SetReadyToDivide; used by the Statechart to control aspects of cell behaviour.
* 
* Basically we only need this class so that StatechartInterface can "test" whether a CellCycleModel
* exposes the right setter methods, without having to faff around referencing a templated class.
*/

class AbstractStatechartCellCycleModel {

public:

	//Two setter methods that expose normally protected variables of AbstractCellCycleModel to the statechart:
    virtual void SetCellCyclePhase(CellCyclePhase_ Phase) = 0;
    virtual void SetReadyToDivide(bool Ready) = 0;  

    virtual ~AbstractStatechartCellCycleModel(){};

};

#endif /*ABSTRACTSTATECHARTCELLCYCLEMODEL_HPP_*/