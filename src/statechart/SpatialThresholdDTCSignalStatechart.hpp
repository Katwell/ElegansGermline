#ifndef SPATIALTHRESHOLDDTCSIGNALSTATECHART_HPP_
#define SPATIALTHRESHOLDDTCSIGNALSTATECHART_HPP_

#include <boost/statechart/event.hpp>
#include <boost/statechart/state_machine.hpp>
#include <boost/statechart/simple_state.hpp>
#include <boost/statechart/state.hpp>
#include <boost/statechart/custom_reaction.hpp>
#include <boost/statechart/transition.hpp>
#include <boost/mpl/list.hpp>
#include <boost/config.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include <bitset>

namespace sc = boost::statechart;
namespace mpl = boost::mpl;

//Predeclare struct for each state
struct Running;
struct GLD2;
struct GLD2_Inactive;
struct GLD2_Active;
struct GLD1;
struct GLD1_Inactive;
struct GLD1_Active;
struct LAG1;
struct LAG1_Inactive;
struct LAG1_Active;
struct GLP1;
struct GLP1_Absent;
struct GLP1_Active;
struct GLP1_Inactive;
struct GLP1_Bound;
struct GLP1_Unbound;
struct CellCycle;
struct CellCycle_Meiosis;
struct CellCycle_Mitosis;
struct CellCycle_Mitosis_M;
struct CellCycle_Mitosis_S;
struct CellCycle_Mitosis_G2;
struct CellCycle_Mitosis_G1;
struct Differentiation;
struct Differentiation_Sperm;
struct Differentiation_Oocyte;
struct Differentiation_EarlyMeiosis;
struct Differentiation_LateMeiosis;
struct Differentiation_Precursor;
struct OocyteEffector;
struct OocyteEffector_Inactive;
struct OocyteEffector_Active;
struct SpermEffector;
struct SpermEffector_Inactive;
struct SpermEffector_Active;

//Declare the update events
struct EvCheckCellData :  sc::event< EvCheckCellData > {};
struct EvCellCycleUpdate : sc::event< EvCellCycleUpdate > {};
struct EvDifferentiationUpdate : sc::event< EvDifferentiationUpdate > {};
struct EvOocyteEffectorUpdate : sc::event< EvOocyteEffectorUpdate > {};
struct EvSpermEffectorUpdate : sc::event< EvSpermEffectorUpdate > {};
struct EvGLD2Update : sc::event< EvGLD2Update > {};
struct EvGLD1Update : sc::event< EvGLD1Update > {};
struct EvLAG1Update : sc::event< EvLAG1Update > {};
struct EvGLP1Update : sc::event< EvGLP1Update > {};

//Declare goto events for each leaf state, to allow copying of a statechart state
struct EvGoToGLD2_Inactive : sc::event< EvGoToGLD2_Inactive > {};
struct EvGoToGLD2_Active : sc::event< EvGoToGLD2_Active > {};
struct EvGoToGLD1_Inactive : sc::event< EvGoToGLD1_Inactive > {};
struct EvGoToGLD1_Active : sc::event< EvGoToGLD1_Active > {};
struct EvGoToLAG1_Inactive : sc::event< EvGoToLAG1_Inactive > {};
struct EvGoToLAG1_Active : sc::event< EvGoToLAG1_Active > {};
struct EvGoToGLP1_Absent : sc::event< EvGoToGLP1_Absent > {};
struct EvGoToGLP1_Active : sc::event< EvGoToGLP1_Active > {};
struct EvGoToGLP1_Inactive : sc::event< EvGoToGLP1_Inactive > {};
struct EvGoToGLP1_Bound : sc::event< EvGoToGLP1_Bound > {};
struct EvGoToGLP1_Unbound : sc::event< EvGoToGLP1_Unbound > {};
struct EvGoToCellCycle_Meiosis : sc::event< EvGoToCellCycle_Meiosis > {};
struct EvGoToCellCycle_Mitosis_M : sc::event< EvGoToCellCycle_Mitosis_M > {};
struct EvGoToCellCycle_Mitosis_S : sc::event< EvGoToCellCycle_Mitosis_S > {};
struct EvGoToCellCycle_Mitosis_G2 : sc::event< EvGoToCellCycle_Mitosis_G2 > {};
struct EvGoToCellCycle_Mitosis_G1 : sc::event< EvGoToCellCycle_Mitosis_G1 > {};
struct EvGoToDifferentiation_Sperm : sc::event< EvGoToDifferentiation_Sperm > {};
struct EvGoToDifferentiation_Oocyte : sc::event< EvGoToDifferentiation_Oocyte > {};
struct EvGoToDifferentiation_EarlyMeiosis : sc::event< EvGoToDifferentiation_EarlyMeiosis > {};
struct EvGoToDifferentiation_LateMeiosis : sc::event< EvGoToDifferentiation_LateMeiosis > {};
struct EvGoToDifferentiation_Precursor : sc::event< EvGoToDifferentiation_Precursor > {};
struct EvGoToOocyteEffector_Inactive : sc::event< EvGoToOocyteEffector_Inactive > {};
struct EvGoToOocyteEffector_Active : sc::event< EvGoToOocyteEffector_Active > {};
struct EvGoToSpermEffector_Inactive : sc::event< EvGoToSpermEffector_Inactive > {};
struct EvGoToSpermEffector_Active : sc::event< EvGoToSpermEffector_Active > {};


//DEFINE THE PARENT STATECHART
struct CellStatechart:  sc::state_machine<CellStatechart,Running>{
  
  //Basic constructor
  CellStatechart();
  
  //Dealing with the associated cell
  CellPtr pCell;
  void SetCell(CellPtr newCell);

  //Deals with copying the StateChart when a cell divides
  boost::shared_ptr<CellStatechart> Copy(boost::shared_ptr<CellStatechart> myNewStatechart);
  std::bitset<100> GetState();
  std::vector<double> GetVariables();
  void SetState(std::bitset<100> state);
  void SetVariables(std::vector<double> variableValues);

  //Declare any state associated variables
  double TimeInPhase;
  double Spermatocyte_divisions;
  double LateMeiosisDelay;
  int Fate;
};


//FIRSTRESPONDER STATE
struct Running:  sc::simple_state<Running,CellStatechart,mpl::list< CellCycle,Differentiation,OocyteEffector,SpermEffector,GLD2,GLD1,LAG1,GLP1 > >{
  typedef sc::custom_reaction< EvCheckCellData > reactions;
  sc::result react( const EvCheckCellData & );
};


//STATES
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD2: sc::state<GLD2,Running::orthogonal<4>,GLD2_Active>{
  GLD2(my_context ctx);

  typedef mpl::list<
  sc::custom_reaction< EvGoToGLD2_Inactive >,
  sc::custom_reaction< EvGoToGLD2_Active > > reactions;

  sc::result react( const EvGoToGLD2_Inactive & ){return transit<GLD2_Inactive>();};
  sc::result react( const EvGoToGLD2_Active & ){return transit<GLD2_Active>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD1: sc::state<GLD1,Running::orthogonal<5>,GLD1_Active>{
  GLD1(my_context ctx);

  typedef mpl::list<
  sc::custom_reaction< EvGoToGLD1_Inactive >,
  sc::custom_reaction< EvGoToGLD1_Active >   > reactions;

  sc::result react( const EvGoToGLD1_Inactive & ){return transit<GLD1_Inactive>();};
  sc::result react( const EvGoToGLD1_Active & ){return transit<GLD1_Active>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct LAG1: sc::state<LAG1,Running::orthogonal<6>,LAG1_Inactive>{
  LAG1(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToLAG1_Inactive >,
   sc::custom_reaction< EvGoToLAG1_Active >   > reactions;

  sc::result react( const EvGoToLAG1_Inactive & ){return transit<LAG1_Inactive>();};
  sc::result react( const EvGoToLAG1_Active & ){return transit<LAG1_Active>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1: sc::state<GLP1,Running::orthogonal<7>,GLP1_Unbound>{
 GLP1(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToGLP1_Absent >,
   sc::custom_reaction< EvGoToGLP1_Active >,
   sc::custom_reaction< EvGoToGLP1_Inactive >,
   sc::custom_reaction< EvGoToGLP1_Bound >,
   sc::custom_reaction< EvGoToGLP1_Unbound >   > reactions;

  sc::result react( const EvGoToGLP1_Absent & ){return transit<GLP1_Absent>();};
  sc::result react( const EvGoToGLP1_Active & ){return transit<GLP1_Active>();};
  sc::result react( const EvGoToGLP1_Inactive & ){return transit<GLP1_Inactive>();};
  sc::result react( const EvGoToGLP1_Bound & ){return transit<GLP1_Bound>();};
  sc::result react( const EvGoToGLP1_Unbound & ){return transit<GLP1_Unbound>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle: sc::state<CellCycle,Running::orthogonal<0>,CellCycle_Mitosis>{
  CellCycle(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToCellCycle_Meiosis >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_M >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_S >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_G2 >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_G1 >   > reactions;

  sc::result react( const EvGoToCellCycle_Meiosis & ){return transit<CellCycle_Meiosis>();};
  sc::result react( const EvGoToCellCycle_Mitosis_M & ){return transit<CellCycle_Mitosis_M>();};
  sc::result react( const EvGoToCellCycle_Mitosis_S & ){return transit<CellCycle_Mitosis_S>();};
  sc::result react( const EvGoToCellCycle_Mitosis_G2 & ){return transit<CellCycle_Mitosis_G2>();};
  sc::result react( const EvGoToCellCycle_Mitosis_G1 & ){return transit<CellCycle_Mitosis_G1>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Mitosis: sc::state<CellCycle_Mitosis,CellCycle,CellCycle_Mitosis_G1>{
  CellCycle_Mitosis(my_context ctx);
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation: sc::state<Differentiation,Running::orthogonal<1>,Differentiation_Precursor>{
  Differentiation(my_context ctx);

  typedef mpl::list< 
   sc::custom_reaction< EvGoToDifferentiation_Sperm >,
   sc::custom_reaction< EvGoToDifferentiation_Oocyte >,
   sc::custom_reaction< EvGoToDifferentiation_EarlyMeiosis >,
   sc::custom_reaction< EvGoToDifferentiation_LateMeiosis >,
   sc::custom_reaction< EvGoToDifferentiation_Precursor >   > reactions;

  sc::result react( const EvGoToDifferentiation_Sperm & ){return transit<Differentiation_Sperm>();};
  sc::result react( const EvGoToDifferentiation_Oocyte & ){return transit<Differentiation_Oocyte>();};
  sc::result react( const EvGoToDifferentiation_EarlyMeiosis & ){return transit<Differentiation_EarlyMeiosis>();};
  sc::result react( const EvGoToDifferentiation_LateMeiosis & ){return transit<Differentiation_LateMeiosis>();};
  sc::result react( const EvGoToDifferentiation_Precursor & ){return transit<Differentiation_Precursor>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct OocyteEffector: sc::state<OocyteEffector,Running::orthogonal<2>,OocyteEffector_Inactive>{
  OocyteEffector(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToOocyteEffector_Inactive >,
   sc::custom_reaction< EvGoToOocyteEffector_Active >   > reactions;

  sc::result react( const EvGoToOocyteEffector_Inactive & ){return transit<OocyteEffector_Inactive>();};
  sc::result react( const EvGoToOocyteEffector_Active & ){return transit<OocyteEffector_Active>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct SpermEffector: sc::state<SpermEffector,Running::orthogonal<3>,SpermEffector_Inactive>{
  SpermEffector(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToSpermEffector_Inactive >,
   sc::custom_reaction< EvGoToSpermEffector_Active >   > reactions;

  sc::result react( const EvGoToSpermEffector_Inactive & ){return transit<SpermEffector_Inactive>();};
  sc::result react( const EvGoToSpermEffector_Active & ){return transit<SpermEffector_Active>();};
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD2_Inactive: sc::state<GLD2_Inactive,GLD2 >{
  GLD2_Inactive(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLD2Update > > reactions;
  sc::result react( const EvGLD2Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD2_Active: sc::state<GLD2_Active,GLD2 >{
  GLD2_Active(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLD2Update > > reactions;
  sc::result react( const EvGLD2Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD1_Inactive: sc::state<GLD1_Inactive,GLD1 >{
  GLD1_Inactive(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLD1Update > > reactions;
  sc::result react( const EvGLD1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD1_Active: sc::state<GLD1_Active,GLD1 >{
  GLD1_Active(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLD1Update > > reactions;
  sc::result react( const EvGLD1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct LAG1_Inactive: sc::state<LAG1_Inactive,LAG1 >{
  LAG1_Inactive(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvLAG1Update > > reactions;
  sc::result react( const EvLAG1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct LAG1_Active: sc::state<LAG1_Active,LAG1 >{
  LAG1_Active(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvLAG1Update > > reactions;
  sc::result react( const EvLAG1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1_Absent: sc::state<GLP1_Absent,GLP1 >{
 GLP1_Absent(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLP1Update > > reactions;
  sc::result react( const EvGLP1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1_Active: sc::state<GLP1_Active,GLP1 >{
 GLP1_Active(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLP1Update > > reactions;
  sc::result react( const EvGLP1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1_Inactive: sc::state<GLP1_Inactive,GLP1 >{
 GLP1_Inactive(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLP1Update > > reactions;
  sc::result react( const EvGLP1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1_Bound: sc::state<GLP1_Bound,GLP1 >{
 GLP1_Bound(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLP1Update > > reactions;
  sc::result react( const EvGLP1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1_Unbound: sc::state<GLP1_Unbound,GLP1 >{
 GLP1_Unbound(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLP1Update > > reactions;
  sc::result react( const EvGLP1Update & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Meiosis: sc::state<CellCycle_Meiosis,CellCycle >{
  double CurrentG1;
  double CurrentS;
  CellCycle_Meiosis(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvCellCycleUpdate > > reactions;
  sc::result react( const EvCellCycleUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Mitosis_M: sc::state<CellCycle_Mitosis_M,CellCycle_Mitosis >{
  double Duration;
  CellCycle_Mitosis_M(my_context ctx);

  typedef sc::custom_reaction< EvCellCycleUpdate> reactions;
  sc::result react( const EvCellCycleUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Mitosis_S: sc::state<CellCycle_Mitosis_S,CellCycle_Mitosis >{
  double Duration;
  CellCycle_Mitosis_S(my_context ctx);

  typedef sc::custom_reaction< EvCellCycleUpdate> reactions;
  sc::result react( const EvCellCycleUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Mitosis_G2: sc::state<CellCycle_Mitosis_G2,CellCycle_Mitosis >{
  double Duration;
  double CompressionThresh;
  CellCycle_Mitosis_G2(my_context ctx);

  typedef sc::custom_reaction< EvCellCycleUpdate> reactions;
  sc::result react( const EvCellCycleUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Mitosis_G1: sc::state<CellCycle_Mitosis_G1,CellCycle_Mitosis >{
  double Duration;
  double CompressionThresh;
  CellCycle_Mitosis_G1(my_context ctx);

  typedef sc::custom_reaction< EvCellCycleUpdate> reactions;
  sc::result react( const EvCellCycleUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation_Sperm: sc::state<Differentiation_Sperm,Differentiation >{
  Differentiation_Sperm(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation_Oocyte: sc::state<Differentiation_Oocyte,Differentiation >{
  Differentiation_Oocyte(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation_EarlyMeiosis: sc::state<Differentiation_EarlyMeiosis,Differentiation >{
  Differentiation_EarlyMeiosis(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation_LateMeiosis: sc::state<Differentiation_LateMeiosis,Differentiation >{
  Differentiation_LateMeiosis(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation_Precursor: sc::state<Differentiation_Precursor,Differentiation >{
  Differentiation_Precursor(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct OocyteEffector_Inactive: sc::state<OocyteEffector_Inactive,OocyteEffector >{
  OocyteEffector_Inactive(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvOocyteEffectorUpdate > > reactions;
  sc::result react( const EvOocyteEffectorUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct OocyteEffector_Active: sc::state<OocyteEffector_Active,OocyteEffector >{
  OocyteEffector_Active(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvOocyteEffectorUpdate > > reactions;
  sc::result react( const EvOocyteEffectorUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct SpermEffector_Inactive: sc::state<SpermEffector_Inactive,SpermEffector >{
  SpermEffector_Inactive(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvSpermEffectorUpdate > > reactions;
  sc::result react( const EvSpermEffectorUpdate & );
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct SpermEffector_Active: sc::state<SpermEffector_Active,SpermEffector >{
  SpermEffector_Active(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvSpermEffectorUpdate > > reactions;
  sc::result react( const EvSpermEffectorUpdate & );
};

#endif