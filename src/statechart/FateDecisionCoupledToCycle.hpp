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

#ifndef FATEDECISIONCOUPLEDTOCYCLE_HPP_
#define FATEDECISIONCOUPLEDTOCYCLE_HPP_

#include <StatechartCellCycleModel.hpp>
#include <ElegansDevStatechartCellCycleModel.hpp>

#include <boost/statechart/event.hpp>
#include <boost/statechart/state_machine.hpp>
#include <boost/statechart/simple_state.hpp>
#include <boost/statechart/state.hpp>
#include <boost/statechart/custom_reaction.hpp>
#include <boost/statechart/transition.hpp>
#include <boost/mpl/list.hpp>
#include <boost/config.hpp>
namespace sc = boost::statechart;
namespace mpl = boost::mpl;

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <bitset>


//Predeclare a struct for each state
struct Running;
struct GLP1;
struct GLP1_Unbound;
struct GLP1_Bound;
struct GLP1_Absent;
struct LAG1;
struct LAG1_Inactive;
struct LAG1_Active;
struct GLD1;
struct GLD1_Inactive;
struct GLD1_Active;
struct GLD2;
struct GLD2_Inactive;
struct GLD2_Active;
struct CellCycle;
struct CellCycle_Mitosis;
struct CellCycle_Mitosis_G1;
struct CellCycle_Mitosis_S;
struct CellCycle_Mitosis_G2;
struct CellCycle_Mitosis_M;
struct CellCycle_ExitedProlif;
struct CellCycle_ExitedProlif_G1;
struct CellCycle_ExitedProlif_MeioticS;
struct CellCycle_ExitedProlif_Meiosis;
struct Differentiation;
struct Differentiation_Precursor;
struct Differentiation_SpermFated;
struct Differentiation_OocyteFated;
struct Differentiation_Sperm;
struct Differentiation_Oocyte;

//Declare an update event for each orthogonal region
struct EvGLP1Update : sc::event< EvGLP1Update > {};
struct EvLAG1Update : sc::event< EvLAG1Update > {};
struct EvGLD1Update : sc::event< EvGLD1Update > {};
struct EvGLD2Update : sc::event< EvGLD2Update > {};
struct EvCellCycleUpdate : sc::event< EvCellCycleUpdate > {};
struct EvDifferentiationUpdate : sc::event< EvDifferentiationUpdate > {};

//Declare goto events for each leaf state, to allow copying of a statechart's state
struct EvGoToGLP1_Unbound : sc::event< EvGoToGLP1_Unbound > {};
struct EvGoToGLP1_Bound : sc::event< EvGoToGLP1_Bound > {};
struct EvGoToGLP1_Absent : sc::event< EvGoToGLP1_Absent > {};
struct EvGoToLAG1_Inactive : sc::event< EvGoToLAG1_Inactive > {};
struct EvGoToLAG1_Active : sc::event< EvGoToLAG1_Active > {};
struct EvGoToGLD1_Inactive : sc::event< EvGoToGLD1_Inactive > {};
struct EvGoToGLD1_Active : sc::event< EvGoToGLD1_Active > {};
struct EvGoToGLD2_Inactive : sc::event< EvGoToGLD2_Inactive > {};
struct EvGoToGLD2_Active : sc::event< EvGoToGLD2_Active > {};
struct EvGoToCellCycle_ExitedProlif_G1 : sc::event< EvGoToCellCycle_ExitedProlif_G1 > {};
struct EvGoToCellCycle_ExitedProlif_MeioticS : sc::event< EvGoToCellCycle_ExitedProlif_MeioticS > {};
struct EvGoToCellCycle_ExitedProlif_Meiosis : sc::event< EvGoToCellCycle_ExitedProlif_Meiosis > {};
struct EvGoToDifferentiation_Precursor : sc::event< EvGoToDifferentiation_Precursor > {};
struct EvGoToDifferentiation_SpermFated : sc::event< EvGoToDifferentiation_SpermFated > {};
struct EvGoToDifferentiation_OocyteFated : sc::event< EvGoToDifferentiation_OocyteFated > {};
struct EvGoToDifferentiation_Sperm : sc::event< EvGoToDifferentiation_Sperm > {};
struct EvGoToDifferentiation_Oocyte : sc::event< EvGoToDifferentiation_Oocyte > {};


//DEFINE THE PARENT STATECHART
struct FateDecisionCoupledToCycle:  sc::state_machine<FateDecisionCoupledToCycle,Running>{
  
  //Basic constructor
  FateDecisionCoupledToCycle();
  
  //Deals with the associated cell
  CellPtr pCell;
  void SetCell(CellPtr newCell);

  //Deals with copying the Statechart when a cell divides
  boost::shared_ptr<FateDecisionCoupledToCycle> CopyInto(boost::shared_ptr<FateDecisionCoupledToCycle> myNewStatechart);
  std::bitset<MAX_STATE_COUNT> GetState();
  std::vector<double> GetVariables();
  void SetState(std::bitset<MAX_STATE_COUNT> state);
  void SetVariables(std::vector<double> variableValues);

  //Declare any chart associated variables
  double TimeInPhase;
  double SpermatocyteDivisions;
  double SpermDevelopmentDelay;
};


//FIRSTRESPONDER STATE
struct Running:  sc::simple_state<Running,FateDecisionCoupledToCycle,mpl::list< GLP1, LAG1, GLD1, GLD2, CellCycle, Differentiation> >{
  typedef sc::custom_reaction< EvCheckCellData > reactions;
  sc::result react( const EvCheckCellData & );
};


//DEFINE STRUCTURE OF NON-LEAF STATES
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1: sc::state<GLP1,Running::orthogonal<0>,GLP1_Unbound>{
 GLP1(my_context ctx);

  typedef mpl::list<
    sc::custom_reaction< EvGoToGLP1_Unbound >, 
    sc::custom_reaction< EvGoToGLP1_Bound >,
    sc::custom_reaction< EvGoToGLP1_Absent >
  > reactions;

  sc::result react( const EvGoToGLP1_Unbound & ){return transit<GLP1_Unbound>();};
  sc::result react( const EvGoToGLP1_Bound & ){return transit<GLP1_Bound>();};
  sc::result react( const EvGoToGLP1_Absent & ){return transit<GLP1_Absent>();};
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct LAG1: sc::state<LAG1,Running::orthogonal<1>,LAG1_Inactive>{
  LAG1(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToLAG1_Inactive >,
   sc::custom_reaction< EvGoToLAG1_Active >   > reactions;

  sc::result react( const EvGoToLAG1_Inactive & ){return transit<LAG1_Inactive>();};
  sc::result react( const EvGoToLAG1_Active & ){return transit<LAG1_Active>();};
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD1: sc::state<GLD1,Running::orthogonal<2>,GLD1_Active>{
  GLD1(my_context ctx);

  typedef mpl::list<
  sc::custom_reaction< EvGoToGLD1_Inactive >,
  sc::custom_reaction< EvGoToGLD1_Active >   > reactions;

  sc::result react( const EvGoToGLD1_Inactive & ){return transit<GLD1_Inactive>();};
  sc::result react( const EvGoToGLD1_Active & ){return transit<GLD1_Active>();};
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLD2: sc::state<GLD2,Running::orthogonal<3>,GLD2_Active>{
  GLD2(my_context ctx);

  typedef mpl::list<
  sc::custom_reaction< EvGoToGLD2_Inactive >,
  sc::custom_reaction< EvGoToGLD2_Active > > reactions;

  sc::result react( const EvGoToGLD2_Inactive & ){return transit<GLD2_Inactive>();};
  sc::result react( const EvGoToGLD2_Active & ){return transit<GLD2_Active>();};
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle: sc::state<CellCycle,Running::orthogonal<4>,CellCycle_Mitosis>{
  CellCycle(my_context ctx);

  typedef mpl::list<
   sc::custom_reaction< EvGoToCellCycle_Mitosis_G1 >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_S >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_G2 >,
   sc::custom_reaction< EvGoToCellCycle_Mitosis_M >,
   sc::custom_reaction< EvGoToCellCycle_ExitedProlif_G1>,
   sc::custom_reaction< EvGoToCellCycle_ExitedProlif_MeioticS>,
   sc::custom_reaction< EvGoToCellCycle_ExitedProlif_Meiosis>
  > reactions;

  sc::result react( const EvGoToCellCycle_Mitosis_G1 & ){return transit<CellCycle_Mitosis_G1>();};
  sc::result react( const EvGoToCellCycle_Mitosis_S & ){return transit<CellCycle_Mitosis_S>();};
  sc::result react( const EvGoToCellCycle_Mitosis_G2 & ){return transit<CellCycle_Mitosis_G2>();};
  sc::result react( const EvGoToCellCycle_Mitosis_M & ){return transit<CellCycle_Mitosis_M>();};
  sc::result react( const EvGoToCellCycle_ExitedProlif_G1 & ){return transit<CellCycle_ExitedProlif_G1>();};
  sc::result react( const EvGoToCellCycle_ExitedProlif_MeioticS & ){return transit<CellCycle_ExitedProlif_MeioticS>();};
  sc::result react( const EvGoToCellCycle_ExitedProlif_Meiosis & ){return transit<CellCycle_ExitedProlif_Meiosis>();};
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_Mitosis: sc::state<CellCycle_Mitosis,CellCycle,CellCycle_Mitosis_G1>{
  CellCycle_Mitosis(my_context ctx);
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_ExitedProlif: sc::state<CellCycle_ExitedProlif,CellCycle,CellCycle_ExitedProlif_G1>{
  CellCycle_ExitedProlif(my_context ctx);
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation: sc::state<Differentiation,Running::orthogonal<5>,Differentiation_Precursor>{
  Differentiation(my_context ctx);

  typedef mpl::list< 
   sc::custom_reaction< EvGoToDifferentiation_Precursor >,
   sc::custom_reaction< EvGoToDifferentiation_SpermFated >,
   sc::custom_reaction< EvGoToDifferentiation_OocyteFated >,
   sc::custom_reaction< EvGoToDifferentiation_Sperm >,
   sc::custom_reaction< EvGoToDifferentiation_Oocyte >
  > reactions;

  sc::result react( const EvGoToDifferentiation_Precursor & ){return transit<Differentiation_Precursor>();};  
  sc::result react( const EvGoToDifferentiation_SpermFated & ){return transit<Differentiation_SpermFated>();};
  sc::result react( const EvGoToDifferentiation_OocyteFated & ){return transit<Differentiation_OocyteFated>();};
  sc::result react( const EvGoToDifferentiation_Sperm & ){return transit<Differentiation_Sperm>();};
  sc::result react( const EvGoToDifferentiation_Oocyte & ){return transit<Differentiation_Oocyte>();};
};


//DEFINE STRUCTURE OF LEAF STATES
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct GLP1_Unbound: sc::state<GLP1_Unbound,GLP1 >{
  GLP1_Unbound(my_context ctx);

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
struct GLP1_Absent: sc::state<GLP1_Absent,GLP1 >{
  GLP1_Absent(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvGLP1Update > > reactions;
  sc::result react( const EvGLP1Update & );
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
struct CellCycle_Mitosis_G1: sc::state<CellCycle_Mitosis_G1,CellCycle_Mitosis >{
  double Duration;
  double CompressionThresh;
  CellCycle_Mitosis_G1(my_context ctx);

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
struct CellCycle_Mitosis_M: sc::state<CellCycle_Mitosis_M,CellCycle_Mitosis >{
  double Duration;
  CellCycle_Mitosis_M(my_context ctx);

  typedef sc::custom_reaction< EvCellCycleUpdate> reactions;
  sc::result react( const EvCellCycleUpdate & );
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_ExitedProlif_G1: sc::state<CellCycle_ExitedProlif_G1,CellCycle_ExitedProlif >{
  double Duration;
  CellCycle_ExitedProlif_G1(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvCellCycleUpdate > > reactions;
  sc::result react( const EvCellCycleUpdate & );
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_ExitedProlif_MeioticS: sc::state<CellCycle_ExitedProlif_MeioticS,CellCycle_ExitedProlif >{
  double Duration;
  CellCycle_ExitedProlif_MeioticS(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvCellCycleUpdate > > reactions;
  sc::result react( const EvCellCycleUpdate & );
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct CellCycle_ExitedProlif_Meiosis: sc::state<CellCycle_ExitedProlif_Meiosis, CellCycle_ExitedProlif >{
  CellCycle_ExitedProlif_Meiosis(my_context ctx);

  typedef mpl::list< sc::custom_reaction< EvCellCycleUpdate > > reactions;
  sc::result react( const EvCellCycleUpdate & );
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
struct Differentiation_SpermFated: sc::state<Differentiation_SpermFated,Differentiation >{
  Differentiation_SpermFated(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
};
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
struct Differentiation_OocyteFated: sc::state<Differentiation_OocyteFated,Differentiation >{
  Differentiation_OocyteFated(my_context ctx);
  typedef mpl::list< sc::custom_reaction< EvDifferentiationUpdate > > reactions;
  sc::result react( const EvDifferentiationUpdate & );
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


// Export StatechartCellCycleModel and ElegansDevStatechartCellCycleModel with this class
// as template parameter
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(StatechartCellCycleModel, FateDecisionCoupledToCycle)
EXPORT_TEMPLATE_CLASS1(ElegansDevStatechartCellCycleModel, FateDecisionCoupledToCycle)

#endif