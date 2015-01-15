/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef LEADERCELLBOUNDARYCONDITION_HPP_
#define LEADERCELLBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "DTCMovementModel.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
* A boundary condition, which takes in a pointer to a leader cell modifier. This class ensures
* that cells stay within a certain distance of that leader cell's path.
*/

template<unsigned DIM>
class LeaderCellBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
    * Serialize the object.
    *
    * @param archive the archive
    * @param version the current version of this class
    */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
    }


    //Pointer to the modifier that generates the leader cell path. Currently DTCMovementModel is the only available choice 
    boost::shared_ptr< DTCMovementModel<DIM> > pLeaderCell;

    //Current maximum distance that cells can be from the leader cell path
    double TubeRadius;

    //Max possible distance a cell can move in a timestep
    double MaxMovementDistance;

public:


    /**
    * Constructor.
    *
    * @param pCellPopulation pointer to the cell population
    * @param pLeaderCellBoundaryModifier pointer to the modifier that handles leader cell movement
    * @param startingRadius initial tube radius
    */
    LeaderCellBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
        boost::shared_ptr<DTCMovementModel<DIM> > pLeaderCellBoundaryModifier,
        double startingRadius);


    /**
    * Overridden ImposeBoundaryCondition() method.
    * Apply the cell population boundary condition.
    *
    * @param rOldLocations the node locations before any boundary conditions are applied
    */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);


    /**
    * Overridden VerifyBoundaryCondition() method.
    * Verify the boundary conditions have been applied.
    * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
    *
    * @return whether the boundary conditions are satisfied.
    */
    bool VerifyBoundaryCondition();


    //Getters for private members
    const boost::shared_ptr< DTCMovementModel<DIM> > GetLeaderCellModifier() const;
    const double GetTubeRadius() const;


    /**
    * Overridden OutputCellPopulationBoundaryConditionParameters() method.
    * Output cell population boundary condition parameters to file.
    *
    * @param rParamsFile the file stream to which the parameters are output
    */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LeaderCellBoundaryCondition)

namespace boost
{
    namespace serialization
    {
        /**
        * Serialize information required to construct a LeaderCellBoundaryCondition.
        */
        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const LeaderCellBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
            const boost::shared_ptr<DTCMovementModel<DIM> > pLeaderCell = t->GetLeaderCellModifier();
            ar << pLeaderCell;
            const double tubeRadius = t->GetTubeRadius();
            ar << tubeRadius;
        }


        /**
        * De-serialize constructor parameters and initialize a LeaderCellBoundaryCondition.
        */
        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, LeaderCellBoundaryCondition<DIM>* t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            AbstractCellPopulation<DIM>* p_cell_population;
            ar >> p_cell_population;
            boost::shared_ptr< DTCMovementModel<DIM> > pLeaderCell;
            ar >> pLeaderCell;
            double startingRadius;
            ar >> startingRadius;

            // Invoke inplace constructor to initialise instance
            ::new(t)LeaderCellBoundaryCondition<DIM>(p_cell_population, pLeaderCell, startingRadius);
        }
    }
} // namespace ...

#endif /*LEADERCELLBOUNDARYCONDITION_HPP_*/