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

#ifndef FERTILISATION_HPP_
#define FERTILISATION_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/**
* A cell killer that loops over all cells in the population and looks for those in the oocyte state. It decides
* whether they can currently be fertilised, by seeing whether the oocyte is in the spermatheca and
* whether a sperm cell is also present in the spermatheca. If so, oocyte and sperm are removed.
*
* This is implemented as a cell killer because both cells involved in fertilisation are removed
* from the simulation.
*/

template<unsigned DIM>
class Fertilisation : public AbstractCellKiller<DIM>
{
private:

    /**
    * Length of the spermatheca, used to determine when a sperm and oocyte are in the position to be fertilised
    */
    double mSpermathecaLength;


    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
    * Archive the object.
    *
    * @param archive the archive
    * @param version the current version of this class
    */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
    }

public:

    /**
    * Constructor.
    *
    * @param pCellPopulation pointer to the cell population
    * @param spermathecaLength length of the spermatheca
    */
    Fertilisation(AbstractCellPopulation<DIM>* pCellPopulation, double spermathecaLength);


    //Destructor
    ~Fertilisation();


    /**
    * @return mSpermathecaLength.
    */
    double GetSpermathecaLength() const;


    /**
    * Once a sperm  and oocyte involved in fertilisation have been labelled for death, this
    * function sends the kill signal.
    *
    * @param pCell the cell to be removed
    */
    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);


    /**
    * Overridden method to test whether a given oocyte should be fertilised.
    * Implemented as a cell killer because a sperm cell and an oocyte are removed when fertilisation
    * occurs.
    *
    * @param pCell the oocyte to test for fertilisation
    */
    void CheckAndLabelCellsForApoptosisOrDeath();


    /**
    * Overridden OutputCellKillerParameters() method.
    *
    * @param rParamsFile the file stream to which the parameters are output
    */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};



#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Fertilisation)

namespace boost
{
    namespace serialization
    {
        /**
        * Serialize information required to construct a Fertilisation cell killer.
        */

        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const Fertilisation<DIM> * t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
            double dist = t->GetSpermathecaLength();
            ar << dist;
        }

        /**
        * De-serialize constructor parameters and initialise a Fertilisation cell killer.
        */
        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, Fertilisation<DIM> * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            AbstractCellPopulation<DIM>* p_cell_population;
            ar >> p_cell_population;
            double dist;
            ar >> dist;

            // Invoke inplace constructor to initialise instance
            ::new(t)Fertilisation<DIM>(p_cell_population, dist);
        }
    }
} // namespace ...

#endif /*FERTILISATION_HPP_*/
