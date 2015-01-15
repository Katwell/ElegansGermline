
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

#ifndef RandomCellKillerByType_HPP_
#define RandomCellKillerByType_HPP_

#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"
#include "GlobalParameterStruct.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

/**
 * A cell killer that loops over all germ cells and for those fated to become oocytes that are NOT in the proximal arm,
 * kills with a certain per hour probability p. 
 */

template<unsigned DIM>
class RandomCellKillerByType : public AbstractCellKiller<DIM>
{
private:

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

    double HourlyProbabilityOfDeath;

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param HourlyProbabilityOfDeath
     */
    RandomCellKillerByType(AbstractCellPopulation<DIM>* pCellPopulation, double HourlyProbabilityOfDeath);
    ~RandomCellKillerByType();


    /**
     * @return HourlyProbabilityOfDeath.
     */
    double GetHourlyProbabilityOfDeath() const;


    /**
     * This function sends the kill signal to cells labelled for death.
     * 
     * @param pCell the cell to kill
     */
    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);


    /**
     * Tests whether each cell should be killed
     *
     * @param pCell the cell to test for RandomCellKillerByType
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomCellKillerByType)


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKillerByType.
 */

template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const RandomCellKillerByType<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    double p = t->GetHourlyProbabilityOfDeath();
    ar << p;
}

/**
 * De-serialize constructor parameters and initialise a RandomCellKillerByType.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, RandomCellKillerByType<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    double p;
    ar >> p;

    // Invoke inplace constructor to initialise instance
    ::new(t)RandomCellKillerByType<DIM>(p_cell_population, p);
}
}
} // namespace ...

#endif /*RandomCellKillerByType_HPP_*/
