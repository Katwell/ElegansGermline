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

#ifndef OOCYTEFATEDCELLAPOPTOSIS_HPP_
#define OOCYTEFATEDCELLAPOPTOSIS_HPP_

#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A cell killer that loops over all germ cells and selects oocyte-fated cells < 250 microns from the DTC
 * (i.e. oocyte-fated cells NOT yet in the proximal arm). These cells are then killed with a probability per
 * hour of p. 
 */

template<unsigned DIM>
class OocyteFatedCellApoptosis : public AbstractCellKiller<DIM>
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

        // Make sure the random number generator is also archived
        SerializableSingleton<RandomNumberGenerator>* p_rng_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_rng_wrapper;
    }


    /*
    * Probability of death for an oocyte-fated cell spending 1 hour outside the proliferative zone
    */
    double mHourlyProbabilityOfDeath;

public:


    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param HourlyProbabilityOfDeath p(death) for each hour spent in target area
     */
    OocyteFatedCellApoptosis(AbstractCellPopulation<DIM>* pCellPopulation, double HourlyProbabilityOfDeath);
    

    /*
    * Destructor
    */
    ~OocyteFatedCellApoptosis();


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
     * Loops over all cells and tests whether they should be killed
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OocyteFatedCellApoptosis)


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OocyteFatedCellApoptosis instance.
 */

template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OocyteFatedCellApoptosis<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    double p = t->GetHourlyProbabilityOfDeath();
    ar << p;
}

/**
 * De-serialize constructor parameters and initialise an OocyteFatedCellApoptosis.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OocyteFatedCellApoptosis<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    
    double p;
    ar >> p;

    // Invoke inplace constructor to initialise instance
    ::new(t)OocyteFatedCellApoptosis<DIM>(p_cell_population, p);
}
}
} // namespace ...

#endif /*OOCYTEFATEDCELLAPOPTOSIS_HPP_*/
