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

#ifndef GONADARMDATAOUTPUT_HPP_
#define GONADARMDATAOUTPUT_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "OutputFileHandler.hpp"
#include "GlobalParameterStruct.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>

/**
 * A modifier class that every mInterval simulation timesteps, saves to a file some
 * C. elegans germline properties. Obviously only relevant for doing C. elegans germ line simulations
 */
template<unsigned DIM>
class GonadArmDataOutput : public AbstractCellBasedSimulationModifier<DIM,DIM>
{

private:
    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    
        SerializableSingleton<GlobalParameterStruct>* p_params_wrapper = GlobalParameterStruct::Instance()->GetSerializationWrapper();
        archive & p_params_wrapper;
    }

    //Output file stream
    out_stream OutputFile;
    
    //Number of timesteps between data recordings
    int mInterval;

public:


    //Constructor 
    GonadArmDataOutput(int interval);


    //Destructor
    virtual ~GonadArmDataOutput();


    //Getters for private members
    int GetInterval() const;


    /**
     * Overriden UpdateAtEndOfTimeStep method
     * Specifies what to do in the simulation at the end of each timestep, in this case record data
     * if the time is appropriate.
     *
     * @param rCellPopulation reference to the cell population
     */
     void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


    /**
     * Overriden SetupSolve method
     * Specifies what to do in the simulation before the start of the time loop. In this case, 
     * open an output file
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
     void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);



     //Output any parameters associated with this class
     void OutputSimulationModifierParameters(out_stream& rParamsFile);

};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GonadArmDataOutput)

namespace boost
{
    namespace serialization
    {
        /**
        * Serialize information required to construct a GonadArmDataOutput.
        */
        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const GonadArmDataOutput<DIM> * t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            int interval = t->GetInterval();
            ar << interval;
        }

        /**
        * De-serialize constructor parameters and initialise a GonadArmDataOutput.
        */
        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, GonadArmDataOutput<DIM> * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            int interval;
            ar >> interval;

            // Invoke inplace constructor to initialise instance
            ::new(t)GonadArmDataOutput<DIM>(interval);
        }
    }
} // namespace ...

#endif /*GONADARMDATAOUTPUT_HPP_*/