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

#ifndef CELLTRACKINGOUTPUT_HPP_
#define CELLTRACKINGOUTPUT_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "OutputFileHandler.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>


/**
 * A modifier that every (samplingInterval) timesteps saves to a file the position of every (cellIdInterval)-th cell
 */

template<unsigned DIM>
class CellTrackingOutput : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

private:

    //Output file stream
    out_stream OutputFile;
    
    //Number of timesteps between position measurements
    int mSamplingInterval;
    
    //How many cells to track (every idInterval-th cell is followed. E.g. for idInterval=5 every 5th cell is tracked)
    int mCellIdInterval;


public:

    /**
     * Default constructor.
     */
    CellTrackingOutput(int samplingInterval, int cellIdInterval);

    /**
     * Destructor.
     */
    virtual ~CellTrackingOutput();


    /*
    * Getters for the private member variables
    */
    int GetSamplingInterval() const;
    int GetCellIdInterval() const;


    /**
     * Overriden SetupSolve method
     * Specifies what to do before the start of the simulation. In this case, open an output file.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);


    /**
     * Overriden UpdateAtEndOfTimeStep method
     * Specifies what to do at the end of each timestep (i.e. check if data output is required,
     * and if it is, save positions to file).
     *
     * @param rCellPopulation reference to the cell population
     */
     void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);



     //Output any associated parameters
     void OutputSimulationModifierParameters(out_stream& rParamsFile);

};



#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellTrackingOutput)


namespace boost
{
    namespace serialization
    {
        /**
        * Serialize information required to construct a CellTrackingOutput.
        */

        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const CellTrackingOutput<DIM> * t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            int interval = t->GetSamplingInterval();
            ar << interval;

            int idInterval = t->GetCellIdInterval();
            ar << idInterval;
        }

        /**
        * De-serialize constructor parameters and initialise a CellTrackingOutput.
        */
        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, CellTrackingOutput<DIM> * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            int interval;
            ar >> interval;

            int idInterval;
            ar >> idInterval;

            // Invoke inplace constructor to initialise instance
            ::new(t)CellTrackingOutput<DIM>(interval, idInterval);
        }
    }
} // namespace ...

#endif /*CELLTRACKINGOUTPUT_HPP_*/