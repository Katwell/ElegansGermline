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

#ifndef DTCMOVEMENTMODEL_HPP_
#define DTCMOVEMENTMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
* An example of a simulation modifier that updates the position of a leader cell at each time step. This class is: 
* 1) Designed to work with a 3D node based cell population only.
* 2) The leader cell is taken to be the first node in the population, with ID=0.
* This modifier can be used together with a LeaderCellBoundaryCondition, to produce a tubular boundary 
* condition enforced allong the path of the moving cell.
*
* This particular modifier is C. elegans germline specific, reflecting the motion of the Distal Tip Cell during 
* gonad development. To model other leader cells, it is recommended to create a new class inheriting from 
* AbstractCellBasedSimulationModifier and use this file as a guide to produce the desired pattern of movement.
*/

template<unsigned DIM>
class DTCMovementModel : public AbstractCellBasedSimulationModifier<DIM, DIM>
{

private:

    /* Needed for serialization. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
    }

    /* Some C. elegans specific variables: */
    bool   Unc5;                    //The gene that is switched on when the DTC starts to turn
    bool   Vab3;                    //The gene that is switched on when the DTC halts
    bool   TurnComplete;            //Set to true when the DTC has reached the middle of the dorsal surface
    double TimeSinceLastUpdate;     //Used in proximal arm stretching to determine when new midline points should be inserted.
    double WormBodyRadius;          //The radius of the gonad turn
    double StretchingRate;          //Rate of growth by stretching during the L4  

    //More generally applicable leader cell variables
    std::vector< c_vector<double, DIM> > PathPointCollection;  //Stores equally spaced points on the leader cell's path
    std::vector< int > PathPointTypes;                         //A flag associated with each point, for additional info
    c_vector < double, DIM > CurrentLocation;                  //Current leader cell position           
    double Spacing;                                            //Separation of points on path
     

public:

    /**
    * Constructor. Takes in an initial collection of path points, and their types (flags that give associated info,
    * such as whether the points form part of a straight or a turn).
    * @param UncInitial initial value of Unc5, determines whether DTC has started turn yet
    * @param VabInitial initial value of Vab3, switches on as the DTC halts
    * @param StartingPointLocations initial collection of points on leader cell path
    * @param StartingPointTypes initial collection of flags for each of those points
    * @param currentLocation = current DTC position
    * @param spacing distance between consecutive points on leader cell path
    */
    DTCMovementModel(bool UncInitial, bool VabInitial,
        double TimeSinceUpdate,
        std::vector< c_vector<double, DIM> > StartingPointLocations,
        std::vector< int > StartingPointTypes, 
        c_vector<double, DIM> currentLocation,
        double spacing);


    /**
    * Destructor.
    */
    virtual ~DTCMovementModel();


    /**
    * Overriden UpdateAtEndOfTimeStep method
    *
    * Specifies what to do in the simulation at the end of each timestep.
    *
    * @param rCellPopulation reference to the cell population
    */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);


    /**
    * Overriden SetupSolve method
    *
    * Specifies what to do in the simulation before the start of the time loop.
    *
    * @param rCellPopulation reference to the cell population
    * @param outputDirectory the output directory, relative to where Chaste output is stored
    */
    void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);



    //Getters for private member variables
    std::vector< c_vector<double, DIM> > getPathPointCollection() const;  
    std::vector< int > getPathPointTypes() const;
    c_vector< double, DIM > getCurrentLocation() const;
    double getSpacing() const;
    double getTimeSinceLastUpdate() const;
    bool getUnc5() const;
    bool getVab3() const;


    /**
    * Overriden OutputSimulationModifierParameters method
    *
    * Outputs parameters relating to this modifier to the simulation parameters file.
    *
    * @param rParamsFile file that stores the simualtion parameters
    */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);


};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DTCMovementModel)

namespace boost
{
    namespace serialization
    {
        /**
        * Serialize information required to construct a DTCMovementModel.
        */
        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const DTCMovementModel<DIM>* t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct an instance
            std::vector< c_vector<double, DIM> > pathPointCollection = t->getPathPointCollection();
            ar << pathPointCollection;
            std::vector< int > pathPointTypes = t->getPathPointTypes();
            ar << pathPointTypes;
            c_vector<double, DIM> currentLocation = t->getCurrentLocation();
            ar << currentLocation;
            // Archive other member variables
            double sp = t->getSpacing();
            ar << sp;
            bool unc5 = t->getUnc5();
            ar << unc5;
            bool vab3 = t->getVab3();
            ar << vab3;
            double timeSinceUpdate = t->getTimeSinceLastUpdate();
            ar << timeSinceUpdate;
        }

        /**
        * De-serialize constructor parameters and initialize a DTCMovementModel.
        */
        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, DTCMovementModel<DIM>* t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector< c_vector<double, DIM> > pathPointCollection;
            ar >> pathPointCollection;
            std::vector< int > pathPointTypes;
            ar >> pathPointTypes;
            c_vector<double, DIM> currentLocation;
            ar >> currentLocation;
            // Retrieve other member variables
            double sp;
            ar >> sp;
            bool unc5;
            ar >> unc5;
            bool vab3;
            ar >> vab3;
            double timeSinceUpdate;
            ar >> timeSinceUpdate;

            // Invoke inplace constructor to initialise instance
            ::new(t)DTCMovementModel<DIM>(unc5, vab3, timeSinceUpdate, pathPointCollection, pathPointTypes, currentLocation, sp);
        }
    }
} // namespace ...

#endif /*DTCMOVEMENTMODEL_HPP_*/