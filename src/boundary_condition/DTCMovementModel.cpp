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

#include "DTCMovementModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GlobalParameterStruct.hpp"

#include <cmath>
#include <vector>

template<unsigned DIM>
DTCMovementModel<DIM>::DTCMovementModel(std::vector< c_vector<double, DIM> > startingLocations,
                                              std::vector< int > startingTypes,
                                              c_vector<double,DIM> currentLocation,
                                              double spacing,
                                              bool Unc5Initial,
                                              bool Vab3Initial,
                                              double TimeSinceUpdate): 
    AbstractCellBasedSimulationModifier<DIM>(),
    PathPointCollection(startingLocations),
    PathPointTypes(startingTypes),
    CurrentLocation(currentLocation),
    Spacing(spacing),
    Unc5(Unc5Initial),
    Vab3(Vab3Initial),
    TimeSinceLastUpdate(TimeSinceUpdate)
{   
    //Get radius of the DTC turn
    WormBodyRadius = GlobalParameterStruct::Instance()->GetParameter(7);
    TurnComplete = false;
    //Get the rate of stretching parameter
    StretchingRate = GlobalParameterStruct::Instance()->GetParameter(1);
}


template<unsigned DIM>
DTCMovementModel<DIM>::~DTCMovementModel()
{
}


template<unsigned DIM>
void DTCMovementModel<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{

    //This block updates the states of the genes Vab3 and Unc5 over time; manually for now. 
    //After a delay, Unc5 switches on and the DTC turns onto the dorsal surface
    //After a longer delay, Vab3 switches on and the DTC halts
    if (SimulationTime::Instance()->GetTime() + 18.5>GlobalParameterStruct::Instance()->GetParameter(6)){
        Vab3 = true;
    }
    if (SimulationTime::Instance()->GetTime() + 18.5>GlobalParameterStruct::Instance()->GetParameter(8)){
        Unc5 = true;
    }

    //This block checks whether or not cells are present close enough behind the DTC to push it 
    bool DTCBeingPushed = false;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End(); ++cell_iter){
        if (cell_iter->GetCellData()->GetItem("IsDTC") == 0.0 &&
            cell_iter->GetCellData()->GetItem("DistanceAwayFromDTC") < 5){
            DTCBeingPushed = true;
            break;
        }
    }

    //Set the current speed of the DTC, based on age of the worm and whether the DTC can move right now
    double currentSpeed;
    double currentTime = SimulationTime::Instance()->GetTime();
    if (currentTime < 3.5){
        currentSpeed = GlobalParameterStruct::Instance()->GetParameter(2);
    }
    else if (currentTime > 3.5 && currentTime < 7.5){
        currentSpeed = GlobalParameterStruct::Instance()->GetParameter(3);
    }
    else if (currentTime > 7.5 && currentTime < 12.5){
        currentSpeed = GlobalParameterStruct::Instance()->GetParameter(4);
    }
    else if (currentTime > 12.5 && currentTime < 17.0){
        currentSpeed = GlobalParameterStruct::Instance()->GetParameter(5);
    }
    if (Vab3 == true || DTCBeingPushed == false){
        currentSpeed = 0;
    }

    //Work out how fast the DTC angle should be changing in a turn, given DTC migration speed  
    double timestep = SimulationTime::Instance()->GetTimeStep();
    double thetaIncrement = (currentSpeed*timestep) / WormBodyRadius;
    

    //Update the DTC location
    int pointType;
    if (currentSpeed > 0){

        if (Unc5 == true){

            //If the turn onto dorsal surface has begun, work out current angle of DTC around worm body
            //To save time, if the PathPointTypes vector indicates that this has aleady been achieved,
            //skip the arcsin calculation and just assume a theta of 0.
            double theta;
            if (TurnComplete==true){
                theta = 0;
            }else{
                double pv = asin(CurrentLocation[0] / WormBodyRadius);
                if (CurrentLocation[1] <= 0 & CurrentLocation[0] >= 0){
                    theta = M_PI - pv;
                }
                else if (CurrentLocation[1] <= 0 & CurrentLocation[0] < 0){
                    theta = -M_PI - pv;
                }
                else{
                    theta = pv;
                }
                if (theta >= -0.1 && theta <= 0.1){
                    TurnComplete = true;
                }
            }

            //If the middle of the dorsal surface has been reached, let the DTC migrate back toward the middle of the worm
            if (theta >= -0.1 && theta <= 0.1){
                CurrentLocation[2] -= currentSpeed*timestep;
                pointType = 2;

            //If the DTC is not in the centre of the dorsal surface, let it migrate toward it, making a turn
            }else if (theta >= 0.1){
                double newTheta = theta - thetaIncrement;
                CurrentLocation[0] = WormBodyRadius * sin(newTheta);
                CurrentLocation[1] = WormBodyRadius * cos(newTheta);
                pointType = 1;
            }else if (theta <= -0.1){
                double newTheta = theta + thetaIncrement;
                CurrentLocation[0] = WormBodyRadius * sin(newTheta);
                CurrentLocation[1] = WormBodyRadius * cos(newTheta);
                pointType = 1;
            }

        }else{
            //If the turn hasn't started, DTC moves away from the worm's centre along the ventral side
            CurrentLocation[2] += currentSpeed*timestep;
            pointType = 0;
        }
        //Update DTC location
        rCellPopulation.GetNode(0)->rGetModifiableLocation() = CurrentLocation;
    }

    //If the DTC has moved further than the spacing interval, put a new point in the PointCollection
    double distanceMoved = norm_2(CurrentLocation - PathPointCollection.back());
    if (distanceMoved > Spacing){
        PathPointCollection.push_back(CurrentLocation);
        PathPointTypes.push_back(pointType);
        if (distanceMoved - Spacing > 0.05){
            printf("Leader cell is moving quickly, and overshooting the spacing separation by > 0.05. Consider reducting the timestep");
        }
    }


    //Deals with gonad growth by stretching, during the L4 when the turn has completed
    if (currentTime > 12.5 && currentTime < 17.0 && PathPointTypes.back() == 2){

        //If enough stretching has happened since the last update, add another point into the collection:
        if (TimeSinceLastUpdate > Spacing / StretchingRate){

            //Work out which points correspond to two ends of the turn 
            int oneStraightEndIndex;
            int otherStraightEndIndex;
            assert(PathPointTypes.size() == PathPointCollection.size());
            for (int i = 1; i<PathPointTypes.size(); i++){
                if (PathPointTypes[i] == 1 && PathPointTypes[i - 1] == 0){
                    oneStraightEndIndex = i - 1;
                }
                if (PathPointTypes[i] == 2 && PathPointTypes[i - 1] == 1){
                    otherStraightEndIndex = i - 1;
                }
                if (PathPointTypes[i] == 1){
                    PathPointCollection[i][2] += Spacing; //Translate the loop points by updating their positions
                }
            }

            //Work out the locations of the two new points
            c_vector<double, DIM> newPoint1 = PathPointCollection[oneStraightEndIndex];
            newPoint1[2] += Spacing;
            c_vector<double, DIM> newPoint2 = PathPointCollection[otherStraightEndIndex + 1];
            newPoint2[2] += Spacing;

            //Make space for new points in the collection vectors
            int currentSize = PathPointTypes.size();
            PathPointCollection.resize(currentSize + 2);
            PathPointTypes.resize(currentSize + 2);

            //Update the point collection vectors, adding 2 new points. All of this will be very slow, but it occurs only 
            //for a short part of the simulation
            for (int i = currentSize + 1; i >= 0; i--){
                if (i>otherStraightEndIndex + 2){
                    PathPointCollection[i] = PathPointCollection[i - 2];
                    PathPointTypes[i] = PathPointTypes[i - 2];
                }
                else if (i == otherStraightEndIndex + 2){
                    PathPointCollection[i] = newPoint2;
                    PathPointTypes[i] = 2;
                }
                else if (i< otherStraightEndIndex + 2 && i>oneStraightEndIndex + 1){
                    PathPointCollection[i] = PathPointCollection[i - 1];
                    PathPointTypes[i] = PathPointTypes[i - 1];
                }
                else if (i == oneStraightEndIndex + 1){
                    PathPointCollection[i] = newPoint1;
                    PathPointTypes[i] = 0;
                }
            }

            //Keeps track of how long it's been since a stretching correction was last applied
            TimeSinceLastUpdate = 0;
        }
        else{
            TimeSinceLastUpdate += timestep;
        }
    }
}


//Getters for private members
template<unsigned DIM>
const std::vector< c_vector<double, DIM> > DTCMovementModel<DIM>::getPathPointCollection() const{
    return PathPointCollection;
}

template<unsigned DIM>
const std::vector< int > DTCMovementModel<DIM>::getPathPointTypes() const{
    return PathPointTypes;
}

template<unsigned DIM>
const c_vector< double,DIM > DTCMovementModel<DIM>::getCurrentLocation() const{
    return CurrentLocation;
}


template<unsigned DIM>
const double DTCMovementModel<DIM>::getSpacing() const{
    return Spacing;
}

template<unsigned DIM>
const double DTCMovementModel<DIM>::getTimeSinceLastUpdate() const{
    return TimeSinceLastUpdate;
}

template<unsigned DIM>
const bool DTCMovementModel<DIM>::getUnc5() const{
    return Unc5;
}

template<unsigned DIM>
const bool DTCMovementModel<DIM>::getVab3() const{
    return Vab3;
}

template<unsigned DIM>
void DTCMovementModel<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    // Nothing to do here for now
}


template<unsigned DIM>
void DTCMovementModel<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Unc5Value>" << getUnc5() << "</Unc5Value>\n";
    *rParamsFile << "\t\t\t<Vab3Value>" << getVab3() << "</Vab3Value>\n";
    *rParamsFile << "\t\t\t<PointOnPathSpacing>" << getSpacing() << "</PointOnPathSpacing>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DTCMovementModel<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DTCMovementModel)