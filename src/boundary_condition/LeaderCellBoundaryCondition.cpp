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

#include "LeaderCellBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GlobalParameterStruct.hpp"


//Constructor
template<unsigned DIM>
LeaderCellBoundaryCondition<DIM>::LeaderCellBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
  boost::shared_ptr< DTCMovementModel<DIM> > pLeaderCellBoundaryModifier,
  double startingRadius)
  : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
  pLeaderCell(pLeaderCellBoundaryModifier),
  TubeRadius(startingRadius) {

  if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
  {
    EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
  }
  if (DIM == 1)
  {
    EXCEPTION("This boundary condition is not implemented in 1D.");
  }
  MaxMovementDistance = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation)->GetAbsoluteMovementThreshold();
}




//Takes cells lying outside the gonad boundary and places them back inside.
template<unsigned DIM>
void LeaderCellBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{

  //Get some relevant information from the leader cell modifier
  std::vector< c_vector<double, DIM> > LeaderCellPointCollection = pLeaderCell->getPathPointCollection();
  std::vector< int > LeaderCellPointTypes = pLeaderCell->getPathPointTypes();
  double Spacing = pLeaderCell->getSpacing();
    
  //!C. ELEGANS SPECIFIC CODE!: Alters the rate of radial gonad growth dependent on the age of worm
  double currentTime = SimulationTime::Instance()->GetTime();
  double timeStep = SimulationTime::Instance()->GetTimeStep();
  if (currentTime < 17.0){
    if (currentTime > 3.5 && currentTime < 7.5){
      TubeRadius = TubeRadius + timeStep*GlobalParameterStruct::Instance()->GetParameter(25);
    }
    else if (currentTime > 7.5 && currentTime < 12.5){
      TubeRadius = TubeRadius + timeStep*GlobalParameterStruct::Instance()->GetParameter(26);
    }
    else if (currentTime > 12.5){
      TubeRadius = TubeRadius + timeStep*GlobalParameterStruct::Instance()->GetParameter(27);
    }
  }

  //Guard against the possibility that the PointCollection may be empty at the start of a simulation
  if ((int)LeaderCellPointCollection.size() > 0){

    //Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
      cell_iter != this->mpCellPopulation->End();
      ++cell_iter)
    {

      //Don't apply the boundary condition to the leader cell itself, assumed to be the first cell in the population
      if (cell_iter != this->mpCellPopulation->Begin()){


        //Read in some properties of this cell
        Node<DIM>* cell_centre_node = this->mpCellPopulation->GetNode(this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter));
        c_vector<double, DIM> cell_location = cell_centre_node->rGetLocation();
        double radius = cell_centre_node->GetRadius();

        //Identify the closest point on leader cell path to this cell, using memoization if possible
        double minDistanceFromPath = DBL_MAX;
        int closestPointIndex = 0;
        //If there is no memoized closest point from the previous timestep, loop over all midline points to find the closest
        if (cell_iter->GetCellData()->GetItem("PreviousClosestPointIndex") == -1){
          double CurrentDistance;
          for (int i = 0; i < (int)LeaderCellPointCollection.size(); i++){
            c_vector<double, DIM> PointOnPathLoc = LeaderCellPointCollection[i];
            CurrentDistance = norm_2(PointOnPathLoc - cell_location);
            if (CurrentDistance < minDistanceFromPath){
              minDistanceFromPath = CurrentDistance;
              closestPointIndex = i;
            }
          }
        //Otherwise search within a few points either side of the previous closest midline point.
        //Decide how far either side to look by how far the cell may have moved since the last timestep
        }else{
          int maxSpheresMoved = (int)(MaxMovementDistance / Spacing) + 5; // safety factor because of stretching
          double previousIndex = cell_iter->GetCellData()->GetItem("PreviousClosestPointIndex");
          double top = fmin(previousIndex + maxSpheresMoved, LeaderCellPointCollection.size() - 1);
          double bottom = fmax(previousIndex - maxSpheresMoved, 0);
          double CurrentDistance;
          for (int i = (int)bottom; i <= (int)top; i++){
            c_vector<double, DIM> PointOnPathLoc = LeaderCellPointCollection[i];
            CurrentDistance = norm_2(PointOnPathLoc - cell_location);
            if (CurrentDistance < minDistanceFromPath){
              minDistanceFromPath = CurrentDistance;
              closestPointIndex = i;
            }
          }
        }


        // The following applies a correction to the cell's position if required... 
        // What to do depends on whether the closest midline point is in one of the 2 endcaps or not.

        //IF CLOSEST POINT IS ONE END OF MIDLINE
        if (closestPointIndex == (int)LeaderCellPointCollection.size() - 1){

          //p1, p2 are the two final points in the leader cell's path. Work out whether the query cell is out past the end of the path
          //i.e. in the endcap, or whether it lies between p1 and p2 i.e. in the final cylindrical portion of the tube 
          c_vector<double, DIM> p1 = LeaderCellPointCollection[closestPointIndex];
          c_vector<double, DIM> p2 = LeaderCellPointCollection[closestPointIndex - 1];
          c_vector<double, DIM> ghostPoint = p1 + (p1 - p2);
          double d1 = norm_2(cell_location - p2);
          double d2 = norm_2(cell_location - ghostPoint);
          if (d2 < d1){
            //Query cell is in the hemispherical endcap - correct if necessary to enforce boundary, or if the cell is in a region with a rachis
            if (minDistanceFromPath > (TubeRadius - radius) || LeaderCellPointTypes[closestPointIndex] == 2 || LeaderCellPointTypes[closestPointIndex] == 1){
              cell_centre_node->rGetModifiableLocation() = p1 + ((TubeRadius - radius) / minDistanceFromPath)*(cell_location - p1);
            }
          }
          else{
            //Otherwise the cell is still in cylindrical part of the tube. 
            //Correct it toward the line between the path end point and previous path point as required
            c_vector<double, DIM> v1 = p1 - p2;
            c_vector<double, DIM> v2 = cell_location - p2;
            double sep = norm_2(p1 - p2);
            double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
            c_vector<double, DIM> path = p2 - p1;
            c_vector<double, DIM> closestPointOnPath = p2 + (dot / (sep*sep))*(p1 - p2);
            double trueMinDistanceFromPath = norm_2(closestPointOnPath - cell_location);
            if (trueMinDistanceFromPath > (TubeRadius - radius) || LeaderCellPointTypes[closestPointIndex] == 2 || LeaderCellPointTypes[closestPointIndex] == 1){
              cell_centre_node->rGetModifiableLocation() = closestPointOnPath + ((TubeRadius - radius) / trueMinDistanceFromPath)
                *(cell_location - closestPointOnPath);
            }
          }

        //ELSE IF CLOSEST POINT IS OTHER END OF MIDLINE (IDENTICAL PROCEDURE)
        }else if (closestPointIndex == 0){

          c_vector<double, DIM> p2 = LeaderCellPointCollection[0];
          c_vector<double, DIM> p1 = LeaderCellPointCollection[1];
          c_vector<double, DIM> ghostPoint = p2 + (p2 - p1);
          double d1 = norm_2(cell_location - p1);
          double d2 = norm_2(cell_location - ghostPoint);
          if (d2 < d1){
            if (minDistanceFromPath > (TubeRadius - radius) || LeaderCellPointTypes[0] == 2){
              cell_centre_node->rGetModifiableLocation() = p2 + ((TubeRadius - radius) / minDistanceFromPath)*(cell_location - p2);
            }
          }
          else{
            c_vector<double, DIM> v1 = p1 - p2;
            c_vector<double, DIM> v2 = cell_location - p2;
            double sep = norm_2(p1 - p2);
            double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
            c_vector<double, DIM> path = p2 - p1;
            c_vector<double, DIM> closestPointOnPath = p2 + (dot / (sep*sep))*(p1 - p2);
            double trueMinDistanceFromPath = norm_2(closestPointOnPath - cell_location);
            if (trueMinDistanceFromPath > (TubeRadius - radius) || LeaderCellPointTypes[closestPointIndex] == 2 || LeaderCellPointTypes[closestPointIndex] == 1){
              cell_centre_node->rGetModifiableLocation() = closestPointOnPath + ((TubeRadius - radius) / trueMinDistanceFromPath)
                *(cell_location - closestPointOnPath);
            }
          }

        //ELSE IF CLOSEST POINT IS ANY OTHER POINT ON PATH 
        }else{

          //Work out which is the second closest point on the path, and therefore which segment of the midline
          //the query cell should be moved toward
          double d1 = norm_2(cell_location - LeaderCellPointCollection[closestPointIndex - 1]);
          double d2 = norm_2(cell_location - LeaderCellPointCollection[closestPointIndex + 1]);
          c_vector<double, DIM> p1;
          c_vector<double, DIM> p2;
          if (d1 < d2){
            p1 = LeaderCellPointCollection[closestPointIndex - 1];
            p2 = LeaderCellPointCollection[closestPointIndex];
          }
          else{
            p1 = LeaderCellPointCollection[closestPointIndex + 1];
            p2 = LeaderCellPointCollection[closestPointIndex];
          }

          //Get query cell perpendicular distance from the closest segment of the path
          c_vector<double, DIM> v1 = p1 - p2;
          c_vector<double, DIM> v2 = cell_location - p2;
          double sep = norm_2(p1 - p2);
          double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
          c_vector<double, DIM> path = p2 - p1;
          c_vector<double, DIM> closestPointOnPath = p2 + (dot / (sep*sep))*(p1 - p2);
          double trueMinDistanceFromPath = norm_2(closestPointOnPath - cell_location);

          //Correct cell position as required
          if (trueMinDistanceFromPath > (TubeRadius - radius) || LeaderCellPointTypes[closestPointIndex] == 2 || LeaderCellPointTypes[closestPointIndex] == 1){
            cell_centre_node->rGetModifiableLocation() = closestPointOnPath + ((TubeRadius - radius) / trueMinDistanceFromPath)
              *(cell_location - closestPointOnPath);
          }

          //Record whether cell is in the proximal arm (useful in some cell cycle models) 
          if(pLeaderCell->getPathPointTypes()[closestPointIndex]==0){
            cell_iter->GetCellData()->SetItem("InProximalArm", 1.0);
          }else{
            cell_iter->GetCellData()->SetItem("InProximalArm", 0.0);
          }

        }

        //Record new closest point on midline path, for memoization purposes
        cell_iter->GetCellData()->SetItem("PreviousClosestPointIndex", (double)closestPointIndex);
        //Also record distance from DTC and max cell radius that fits in the gonad, for use by other classes
        cell_iter->GetCellData()->SetItem("DistanceAwayFromDTC", Spacing*(LeaderCellPointCollection.size() - 1 - closestPointIndex));
        cell_iter->GetCellData()->SetItem("MaxRadius", TubeRadius);

      }

    }
  }
}




//Getter methods for private members
template<unsigned DIM>
  boost::shared_ptr<DTCMovementModel<DIM> > LeaderCellBoundaryCondition<DIM>::GetLeaderCellModifier() const{
  return pLeaderCell;
};
template<unsigned DIM>
  double LeaderCellBoundaryCondition<DIM>::GetTubeRadius() const{
  return TubeRadius;
};



//Boundary condition verification. Required when a BC is used in combination with other, possibly 
//incompatible BCs.
template<unsigned DIM>
bool LeaderCellBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
  //I'm not verifying the BC for speed in this application.
  //There's only one boundary condition in the C. elegans project, so there is no 
  //reason why any incompatible corrections should be being applied. The form of the code would be similar to
  //doing ImposeBoundaryCondition all over again, for future reference
  bool condition_satisfied = true;
  return condition_satisfied;
}



//Parameter output to log file
template<unsigned DIM>
void LeaderCellBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
  // Call method on parent class
  AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class LeaderCellBoundaryCondition<1>;
template class LeaderCellBoundaryCondition<2>;
template class LeaderCellBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LeaderCellBoundaryCondition)
