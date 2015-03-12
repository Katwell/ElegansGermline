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

#ifndef REPULSIONFORCESIZECORRECTED_HPP_
#define REPULSIONFORCESIZECORRECTED_HPP_

#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
* A two-body repulsion force law, designed for use in node-based simulations.
* Creates a repulsion between pairs of cells that overlap. 
*
* This is a modification of the GeneralisedLinearSpringForce class. The main change is that the force 
* experienced by a cell is scaled linearly downward proportional to its radius, to account for the 
* different drag forces experienced by cells of different sizes.
*/

template<unsigned DIM>
class RepulsionForceSizeCorrected : public GeneralisedLinearSpringForce<DIM>
{
private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
    }

public :

    /**
     * Constructor.
     */
    RepulsionForceSizeCorrected();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the CellPopulation
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Outputs force Parameters to file
     * Nothing to add here over and above the parameters that exist in the parent class.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionForceSizeCorrected)

#endif /*REPULSIONFORCESIZECORRECTED_HPP_*/
