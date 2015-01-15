/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef GLOBALPARAMETERSTRUCT_HPP_
#define GLOBALPARAMETERSTRUCT_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include "SerializableSingleton.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


//Creates a single, globally available instance of a parameter struct. Allows parameter values 
//across your code to be replaced with lookup calls to this object, useful for sweeps
//and for setting all parameter values in one place. Reads in the values from a config file specified
//the first time an instance is requested.

class GlobalParameterStruct : public SerializableSingleton<GlobalParameterStruct>
{
public:

    /**
    * @set a pointer to the parameter struct. Requires the name of a config file.
    */
    static GlobalParameterStruct* Instance(std::string configFilename, std::string parametersDirectory);

    /**
    * @return a pointer to the parameters object. Used to retreive instance once config file has been set.
    */
    static GlobalParameterStruct* Instance();

    /**
    * @destroy parameters object
    */
    static void Destroy();

    /**
    * @return a particular parameter value
    */
    double GetParameter(int index);

    /**
    * @get name of the results output directory
    */
    std::string GetDirectory();

    /**
    * @reset a parameter value - sometimes useful when conducting a parameter sweep
    */
    void ResetParameter(int index, double newValue);

    /**
    * @reset results directory name - sometimes useful when conducting a parameter sweep
    */
    void ResetDirectoryName(std::string newName);


protected:

    GlobalParameterStruct(std::string filename, std::string parametersDirectory);

private: 

    /**
     * A pointer to the singleton instance of this class.
     */
    static GlobalParameterStruct* mpInstance;

    /**
    * Results directory name
    **/
    std::string Directory;

    /**
    *  A vector that stores actual parameter values
    **/
    std::vector<double> Params;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialization must be done with care.
     * Do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & Directory;
        archive & Params;
    }
};

#endif /*GLOBALPARAMETERSTRUCT_HPP_*/