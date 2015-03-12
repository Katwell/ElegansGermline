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

#ifndef GLOBALPARAMETERSTRUCT_HPP_
#define GLOBALPARAMETERSTRUCT_HPP_

#include "SerializableSingleton.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "ChasteSerialization.hpp"

/*
* Creates a single, globally available instance of a parameter vector. Parameter values 
* across the code can then be replaced with lookup calls to this object, useful for 
* setting all parameter values in one place. At the start of the simulation, in the test file,
* this class reads in its parameter values from the config file specified in the command line 
* arguments.
*/

class GlobalParameterStruct : public SerializableSingleton<GlobalParameterStruct>
{
private: 

    /*
     * A pointer to the singleton instance of this class.
     */
    static GlobalParameterStruct* mpInstance;

    /*
    * Results directory name for this simulation
    */
    std::string Directory;

    /*
    *  A vector that stores the actual parameter values
    */
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
    void save(Archive & archive, const unsigned int version) const
    {
        archive & Params;
        archive & Directory;
    }

    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & Params;
        archive & Directory;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()


protected:

    /*
    * Constructor. Protected, should only be called by the method Instance()
    */
    GlobalParameterStruct();


public:

    /*
    * @return a pointer to this parameters object.
    */
    static GlobalParameterStruct* Instance();


    /*
    * @set the directory name and parameter values for this object, using a config file.
    *
    * configFilename = name of file to use, parametersDirectory = the file's location.
    */
    void ConfigureFromFile(std::string configFilename, std::string parametersDirectory);


    /**
    * @destroy parameters object
    */
    static void Destroy();

    /**
    * @return a particular parameter value by number
    */
    double GetParameter(int index);


    /**
    * @get the name of the results directory
    */
    std::string GetDirectory();


    /**
    * @reset a parameter value - occasionaly useful when conducting a parameter sweep
    */
    void ResetParameter(int index, double newValue);


    /**
    * @reset results directory name - occasionaly useful when conducting a parameter sweep
    */
    void ResetDirectoryName(std::string newName);

};

#endif /*GLOBALPARAMETERSTRUCT_HPP_*/