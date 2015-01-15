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

#include "GlobalParameterStruct.hpp"
#include "Exception.hpp"

//A pointer to the single parameter struct instance
GlobalParameterStruct* GlobalParameterStruct::mpInstance = NULL;


//For initialising the parameters object from a config file the first time it is used.
GlobalParameterStruct* GlobalParameterStruct::Instance(std::string filename, std::string parametersDirectory)
{
  if (mpInstance == NULL)
  {
    mpInstance = new GlobalParameterStruct(filename, parametersDirectory);
    std::atexit(Destroy);
    return mpInstance;
  }
  else{
    EXCEPTION("A parameter struct has already been created. To retrieve parameters, call GlobalParameterStruct::Instance()->GetParameter(*int*)");
  }
}


//For retrieving a parameter struct that has already been created 
GlobalParameterStruct* GlobalParameterStruct::Instance()
{
  if (mpInstance == NULL)
  {
    EXCEPTION("Please specify a parameter file name when calling GlobalParameterStruct for the first time.");
  }
  else{
    return mpInstance;
  }
}


//Reads in parameter data from a config file located in the directory given below
GlobalParameterStruct::GlobalParameterStruct(std::string filename, std::string parametersDirectory){

  std::string filepath;
  filepath.append(parametersDirectory);
  filepath.append(filename);
  std::cout << filepath.c_str() << std::endl;
  std::ifstream CONFIG(filepath.c_str());

  //read in line by line      
  double param;
  char temp[256];
  
  if (CONFIG.is_open())
  {
    std::cout << "parameters file opened" << std::endl;
    std::cout << "setting results directory name" << std::endl;
    CONFIG.getline(temp, 256);
    Directory = temp;
    std::cout << Directory << std::endl;

    while (!CONFIG.getline(temp, 256, '\t').eof())
    {
      param = strtod(temp, 0);
      std::cout << "reading parameter "<< Params.size() << "  "  << param << std::endl;
      Params.push_back(param);  
      CONFIG.getline(temp, 256);
    }
    std::cout << "parameters file closing" << std::endl;
    CONFIG.close();
  }else{
    EXCEPTION("Failed to open parameters file");
  }
}


//Resets the results directory name, can be useful for parameter sweeps
void GlobalParameterStruct::ResetDirectoryName(std::string newName){
  Directory = newName;
}


//Resets a single parameter value, can be useful in parameter sweeps
void GlobalParameterStruct::ResetParameter(int index, double newValue){
  Params[index] = newValue;
};


//Destroys the parameter struct object
void GlobalParameterStruct::Destroy()
{
  if (mpInstance)
  {
    delete mpInstance;
    mpInstance = NULL;
  }
}


//Retreives a parameter value by index
double GlobalParameterStruct::GetParameter(int index){
  return Params[index];
}


//Retreives results directory name
std::string GlobalParameterStruct::GetDirectory(){
  return Directory;
}