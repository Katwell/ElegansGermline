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

#include "GlobalParameterStruct.hpp"
#include "Exception.hpp"


//A pointer to the single parameter struct instance. Initially null.
GlobalParameterStruct* GlobalParameterStruct::mpInstance = NULL;


//For retrieving a pointer to the current GlobalParameterStruct struct 
GlobalParameterStruct* GlobalParameterStruct::Instance()
{
  if (mpInstance == NULL)
  {
    mpInstance = new GlobalParameterStruct();
  }
  return mpInstance;
}


//Protected constructor. Leaves directory and parameters vector blank for now
GlobalParameterStruct::GlobalParameterStruct()
{
    Directory =  std::string();
    Params = std::vector<double>();
    assert(mpInstance == NULL); 
}



//Initialises the parameters vector and output directory from a config file.
void GlobalParameterStruct::ConfigureFromFile(std::string filename, std::string parametersDirectory)
{

  //Open the file specified in the arguments
  std::string filepath;
  filepath.append(parametersDirectory);
  filepath.append(filename);
  std::ifstream CONFIG(filepath.c_str());
    
  double param;   // <- temp parameter storage
  char temp[256]; // <- temp string buffer
  if (CONFIG.is_open())
  {
    //If file opened fine, read it line by line  

    //Set output directory
    //std::cout << "Parameters file open" << std::endl;
    //std::cout << "Setting output directory to ";
    CONFIG.getline(temp, 256);
    Directory = temp;
    //std::cout << Directory << std::endl;

    while (!CONFIG.getline(temp, 256, '\t').eof())
    {
      //Get each param and ftor in Params vector.
      param = strtod(temp, 0);
      //std::cout << "Parameter " << Params.size() << " = " << param << std::endl;
      Params.push_back(param);  
      CONFIG.getline(temp, 256);
    }
    CONFIG.close();

  }else{
    EXCEPTION("Failed to open parameters file");
  }
}



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
  if(index > (int)Params.size()-1){
    EXCEPTION("Parameter has yet to be initialised. Check that ConfigureFromFile was called and that you are using the correct input file.");
  }
  return Params[index];
}


//Retreives the results directory name
std::string GlobalParameterStruct::GetDirectory(){
  return Directory;
}


//Resets the results directory name. Can be useful for parameter sweeps
void GlobalParameterStruct::ResetDirectoryName(std::string newName){
  Directory = newName;
}


//Resets a single parameter value. Can be useful in parameter sweeps
void GlobalParameterStruct::ResetParameter(int index, double newValue){
  Params[index] = newValue;
};

