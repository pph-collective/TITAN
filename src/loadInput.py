#!/usr/bin/env python

"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University


Description:
    Module responsible for loading input parameters from csv files.


Copyright (c) 2016, Maximilian King
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
"""


from __future__ import print_function
from six.moves import range
def read_input_parameter(paramFile=None):
    # Read parameters from file, put them into the
    # dictionary parameter_dict and return parameter_dict

    # Read scalars
    if paramFile:
        infile = open(('input/' + paramFile + '.csv'),'r')
    else:
        infile = open('input/InputParameters.csv','r')
    lines = infile.readlines()
    infile.close()
    data_dict = {} # data[run#][parameter]
    first_line = lines[0]
    properties = first_line.split(',')
    properties.remove('Name')
    properties.remove('Description')
    #NumSimulations = len(properties)
    NumSimulations = num_Simulations
    for num_sim in range(NumSimulations):
        data_dict[num_sim] = {}
    for line in lines[1:]:
        words = line.split(',')
        name_key = words[0]
        values = words[2:]
        for i, value in enumerate(values):
            data_dict[i][name_key] = float(value.strip())

    # Read vector parameters (function of time t)
    for vector_value in ['NSP_SAT','NSP_NoSAT']:
        text = open(('input/' + vector_value + '.csv'),'r').read()
        if '\r\n' in text: lines = text.split('\r\n')
        else: lines = text.split('\r')
        first_line = lines[0]
        properties = first_line.split(',')
        properties.remove('Time')
        properties.remove('Description')
        #assert len(properties)==NumSimulations,('Inconsistent parameter files!'+
        #'\nInputParameter.csv, NSP_SAT.csv, and NSP_NoSAT.csv must have '+
        #'the same number of simulations! Each column contains the parameter '+
        #'set for one simulation.\nlen(properties) = %d\nNumSimulations = %d\n'%(
        #len(properties),NumSimulations))
        for num_run in range(NumSimulations):
            data_dict[num_run].update({vector_value:{}})
        for line in lines[1:]:
            words = line.split(',')
            t = int(words[0].strip())
            values = words[2:] # Time and Description column offset
            for i, value in enumerate(values):
                data_dict[i][vector_value][t] = float(value.strip())

    return data_dict

def read_classifier_dict(input_params=None):
    # Read parameters from file, put them into the
    # dictionary parameter_dict and return parameter_dict

    # Read scalars
    num_Races = 1#input_params['numRaceClassifiers']
    infile = open('input/DemographicParameters.csv','r')
    lines = infile.readlines() # name - race - so - variable - description
    infile.close()
    data_dict = {} # data[race][so][variable] = value
    for line in lines[1:]:
        words = line.split(',')
        name_key = words[0]
        race_key = words[1]
        so_key = words[2]
        values = words[3]
        print(name_key)
        print(race_key)
        print(so_key)
        print(values)
        for value in enumerate(values):
            data_dict[race_key][so_key][name_key] = value

    print(data_dict)