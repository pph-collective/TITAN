#!/usr/bin/env python
# encoding: utf-8

"""
*****************************************************************************
Author(s):  Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University

Description:
    Library of functions needed for multiprocessing or MPI simulations.


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


import os
import time as time_mod
import itertools
import pprint
# import random
import matplotlib.pyplot as plt
import params
import numpy as np

from ABM_core import HIVModel

try:
    from HIVABM_Population import PopulationClass, print_population
except ImportError:
    raise ImportError("Can't import PopulationClass")


# from Evolution import HIVModel

def read_parameter_dict(num_Simulations):
    # Read parameters from file, put them into the
    # dictionary parameter_dict and return parameter_dict

    # Read scalars
    infile = open('input/InputParameters.csv', 'r')
    lines = infile.readlines()
    infile.close()
    data_dict = {}  # data[run#][parameter]
    first_line = lines[0]
    properties = first_line.split(',')
    properties.remove('Name')
    properties.remove('Description')
    # NumSimulations = len(properties)
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
    for vector_value in ['NSP_SAT', 'NSP_NoSAT']:
        text = open(('input/' + vector_value + '.csv'), 'r').read()
        if '\r\n' in text:
            lines = text.split('\r\n')
        else:
            lines = text.split('\r')
        first_line = lines[0]
        properties = first_line.split(',')
        properties.remove('Time')
        properties.remove('Description')
        # assert len(properties)==NumSimulations,('Inconsistent parameter files!'+
        # '\nInputParameter.csv, NSP_SAT.csv, and NSP_NoSAT.csv must have '+
        # 'the same number of simulations! Each column contains the parameter '+
        # 'set for one simulation.\nlen(properties) = %d\nNumSimulations = %d\n'%(
        # len(properties),NumSimulations))
        for num_run in range(NumSimulations):
            data_dict[num_run].update({vector_value: {}})
        for line in lines[1:]:
            words = line.split(',')
            t = int(words[0].strip())
            values = words[2:]  # Time and Description column offset
            for i, value in enumerate(values):
                data_dict[i][vector_value][t] = float(value.strip())

    return data_dict


def simulation(nreps, save_adjlist_flag, time_range,
               N_pop, outfile_dir, parameters, rSeed, uniqueSeed=False, model=None):
    # Check input
    if save_adjlist_flag not in [0, 1]:
        raise ValueError('Invalid input! save_adjlist_flag = %s' %
                         str(save_adjlist_flag))
    # if time_range != 30:
    #   raise Warning('time_range=%d'%time_range)
    # Run nreps simulations using the given parameters.
    # Information are printed to outfile_dir directory.
    pid = os.getpid()
    # print "Process %5s runs %d simulations:"%(pid, nreps)
    # pprint.pprint(parameters)
    result_dict = {}

    for num_sim in range(nreps):
        # random.seed(rSeed)
        inputSeed = rSeed
        if rSeed == -1:
            inputSeed = num_sim + 1
        #print "\n\n------------------------------------------------------------------------------------------------------------------------------------------"
        print "\tProcess %5s runs simulation %d/%d\t.:.\tInput rSeed: %d" % (pid, num_sim + 1, nreps, inputSeed)
        # fixedPop = PopulationClass(n=N_pop, rSeed = rSeed, model=model)._return_agent_set()
        fixedPop = None

        MyModel = HIVModel(N=N_pop, tmax=time_range, parameter_dict=parameters, rseed=inputSeed, runtime_diffseed=uniqueSeed, model=model,
                           network_type=fixedPop)

        if save_adjlist_flag == 1 and num_sim == 0:
            MyModel.run(save_adjlist_flag=1, dir_prefix=outfile_dir)
        else:
            MyModel.run(save_adjlist_flag=0, dir_prefix=outfile_dir)

        #print MyModel
        result_dict_tmp = MyModel.return_results()
        for (key, x_v) in result_dict_tmp.iteritems():
            if key not in result_dict:
                result_dict.update({key: {}})
            for t, x in x_v.iteritems():
                if t not in result_dict[key]:
                    result_dict[key].update({t: []})
                if np.isnan(x):
                    # mssg = '\tWARNING:\
                    # NaN encountered in pid %5s run %d for %s at time %d'
                    # print mssg%(pid, num_sim,key, t)
                    result_dict[key][t].append(0)
                else:
                    result_dict[key][t].append(x)
    return result_dict


def simulation_star(zipped_input):
    """
    This helper function is needed because multithreading pool only
    accepts one input argument. This helper function converts combined
    arguments and calls the simulation function.
    """
    try:
        return simulation(*zipped_input)
    except TypeError:
        for input_info in zipped_input:
            print input_info
        raise TypeError("Wrong input for simulation_star()!")


def save_results(N_MC, time_range, rslts, outfile_dir, num_sim):
    """
    Save the result dictionary.
    Input:   N_MC : Number of Monte Carlo runs for each simulation.
             time_range : Number of time steps in each simulation.
             rslts : Result dictionary.
             outfile_dir : Directory for output file.
             num_sim : Number of simulation. IMPORTANT to identify
                       the applied parameters. Must be consistent with
                       parameter input file columns.
    """
    # convert list of dicts into one nested result dict
    result_dict = {}
    """for tmp_dict in rslts:
       for (key, tmp_dict_L2) in tmp_dict.iteritems():
          if key not in result_dict:
             result_dict.update({key:{}})
          for (t,x) in tmp_dict_L2.iteritems():
             if t not in result_dict[key]:
                result_dict[key].update({t:[]})
             result_dict[key][t].extend(x)"""

    # save results
    if not os.path.isdir(outfile_dir):
        os.mkdir(outfile_dir)
    OutFileName = os.path.join(outfile_dir, ('Result_simulation_%d.txt' % num_sim))
    print "\n\tSaving results to:\n\t%s\n" % str(OutFileName)
    if os.path.isfile(OutFileName):
        os.remove(OutFileName)
    outfile = open(OutFileName, 'w')
    outfile.write('t,model,coverage')
    for result_property in sorted(rslts):
        outfile.write(',%s_mean,%s_std,%s_10th,%s_90th'% (result_property,result_property,result_property,result_property))
    outfile.write('\n')
    #for result_property in sorted(rslts):  # result_dict.keys()):
    for t in range(0, time_range + 1):
        #outfile.write('%s\t' % result_property)
        #print result_property
        """
        # print SUM
        outfile.write('%s\tSum\t'%result_property)
        for t in sorted(rslts[result_property]):  # result_dict[result_property].keys()):
            x_v = np.array(rslts[result_property][t])  # result_dict[result_property][t])
            x_v = x_v[np.logical_not(np.isnan(x_v))]
            if len(x_v) > 0:
                mean = np.sum(x_v)
                outfile.write('%4.2f\t' % mean)
            else:
                outfile.write('%4.5f\t' % np.NaN)
        outfile.write('\n')
        """

        outfile.write('%d,%s,%0.2f'%(t,params.PrEP_type,params.PrEP_Target))
        for result_property in sorted(rslts):  # result_dict.keys()):
            #outfile.write('%s\tMean\t'%result_property)
            x_v = []

            try:
                x_v = np.array(rslts[result_property][t])  # result_dict[result_property][t])
                x_v = x_v[np.logical_not(np.isnan(x_v))]
            except:pass


            # print sum
            # if len(x_v) > 0:
            #     sum = np.sum(x_v)
            #     outfile.write(',%4.5f' % sum)
            # else:
            #     outfile.write(',%4.5f' % np.NaN)

            # print mean
            if len(x_v) > 0:
                mean = np.mean(x_v)
                outfile.write(',%4.5f' % mean)
            else:
                outfile.write(',%4.5f' % np.NaN)
            # print std
            if len(x_v)>0:
                std_dev=np.std(x_v)
                outfile.write(',%4.5f'%std_dev)
            else:
                outfile.write(',%4.5f'%np.NaN)

            # print 10th
            if len(x_v)>0:
              p10=np.percentile(x_v,10)
            else:
                p10=np.NaN
            outfile.write(',%4.5f'%p10)

            # print 90th
            if len(x_v)>0:
              p90=np.percentile(x_v,90)
            else:
                p90=np.NaN
            outfile.write(',%4.5f'%p90)
            #outfile.write('\n')

            # print std. deviation of mean
            #outfile.write('%s\tStd.Dev\t'%result_property)
        outfile.write('\n')

    outfile.close()


if __name__ == '__main__':
    print 'Use MP_simulation and MPI_simulation!\n'
