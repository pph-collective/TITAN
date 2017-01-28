#!/usr/bin/env python

"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University


Description:
    Module responsible for simulation handling and total loops. Manages input for
    simulation population size, number of total jobs, and job duration.


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

#import os
import time as time_mod
#import itertools
#import multiprocessing as mp
#import matplotlib.pyplot as plt

from simulation_lib import *

open('Results/dynnetworkReport.txt','w').write('')

# Global parameters:
PROCESSES = 1       # number of processes in parallel (quadcore)
N_MC = 1            # total number of iterations (Monte Carlo runs)
N_POP = 23815         # population size
TIME_RANGE = 36    # time range for simulation

def main():
    parameter_dict = read_parameter_dict(1)   # get parameters
    wct = []                                 # wall clock times

    whiteReport = open('Results/W_pop_report.txt', 'w')
    whiteReport.write("t\tTotal-HIV\tMSM\tTested+\tHAART\n")
    #whiteReport.write("0\t0\t0\t0\t0\n")
    whiteReport.close()

    blackReport = open('Results/B_pop_report.txt', 'w')
    blackReport.write("t\tTotal-HIV\tMSM\tTested+\tHAART\n")
    #blackReport.write("0\t0\t0\t0\t0\n")
    blackReport.close()

    incidenceReport = open('Results/IncidenceReport.txt', 'w')
    incidenceReport.write("t\tTotal\tHM\tHF\n")
    #incidenceReport.write("0\t0\t0\t0\n")
    incidenceReport.close()

    prevalenceReport = open('Results/PrevalenceReport.txt', 'w')
    prevalenceReport.write("t\tTotal\tHM\tHF\tHIV_tot\tHIV_HM\tHIV_HF\n")
    #prevalenceReport.write("0\t0\t0\t0\t0\t0\n")
    prevalenceReport.close()

    incarReport = open('Results/IncarReport.txt', 'w')
    incarReport.write("t\tTotal\tHM\tHF\tHIV_tot\tHIV_HM\tHIV_HF\n")
    # prevalenceReport.write("0\t0\t0\t0\t0\t0\n")
    incarReport.close()

    highriskReport = open('Results/highriskReport.txt', 'w')
    highriskReport.write("t\tTot_Ev\tHM_Ev\tHF_Ev\tTot_6mo\tHM_6mo\tHF_6mo\n")
    # prevalenceReport.write("0\t0\t0\t0\t0\t0\n")
    highriskReport.close()

    deathReport = open('Results/DeathReport.txt', 'w')
    deathReport.write("t\tTotal\tHM\tHF\tHIV_tot\tHIV_HM\tHIV_HF\n")
    # prevalenceReport.write("0\t0\t0\t0\t0\t0\n")
    deathReport.close()



    iduReport = open('Results/iduReport.txt', 'w')
    iduReport.write("t\tTotal-IDU\tIDU-HIV\tIDU-AIDS\tIDU-HAART\tIDU-tested\n")
    #iduReport.write("0\t0\t0\t0\t0\t0\n")
    iduReport.close()
    #plt.ion()

    for single_sim in parameter_dict:
        outfile_dir = os.path.join(os.getcwd(),
                                   'Results/Results_simulation_MP_%d'%single_sim)
        if not os.path.isdir(outfile_dir):
            os.mkdir(outfile_dir)
        tic = time_mod.time()
        parameters = parameter_dict[0]
        rSeed = 0
        # distribute simulations manually
        if N_MC%PROCESSES == 0:
            nreps = N_MC/PROCESSES
            num_runs = [nreps]*PROCESSES
        else:
            num_runs = [int(N_MC/PROCESSES)]*PROCESSES
            for i in range(N_MC%PROCESSES):
                num_runs[i] += 1

        # runs simulations
        """pool = mp.Pool(PROCESSES)
        save_adj_list_flags = [1] + [0]*(len(num_runs)-1)
        combined_input = itertools.izip(num_runs,
                                        save_adj_list_flags,
                                        itertools.repeat(TIME_RANGE),
                                        itertools.repeat(N_POP),
                                        itertools.repeat(outfile_dir),
                                        itertools.repeat(parameters))
        """
        #rslts = pool.map(simulation_star,combined_input)
        #rslts = simulation(combined_input)
        rslts = simulation(nreps, 1, TIME_RANGE, N_POP, outfile_dir, parameters, rSeed)
        wct.append(time_mod.time() - tic)
        save_results(N_MC, TIME_RANGE, rslts, outfile_dir, single_sim)
        #if single_sim == 1:
        #print rslts

    for task, time_t in enumerate(wct):
        print('wall clock time on for simulation %d: %8.4f seconds' %
              (task, time_t))
    def mean(seq): return sum(seq)/len(seq)
    print    ('\nSUMMARY:\nall tasks - mean: %8.4f seconds' % mean(wct))
    print    ('all tasks - min:  %8.4f seconds' % min(wct))
    print    ('all tasks - max:  %8.4f seconds' % max(wct))
    print    ('all tasks - sum:  %8.4f seconds' % sum(wct))

    #plt.show(block=True)

if __name__ == '__main__':
    main()
