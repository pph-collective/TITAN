#!/usr/bin/env python3
# encoding: utf-8

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

# import os
import time as time_mod

# import itertools
# import multiprocessing as mp
# import matplotlib.pyplot as plt

from simulation_lib import *
import params
from loadInput import *
import sys, os

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def main():
    # parameter_dict = read_parameter_dict(1)   # get parameters
    wct = []  # wall clock times
    open_outputs()

    # read_classifier_dict()
    for single_sim in range(params.N_MC):
        outfile_dir = os.path.join(
            os.getcwd(), "results/results_simulation_MP_%d" % single_sim
        )
        if not os.path.isdir(outfile_dir):
            os.mkdir(outfile_dir)
        tic = time_mod.time()
        parameters = None

        inputPopSeed = params.rSeed_pop
        inputNetSeed = params.rSeed_net
        inputRunSeed = params.rSeed_run

        if inputPopSeed == -1:
            inputPopSeed = single_sim + 1

        if inputNetSeed == -1:
            inputNetSeed = single_sim + 1

        # distribute simulations manually
        if params.N_MC % params.PROCESSES == 0:
            nreps = params.N_MC / params.PROCESSES
            num_runs = [nreps] * params.PROCESSES
        else:
            num_runs = [int(params.N_MC / params.PROCESSES)] * params.PROCESSES
            for i in range(params.N_MC % params.PROCESSES):
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
        # rslts = pool.map(simulation_star,combined_input)
        # rslts = simulation(combined_input)
        rslts = simulation(
            params.N_REPS,
            1,
            params.TIME_RANGE,
            params.N_POP,
            outfile_dir,
            parameters,
            runSeed=inputRunSeed,
            popSeed=inputPopSeed,
            netSeed=inputNetSeed,
            model=params.model,
        )
        wct.append(time_mod.time() - tic)
        save_results(params.N_MC, params.TIME_RANGE, rslts, outfile_dir, single_sim)
        # print rslts

    for task, time_t in enumerate(wct):
        print(("wall clock time on for simulation %d: %8.4f seconds" % (task, time_t)))

    def mean(seq):
        return sum(seq) / len(seq)

    print(("\nSUMMARY:\nall tasks - mean: %8.4f seconds" % mean(wct)))
    print(("all tasks - min:  %8.4f seconds" % min(wct)))
    print(("all tasks - max:  %8.4f seconds" % max(wct)))
    print(("all tasks - sum:  %8.4f seconds" % sum(wct)))

    # plt.show(block=True)


def open_outputs():
    off = 0
    kickoff = 0
    # dynReport = open('results/dynnetworkReport.txt','w')
    # dynReport.write('')
    # dynReport.close()
    for agentRace in ["WHITE", "BLACK", "ALL"]:
        for agentTypes in params.agentPopulations:
            name = "basicReport_" + agentTypes + "_" + agentRace
            # print name
            tmpReport = open("results/" + name + ".txt", "w")
            tmpReport.write(
                "rseed\tpseed\tnseed\tt\tTotal\tHIV\tAIDS\tTstd\tART\tnHR\tIncid\tHR_6mo\tHR_Ev\tNewDiag\tDeaths\tPrEP\n"
            )
            tmpReport.close()

    for demographicTypes in list(params.DemographicParams.keys()):
        name = "basicReport_" + demographicTypes
        tmpReport = open("results/" + name + ".txt", "w")
        tmpReport.write(
            "rseed\tpseed\tnseed\tt\tTotal\tHIV\tAIDS\tTstd\tART\tnHR\tIncid\tHR_6mo\tHR_Ev\tNewDiag\tDeaths\tPrEP\n"
        )
        # whiteReport.write("0,0,0,0,0\n")
        tmpReport.close()

    # component report file creation
    name = "componentReport_ALL"
    tmpReport = open("results/" + name + ".txt", "w")
    tmpReport.write(
        "rseed\tpseed\tnseed\tt\tcompID\ttotalN\tNhiv\tNprepElig\tNprep\tNnewinf\n"
    )
    tmpReport.close()

    whiteReport = open("results/W_pop_report.txt", "w")
    whiteReport.write("seed\tt\tTotal-HIV\tMSM\tTested+\tHAART\n")
    # whiteReport.write("0,0,0,0,0\n")
    whiteReport.close()

    blackReport = open("results/B_pop_report.txt", "w")
    blackReport.write("seed\tt\tTotal-HIV\tMSM\tTested+\tHAART\n")
    # blackReport.write("0,0,0,0,0\n")
    blackReport.close()

    incidenceReport = open("results/IncidenceReport.txt", "w")
    incidenceReport.write(
        "seed\tt\tTotal\tW_HM\tB_HM\tHM\tW_HF\tB_HF\tHF\tW_MSM\tB_MSM\tMSM\n"
    )
    # incidenceReport.write("0,0,0,0\n")
    incidenceReport.close()

    prevalenceReport = open("results/PrevalenceReport.txt", "w")
    prevalenceReport.write("seed\tt\tTotal\tHM\tHF\tHIV_tot\tHIV_HM\tHIV_HF\n")
    # prevalenceReport.write("0,0,0,0,0,0\n")
    prevalenceReport.close()

    incarReport = open("results/IncarReport.txt", "w")
    incarReport.write(
        "seed\tt\tTotal\tW_HM\tB_HM\tW_HF\tB_HF\tW_MSM\tB_MSM\tW_HIV\tB_HIV\tW_rlsd_HM\tW_rlsd_HF\tB_rlsd_HM\tB_rlsd_HF\tW_rlsdHIV\tB_rlsdHIV\n"
    )
    # prevalenceReport.write("0,0,0,0,0,0\n")
    incarReport.close()

    newlyhighriskReport = open("results/newlyHR_Report.txt", "w")
    newlyhighriskReport.write(
        "seed\tt\tnewHR_HM\tnewHR_HIV_HM\tnewHR_AIDS_HM\tnewHR_Tested_HM\tnewHR_ART_HM\tnewHR_HF\tnewHR_HIV_HF\tnewHR_AIDS_HF\tnewHR_Tested_HF\tnewHR_ART_HF\n"
    )
    newlyhighriskReport.close()

    highriskReport = open("results/HR_incidenceReport.txt", "w")
    highriskReport.write("seed\tt\tTot_Ev\tHM_Ev\tHF_Ev\tTot_6mo\tHM_6mo\tHF_6mo\n")
    # prevalenceReport.write("0,0,0,0,0,0\n")
    highriskReport.close()

    deathReport = open("results/DeathReport.txt", "w")
    deathReport.write("seed\tt\tTotal\tHM\tMSM\tHF\tHIV_tot\tHIV_HM\tHIV_MSM\tHIV_HF\n")
    # prevalenceReport.write("0,0,0,0,0,0\n")
    deathReport.close()

    iduReport = open("results/iduReport.txt", "w")
    iduReport.write("seed\tt\tTotal-IDU\tIDU-HIV\tIDU-AIDS\tIDU-HAART\tIDU-tested\n")
    # iduReport.write("0,0,0,0,0,0\n")
    iduReport.close()
    # plt.ion()

    femaleReport = open("results/FemaleReport.txt", "w")
    femaleReport.write(
        "seed\tt\tTotal\tHIV\tAIDS\tTested\tART\tIncidence\tHRInc_6mo\tHRInc_Ev\tNewlyDiag\tDeaths\n"
    )
    # prevalenceReport.write("0,0,0,0,0,0\n")
    femaleReport.close()

    maleReport = open("results/MaleReport.txt", "w")
    maleReport.write(
        "seed\tt\tTotal\tHIV\tAIDS\tTested\tART\tIncidence\tHRInc_6mo\tHRInc_Ev\tNewlyDiag\tDeaths\tPrEP\n"
    )
    # prevalenceReport.write("0,0,0,0,0,0\n")
    maleReport.close()

    msmReport = open("results/MSMReport.txt", "w")
    msmReport.write(
        "seed\tt\tTotal\tHIV\tAIDS\tTested\tART\tIncidence\tHRInc_6mo\tHRInc_Ev\tNewlyDiag\tDeaths\tPrEP\n"
    )
    # prevalenceReport.write("0,0,0,0,0,0\n")
    msmReport.close()

    # for i in range(1,6):
    #     ageNReport = open('results/MSMReport_a%d.txt' %i, 'w')
    #     ageNReport.write('seed,t,total_N,HIV,Tested,ART,PrEP\n')
    #     ageNReport.close()


if __name__ == "__main__":
    main()
