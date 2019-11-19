#!/usr/bin/env python3
# encoding: utf-8

import time as time_mod
import sys
import os
import shutil

from titan.simulation_lib import simulation, save_results
from titan import params


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def main():
    wct = []  # wall clock times

    # delete old results before overwriting with new results
    outfile_dir = os.path.join(os.getcwd(), "results")
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)

    for single_sim in range(params.N_MC):
        outfile_dir = os.path.join(
            os.getcwd(), "results/results_simulation_MP_%d" % single_sim
        )
        if not os.path.isdir(outfile_dir):
            os.mkdir(outfile_dir)
        tic = time_mod.time()

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
        rslts = simulation(
            params.N_REPS,
            params.TIME_RANGE,
            params.N_POP,
            outfile_dir,
            runSeed=inputRunSeed,
            popSeed=inputPopSeed,
            netSeed=inputNetSeed,
        )
        wct.append(time_mod.time() - tic)
        save_results(params.TIME_RANGE, rslts, outfile_dir, single_sim)

    for task, time_t in enumerate(wct):
        print(("wall clock time on for simulation %d: %8.4f seconds" % (task, time_t)))

    def mean(seq):
        return sum(seq) / len(seq)

    print(("\nSUMMARY:\nall tasks - mean: %8.4f seconds" % mean(wct)))
    print(("all tasks - min:  %8.4f seconds" % min(wct)))
    print(("all tasks - max:  %8.4f seconds" % max(wct)))
    print(("all tasks - sum:  %8.4f seconds" % sum(wct)))


if __name__ == "__main__":
    main()
