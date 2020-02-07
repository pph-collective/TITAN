#!/usr/bin/env python3
# encoding: utf-8

import time as time_mod
import sys
import os
import shutil
import argparse

from titan.simulation_lib import simulation, save_results
from titan.params_parse import create_params

# set up args parsing
parser = argparse.ArgumentParser(description="Run TITAN model")
parser.add_argument(
    "-n",
    "--nMC",
    type=int,
    nargs="?",
    default=1,
    help="number of monte carlo runs to complete",
)
parser.add_argument(
    "-s", "--setting", nargs="?", default="custom", help="setting directory to use"
)
parser.add_argument(
    "-p", "--params", required=True, help="directory or file with params yaml(s)"
)

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def main(setting, paramsPath, nMC):
    wct = []  # wall clock times

    # delete old results before overwriting with new results
    outfile_dir = os.path.join(os.getcwd(), "results")
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)
    if not os.path.exists("results/network"):
        os.makedirs("results/network")

    # generate params - if no setting, set to null
    if setting == "custom":
        setting = {}
    # FIGURE OUT ELSE

    params = create_params(setting, paramsPath)

    for single_sim in range(nMC):
        outfile_dir = os.path.join(
            os.getcwd(), f"results/results_simulation_MP_{single_sim}"
        )
        if not os.path.isdir(outfile_dir):
            os.mkdir(outfile_dir)
        tic = time_mod.time()

        inputPopSeed = params.model.seed.ppl
        inputNetSeed = params.model.seed.net
        inputRunSeed = params.model.seed.run

        # get rid of num_reps and -1 seed shenanigans
        if inputPopSeed == -1:
            inputPopSeed = single_sim + 1

        if inputNetSeed == -1:
            inputNetSeed = single_sim + 1

        # runs simulations
        rslts = simulation(params)

        wct.append(time_mod.time() - tic)
        save_results(params, rslts, outfile_dir, single_sim)

    for task, time_t in enumerate(wct):
        print(("wall clock time on for simulation %d: %8.4f seconds" % (task, time_t)))

    def mean(seq):
        return sum(seq) / len(seq)

    print(("\nSUMMARY:\nall tasks - mean: %8.4f seconds" % mean(wct)))
    print(("all tasks - min:  %8.4f seconds" % min(wct)))
    print(("all tasks - max:  %8.4f seconds" % max(wct)))
    print(("all tasks - sum:  %8.4f seconds" % sum(wct)))


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.setting, args.params, args.nMC)
