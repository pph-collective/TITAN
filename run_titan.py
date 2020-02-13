#!/usr/bin/env python3
# encoding: utf-8

import time as time_mod
import sys
import os
import shutil
import argparse

from titan.ABM_core import HIVModel
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
        tic = time_mod.time()

        # runs simulations
        model = HIVModel(params)
        stats = model.run()

        wct.append(time_mod.time() - tic)

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
