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
    "-S", "--setting", nargs="?", default="custom", help="setting directory to use"
)
parser.add_argument(
    "-p", "--params", required=True, help="directory or file with params yaml(s)"
)
parser.add_argument(
    "-o",
    "--outdir",
    nargs="?",
    default="results",
    help="directory name to save results to",
)
parser.add_argument(
    "-b",
    "--base",
    nargs="?",
    type=bool,
    default=True,
    help="whether to use base setting",
)

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def main(setting, paramsPath, nMC, outdir, use_base):
    wct = []  # wall clock times

    # delete old results before overwriting with new results
    outfile_dir = os.path.join(os.getcwd(), outdir)
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)
    os.mkdir(os.path.join(outfile_dir, "network"))

    # generate params - if no setting, set to null
    setting = setting.lower()
    if setting == "custom":
        setting = None
    else:
        setting = os.path.join("settings", setting)
        assert os.path.isdir(setting)

    params = create_params(setting, paramsPath, outfile_dir, use_base=use_base)

    for single_sim in range(nMC):
        tic = time_mod.time()

        # runs simulations
        model = HIVModel(params)
        stats = model.run(outfile_dir)

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
    main(
        args.setting.strip(),
        args.params.strip(),
        args.nMC,
        args.outdir.strip(),
        args.base,
    )
