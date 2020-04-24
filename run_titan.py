#!/usr/bin/env pypy3
# encoding: utf-8

import time as time_mod
from copy import copy
import sys
import os
import shutil
import argparse
import itertools
import json
from multiprocessing import Pool, cpu_count

from titan.model import HIVModel
from titan.parse_params import create_params

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

# how many cores can we use
NCORES = os.environ.get("SLURM_CPUS_PER_TASK", cpu_count())
NCORES = int(NCORES)  # environment variable returns string


def sweep_range(string):
    error_msg = "Sweep range must have format param:start:stop[:step]"
    parts = string.split(":")
    assert len(parts) in (3, 4), error_msg
    if len(parts) == 4:
        step = parts[3]
    else:
        step = "1"

    try:
        start = int(parts[1])
        stop = int(parts[2])
        step = int(step)
    except:
        try:
            start = float(parts[1])
            stop = float(parts[2])
            step = float(step)
        except:
            raise ValueError("start, stop, and step must have same type (int or float)")

    return {"param": parts[0], "start": start, "stop": stop, "step": step}


parser.add_argument(
    "-w",
    "--sweep",
    nargs="*",
    type=sweep_range,
    default=[],
    help="Optional and repeatable definitions of numeric params to sweep. Expected format is param:start:stop[:step]",
)

parser.add_argument(
    "-F",
    "--force",
    action="store_true",
    help="Run model even if number of sweeps excedes 100",
)


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step
        r = round(r, 3)


def setup_sweeps(sweeps):
    sweep_params = [s["param"] for s in sweeps]
    sweep_vals = list(
        itertools.product(
            *[list(drange(s["start"], s["stop"], s["step"])) for s in sweeps]
        )
    )

    defs = []
    for val in sweep_vals:
        defn = {}
        for i, param in enumerate(sweep_params):
            defn[param] = val[i]
        defs.append(defn)

    return defs


def consolidate_files(outdir):
    for item in os.listdir(outdir):
        subdir = os.path.join(outdir, item)
        if os.path.isdir(subdir) and item != "network":
            for report in os.listdir(subdir):
                # network folder
                if report == "network":
                    for file in os.listdir(os.path.join(subdir, report)):
                        shutil.move(
                            os.path.join(subdir, report, file),
                            os.path.join(outdir, network),
                        )
                else:
                    # copy data to existing file
                    if os.path.isfile(os.path.join(outdir, report)):
                        if report == "SweepVals.json":
                            header_skipped = True
                        else:
                            header_skipped = False
                        with open(os.path.join(outdir, report), "a") as tgt, open(
                            os.path.join(subdir, report), "r"
                        ) as src:
                            for line in src:
                                if not header_skipped:
                                    header_skipped = True
                                else:
                                    tgt.write(line)

                    else:
                        shutil.move(os.path.join(subdir, report), os.path.join(outdir))

            # after copying remove directory
            shutil.rmtree(subdir)


def update_sweep_file(run_id, defn, outdir):
    f = open(os.path.join(outdir, "SweepVals.json"), "a")
    res = copy(defn)
    res["run_id"] = str(run_id)
    f.write(json.dumps(res))
    f.write("\n")
    f.close()


def single_run(sweep, outfile_dir, params):
    pid = str(os.getpid())
    pid_outfile_dir = os.path.join(outfile_dir, pid)
    if not os.path.isdir(pid_outfile_dir):
        os.mkdir(pid_outfile_dir)
        os.mkdir(os.path.join(pid_outfile_dir, "network"))

    for param, val in sweep.items():
        print(f"\t{param}: {val}")
        path = param.split(".")
        sweep_item = params
        for p in path[:-1]:
            sweep_item = sweep_item[p]
        sweep_item[path[-1]] = val

    tic = time_mod.time()

    # runs simulations
    model = HIVModel(params)
    run_id = model.run(pid_outfile_dir)

    update_sweep_file(run_id, sweep, pid_outfile_dir)

    return time_mod.time() - tic


def main(setting, params_path, num_reps, outdir, use_base, sweeps, force):
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
        setting = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "settings", setting
        )
        assert os.path.isdir(setting), f"{setting} is not a directory"

    params = create_params(setting, params_path, outfile_dir, use_base=use_base)

    # set up sweeps
    if sweeps:
        sweep_defs = setup_sweeps(sweeps)
    else:
        sweep_defs = [{}]  # make sure runs once

    if len(sweep_defs) > 100 and not force:
        raise ValueError(
            "Sweeping more than 100 models. Set `-F` (force) flag if you really want to do this."
        )

    sweep_defs *= num_reps

    tic = time_mod.time()
    wct = []  # wall clock times

    with Pool(processes=NCORES) as pool:
        results = [
            pool.apply_async(single_run, (sweep_def, outfile_dir, params))
            for sweep_def in sweep_defs
        ]
        while True:
            if all([r.ready() for r in results]):
                break

        for r in results:
            t = r.get()
            wct.append(t)

    toc = time_mod.time() - tic

    consolidate_files(outfile_dir)

    for task, time_t in enumerate(wct):
        print(("wall clock time on for simulation %d: %8.4f seconds" % (task, time_t)))

    def mean(seq):
        return sum(seq) / len(seq)

    print(("\nSUMMARY:\nall tasks - mean: %8.4f seconds" % mean(wct)))
    print(("all tasks - min:  %8.4f seconds" % min(wct)))
    print(("all tasks - max:  %8.4f seconds" % max(wct)))
    print(("all tasks - sum:  %8.4f seconds" % sum(wct)))
    print(f"all tasks - total: {toc} seconds")


if __name__ == "__main__":
    args = parser.parse_args()
    main(
        args.setting.strip(),
        args.params.strip(),
        args.nMC,
        args.outdir.strip(),
        args.base,
        args.sweep,
        args.force,
    )
