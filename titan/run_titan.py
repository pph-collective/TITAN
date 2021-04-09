import time as time_mod
from copy import copy
import sys
import os
import shutil
import argparse
import itertools
import json
from multiprocessing import Pool, cpu_count
import csv
import traceback
from typing import List, Optional
import logging

# allow imports to work if running it as a script for development locally
if __name__ == "__main__":
    PACKAGE_PARENT = ".."
    SCRIPT_DIR = os.path.dirname(
        os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__)))
    )
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from titan.model import TITAN
from titan.population import Population
import titan.population_io as pop_io
from titan.parse_params import create_params
from titan import utils

# how many cores can we use, environment variable returns string
NCORES = int(os.environ.get("SLURM_CPUS_PER_TASK", cpu_count()))

# set up args parsing
parser = argparse.ArgumentParser(description="Run TITAN model")
parser.add_argument(
    "-p", "--params", required=True, help="directory or file with params yaml(s)"
)
parser.add_argument(
    "-S", "--setting", default="custom", help="setting directory to use"
)
parser.add_argument(
    "-n",
    "--nMC",
    type=int,
    nargs="?",
    default=1,
    help="number of monte carlo runs to complete",
)
parser.add_argument(
    "-o", "--outdir", default="results", help="directory name to save results to"
)

parser.add_argument(
    "-e",
    "--error",
    action="store_true",
    help="Error on unused parameters instead of warning",
)

parser.add_argument(
    "--savepop",
    action="store_true",
    help="Save population after creation, but before model run.",
)

parser.add_argument(
    "--poppath",
    type=str,
    default=None,
    help="Path to saved population (directory or .tar.gz file)",
)


def sweep_range(string):
    """
    parsing and checking of param:start:stop:step into dictionary
    """
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
    except ValueError:
        try:
            start = float(parts[1])
            stop = float(parts[2])
            step = float(step)
        except ValueError:
            raise ValueError("start, stop, and step must have same type (int or float)")

    return {"param": parts[0].strip(), "start": start, "stop": stop, "step": step}


parser.add_argument(
    "-w",
    "--sweep",
    nargs="+",
    type=sweep_range,
    default=[],
    help="Optional and repeatable definitions of numeric params to sweep. Expected format is param:start:stop[:step]",
)

parser.add_argument(
    "-W",
    "--sweepfile",
    type=str,
    default=None,
    help="Optional. CSV file with param sweep definitions. Header row must contain param paths, with data rows containing values. If this is passed, any `-w` args will be ignored.",
)

parser.add_argument(
    "-r",
    "--rows",
    type=str,
    default=None,
    help="Optional. Which data rows of sweepfile to use in format start:stop.",
)


parser.add_argument(
    "-F",
    "--force",
    action="store_true",
    help="Run model even if number of sweeps exceeds 100",
)


def drange(start, stop, step):
    """
    Like range, but allows decimal steps
    """
    r = start
    while r < stop:
        yield r
        r += step
        r = round(r, 6)


def setup_sweeps(sweeps):
    """
    Set up sweeps definitions from command line definitions
    """
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


def setup_sweeps_file(sweepfile, rows):
    """
    Set up sweeps definitions (params to values) from csv file
    """
    if rows is not None:
        start, stop = (int(item) for item in rows.split(":"))
    else:
        start = 1
        stop = sys.maxsize  # very very big number

    defs = []
    with open(sweepfile, newline="") as f:
        reader = csv.DictReader(f)
        n = 1
        for row in reader:
            if start <= n <= stop:
                defn = {}
                for k, v in row.items():
                    try:
                        val = int(v)
                    except ValueError:
                        try:
                            val = round(float(v), 6)
                        except ValueError:
                            raise ValueError(
                                "sweep values must be numbers (int or float)"
                            )
                    defn[k] = val
                defs.append(defn)
            n += 1

    return defs


def consolidate_files(outdir):
    """
    After running multiple processes, consolidate all of the different processes results
    into the root result directory (outdir)
    """
    for item in os.listdir(outdir):
        subdir = os.path.join(outdir, item)
        if os.path.isdir(subdir) and item not in ("network", "pop"):
            for report in os.listdir(subdir):
                # network folder
                if report == "network":
                    for file in os.listdir(os.path.join(subdir, report)):
                        shutil.move(
                            os.path.join(subdir, report, file),
                            os.path.join(outdir, "network"),
                        )
                elif report == "pop":
                    for file in os.listdir(os.path.join(subdir, report)):
                        shutil.move(
                            os.path.join(subdir, report, file),
                            os.path.join(outdir, "pop"),
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


def update_sweep_file(run_id, pop_id, defn, outdir):
    """
    Add this run to the sweep file json
    """
    f = open(os.path.join(outdir, "SweepVals.json"), "a")
    res = copy(defn)
    res["run_id"] = run_id
    res["pop_id"] = pop_id
    f.write(json.dumps(res))
    f.write("\n")
    f.close()


def single_run(sweep, outfile_dir, params, save_pop, pop_path):
    """
    A single run of titan.  Dispatched from main using parallel processes.
    """
    utils.set_up_logging(params)

    pid = str(os.getpid())
    pid_outfile_dir = os.path.join(outfile_dir, pid)
    save_pop_dir = (
        os.path.join(pid_outfile_dir, "pop") if save_pop is not None else None
    )
    if not os.path.isdir(pid_outfile_dir):
        os.mkdir(pid_outfile_dir)
        os.mkdir(os.path.join(pid_outfile_dir, "network"))
        if save_pop:
            os.mkdir(save_pop_dir)
        else:
            save_pop_dir = None

    # apply params from sweep for this run
    for param, val in sweep.items():
        utils.override_param(params, param, val, delimiter=".")

    tic = time_mod.time()

    # runs simulations
    if pop_path is None:
        pop = Population(params)
    else:
        pop = pop_io.read(params, pop_path)

    if save_pop_dir is not None:
        pop_io.write(pop, save_pop_dir)
        logging.info(f"Population saved to: {save_pop_dir}")

    try:
        model = TITAN(params, pop=pop)
    except Exception as e:
        raise Exception(f"Model creation failed: {e}")

    try:
        model.run(pid_outfile_dir)
    except Exception as e:
        raise Exception(f"Model run failed for run {model.id}: {e}")

    update_sweep_file(model.id, model.pop.id, sweep, pid_outfile_dir)

    return time_mod.time() - tic


def setup_outdir(outdir, save_pop):
    """
    Set up the results folder - will delete any files already present.
    """
    # delete old results before overwriting with new results
    outfile_dir = os.path.join(os.getcwd(), outdir)
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)
    os.mkdir(os.path.join(outfile_dir, "network"))
    if save_pop:
        os.mkdir(os.path.join(outfile_dir, "pop"))

    return outfile_dir


def get_sweep_defs(sweepfile, rows, sweeps, num_reps, force):
    """
    Get the sweep definitions from file if available, then from the command line definitions, otherwise return an empty definition to ensure the model runs.
    """
    if sweepfile is not None:
        sweep_defs = setup_sweeps_file(sweepfile, rows)
    elif sweeps:
        sweep_defs = setup_sweeps(sweeps)
    else:
        sweep_defs = [{}]  # make sure runs once

    if len(sweep_defs) > 100 and not force:
        raise ValueError(
            "Sweeping more than 100 models. Set `-F` (force) flag if you really want to do this."
        )

    sweep_defs *= num_reps

    return sweep_defs


def main(
    setting: str,
    params_path: str,
    num_reps: int,
    outdir: str,
    sweeps: List[str],
    force: bool,
    sweepfile: Optional[str] = None,
    rows: Optional[str] = None,
    error_on_unused: bool = False,
    save_pop: bool = False,
    pop_path: Optional[str] = None,
):
    """
    Run TITAN!

    args:
        setting: setting name to use, matches a folder name in `settings/`
        params_path: path to params file or directory
        num_reps: number of time to repeat each sweep
        outdir: directory where results are to be saved
        sweeps: array of strings in param:start:stop:step format
        force: if true, will run even if combination of sweeps results in greater than 100 runs
        sweepfile: path to csv file of sweep definitions
        rows: which rows of the csv to load to create sweeps in start:stop format
        error_on_unused: error if there are parameters that are unused by the model
        save_pop: if true, will save the population to file after creation
        pop_path: path to a population to load instead of creating a new population for each run
    """
    outfile_dir = setup_outdir(outdir, save_pop)

    # generate params - if no setting, set to none
    setting = setting.lower()
    setting_parsed = None if setting == "custom" else setting

    params = create_params(
        setting_parsed,
        params_path,
        outfile_dir,
        error_on_unused=error_on_unused,
    )

    # set up sweeps
    sweep_defs = get_sweep_defs(sweepfile, rows, sweeps, num_reps, force)

    tic = time_mod.time()
    wct = []  # wall clock times

    with Pool(
        processes=NCORES, maxtasksperchild=1
    ) as pool:  # set max tasks/child to prevent processor drift
        results = [
            pool.apply_async(
                single_run, (sweep_def, outfile_dir, params, save_pop, pop_path)
            )
            for sweep_def in sweep_defs
        ]
        while True:
            if all([r.ready() for r in results]):
                break
            else:
                time_mod.sleep(1)

        for r in results:
            try:
                t = r.get()
                wct.append(t)
            except Exception:
                traceback.print_exc()

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


def script_init():
    args = parser.parse_args()
    rows = args.rows.strip() if args.rows is not None else None
    sweepfile = args.sweepfile.strip() if args.sweepfile is not None else None
    poppath = args.poppath.strip() if args.poppath is not None else None
    main(
        args.setting.strip(),
        args.params.strip(),
        args.nMC,
        args.outdir.strip(),
        args.sweep,
        args.force,
        sweepfile=sweepfile,
        rows=rows,
        error_on_unused=args.error,
        save_pop=args.savepop,
        pop_path=poppath,
    )


if __name__ == "__main__":
    script_init()
