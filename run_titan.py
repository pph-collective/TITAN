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
import csv

from titan.model import HIVModel
from titan.population import Population
import titan.population_io as pop_io
from titan.parse_params import create_params

# how many cores can we use
NCORES = os.environ.get("SLURM_CPUS_PER_TASK", cpu_count())
NCORES = int(NCORES)  # environment variable returns string

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
    "-S", "--setting", default="custom", help="setting directory to use"
)
parser.add_argument(
    "-p", "--params", required=True, help="directory or file with params yaml(s)"
)
parser.add_argument(
    "-o", "--outdir", default="results", help="directory name to save results to",
)
parser.add_argument(
    "-b", "--base", type=bool, default=True, help="whether to use base setting",
)

parser.add_argument(
    "-e",
    "--error",
    action="store_true",
    help="Error on unused parameters instead of warning",
)

parser.add_argument(
    "--savepop",
    type=str,
    default=None,
    help="Save population after creation, but before model run. 'all' = save all atributes, 'core' = save core (non-intervention) attributes.",
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

    return {"param": parts[0], "start": start, "stop": stop, "step": step}


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
    help="Run model even if number of sweeps excedes 100",
)


def drange(start, stop, step):
    """
    Like range, but allows decimal steps
    """
    r = start
    while r < stop:
        yield r
        r += step
        r = round(r, 3)


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
                            val = float(v)
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
    pid = str(os.getpid())
    pid_outfile_dir = os.path.join(outfile_dir, pid)
    save_pop_dir = (
        os.path.join(pid_outfile_dir, "pop") if save_pop is not None else None
    )
    if not os.path.isdir(pid_outfile_dir):
        os.mkdir(pid_outfile_dir)
        os.mkdir(os.path.join(pid_outfile_dir, "network"))
        if save_pop in ("all", "core"):
            os.mkdir(save_pop_dir)
        else:
            save_pop_dir = None

    # apply params from sweep for this run
    for param, val in sweep.items():
        path = param.split(".")
        sweep_item = params
        for p in path[:-1]:
            sweep_item = sweep_item[p]
        sweep_item[path[-1]] = val

    tic = time_mod.time()

    # runs simulations
    if pop_path is None:
        pop = Population(params)
    else:
        pop = pop_io.read(params, pop_path)

    if save_pop_dir is not None:
        intervention_attrs = True if save_pop == "all" else False
        pop_io.write(pop, save_pop_dir, intervention_attrs=intervention_attrs)
        print(save_pop_dir)

    try:
        model = HIVModel(params, population=pop)
    except Exception as e:
        raise Exception(f"Model creation failed: {e}")

    try:
        model.run(pid_outfile_dir)
    except Exception as e:
        raise Exception(f"Model run failed for run {model.id}: {e}")

    update_sweep_file(model.id, model.pop.id, sweep, pid_outfile_dir)

    return time_mod.time() - tic


def main(
    setting,
    params_path,
    num_reps,
    outdir,
    use_base,
    sweeps,
    force,
    sweepfile=None,
    rows=None,
    error_on_unused=False,
    save_pop=None,
    pop_path=None,
):
    """
    Run TITAN!

    :Input:
        setting : str
            setting name to use, matches a folder name in `settings/`
        params_path : str
            path to params file or directory
        num_reps : int
            number of time to repeat each sweep
        outdir : str
            directory where results are to be saved
        use_base : boolean
            whether to use the "base" setting (includes some more complicated defaults)
        sweeps : array[str]
            array of strings in param:start:stop:step format
        force : boolean
            if true, will run even if combination of sweeps results in greater than 100 runs
        sweepfile : str [default: None]
            path to csv file of sweep definitions
        rows: str [default: None]
            which rows of the csv to load to create sweeps in start:stop format
        error_on_unused : boolean [default: False]
            error if there are parameters that are unused by the model
        save_pop : str [default: None]
            'all' or 'core' will save the population with agents have all or core attributes saved
        pop_path : str [default: None]
            path to a population to load instead of creating a new population for each run
    """
    # delete old results before overwriting with new results
    outfile_dir = os.path.join(os.getcwd(), outdir)
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)
    os.mkdir(os.path.join(outfile_dir, "network"))
    if save_pop in ("all", "core"):
        os.mkdir(os.path.join(outfile_dir, "pop"))

    # generate params - if no setting, set to none
    setting = setting.lower()
    if setting == "custom":
        setting = None
    else:
        setting = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "settings", setting
        )
        assert os.path.isdir(setting), f"{setting} is not a directory"

    params = create_params(
        setting,
        params_path,
        outfile_dir,
        use_base=use_base,
        error_on_unused=error_on_unused,
    )

    # set up sweeps
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

    tic = time_mod.time()
    wct = []  # wall clock times

    with Pool(processes=NCORES) as pool:
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
                print([r.ready() for r in results])
                time_mod.sleep(1)

        for r in results:
            try:
                t = r.get()
                wct.append(t)
            except Exception as e:
                print(e)

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
    rows = args.rows.strip() if args.rows is not None else None
    sweepfile = args.sweepfile.strip() if args.sweepfile is not None else None
    savepop = args.savepop.strip() if args.savepop is not None else None
    poppath = args.poppath.strip() if args.poppath is not None else None
    main(
        args.setting.strip(),
        args.params.strip(),
        args.nMC,
        args.outdir.strip(),
        args.base,
        args.sweep,
        args.force,
        sweepfile=sweepfile,
        rows=rows,
        error_on_unused=args.error,
        save_pop=savepop,
        pop_path=poppath,
    )
