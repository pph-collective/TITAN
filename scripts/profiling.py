import cProfile, pstats, io
from sys import path
import os

path_to_this_dir = os.path.dirname(os.path.realpath(__file__))
path.append(os.path.join(path_to_this_dir, ".."))

from titan import run_titan
from titan.parse_params import create_params

outdir = "results"
save_pop = False
setting = "nyc-msm"
params_path = "tests/params/integration_base.yml"
sweepfile = None
rows = None
sweeps = []
num_reps = 1
force = True
pop_path = None
error_on_unused = False

pr = cProfile.Profile()
pr.enable()

# running the normal way only profiles the main process, get a single run going
outfile_dir = run_titan.setup_outdir(outdir, save_pop)

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
sweep_defs = run_titan.get_sweep_defs(sweepfile, rows, sweeps, num_reps, force)
run_titan.single_run(sweep_defs[0], outfile_dir, params, save_pop, pop_path)

pr.disable()
s = io.StringIO()
sortby = "tottime"
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
pr.dump_stats("profile.prof")

print(s.getvalue())

# run `snakeviz profile.prof` to visualize results
