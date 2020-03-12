import cProfile, pstats, io
from pstats import SortKey
import run_titan
pr = cProfile.Profile()
pr.enable()
run_titan.main("custom", "tests/params/basic_seeded.yml", 1, "results", True)
pr.disable()
s = io.StringIO()
sortby = 'tottime'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())


from titan import ObjMap
import oyaml as yaml
param_path = "tests/params/basic_seeded.yml"
with open(param_path) as f:
    yml = yaml.safe_load(f)

a = ObjMap.ObjMap(yml)
