import cProfile, pstats, io
import run_titan

pr = cProfile.Profile()
pr.enable()
run_titan.main("custom", "tests/params/basic_seeded.yml", 1, "results", True, [], False)

pr.disable()
s = io.StringIO()
sortby = "tottime"
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
pr.dump_stats("profile.prof")

print(s.getvalue())
