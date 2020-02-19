import oyaml as yaml
import sys
import csv
from dotmap import DotMap

param_file = sys.argv[1]  # format = titan_new.<name>
print(param_file)
param = __import__(param_file, globals(), locals(), ["params"], 0)

params = DotMap(param.params.toDict())  # make sure it's dotmaps all the way down

parts = param_file.split(".")
base_name = parts[2]
dir_name = parts[1]
root_name = parts[0]
yaml_name = base_name + ".yml"

with open("post_process_map.txt", "r") as fmap:
    freader = csv.DictReader(fmap)
    for row in freader:
        if row["old"] != "":
            for race in ["WHITE", "BLACK"]:
                for st in params.classes.sex_types:
                    if row["old"] in params.demographics[race][st]:
                        val = params.demographics[race][st][row["old"]]
                        expr = f"params.demographics[race][st].{row['dot']} = val"
                        params.demographics[race][st].__delattr__(row["old"])
                        exec(expr, globals(), locals())

                # and also PWID
                if row["old"] in params.demographics[race]["PWID"]:
                    val = params.demographics[race]["PWID"][row["old"]]
                    expr = f"params.demographics[race]['PWID'].{row['dot']} = val"
                    params.demographics[race]["PWID"].__delattr__(row["old"])
                    exec(expr, globals(), locals())


with open(f"{root_name}/{dir_name}/{yaml_name}", "w") as f:
    yaml.dump(params.toDict(), f)
