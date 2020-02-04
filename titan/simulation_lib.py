#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np # type: ignore
from typing import Sequence, Dict, Any

from titan import params
from .ABM_core import HIVModel


def initiate_result_dict(tmax: int) -> Dict[str, Dict[int, Sequence]]:
    # nested dictionary for results (inner dictionary has the form: time:result)
    d: Dict[str, Any] = {
        "Prv_HIV": {},
        "Prv_AIDS": {},
        "Prv_Test": {},
        "Prv_ART": {},
        "Prv_PrEP": {},
        "n_Relations": {},
        "Inc_t_HM": {},
        "Inc_t_HF": {},
    }

    for key in d:
        for t in range(1, tmax + 1):
            d[key].update({t: []})

    return d


def safe_divide(numerator: int, denominator: int):
    """
    Default 0 if denominator is 0, otherwise divide as normal
    """
    if denominator == 0:
        return 0.0
    else:
        return 1.0 * numerator / denominator


def stats_to_results(stats: Dict[str, Any], results: Dict[str, Any]):
    """
    Update results dict with stats from this simulation
    """
    for t in stats.keys():
        stat = stats[t]["ALL"]["ALL"]
        results["Prv_HIV"][t].append(safe_divide(stat["numHIV"], stat["numAgents"]))
        results["Prv_AIDS"][t].append(safe_divide(stat["numAIDS"], stat["numHIV"]))
        results["Prv_Test"][t].append(safe_divide(stat["numTested"], stat["numHIV"]))
        results["Prv_ART"][t].append(safe_divide(stat["numART"], stat["numTested"]))
        results["Prv_PrEP"][t].append(safe_divide(stat["numPrEP"], stat["numAgents"]))

        results["n_Relations"][t].append(stat["numRels"])
        if "HM" in params.agentSexTypes:
            results["Inc_t_HM"][t].append(stats[t]["WHITE"]["HM"]["inf_newInf"])
        if "HF" in params.agentSexTypes:
            results["Inc_t_HF"][t].append(stats[t]["WHITE"]["HF"]["inf_newInf"])


def simulation(
    nreps: int, time_range: int, N_pop: int, runSeed: int, popSeed: int, netSeed: int,
):

    # Run nreps simulations using the given parameters.
    pid = os.getpid()
    result_dict: Dict[str, Any] = initiate_result_dict(time_range)

    for rep in range(nreps):
        inputSeed = runSeed
        if runSeed == -1:
            inputSeed = rep + 1
            print(inputSeed)

        print(
            "\tProcess %5s runs simulation %d/%d\t.:.\tInput rSeed: %d, pSeed: %d, nSeed: %d"
            % (pid, rep + 1, nreps, inputSeed, popSeed, netSeed)
        )

        MyModel = HIVModel(
            N=N_pop,
            tmax=time_range,
            runseed=inputSeed,
            popseed=popSeed,
            netseed=netSeed,
            network_type=params.network_type,
        )

        stats = MyModel.run()
        stats_to_results(stats, result_dict)

    return result_dict


# REVIEWED if time_range is number of timesteps, why not get this from the rslts? - get it from the dict
def save_results(
    time_range: int,
    rslts: Dict[str, Dict[int, Sequence]],
    outfile_dir: str,
    num_sim: int,
):
    """
    Save the result dictionary.
    Input:
             time_range : Number of time steps in each simulation.
             rslts : Result dictionary.
             outfile_dir : Directory for output file.
             num_sim : Number of simulation. IMPORTANT to identify
                       the applied parameters. Must be consistent with
                       parameter input file columns.
    """
    OutFileName = os.path.join(outfile_dir, ("Result_simulation_%d.txt" % num_sim))
    print("\n\tSaving results to:\n\t%s\n" % str(OutFileName))

    outfile = open(OutFileName, "w")
    outfile.write("t,model,coverage")

    for result_property in sorted(rslts):
        outfile.write(
            ",%s_mean,%s_std,%s_10th,%s_90th"
            % (result_property, result_property, result_property, result_property)
        )
    outfile.write("\n")

    for t in range(1, time_range + 1):

        # write timestep, model, coverage
        prep_type = "+".join(params.PrEP_type)

        outfile.write("%d,%s,%0.2f" % (t, prep_type, params.PrEP_Target))

        # for each property in the result dict, write the mean, std dev, 10th % and 90th % over the mante carlo
        # iterations (params.N_MC) in the simulation
        for result_property in sorted(rslts):
            x_v = np.array(rslts[result_property][t])

            if len(x_v) != 0:
                outfile.write(",%4.5f" % np.mean(x_v))
                outfile.write(",%4.5f" % np.std(x_v))
                outfile.write(",%4.5f" % np.percentile(x_v, 10))
                outfile.write(",%4.5f" % np.percentile(x_v, 90))
            else:
                for i in range(4):
                    outfile.write(",%4.5f" % 0)
        outfile.write("\n")

    outfile.close()


if __name__ == "__main__":
    print("Use run_titan!\n")
