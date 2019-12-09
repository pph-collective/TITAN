#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np  # type: ignore
from typing import Sequence, List, Dict, Optional, Any

from . import params  # type: ignore
from .ABM_core import HIVModel


def initiate_ResultDict(tmax: int) -> Dict[str, Any]:
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


def statsToResults(stats: Dict[str, Any], results: Dict[str, Any]):
    """
    Update results dict with stats from this simulation
    """
    for t in stats.keys():
        results["Prv_HIV"][t].append(
            1.0 * stats[t]["ALL"]["ALL"]["numHIV"] / stats[t]["ALL"]["ALL"]["numAgents"]
        )

        results["Prv_AIDS"][t].append(
            1.0 * stats[t]["ALL"]["ALL"]["numAIDS"] / stats[t]["ALL"]["ALL"]["numHIV"]
        )

        results["Prv_Test"][t].append(
            1.0
            * stats[t]["ALL"]["ALL"]["numTested"]
            / max(stats[t]["ALL"]["ALL"]["numHIV"], 1)
        )

        results["Prv_ART"][t].append(
            1.0 * stats[t]["ALL"]["ALL"]["numART"] / stats[t]["ALL"]["ALL"]["numTested"]
        )

        results["Prv_PrEP"][t].append(
            1.0
            * stats[t]["ALL"]["ALL"]["numPrEP"]
            / stats[t]["ALL"]["ALL"]["numAgents"]
        )

        results["n_Relations"][t].append(stats[t]["ALL"]["ALL"]["numRels"])

        results["Inc_t_HM"][t].append(stats[t]["WHITE"]["HM"]["inf_newInf"])
        results["Inc_t_HF"][t].append(stats[t]["WHITE"]["HF"]["inf_newInf"])


def simulation(
    nreps: int,
    time_range: int,
    N_pop: int,
    outfile_dir: str,
    runSeed: int,
    popSeed: int,
    netSeed: int,
):

    # Run nreps simulations using the given parameters.
    # Information are printed to outfile_dir directory.
    pid = os.getpid()
    result_dict: Dict[str, Any] = initiate_ResultDict(time_range)

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
        statsToResults(stats, result_dict)

    return result_dict


def save_results(time_range, rslts, outfile_dir, num_sim):
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
        outfile.write("%d,%s,%0.2f" % (t, params.PrEP_type, params.PrEP_Target))

        # for each property in the result dict, write the mean, std dev, 10th % and 90th % over the mante carlo iterations (params.N_MC) in the simulation
        for result_property in sorted(rslts):
            x_v = np.array(rslts[result_property][t])

            outfile.write(",%4.5f" % np.mean(x_v))
            outfile.write(",%4.5f" % np.std(x_v))
            outfile.write(",%4.5f" % np.percentile(x_v, 10))
            outfile.write(",%4.5f" % np.percentile(x_v, 90))

        outfile.write("\n")

    outfile.close()


if __name__ == "__main__":
    print("Use run_titan!\n")
