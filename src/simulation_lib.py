#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np

from . import params
from .ABM_core import HIVModel


def simulation(
    nreps,
    save_adjlist_flag,
    time_range,
    N_pop,
    outfile_dir,
    parameters,
    runSeed,
    popSeed,
    netSeed,
    uniqueSeed=False,
    model=None,
):
    # Check input
    if save_adjlist_flag not in [0, 1]:
        raise ValueError(
            "Invalid input! save_adjlist_flag = %s" % str(save_adjlist_flag)
        )

    # Run nreps simulations using the given parameters.
    # Information are printed to outfile_dir directory.
    pid = os.getpid()
    result_dict = {}

    for num_sim in range(nreps):
        inputSeed = runSeed
        if runSeed == -1:
            inputSeed = num_sim + 1
            print(inputSeed)

        print(
            "\tProcess %5s runs simulation %d/%d\t.:.\tInput rSeed: %d, pSeed: %d, nSeed: %d"
            % (pid, num_sim + 1, nreps, inputSeed, popSeed, netSeed)
        )

        MyModel = HIVModel(
            N=N_pop,
            tmax=time_range,
            parameter_dict=parameters,
            runseed=inputSeed,
            popseed=popSeed,
            netseed=netSeed,
            runtime_diffseed=uniqueSeed,
            model=model,
            network_type=params.network_type,
        )

        if save_adjlist_flag == 1 and num_sim == 0:
            MyModel.run(save_adjlist_flag=1, dir_prefix=outfile_dir)
        else:
            MyModel.run(save_adjlist_flag=0, dir_prefix=outfile_dir)

        # print MyModel
        result_dict_tmp = MyModel.return_results()
        for (key, x_v) in result_dict_tmp.items():
            if key not in result_dict:
                result_dict.update({key: {}})
            for t, x in x_v.items():
                if t not in result_dict[key]:
                    result_dict[key].update({t: []})
                if np.isnan(x):
                    # mssg = '\tWARNING:\
                    # NaN encountered in pid %5s run %d for %s at time %d'
                    # print mssg%(pid, num_sim,key, t)
                    result_dict[key][t].append(0)
                else:
                    result_dict[key][t].append(x)
    return result_dict


def simulation_star(zipped_input):
    """
    This helper function is needed because multithreading pool only
    accepts one input argument. This helper function converts combined
    arguments and calls the simulation function.
    """
    try:
        return simulation(*zipped_input)
    except TypeError:
        for input_info in zipped_input:
            print(input_info)
        raise TypeError("Wrong input for simulation_star()!")


def save_results(N_MC, time_range, rslts, outfile_dir, num_sim):
    """
    Save the result dictionary.
    Input:   N_MC : Number of Monte Carlo runs for each simulation.
             time_range : Number of time steps in each simulation.
             rslts : Result dictionary.
             outfile_dir : Directory for output file.
             num_sim : Number of simulation. IMPORTANT to identify
                       the applied parameters. Must be consistent with
                       parameter input file columns.
    """
    # convert list of dicts into one nested result dict
    result_dict = {}

    # save results
    if not os.path.isdir(outfile_dir):
        os.mkdir(outfile_dir)
    OutFileName = os.path.join(outfile_dir, ("Result_simulation_%d.txt" % num_sim))
    print("\n\tSaving results to:\n\t%s\n" % str(OutFileName))
    if os.path.isfile(OutFileName):
        os.remove(OutFileName)
    outfile = open(OutFileName, "w")
    outfile.write("t,model,coverage")
    for result_property in sorted(rslts):
        outfile.write(
            ",%s_mean,%s_std,%s_10th,%s_90th"
            % (result_property, result_property, result_property, result_property)
        )
    outfile.write("\n")

    for t in range(0, time_range + 1):
        outfile.write("%d,%s,%0.2f" % (t, params.PrEP_type, params.PrEP_Target))
        for result_property in sorted(rslts):
            x_v = []

            try:
                x_v = np.array(rslts[result_property][t])
                x_v = x_v[np.logical_not(np.isnan(x_v))]
            except:
                pass

            if len(x_v) > 0:
                mean = np.mean(x_v)
                outfile.write(",%4.5f" % mean)
            else:
                outfile.write(",%4.5f" % np.NaN)

            # print std
            if len(x_v) > 0:
                std_dev = np.std(x_v)
                outfile.write(",%4.5f" % std_dev)
            else:
                outfile.write(",%4.5f" % np.NaN)

            # print 10th
            if len(x_v) > 0:
                p10 = np.percentile(x_v, 10)
            else:
                p10 = np.NaN
            outfile.write(",%4.5f" % p10)

            # print 90th
            if len(x_v) > 0:
                p90 = np.percentile(x_v, 90)
            else:
                p90 = np.NaN
            outfile.write(",%4.5f" % p90)

        outfile.write("\n")

    outfile.close()


if __name__ == "__main__":
    print("Use run_titan!\n")
