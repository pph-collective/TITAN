#!/usr/bin/env python
# encoding: utf-8

from . import params  # type: ignore
from typing import Dict, Any, List, Sequence, Optional
from .agent import Agent_set, Relationship, Agent
from copy import deepcopy
from uuid import UUID

MAIN_CAT = list(params.DemographicParams.keys())
MAIN_CAT.append("ALL")
SUB_CAT = params.agentPopulations
SUB_CAT.append("ALL")


def get_stats(
    totalAgents: Agent_set,
    HIVAgents: Agent_set,
    IncarAgents: Agent_set,
    PrEPAgents: Agent_set,
    newPrEPAgents: Agent_set,
    NewInfections: Agent_set,
    NewDiagnosis: Agent_set,
    Relationships: List[Relationship],
    newHR: Agent_set,
    newIncarRelease: Agent_set,
    deathSet: List[Agent],
    LAIagents: Agent_set,
    oralAgents: Agent_set,
):

    stats_template = {
        "numAgents": 0,
        "inf_HR6m": 0,
        "inf_HRever": 0,
        "inf_newInf": 0,
        "newHR": 0,
        "newHR_HIV": 0,
        "newHR_AIDS": 0,
        "newHR_tested": 0,
        "newHR_ART": 0,
        "newRelease": 0,
        "newReleaseHIV": 0,
        "numHIV": 0,
        "numTested": 0,
        "numAIDS": 0,
        "numART": 0,
        "numHR": 0,
        "newlyTested": 0,
        "deaths": 0,
        "deaths_HIV": 0,
        "incar": 0,
        "incarHIV": 0,
        "numPrEP": 0,
        "newNumPrEP": 0,
        "iduPartPrep": 0,
        "msmwPartPrep": 0,
        "testedPartPrep": 0,
        "numVaccinated": 0,
        "LAIagents": 0,
        "oralAgents": 0,
        "awareAgents": 0,
    }

    stats = {}
    for cat in MAIN_CAT:
        stats[cat] = {sc: dict(stats_template) for sc in SUB_CAT}

    # Incarceration metrics
    for tmpA in IncarAgents.iter_agents():
        stats[tmpA._race][tmpA._SO]["incar"] += 1
        if tmpA._HIV_bool:
            stats[tmpA._race][tmpA._SO]["incarHIV"] += 1

    for tmpA in newIncarRelease.iter_agents():
        stats[tmpA._race][tmpA._SO]["newRelease"] += 1
        if tmpA._HIV_bool:
            stats[tmpA._race][tmpA._SO]["newReleaseHIV"] += 1

    # Newly infected tracker statistics (with HR within 6mo and HR ever bool check)
    for tmpA in NewInfections.iter_agents():
        stats[tmpA._race][tmpA._SO]["inf_newInf"] += 1
        if tmpA._everhighrisk_bool:
            stats[tmpA._race][tmpA._SO]["inf_HRever"] += 1
        if tmpA._highrisk_bool:
            stats[tmpA._race][tmpA._SO]["inf_HR6m"] += 1

    # PrEP reason tracker
    for tmpA in totalAgents.iter_agents():
        if tmpA._PrEP_bool:
            stats[tmpA._race][tmpA._SO]["numPrEP"] += 1
            if "IDU" in tmpA._PrEP_reason:
                stats[tmpA._race][tmpA._SO]["iduPartPrep"] += 1
            if "MSMW" in tmpA._PrEP_reason:
                stats[tmpA._race][tmpA._SO]["msmwPartPrep"] += 1
            if "HIV test" in tmpA._PrEP_reason:
                stats[tmpA._race][tmpA._SO]["testedPartPrep"] += 1
        if tmpA.awareness:
            stats[tmpA._race][tmpA._SO]["awareAgents"] += 1

    # Newly PrEP tracker statistics
    for tmpA in newPrEPAgents.iter_agents():
        stats[tmpA._race][tmpA._SO]["newNumPrEP"] += 1

    # Newly diagnosed tracker statistics
    for tmpA in NewDiagnosis.iter_agents():
        stats[tmpA._race][tmpA._SO]["newlyTested"] += 1

    # Newly HR agents
    for tmpA in newHR.iter_agents():
        stats[tmpA._race][tmpA._SO]["newHR"] += 1
        if tmpA._HIV_bool:
            stats[tmpA._race][tmpA._SO]["newHR_HIV"] += 1
            if tmpA._AIDS_bool:
                stats[tmpA._race][tmpA._SO]["newHR_AIDS"] += 1
            if tmpA._tested:
                stats[tmpA._race][tmpA._SO]["newHR_tested"] += 1
                if tmpA._HAART_bool:
                    stats[tmpA._race][tmpA._SO]["newHR_ART"] += 1

    # Total HIV summary snapshot for timestep
    for tmpA in HIVAgents.iter_agents():
        stats[tmpA._race][tmpA._SO]["numHIV"] += 1
        if tmpA._AIDS_bool:
            stats[tmpA._race][tmpA._SO]["numAIDS"] += 1
        if tmpA._tested:
            stats[tmpA._race][tmpA._SO]["numTested"] += 1
        if tmpA._HAART_bool:
            stats[tmpA._race][tmpA._SO]["numART"] += 1

    for tmpA in totalAgents._subset["DU"]._subset["IDU"].iter_agents():
        stats[tmpA._race]["IDU"]["numAgents"] += 1
        if tmpA._HIV_bool:
            stats[tmpA._race]["IDU"]["numHIV"] += 1
        if tmpA._AIDS_bool:
            stats[tmpA._race]["IDU"]["numAIDS"] += 1
        if tmpA._tested:
            stats[tmpA._race]["IDU"]["numTested"] += 1
        if tmpA._HAART_bool:
            stats[tmpA._race]["IDU"]["numART"] += 1

    for tmpA in totalAgents.iter_agents():
        stats[tmpA._race][tmpA._SO]["numAgents"] += 1
        if tmpA.vaccine_bool:
            stats[tmpA._race][tmpA._SO]["numVaccinated"] += 1

    for tmpA in deathSet:
        stats[tmpA._race][tmpA._SO]["deaths"] += 1
        if tmpA._HIV_bool:
            stats[tmpA._race][tmpA._SO]["deaths_HIV"] += 1

    for (
        tmpA
    ) in (
        LAIagents.iter_agents()
    ):  # TODO this should have a way of doing a for loop of PrEP types
        # TODO: use num_members and subsets? Or other non-loop method
        stats[tmpA._race][tmpA._SO]["LAIagents"] += 1
    for tmpA in oralAgents.iter_agents():
        stats[tmpA._race][tmpA._SO]["oralAgents"] += 1

    # Sum 'ALL' categories for race/SO bins

    for race in stats:
        for param in stats_template:
            for sc in SUB_CAT:
                if sc in params.agentSexTypes:
                    stats[race]["ALL"][param] += stats[race][sc][param]
            for sc in SUB_CAT:
                stats["ALL"][sc][param] += stats[race][sc][param]
    # add relationship count (only makes sense at the all level)
    stats["ALL"]["ALL"]["numRels"] = len(Relationships)

    return stats


# ================== Printer Functions =========================
# Each of the following functions takes in the time, seeds, and stats dict for that time
# and prints the appropriate stats to file


def deathReport(
    run_id: UUID,
    t: int,
    runseed: int,
    popseed: int,
    netseed: int,
    stats: Dict[str, Any],
):
    f = open("results/DeathReport.txt", "a")
    sex_types = deepcopy(params.agentSexTypes)
    sex_types.append("ALL")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\ttot_{st}\tHIV_{st}"
        for sex_type in sex_types:
            f.write(template.format(st=sex_type))

        f.write("\n")

    f.write("%d\t%d\t%d" % (run_id, runseed, t))  # start row

    for sex_type in sex_types:
        f.write(
            "\t%d\t%d"
            % (stats["ALL"][sex_type]["deaths"], stats["ALL"][sex_type]["deaths_HIV"])
        )

    f.write("\n")
    f.close()


def incarReport(
    run_id: UUID,
    t: int,
    runseed: int,
    popseed: int,
    netseed: int,
    stats: Dict[str, Any],
):
    f = open("results/IncarReport.txt", "a")

    name_map = {
        "incar": "tot",
        "incarHIV": "HIV",
        "newRelease": "rlsd",
        "newReleaseHIV": "rlsdHIV",
    }

    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\t{mc}_{st}_{p}"
        for p in name_map.values():
            for mc in MAIN_CAT:
                for sex_type in params.agentSexTypes:
                    f.write(template.format(mc=mc, st=sex_type, p=p))

        f.write("\n")

    f.write("%d\t%d\t%d" % (run_id, runseed, t))

    for p in name_map:
        for mc in MAIN_CAT:
            for st in params.agentSexTypes:
                f.write("\t")
                f.write(str(stats[mc][st][p]))

    f.write("\n")
    f.close()


def newlyhighriskReport(
    run_id: UUID,
    t: int,
    runseed: int,
    popseed: int,
    netseed: int,
    stats: Dict[str, Any],
):
    f = open("results/newlyHR_Report.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\tnewHR_{st}\tnewHR_HIV_{st}\tnewHR_AIDS_{st}\tnewHR_Tested_{st}\tnewHR_ART_{st}"
        for sex_type in params.agentSexTypes:
            f.write(template.format(st=sex_type))

        f.write("\n")

    f.write("%d\t%d\t%d" % (run_id, runseed, t))  # start row

    for sex_type in params.agentSexTypes:
        f.write(
            "\t%d\t%d\t%d\t%d\t%d"
            % (
                stats["ALL"][sex_type]["newHR"],
                stats["ALL"][sex_type]["newHR_HIV"],
                stats["ALL"][sex_type]["newHR_AIDS"],
                stats["ALL"][sex_type]["newHR_tested"],
                stats["ALL"][sex_type]["newHR_ART"],
            )
        )
        f.write("\n")
    f.close()


def prepReport(
    run_id: UUID,
    t: int,
    runseed: int,
    popseed: int,
    netseed: int,
    stats: Dict[str, Any],
):
    f = open("results/PrEPReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt\tNewEnroll\tIDUpartner\tTestedPartner\tMSMWpartner\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            run_id,
            runseed,
            t,
            stats["ALL"]["ALL"]["newNumPrEP"],
            stats["ALL"]["ALL"]["iduPartPrep"],
            stats["ALL"]["ALL"]["testedPartPrep"],
            stats["ALL"]["ALL"]["msmwPartPrep"],
        )
    )
    f.close()


def basicReport(
    run_id: UUID,
    t: int,
    runseed: int,
    popseed: int,
    netseed: int,
    stats: Dict[str, Any],
):
    for agentRace in MAIN_CAT:
        for agentTypes in SUB_CAT:
            name = "basicReport_" + agentTypes + "_" + agentRace
            tmpReport = open("results/" + name + ".txt", "a")

            # if this is a new file, write the header info
            if tmpReport.tell() == 0:
                tmpReport.write(
                    "run_id\trseed\tpseed\tnseed\tt\tTotal\tHIV\tAIDS\tTstd\tART\tnHR\tIncid\tHR_6mo\tHR_Ev\tNewDiag\tDeaths\tPrEP\tIDUpart_PrEP\tMSMWpart_PrEP\ttestedPart_PrEP\tVaccinated\tLAI\tOral\tAware\n"
                )

            tmpReport.write(
                (
                    "{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                        str(run_id),
                        runseed,
                        popseed,
                        netseed,
                        t,
                        stats[agentRace][agentTypes]["numAgents"],
                        stats[agentRace][agentTypes]["numHIV"],
                        stats[agentRace][agentTypes]["numAIDS"],
                        stats[agentRace][agentTypes]["numTested"],
                        stats[agentRace][agentTypes]["numART"],
                        stats[agentRace][agentTypes]["numHR"],
                        stats[agentRace][agentTypes]["inf_newInf"],
                        stats[agentRace][agentTypes]["inf_HR6m"],
                        stats[agentRace][agentTypes]["inf_HRever"],
                        stats[agentRace][agentTypes]["newlyTested"],
                        stats[agentRace][agentTypes]["deaths"],
                        stats[agentRace][agentTypes]["numPrEP"],
                        stats[agentRace][agentTypes]["iduPartPrep"],
                        stats[agentRace][agentTypes]["msmwPartPrep"],
                        stats[agentRace][agentTypes]["testedPartPrep"],
                        stats[agentRace][agentTypes]["numVaccinated"],
                        stats[agentRace][agentTypes]["LAIagents"],
                        stats[agentRace][agentTypes]["oralAgents"],
                        stats[agentRace][agentTypes]["awareAgents"]
                        # np.mean(stats[agentRace][agentTypes]["agentOpinions"])
                    )
                )
            )
            tmpReport.close()


# ========================== Other Print Functions =============================


def print_components(
    run_id: UUID, t: int, runseed: int, popseed: int, netseed: int, components
):
    """
    Write stats describing the components (sub-graphs) in a graph to file
    """
    f = open("results/componentReport_ALL.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "run_id\trunseed\tpopseed\tnetseed\tt\tcompID\ttotalN\tNhiv\tNprep\tNtrtHIV\tNprepHIV\tTrtComponent\tOral\tLAI\tAware\n"
        )

    compID = 0
    for comp in components:
        totN = nhiv = ntrthiv = nprep = PrEP_ever_HIV = trtbool = LAI = oral = aware = 0
        for agent in comp.nodes():
            totN += 1
            if agent._HIV_bool:
                nhiv += 1
                if agent._treatment_bool:
                    ntrthiv += 1
                if agent._PrEP_ever_bool:
                    PrEP_ever_HIV += 1
            elif agent._treatment_bool:
                if agent._PrEP_bool:
                    nprep += 1
                    if agent.PrEP_type == "LAI":
                        LAI += 1
                    elif agent.PrEP_type == "Oral":
                        oral += 1
            trtbool += agent._PCA
            if agent.awareness:
                aware += 1

        f.write(
            "{run_id}\t{runseed}\t{pseed}\t{nseed}\t{t}\t{compID}\t{totalN}\t{Nhiv}\t{Nprep}\t{NtrtHIV}"
            "\t{NprepHIV}\t{trtbool}\t{Oral}\t{LAI}\t{aware}\n".format(
                run_id=run_id,
                runseed=runseed,
                pseed=popseed,
                nseed=netseed,
                t=t,
                compID=compID,
                totalN=totN,
                Nhiv=nhiv,
                Nprep=nprep,
                NtrtHIV=ntrthiv,
                NprepHIV=PrEP_ever_HIV,
                trtbool=trtbool,
                Oral=oral,
                LAI=LAI,
                aware=aware,
            )
        )

        compID += 1
    f.close()
