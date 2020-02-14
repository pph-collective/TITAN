#!/usr/bin/env python
# encoding: utf-8

from typing import Dict, Any, List, Sequence, Optional
from .agent import AgentSet, Relationship, Agent
from copy import deepcopy
from uuid import UUID
import os

from dotmap import DotMap  # type: ignore


def get_stats(
    totalAgents: AgentSet,
    HIVAgents: AgentSet,
    IncarAgents: AgentSet,
    PrEPAgents: AgentSet,
    newPrEPAgents: AgentSet,
    NewInfections: AgentSet,
    NewDiagnosis: AgentSet,
    Relationships: List[Relationship],
    newHR: AgentSet,
    newIncarRelease: AgentSet,
    deathSet: List[Agent],
    params: DotMap,
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
        "Vaccinated": 0,
        "injectable_prep": 0,
        "oral_prep": 0,
        "prep_aware": 0,
    }

    MAIN_CAT = deepcopy(params.classes.races)
    MAIN_CAT.append("ALL")
    SUB_CAT = deepcopy(params.classes.populations)
    SUB_CAT.append("ALL")

    stats = {}
    for cat in MAIN_CAT:
        stats[cat] = {sc: dict(stats_template) for sc in SUB_CAT}

    # Incarceration metrics
    for tmpA in IncarAgents:
        stats[tmpA.race][tmpA.so]["incar"] += 1
        if tmpA.hiv:
            stats[tmpA.race][tmpA.so]["incarHIV"] += 1

    for tmpA in newIncarRelease:
        stats[tmpA.race][tmpA.so]["newRelease"] += 1
        if tmpA.hiv:
            stats[tmpA.race][tmpA.so]["newReleaseHIV"] += 1

    # Newly infected tracker statistics (with HR within 6mo and HR ever bool check)
    for tmpA in NewInfections:
        stats[tmpA.race][tmpA.so]["inf_newInf"] += 1
        if tmpA.high_risk_ever:
            stats[tmpA.race][tmpA.so]["inf_HRever"] += 1
        if tmpA.high_risk:
            stats[tmpA.race][tmpA.so]["inf_HR6m"] += 1

    # PrEP reason tracker
    for tmpA in totalAgents:
        if tmpA.prep:
            stats[tmpA.race][tmpA.so]["numPrEP"] += 1
            if "PWID" in tmpA.prep_reason:
                stats[tmpA.race][tmpA.so]["iduPartPrep"] += 1
            if "MSMW" in tmpA.prep_reason:
                stats[tmpA.race][tmpA.so]["msmwPartPrep"] += 1
            if "HIV test" in tmpA.prep_reason:
                stats[tmpA.race][tmpA.so]["testedPartPrep"] += 1
            if tmpA.prep_type == "Inj":
                stats[tmpA.race][tmpA.so]["injectable_prep"] += 1
            elif tmpA.prep_type == "Oral":
                stats[tmpA.race][tmpA.so]["oral_prep"] += 1

    # Newly PrEP tracker statistics
    for tmpA in newPrEPAgents:
        stats[tmpA.race][tmpA.so]["newNumPrEP"] += 1

    # Newly diagnosed tracker statistics
    for tmpA in NewDiagnosis:
        stats[tmpA.race][tmpA.so]["newlyTested"] += 1

    # Newly HR agents
    for tmpA in newHR:
        stats[tmpA.race][tmpA.so]["newHR"] += 1
        if tmpA.hiv:
            stats[tmpA.race][tmpA.so]["newHR_HIV"] += 1
            if tmpA.aids:
                stats[tmpA.race][tmpA.so]["newHR_AIDS"] += 1
            if tmpA.hiv_dx:
                stats[tmpA.race][tmpA.so]["newHR_tested"] += 1
                if tmpA.haart:
                    stats[tmpA.race][tmpA.so]["newHR_ART"] += 1

    # Total HIV summary snapshot for timestep
    for tmpA in HIVAgents:
        stats[tmpA.race][tmpA.so]["numHIV"] += 1
        if tmpA.aids:
            stats[tmpA.race][tmpA.so]["numAIDS"] += 1
        if tmpA.hiv_dx:
            stats[tmpA.race][tmpA.so]["numTested"] += 1
        if tmpA.haart:
            stats[tmpA.race][tmpA.so]["numART"] += 1

    # PWID agent summary
    for tmpA in totalAgents.subset["DU"].subset["Inj"]:
        stats[tmpA.race]["PWID"]["numAgents"] += 1
        if tmpA.hiv:
            stats[tmpA.race]["PWID"]["numHIV"] += 1
        if tmpA.aids:
            stats[tmpA.race]["PWID"]["numAIDS"] += 1
        if tmpA.hiv_dx:
            stats[tmpA.race]["PWID"]["numTested"] += 1
        if tmpA.haart:
            stats[tmpA.race]["PWID"]["numART"] += 1

    # total number of agents
    for tmpA in totalAgents:
        stats[tmpA.race][tmpA.so]["numAgents"] += 1
        if tmpA.vaccine:
            stats[tmpA.race][tmpA.so]["Vaccinated"] += 1

    for tmpA in deathSet:
        stats[tmpA.race][tmpA.so]["deaths"] += 1
        if tmpA.hiv:
            stats[tmpA.race][tmpA.so]["deaths_HIV"] += 1

    # Sum 'ALL' categories for race/SO bins
    for race in stats:
        if race != "ALL":
            for param in stats_template:
                for sc in SUB_CAT:
                    if sc in params.classes.sex_types:
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
    params: DotMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "DeathReport.txt"), "a")
    sex_types = deepcopy(params.classes.sex_types)
    sex_types.append("ALL")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\ttot_{st}\tHIV_{st}"
        for sex_type in sex_types:
            f.write(template.format(st=sex_type))

        f.write("\n")

    f.write("%s\t%d\t%d" % (run_id, runseed, t))  # start row

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
    params: DotMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "IncarReport.txt"), "a")

    name_map = {
        "incar": "tot",
        "incarHIV": "HIV",
        "newRelease": "rlsd",
        "newReleaseHIV": "rlsdHIV",
    }

    MAIN_CAT = deepcopy(params.classes.races)
    MAIN_CAT.append("ALL")

    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\t{mc}_{st}_{p}"
        for p in name_map.values():
            for mc in MAIN_CAT:
                for sex_type in params.classes.sex_types:
                    f.write(template.format(mc=mc, st=sex_type, p=p))

        f.write("\n")

    f.write("%s\t%d\t%d" % (run_id, runseed, t))

    for p in name_map:
        for mc in MAIN_CAT:
            for st in params.classes.sex_types:
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
    params: DotMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "newlyHR_Report.txt"), "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\tnewHR_{st}\tnewHR_HIV_{st}\tnewHR_AIDS_{st}\tnewHR_Tested_{st}\tnewHR_ART_{st}"
        for sex_type in params.classes.sex_types:
            f.write(template.format(st=sex_type))

        f.write("\n")

    f.write("%s\t%d\t%d" % (run_id, runseed, t))  # start row

    for sex_type in params.classes.sex_types:
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
    params: DotMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "PrEPReport.txt"), "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt\tNewEnroll\tPWIDpartner\tTestedPartner\tMSMWpartner\n")

    f.write(
        "%s\t%d\t%d\t%d\t%d\t%d\t%d\n"
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
    params: DotMap,
    outdir: str,
):
    MAIN_CAT = deepcopy(params.classes.races)
    MAIN_CAT.append("ALL")
    SUB_CAT = deepcopy(params.classes.populations)
    SUB_CAT.append("ALL")

    for agentRace in MAIN_CAT:
        for agentTypes in SUB_CAT:
            name = "basicReport_" + agentTypes + "_" + agentRace + ".txt"
            tmpReport = open(os.path.join(outdir, name), "a")

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
                        stats[agentRace][agentTypes]["Vaccinated"],
                        stats[agentRace][agentTypes]["injectable_prep"],
                        stats[agentRace][agentTypes]["oral_prep"],
                        stats[agentRace][agentTypes]["prep_aware"],
                    )
                )
            )
            tmpReport.close()


# ========================== Other Print Functions =============================


def print_components(
    run_id: UUID,
    t: int,
    runseed: int,
    popseed: int,
    netseed: int,
    components,
    outdir: str,
):
    """
    Write stats describing the components (sub-graphs) in a graph to file
    """
    f = open(os.path.join(outdir, "componentReport_ALL.txt"), "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "run_id\trunseed\tpopseed\tnetseed\tt\tcompID\ttotalN\tNhiv\tNprep\tNtrtHIV\tTrtComponent\tPCA\tOral\tLAI\tAware\n"
        )

    comp_id = 0
    for comp in components:
        tot_agents = (
            nhiv
        ) = ntrthiv = nprep = trtbool = injectable_prep = oral = aware = pca = 0
        for agent in comp.nodes():
            tot_agents += 1
            if agent.hiv:
                nhiv += 1
                if agent.intervention_ever:
                    ntrthiv += 1
            if agent.prep:
                nprep += 1
                if agent.prep_type == "Inj":
                    injectable_prep += 1
                elif agent.prep_type == "Oral":
                    oral += 1
            if agent.pca:
                trtbool += 1
                if agent.pca_suitable:
                    pca += 1
            if agent.awareness:
                aware += 1

        f.write(
            f"{run_id}\t{runseed}\t{popseed}\t{netseed}\t{t}\t{comp_id}\t{tot_agents}\t{nhiv}\t{nprep}\t{ntrthiv}"
            f"\t{trtbool}\t{pca}\t{oral}\t{injectable_prep}\t{aware}\n"
        )

        comp_id += 1

    f.close()
