#!/usr/bin/env python
# encoding: utf-8

from typing import Dict, Any, List
from .agent import AgentSet, Agent
from copy import deepcopy
import os

from networkx import betweenness_centrality, effective_size, density  # type: ignore
from .parse_params import ObjMap


def get_stats(
    all_agents: AgentSet,
    new_prep_agents: AgentSet,
    new_hiv: AgentSet,
    new_hiv_dx: AgentSet,
    new_high_risk: AgentSet,
    new_incar_release: AgentSet,
    deaths: List[Agent],
    params: ObjMap,
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

    MAIN_CAT = [race for race in params.classes.races.keys()]
    MAIN_CAT.append("ALL")
    SUB_CAT = deepcopy(params.classes.populations)
    SUB_CAT.append("ALL")

    stats = {}
    for cat in MAIN_CAT:
        stats[cat] = {sc: dict(stats_template) for sc in SUB_CAT}

    # Incarceration metrics
    for a in new_incar_release:
        stats[a.race][a.so]["newRelease"] += 1
        if a.hiv:
            stats[a.race][a.so]["newReleaseHIV"] += 1

    # Newly infected tracker statistics (with HR within 6mo and HR ever bool check)
    for a in new_hiv:
        stats[a.race][a.so]["inf_newInf"] += 1
        if a.high_risk_ever:
            stats[a.race][a.so]["inf_HRever"] += 1
        if a.high_risk:
            stats[a.race][a.so]["inf_HR6m"] += 1

    for a in all_agents:
        stats[a.race][a.so]["numAgents"] += 1

        if a.prep:
            stats[a.race][a.so]["numPrEP"] += 1
            if "PWID" in a.prep_reason:
                stats[a.race][a.so]["iduPartPrep"] += 1
            if "MSMW" in a.prep_reason:
                stats[a.race][a.so]["msmwPartPrep"] += 1
            if "HIV test" in a.prep_reason:
                stats[a.race][a.so]["testedPartPrep"] += 1
            if a.prep_type == "Inj":
                stats[a.race][a.so]["injectable_prep"] += 1
            elif a.prep_type == "Oral":
                stats[a.race][a.so]["oral_prep"] += 1

        if a.incar:
            stats[a.race][a.so]["incar"] += 1
            if a.hiv:
                stats[a.race][a.so]["incarHIV"] += 1

        if a.hiv:
            stats[a.race][a.so]["numHIV"] += 1
            if a.aids:
                stats[a.race][a.so]["numAIDS"] += 1
            if a.hiv_dx:
                stats[a.race][a.so]["numTested"] += 1
            if a.haart:
                stats[a.race][a.so]["numART"] += 1

        if a.drug_use == "Inj":
            stats[a.race]["PWID"]["numAgents"] += 1
            if a.hiv:
                stats[a.race]["PWID"]["numHIV"] += 1
            if a.aids:
                stats[a.race]["PWID"]["numAIDS"] += 1
            if a.hiv_dx:
                stats[a.race]["PWID"]["numTested"] += 1
            if a.haart:
                stats[a.race]["PWID"]["numART"] += 1

        if a.vaccine:
            stats[a.race][a.so]["Vaccinated"] += 1

    # Newly PrEP tracker statistics
    for a in new_prep_agents:
        stats[a.race][a.so]["newNumPrEP"] += 1

    # Newly diagnosed tracker statistics
    for a in new_hiv_dx:
        stats[a.race][a.so]["newlyTested"] += 1

    # Newly HR agents
    for a in new_high_risk:
        stats[a.race][a.so]["newHR"] += 1
        if a.hiv:
            stats[a.race][a.so]["newHR_HIV"] += 1
            if a.aids:
                stats[a.race][a.so]["newHR_AIDS"] += 1
            if a.hiv_dx:
                stats[a.race][a.so]["newHR_tested"] += 1
                if a.haart:
                    stats[a.race][a.so]["newHR_ART"] += 1

    for a in deaths:
        stats[a.race][a.so]["deaths"] += 1
        if a.hiv:
            stats[a.race][a.so]["deaths_HIV"] += 1

    # Sum 'ALL' categories for race/SO bins
    for race in stats:
        if race != "ALL":
            for param in stats_template:
                for sc in SUB_CAT:
                    if sc in params.classes.sex_types:
                        stats[race]["ALL"][param] += stats[race][sc][param]
                for sc in SUB_CAT:
                    stats["ALL"][sc][param] += stats[race][sc][param]

    return stats


# ================== Printer Functions =========================
# Each of the following functions takes in the time, seeds, and stats dict for that time
# and prints the appropriate stats to file


def deathReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "DeathReport.txt"), "a")
    sex_types = deepcopy(list(params.classes.sex_types.keys()))
    sex_types.append("ALL")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\ttot_{st}\tHIV_{st}"
        for sex_type in sex_types:
            f.write(template.format(st=sex_type))

        f.write("\n")

    f.write(f"{run_id}\t{runseed}\t{t}")  # start row

    for sex_type in sex_types:
        f.write(
            "\t{}\t{}".format(
                stats["ALL"][sex_type]["deaths"], stats["ALL"][sex_type]["deaths_HIV"]
            )
        )

    f.write("\n")
    f.close()


def incarReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "IncarReport.txt"), "a")

    name_map = {
        "incar": "tot",
        "incarHIV": "HIV",
        "newRelease": "rlsd",
        "newReleaseHIV": "rlsdHIV",
    }

    MAIN_CAT = [race for race in params.classes.races.keys()]
    MAIN_CAT.append("ALL")

    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = "\t{mc}_{st}_{p}"
        for p in name_map.values():
            for mc in MAIN_CAT:
                for sex_type in params.classes.sex_types:
                    f.write(template.format(mc=mc, st=sex_type, p=p))

        f.write("\n")

    f.write(f"{run_id}\t{runseed}\t{t}")

    for p in name_map:
        for mc in MAIN_CAT:
            for st in params.classes.sex_types:
                f.write("\t")
                f.write(str(stats[mc][st][p]))

    f.write("\n")
    f.close()


def newlyhighriskReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "newlyHR_Report.txt"), "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt")  # start header

        template = (
            "\tnewHR_{st}\tnewHR_HIV_{st}\tnewHR_AIDS_{st}\t"
            "newHR_Tested_{st}\tnewHR_ART_{st}"
        )
        for sex_type in params.classes.sex_types:
            f.write(template.format(st=sex_type))

        f.write("\n")

    f.write(f"{run_id}\t{runseed}\t{t}")  # start row

    for sex_type in params.classes.sex_types:
        f.write(
            "\t{}\t{}\t{}\t{}\t{}".format(
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
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    f = open(os.path.join(outdir, "PrEPReport.txt"), "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("run_id\tseed\tt\tNewEnroll\tPWIDpartner\tTestedPartner\tMSMWpartner\n")

    f.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
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
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    MAIN_CAT = [race for race in params.classes.races.keys()]
    MAIN_CAT.append("ALL")
    SUB_CAT = deepcopy(params.classes.populations)
    SUB_CAT.append("ALL")

    for race in MAIN_CAT:
        for population in SUB_CAT:
            name = "basicReport_" + population + "_" + race + ".txt"
            f = open(os.path.join(outdir, name), "a")

            # if this is a new file, write the header info
            if f.tell() == 0:
                f.write(
                    "run_id\trseed\tpseed\tt\tTotal\tHIV\tAIDS\tTstd\tART\tnHR\tIncid\t"
                    "HR_6mo\tHR_Ev\tNewDiag\tDeaths\tPrEP\tIDUpart_PrEP\t"
                    "MSMWpart_PrEP\ttestedPart_PrEP\tVaccinated\tLAI\tOral\tAware\n"
                )

            f.write(
                (
                    "{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t"
                    "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t"
                    "{:d}\n".format(
                        str(run_id),
                        runseed,
                        popseed,
                        t,
                        stats[race][population]["numAgents"],
                        stats[race][population]["numHIV"],
                        stats[race][population]["numAIDS"],
                        stats[race][population]["numTested"],
                        stats[race][population]["numART"],
                        stats[race][population]["numHR"],
                        stats[race][population]["inf_newInf"],
                        stats[race][population]["inf_HR6m"],
                        stats[race][population]["inf_HRever"],
                        stats[race][population]["newlyTested"],
                        stats[race][population]["deaths"],
                        stats[race][population]["numPrEP"],
                        stats[race][population]["iduPartPrep"],
                        stats[race][population]["msmwPartPrep"],
                        stats[race][population]["testedPartPrep"],
                        stats[race][population]["Vaccinated"],
                        stats[race][population]["injectable_prep"],
                        stats[race][population]["oral_prep"],
                        stats[race][population]["prep_aware"],
                    )
                )
            )
            f.close()


# ========================== Other Print Functions =============================


def print_components(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    components,
    outdir: str,
    races: list,
    pca_bool: bool,
):
    """
    Write stats describing the components (sub-graphs) in a graph to file
    """
    f = open(os.path.join(outdir, f"{run_id}_componentReport_ALL.txt"), "a")

    race_count: Dict[str, int] = {}
    race_list = []
    header = ""
    for race in races:
        race_list.append(race)
        header += "\t" + race

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "run_id\trunseed\tpopseed\tt\tcompID\ttotalN\tNhiv\tNprep\tNtrtHIV"
            "\tTrtComponent\tPCA\tOral\tInjectable\tAware\tnidu\tcentrality\tDensity\tEffectiveSize"
            + header
            + "\n"
        )

    comp_id = 0
    for comp in components:
        for race in races:
            race_count[race] = 0
        assert comp.number_of_nodes() >= 0
        tot_agents = (
            nhiv
        ) = ntrthiv = nprep = trtbool = injectable_prep = oral = aware = pca = nidu = 0

        for agent in comp.nodes():
            tot_agents += 1
            race_count[agent.race] += 1
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

            if pca_bool:  # pca intervention characteristics
                if agent.pca:
                    trtbool += 1
                    if agent.pca_suitable:
                        pca += 1
                elif agent.intervention_ever:
                    trtbool = 1
            else:
                if agent.intervention_ever:
                    trtbool = 1

            if agent.prep_awareness:
                aware += 1

            if agent.drug_use == "NonInj":
                nidu += 1

        comp_centrality = (
            sum(betweenness_centrality(comp).values()) / comp.number_of_nodes()
        )
        average_size = sum(effective_size(comp).values()) / comp.number_of_nodes()
        comp_density = density(comp)

        race_str = ""
        for race in race_list:
            race_str += "\t" + str(race_count[race])

        f.write(
            f"{run_id}\t{runseed}\t{popseed}\t{t}\t{comp_id}\t{tot_agents}\t"
            f"{nhiv}\t"
            f"{nprep}\t{ntrthiv}\t{trtbool}\t{pca}\t{oral}\t{injectable_prep}"
            f"\t{aware}\t{nidu}\t{comp_centrality:.4f}\t{comp_density:.4f}"
            f"\t{average_size:.4f}{race_str}\n"
        )

        comp_id += 1

    f.close()
