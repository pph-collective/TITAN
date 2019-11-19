#!/usr/bin/env python
# encoding: utf-8

from . import params  # type: ignore
from typing import Dict, Any, List, Sequence, Optional
from .agent import Agent_set, Relationship, Agent

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
):

    stats_template = {
        "numAgents": 0,
        "inf_HR6m": 0,
        "inf_HRever": 0,
        "inf_newInf": 0,
        "newHR": 0,
        "newHR_HIV": 0,
        "newHR_AIDS": 0,
        "newHR_Tested": 0,
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

    for tmpA in deathSet:
        stats[tmpA._race][tmpA._SO]["deaths"] += 1
        if tmpA._HIV_bool:
            stats[tmpA._race][tmpA._SO]["deaths_HIV"] += 1

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


def incidenceReport(
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    f = open("results/IncidenceReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("seed\tt\tTotal\tW_HM\tB_HM\tHM\tW_HF\tB_HF\tHF\tW_MSM\tB_MSM\tMSM\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["ALL"]["inf_newInf"],
            stats["WHITE"]["HM"]["inf_newInf"],
            stats["BLACK"]["HM"]["inf_newInf"],
            stats["ALL"]["HM"]["inf_newInf"],
            stats["WHITE"]["HF"]["inf_newInf"],
            stats["BLACK"]["HF"]["inf_newInf"],
            stats["ALL"]["HF"]["inf_newInf"],
            stats["WHITE"]["MSM"]["inf_newInf"],
            stats["BLACK"]["MSM"]["inf_newInf"],
            stats["ALL"]["MSM"]["inf_newInf"],
        )
    )
    f.close()


def prevalenceReport(
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    f = open("results/PrevalenceReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("seed\tt\tTotal\tHM\tHF\tHIV_tot\tHIV_HM\tHIV_HF\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["ALL"]["numAgents"],
            stats["ALL"]["HM"]["numAgents"],
            stats["ALL"]["HF"]["numAgents"],
            stats["ALL"]["ALL"]["numHIV"],
            stats["ALL"]["HM"]["numHIV"],
            stats["ALL"]["HF"]["numHIV"],
        )
    )
    f.close()


def deathReport(
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    f = open("results/DeathReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("seed\tt\tTotal\tHM\tMSM\tHF\tHIV_tot\tHIV_HM\tHIV_MSM\tHIV_HF\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["ALL"]["deaths"],
            stats["ALL"]["HM"]["deaths"],
            stats["ALL"]["MSM"]["deaths"],
            stats["ALL"]["HF"]["deaths"],
            stats["ALL"]["ALL"]["deaths_HIV"],
            stats["ALL"]["HM"]["deaths_HIV"],
            stats["ALL"]["MSM"]["deaths_HIV"],
            stats["ALL"]["HF"]["deaths_HIV"],
        )
    )
    f.close()


def incarReport(
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    f = open("results/IncarReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "seed\tt\tTotal\tW_HM\tB_HM\tW_HF\tB_HF\tW_MSM\tB_MSM\tW_HIV\tB_HIV\tW_rlsd_HM\tW_rlsd_HF\tB_rlsd_HM\tB_rlsd_HF\tW_rlsdHIV\tB_rlsdHIV\n"
        )

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["ALL"]["incar"],
            stats["WHITE"]["HM"]["incar"],
            stats["BLACK"]["HM"]["incar"],
            stats["WHITE"]["HF"]["incar"],
            stats["BLACK"]["HF"]["incar"],
            stats["WHITE"]["MSM"]["incar"],
            stats["BLACK"]["MSM"]["incar"],
            stats["WHITE"]["ALL"]["incarHIV"],
            stats["BLACK"]["ALL"]["incarHIV"],
            stats["WHITE"]["HM"]["newRelease"],
            stats["WHITE"]["HF"]["newRelease"],
            stats["BLACK"]["HM"]["newRelease"],
            stats["BLACK"]["HF"]["newRelease"],
            stats["WHITE"]["ALL"]["newReleaseHIV"],
            stats["BLACK"]["ALL"]["newReleaseHIV"],
        )
    )
    f.close()


def iduReport(t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]):
    f = open("results/iduReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("seed\tt\tTotal-IDU\tIDU-HIV\tIDU-AIDS\tIDU-HAART\tIDU-tested\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["IDU"]["numAgents"],
            stats["ALL"]["IDU"]["numHIV"],
            stats["ALL"]["IDU"]["numAIDS"],
            stats["ALL"]["IDU"]["numART"],
            stats["ALL"]["IDU"]["numTested"],
        )
    )
    f.close()


def highriskReport(
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    f = open("results/HR_incidenceReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("seed\tt\tTot_Ev\tHM_Ev\tHF_Ev\tTot_6mo\tHM_6mo\tHF_6mo\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["ALL"]["inf_HRever"],
            stats["ALL"]["HM"]["inf_HRever"],
            stats["ALL"]["HF"]["inf_HRever"],
            stats["ALL"]["ALL"]["inf_HR6m"],
            stats["ALL"]["HM"]["inf_HR6m"],
            stats["ALL"]["HF"]["inf_HR6m"],
        )
    )
    f.close()


def newlyhighriskReport(
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    f = open("results/newlyHR_Report.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "seed\tt\tnewHR_HM\tnewHR_HIV_HM\tnewHR_AIDS_HM\tnewHR_Tested_HM\tnewHR_ART_HM\tnewHR_HF\tnewHR_HIV_HF\tnewHR_AIDS_HF\tnewHR_Tested_HF\tnewHR_ART_HF\n"
        )

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            runseed,
            t,
            stats["ALL"]["HM"]["newHR"],
            stats["ALL"]["HM"]["newHR_HIV"],
            stats["ALL"]["HM"]["newHR_AIDS"],
            stats["ALL"]["HM"]["newHR_Tested"],
            stats["ALL"]["HM"]["newHR_ART"],
            stats["ALL"]["HF"]["newHR"],
            stats["ALL"]["HF"]["newHR_HIV"],
            stats["ALL"]["HF"]["newHR_AIDS"],
            stats["ALL"]["HF"]["newHR_Tested"],
            stats["ALL"]["HF"]["newHR_ART"],
        )
    )
    f.close()


# TO_REVIEW this is really the HM vs HF vs MSM report - should it be more general?
def sexReport(t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]):
    for sexType in ["HF", "HM", "MSM"]:
        if sexType == "HF":
            longName = "Female"
        elif sexType == "HM":
            longName = "Male"
        else:
            longName = sexType

        f = open(f"results/{longName}Report.txt", "a")

        # if this is a new file, write the header info
        if f.tell() == 0:
            f.write(
                "seed\tt\tTotal\tHIV\tAIDS\tTested\tART\tIncidence\tHRInc_6mo\tHRInc_Ev\tNewlyDiag\tDeaths\n"
            )

        f.write(
            (
                "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
                % (
                    runseed,
                    t,
                    stats["ALL"][sexType]["numAgents"],
                    stats["ALL"][sexType]["numHIV"],
                    stats["ALL"][sexType]["numAIDS"],
                    stats["ALL"][sexType]["numTested"],
                    stats["ALL"][sexType]["numART"],
                    stats["ALL"][sexType]["inf_newInf"],
                    stats["ALL"][sexType]["inf_HR6m"],
                    stats["ALL"][sexType]["inf_HRever"],
                    stats["ALL"][sexType]["newlyTested"],
                    stats["ALL"][sexType]["deaths"],
                    stats["ALL"][sexType]["numPrEP"],
                )
            )
        )
        f.close()


# TO_REVIEW data doesn't match headers
def raceReport(t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]):
    for race in ["WHITE", "BLACK"]:
        r_init = race[0]
        f = open(f"results/{r_init}_pop_report.txt", "a")

        # if this is a new file, write the header info
        if f.tell() == 0:
            f.write("seed\tt\tTotal\tHM\tHF\tHIV_tot\tHIV_HM\tHIV_HF\n")

        f.write(
            "%d\t%d\t%d\t%d\t%d\t%d\n"
            % (
                runseed,
                t,
                stats[race]["ALL"]["numHIV"],
                stats[race]["MSM"]["numHIV"],
                stats[race]["ALL"]["numTested"],
                stats[race]["ALL"]["numART"],
            )
        )
        f.close()


def prepReport(t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]):
    f = open("results/PrEPReport.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write("seed\tt\tNewEnroll\tIDUpartner\tTestedPartner\tMSMWpartner\n")

    f.write(
        "%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
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
    t: int, runseed: int, popseed: int, netseed: int, stats: Dict[str, Any]
):
    for agentRace in MAIN_CAT:
        for agentTypes in SUB_CAT:
            name = "basicReport_" + agentTypes + "_" + agentRace
            tmpReport = open("results/" + name + ".txt", "a")

            # if this is a new file, write the header info
            if tmpReport.tell() == 0:
                tmpReport.write(
                    "rseed\tpseed\tnseed\tt\tTotal\tHIV\tAIDS\tTstd\tART\tnHR\tIncid\tHR_6mo\tHR_Ev\tNewDiag\tDeaths\tPrEP\tIDUpart_PrEP\tMSMWpart_PrEP\ttestedPart_PrEP\n"
                )

            tmpReport.write(
                (
                    "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
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
                    )
                )
            )
            tmpReport.close()


# ========================== Other Print Functions ============================


def print_components(t: int, runseed: int, popseed: int, netseed: int, components):
    """
    Write stats describing the components (sub-graphs) in a graph to file
    """
    f = open("results/componentReport_ALL.txt", "a")

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "runseed\tpopseed\tnetseed\tt\tcompID\ttotalN\tTestedPartner\tMSMWpartner\n"
        )

    compID = 0
    for comp in components:
        totN = nhiv = ntrtmt = ntrthiv = nprep = PrEP_ever_HIV = 0
        for ag in comp.nodes():
            totN += 1
            if ag._HIV_bool:
                nhiv += 1
                if ag._treatment_bool:
                    ntrthiv += 1
                if ag._PrEP_ever_bool:
                    PrEP_ever_HIV += 1
            elif ag._treatment_bool:
                ntrtmt += 1
                if ag._PrEP_bool:
                    nprep += 1
        f.write(
            "{runseed}\t{pseed}\t{nseed}\t{t}\t{compID}\t{totalN}\t{Nhiv}\t{Ntrtmt}\t{Nprep}\t{NtrtHIV}\t{NprepHIV}\n".format(
                runseed=runseed,
                pseed=popseed,
                nseed=netseed,
                t=t,
                compID=compID,
                totalN=totN,
                Nhiv=nhiv,
                Ntrtmt=ntrtmt,
                Nprep=nprep,
                NtrtHIV=ntrthiv,
                NprepHIV=PrEP_ever_HIV,
            )
        )

        compID += 1
    f.close()
