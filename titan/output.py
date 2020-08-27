#!/usr/bin/env python
# encoding: utf-8

from typing import Dict, Any, List, Iterator
from .agent import AgentSet, Agent
import itertools
import os

from networkx import betweenness_centrality, effective_size, density  # type: ignore
from .parse_params import ObjMap


def setup_aggregates(params: ObjMap, classes: List[str]) -> Dict:
    """
    Recursively create a nested dictionary of attribute values to items to count.

    Attributes are classes defined in params, the items counted are:

    * "numAgents"
    * "inf_HR6m"
    * "inf_HRever"
    * "inf_newInf"
    * "newHR"
    * "newHR_HIV"
    * "newHR_AIDS"
    * "newHR_dx"
    * "newHR_ART"
    * "newRelease"
    * "newReleaseHIV"
    * "numHIV"
    * "numTested"
    * "numAIDS"
    * "numART"
    * "numHR"
    * "newlyTested"
    * "deaths"
    * "deaths_HIV"
    * "incar"
    * "incarHIV"
    * "numPrEP"
    * "newNumPrEP"
    * "iduPartPrep"
    * "msmwPartPrep"
    * "testedPartPrep"
    * "vaccinated"
    * "injectable_prep"
    * "oral_prep"

    args:
        params: model parameters
        classes: which classes to aggregate by [params.outputs.classes]

    returns:
        dictionary of class values to counts
    """
    if classes == []:
        return {
            "numAgents": 0,
            "inf_HR6m": 0,
            "inf_HRever": 0,
            "inf_newInf": 0,
            "newHR": 0,
            "newHR_HIV": 0,
            "newHR_AIDS": 0,
            "newHR_dx": 0,
            "newHR_ART": 0,
            "newRelease": 0,
            "newReleaseHIV": 0,
            "numHIV": 0,
            "numDiagnosed": 0,
            "numAIDS": 0,
            "numART": 0,
            "numHR": 0,
            "newlyDiagnosed": 0,
            "deaths": 0,
            "deaths_HIV": 0,
            "incar": 0,
            "incarHIV": 0,
            "numPrEP": 0,
            "newNumPrEP": 0,
            "iduPartPrep": 0,
            "msmwPartPrep": 0,
            "testedPartPrep": 0,
            "vaccinated": 0,
            "injectable_prep": 0,
            "oral_prep": 0,
        }

    stats = {}
    clss, *rem_clss = classes  # head, tail
    keys = [k for k in params.classes[clss]]
    for key in keys:
        stats[key] = setup_aggregates(params, rem_clss)

    return stats


def get_aggregates(params: ObjMap) -> Iterator:
    """
    Get iterator over all attribute combinations for output classes

    args:
        params: model parameters

    returns:
        iterator over attribute combinations
    """
    return itertools.product(
        *[list(k for k in params.classes[clss]) for clss in params.outputs.classes]
    )


def get_agg_val(stats: Dict, attrs: List, key: str) -> int:
    """
    Get the value of a key in stats given the attribute values

    args:
        stats: a nested dictionariy of attributes to counts
        attrs: a list of attribute values to find the count for
        key: the type of count to get the value of

    returns:
        the count of key for the given attributes
    """
    stats_item = stats
    for attr in attrs:
        stats_item = stats_item[attr]

    return stats_item[key]


def add_agent_to_stats(stats: Dict[str, Any], attrs: List[str], agent: Agent, key: str):
    """
    Update the stats dictionary counts for the key given the agent's attributes

    args:
        stats: a nested dictionary of attributes to counts
        attrs: a list of attribute types (e.g. "race")
        agent: the agent whose attribute values will be evaluated
        key: the type of count to increment
    """
    stats_item = stats
    for attr in attrs:
        stats_item = stats_item[str(getattr(agent, attr))]

    stats_item[key] += 1


def get_stats(
    all_agents: AgentSet,
    new_prep_agents: AgentSet,
    new_hiv: AgentSet,
    new_hiv_dx: AgentSet,
    new_high_risk: AgentSet,
    new_incar_release: AgentSet,
    deaths: List[Agent],
    params: ObjMap,
) -> Dict:
    """
    Get the current statistics for a model based on the population, and tracking agent sets from the model.

    args:
        all_agents: all of the agents in the population
        new_prep_agents: agents who are newly on prep this timestep
        new_hiv: agents are newly hiv this timestep
        new_hiv_dx: agents who are newly diagnosed with hiv this timestep
        new_high_risk: agents are newly high risk this timestep
        new_incar_release: agents are released from incarceration this timestep
        deaths: agents who died this timestep
        params: model parameters

    returns:
        nested dictionary of agent attributes to counts of various items
    """
    stats = setup_aggregates(params, params.outputs.classes)
    attrs = [
        clss[:-1] for clss in params.outputs.classes
    ]  # attribute version (non-plural)

    # Incarceration metrics
    for a in new_incar_release:
        add_agent_to_stats(stats, attrs, a, "newRelease")
        if a.hiv:
            add_agent_to_stats(stats, attrs, a, "newReleaseHIV")

    # Newly infected tracker statistics (with HR within 6mo and HR ever bool check)
    for a in new_hiv:
        add_agent_to_stats(stats, attrs, a, "inf_newInf")
        if a.high_risk_ever:
            add_agent_to_stats(stats, attrs, a, "inf_HRever")
        if a.high_risk:
            add_agent_to_stats(stats, attrs, a, "inf_HR6m")

    for a in all_agents:
        add_agent_to_stats(stats, attrs, a, "numAgents")

        if a.prep:
            add_agent_to_stats(stats, attrs, a, "numPrEP")
            if "PWID" in a.prep_reason:
                add_agent_to_stats(stats, attrs, a, "iduPartPrep")
            if "MSMW" in a.prep_reason:
                add_agent_to_stats(stats, attrs, a, "msmwPartPrep")
            if "HIV test" in a.prep_reason:
                add_agent_to_stats(stats, attrs, a, "testedPartPrep")
            if a.prep_type == "Inj":
                add_agent_to_stats(stats, attrs, a, "injectable_prep")
            elif a.prep_type == "Oral":
                add_agent_to_stats(stats, attrs, a, "oral_prep")

        if a.incar:
            add_agent_to_stats(stats, attrs, a, "incar")
            if a.hiv:
                add_agent_to_stats(stats, attrs, a, "incarHIV")

        if a.hiv:
            add_agent_to_stats(stats, attrs, a, "numHIV")
            if a.aids:
                add_agent_to_stats(stats, attrs, a, "numAIDS")
            if a.hiv_dx:
                add_agent_to_stats(stats, attrs, a, "numDiagnosed")
            if a.haart:
                add_agent_to_stats(stats, attrs, a, "numART")

        if a.vaccine:
            add_agent_to_stats(stats, attrs, a, "vaccinated")

    # Newly PrEP tracker statistics
    for a in new_prep_agents:
        add_agent_to_stats(stats, attrs, a, "newNumPrEP")

    # Newly diagnosed tracker statistics
    for a in new_hiv_dx:
        add_agent_to_stats(stats, attrs, a, "newlyDiagnosed")

    # Newly HR agents
    for a in new_high_risk:
        add_agent_to_stats(stats, attrs, a, "newHR")
        if a.hiv:
            add_agent_to_stats(stats, attrs, a, "newHR_HIV")
            if a.aids:
                add_agent_to_stats(stats, attrs, a, "newHR_AIDS")
            if a.hiv_dx:
                add_agent_to_stats(stats, attrs, a, "newHR_dx")
                if a.haart:
                    add_agent_to_stats(stats, attrs, a, "newHR_ART")

    for a in deaths:
        add_agent_to_stats(stats, attrs, a, "deaths")
        if a.hiv:
            add_agent_to_stats(stats, attrs, a, "deaths_HIV")

    return stats


# ================== Printer Functions =========================
# Each of the following functions takes in the time, seeds, and stats dict for that time
# and prints the appropriate stats to file


def write_report(
    file_name: str,
    name_map: Dict[str, str],
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict,
    params: ObjMap,
    outdir: str,
):
    """
    Core function for writing reports, writes header if file is new, then data based
    on the `params` and `name_map`

    args:
        file_name: Name of the file to write, including the extension (e.g. `MyReport.txt`)
        name_map: Map from keys in the stats dictionary to column headers in this report
        run_id: unique identifier for this model
        t: current timestep
        runseed: integer used to seed the random number generator for the model
        popseed: integer used to seed the random number generator for the population
        stats: nested dictionary of agent attributes to counts
        params: model parameters
        outdir: path of where to save this file
    """
    f = open(os.path.join(outdir, file_name), "a")
    attrs = [clss[:-1] for clss in params.outputs.classes]

    if f.tell() == 0:
        f.write("run_id\trseed\tpseed\tt\t")  # start header

        # attributes in stats
        f.write("\t".join(attrs))

        # report specific fields
        for name in name_map.values():
            f.write(f"\t{name}")

        f.write("\n")

    for agg in get_aggregates(params):
        f.write(f"{run_id}\t{runseed}\t{popseed}\t{t}\t")

        f.write("\t".join(agg))  # write attribute values

        for key, name in name_map.items():
            f.write(f"\t{(get_agg_val(stats, agg, key))}")

        f.write("\n")

    f.close()


def deathReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    """
    Standard report writer for agent deaths, columns include:

    * `tot`: total number of deaths at this timestep
    * `HIV`: number of deaths of agents with HIV at this timestep
    """
    name_map = {
        "deaths": "tot",
        "deaths_HIV": "HIV",
    }
    write_report(
        "DeathReport.txt", name_map, run_id, t, runseed, popseed, stats, params, outdir
    )


def incarReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    """
    Standard report writer for agent incarcerations, columns include:

    * `tot`: total number of incarcerated agents at this timestep
    * `HIV`: number of incarcerated agents with HIV at this timestep
    * `rlsd`: number of agents released at this timestep
    * `rlsdHIV`: number of agents with HIV released at this timestep
    """
    name_map = {
        "incar": "tot",
        "incarHIV": "HIV",
        "newRelease": "rlsd",
        "newReleaseHIV": "rlsdHIV",
    }
    write_report(
        "IncarReport.txt", name_map, run_id, t, runseed, popseed, stats, params, outdir
    )


def newlyhighriskReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    """
    Standard report writer for newly high risk agents, columns include:

    * `newHR`: number of agents that became high risk at this timestep
    * `newHR_HIV`: number of agents with HIV that became high risk at this timestep
    * `newHR_AIDS`: number of agents with AIDS that became high risk at this timestep
    * `newHR_Tested`: number of agents with diagnosed HIV that became high risk at this timestep
    * `newHR_ART`: number of agents on HAART that became high risk at this timestep
    """
    name_map = {
        "newHR": "newHR",
        "newHR_HIV": "newHR_HIV",
        "newHR_AIDS": "newHR_AIDS",
        "newHR_dx": "newHR_Diagnosed",
        "newHR_ART": "newHR_ART",
    }
    write_report(
        "newlyHR_Report.txt",
        name_map,
        run_id,
        t,
        runseed,
        popseed,
        stats,
        params,
        outdir,
    )


def prepReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    """
    Standard report writer for prep agents, columns include:

    * `newEnroll`: number of agents newly enrolled in PrEP at this timestep
    * `PWIDpartner`: number of prep agents with a PWID partner at this timestep
    * `TestedPartner`: number of prep agents with a partner with diagnosed HIV at this timestep
    * `MSMWpartner`: number of prep agents with an MSMW partner at this timestep
    """
    name_map = {
        "newNumPrEP": "NewEnroll",
        "iduPartPrep": "PWIDpartner",
        "testedPartPrep": "DiagnosedPartner",
        "msmwPartPrep": "MSMWpartner",
    }
    write_report(
        "PrEPReport.txt", name_map, run_id, t, runseed, popseed, stats, params, outdir
    )


def basicReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    """
    Standard report writer for basic agent statistics, columns include:

    * `Total`: number of agents in the population
    * `HIV`: number of agents with HIV
    * `AIDS`: number of agents with AIDS
    * `Dx`: number of agents with HIV who are diagnosed
    * `ART`: number of agents on HAART
    * `nHR`: number of agents who are high risk
    * `Incid`: number of agents who HIV converted this time period
    * `HR_6mo`: number of agents with HIV converted this time period who are high risk
    * `HR_Ev`: number of agents with HIV converted this time period who have ever been high risk
    * `NewDiag`: number of agents with HIV who were diagnosed this time period
    * `Deaths`: number of agents who died this time period
    * `PrEP`: number of agents enrolled in PrEP
    * `IDUpart_PrEP`: number of agents enrolled in PrEP who have a PWID partner
    * `MSMWpart_PrEP`: number of agents enrolled in PrEP who have an MSMW partner
    * `testedPart_PrEP`: number of agents enrolled in PrEP who have an HIV diagnosed partner
    * `Vaccinated`: number of agents who have been vaccinated
    * `LAI`: number of agents enrolled in PrEP with a type of LAI
    * `Oral`: number of agents enrolled in PrEP with a type of Oral
    """
    name_map = {
        "numAgents": "Total",
        "numHIV": "HIV",
        "numAIDS": "AIDS",
        "numDiagnosed": "Dx",
        "numART": "ART",
        "numHR": "nHR",
        "inf_newInf": "Incid",
        "inf_HR6m": "HR_6mo",
        "inf_HRever": "HR_Ev",
        "newlyDiagnosed": "NewDx",
        "deaths": "Deaths",
        "numPrEP": "PrEP",
        "iduPartPrep": "IDUpart_PrEP",
        "msmwPartPrep": "MSMWpart_PrEP",
        "testedPartPrep": "testedPart_PrEP",
        "vaccinated": "Vaccinated",
        "injectable_prep": "LAI",
        "oral_prep": "Oral",
    }
    write_report(
        "basicReport.txt", name_map, run_id, t, runseed, popseed, stats, params, outdir
    )


# ========================== Other Print Functions =============================


def print_components(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    components: List,
    outdir: str,
    races: list,
):
    """
    Write stats describing the components (sub-graphs) in a graph to file

    args:
        run_id: unique identifer for this run of the model
        t: current timestep
        runseed: integer used to seed the model's random number generator
        popseed: integer used to seed the population's random number generator
        components: a list of graph components
        outdir: path where the file should be saved
        races: the races in the population
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
            "\tComp_Status\tPCA\tOral\tInjectable\tAware\tnidu\tcentrality\tDensity"
            "\tEffectiveSize" + header + "\n"
        )

    comp_id = 0
    for comp in components:
        for race in races:
            race_count[race] = 0
        assert comp.number_of_nodes() >= 0
        tot_agents = (
            nhiv
        ) = ntrthiv = nprep = injectable_prep = oral = aware = pca = nidu = 0
        component_status = "control"
        trt_comp = trt_agent = False

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

            if agent.pca_suitable and agent.pca:
                pca += 1

            if agent.intervention_ever:  # treatment component
                trt_agent = True

            if agent.random_trial_enrolled:
                trt_comp = True

            if agent.prep_awareness:
                aware += 1

            if agent.drug_type == "NonInj":
                nidu += 1

        comp_centrality = (
            sum(betweenness_centrality(comp).values()) / comp.number_of_nodes()
        )
        average_size = sum(effective_size(comp).values()) / comp.number_of_nodes()
        comp_density = density(comp)

        if trt_comp:
            if trt_agent:
                component_status = "treatment"
            else:
                component_status = "treatment_no_eligible"

        race_str = ""
        for race in race_list:
            race_str += "\t" + str(race_count[race])

        f.write(
            f"{run_id}\t{runseed}\t{popseed}\t{t}\t{comp_id}\t{tot_agents}\t"
            f"{nhiv}\t"
            f"{nprep}\t{ntrthiv}\t{component_status}\t{pca}\t{oral}\t{injectable_prep}"
            f"\t{aware}\t{nidu}\t{comp_centrality:.4f}\t{comp_density:.4f}"
            f"\t{average_size:.4f}{race_str}\n"
        )

        comp_id += 1

    f.close()
