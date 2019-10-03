#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt

from . import params

try:
    from .agent import *
except ImportError:
    raise ImportError("Can't import Agent class")


def initiate_ResultDict():
    # nested dictionary for results (inner dictionary has the form: time:result)
    ResultDict = {
        "Prv_HIV": {},
        "Prv_AIDS": {},
        "Prv_Test": {},
        "Prv_ART": {},
        "Prv_PrEP": {},
        "n_Relations": {},
        "Inc_c_Tot": {},
        "Inc_c_HM": {},
        "Inc_c_HF": {},
        "Inc_t_HM": {},
        "Inc_t_HF": {},
    }

    return ResultDict


def print_stats(
    self,
    rseed,
    t,
    totalAgents,
    HIVAgents,
    IncarAgents,
    PrEPAgents,
    NewInfections,
    NewDiagnosis,
    deaths,
    ResultDict,
    Relationships,
    newHR,
    newIncarRelease,
    deathSet,
    outifle=None,
):
    incidenceReport = open("results/IncidenceReport.txt", "a")
    prevalenceReport = open("results/PrevalenceReport.txt", "a")
    deathReport = open("results/DeathReport.txt", "a")
    incarReport = open("results/IncarReport.txt", "a")
    iduReport = open("results/iduReport.txt", "a")
    highriskReport = open("results/HR_incidenceReport.txt", "a")
    newlyhighriskReport = open("results/newlyHR_Report.txt", "a")
    femaleReport = open("results/FemaleReport.txt", "a")
    maleReport = open("results/MaleReport.txt", "a")
    msmReport = open("results/MSMReport.txt", "a")
    nalReport = open("results/nalReport.txt", "a")
    oatReport = open("results/oatReport.txt", "a")
    peopleOff = open("results/peopleOff.txt", "a")
    whiteReport = open("results/W_pop_report.txt", "a")
    blackReport = open("results/B_pop_report.txt", "a")
    num_SEP = 0

    OAT_IDU_F = 0
    OAT_IDU_M = 0
    OAT_NIDU_F = 0
    OAT_NIDU_M = 0
    Prior_Year_OAT = 0
    Naltrex_NIDU_M = 0
    Naltrex_NIDU_F = 0
    Naltrex_IDU_M = 0
    Naltrex_IDU_F = 0
    Prior_Year_Naltrex = 0
    DOC_OAT_M = 0
    DOC_OAT_F = 0
    DOC_Naltrex_M = 0
    DOC_Naltrex_F = 0

    MSMW_part = 0
    IDU_part = 0
    Test_part = 0

    newHR_HM = 0
    newHR_HIV_HM = 0
    newHR_AIDS_HM = 0
    newHR_Tested_HM = 0
    newHR_ART_HM = 0
    off_HF = 0
    off_HM = 0
    newHR_HF = 0
    newHR_HIV_HF = 0
    newHR_AIDS_HF = 0
    newHR_Tested_HF = 0
    newHR_ART_HF = 0

    rc_template = {
        "numAgents": 0,
        "inf_HR6m": 0,
        "inf_HRever": 0,
        "inf_newInf": 0,
        "newHighRisk": 0,
        "newRelease": 0,
        "newReleaseHIV": 0,
        "numHIV": 0,
        "numTested": 0,
        "numAIDS": 0,
        "numART": 0,
        "numHR": 0,
        "newlyTested": 0,
        "deaths": 0,
        "incar": 0,
        "incarHIV": 0,
        "numPrEP": 0,
        "iduPartPrep": 0,
        "msmwPartPrep": 0,
        "testedPartPrep": 0,
    }

    rc1_infections = {
        "MTF": dict(rc_template),
        "MSM": dict(rc_template),
        "HM": dict(rc_template),
        "HF": dict(rc_template),
        "IDU": dict(rc_template),
        "ALL": dict(rc_template),
    }

    rc2_infections = {
        "MTF": dict(rc_template),
        "MSM": dict(rc_template),
        "HM": dict(rc_template),
        "HF": dict(rc_template),
        "IDU": dict(rc_template),
        "ALL": dict(rc_template),
    }
    all_infections = {
        "MTF": dict(rc_template),
        "MSM": dict(rc_template),
        "HM": dict(rc_template),
        "HF": dict(rc_template),
        "IDU": dict(rc_template),
        "ALL": dict(rc_template),
    }
    rsltdic = {"WHITE": rc1_infections, "BLACK": rc2_infections, "ALL": all_infections}
    tot_rsltdic = {"ALL": all_infections}

    # Incarceration metrics
    for tmpA in IncarAgents.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["incar"] += 1
        if tmpA._HIV_bool:
            rsltdic[tmpA._race][tmpA._SO]["incarHIV"] += 1

    for tmpA in newIncarRelease.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["newRelease"] += 1
        if tmpA._HIV_bool:
            rsltdic[tmpA._race][tmpA._SO]["newReleaseHIV"] += 1

    # Newly infected tracker statistics (with HR within 6mo and HR ever bool check)

    for tmpA in NewInfections.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["inf_newInf"] += 1
        if tmpA._everhighrisk_bool:
            rsltdic[tmpA._race][tmpA._SO]["inf_HRever"] += 1
        if tmpA._highrisk_bool:
            rsltdic[tmpA._race][tmpA._SO]["inf_HR6m"] += 1

    # MAT statistics

    # PrEP reason tracker
    for tmpA in totalAgents.iter_agents():
        if tmpA._PrEP_bool:
            rsltdic[tmpA._race][tmpA._SO]["numPrEP"] += 1
            if "IDU" in tmpA._PrEP_reason:
                rsltdic[tmpA._race][tmpA._SO]["iduPartPrep"] += 1
            if "MSMW" in tmpA._PrEP_reason:
                rsltdic[tmpA._race][tmpA._SO]["msmwPartPrep"] += 1
            if "HIV test" in tmpA._PrEP_reason:
                rsltdic[tmpA._race][tmpA._SO]["testedPartPrep"] += 1

    # Newly diagnosed tracker statistics
    for tmpA in NewDiagnosis.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["newlyTested"] += 1

    # Newly HR agents
    for tmpA in newHR.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["newHighRisk"] += 1
        if tmpA._SO == "HM":
            newHR_HM += 1
            if tmpA._HIV_bool:
                newHR_HIV_HM += 1
                if tmpA._AIDS_bool:
                    newHR_AIDS_HM += 1
                if tmpA._tested:
                    newHR_Tested_HM += 1
                    if tmpA._HAART_bool:
                        newHR_ART_HM += 1
        elif tmpA._SO == "HF":
            newHR_HF += 1
            if tmpA._HIV_bool:
                newHR_HIV_HF += 1
                if tmpA._AIDS_bool:
                    newHR_AIDS_HF += 1
                if tmpA._tested:
                    newHR_Tested_HF += 1
                    if tmpA._HAART_bool:
                        newHR_ART_HF += 1

    # Total HIV summary snapshot for timestep
    for tmpA in self.HIV_agentSet.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["numHIV"] += 1
        if tmpA._AIDS_bool:
            rsltdic[tmpA._race][tmpA._SO]["numAIDS"] += 1
        if tmpA._tested:
            rsltdic[tmpA._race][tmpA._SO]["numTested"] += 1
        if tmpA._HAART_bool:
            rsltdic[tmpA._race][tmpA._SO]["numART"] += 1

    for tmpA in totalAgents._subset["DU"]._subset["IDU"].iter_agents():
        if tmpA._HIV_bool:
            rsltdic[tmpA._race]["IDU"]["numHIV"] += 1
        if tmpA._AIDS_bool:
            rsltdic[tmpA._race]["IDU"]["numAIDS"] += 1
        if tmpA._tested:
            rsltdic[tmpA._race]["IDU"]["numTested"] += 1
        if tmpA._HAART_bool:
            rsltdic[tmpA._race]["IDU"]["numART"] += 1

    for tmpA in totalAgents.iter_agents():
        rsltdic[tmpA._race][tmpA._SO]["numAgents"] += 1
        if tmpA._highrisk_type == "postIncar":
            rsltdic[tmpA._race][tmpA._SO]["numHR"] += 1

    for tmpA in deathSet:
        rsltdic[tmpA._race][tmpA._SO]["deaths"] += 1

    deaths_total = (
        deaths["Total"]["HM"] + deaths["Total"]["HF"] + deaths["Total"]["MSM"]
    )
    deaths_HM = deaths["Total"]["HM"]
    deaths_MSM = deaths["Total"]["MSM"]
    deaths_HF = deaths["Total"]["HF"]
    deaths_HIV_total = (
        deaths["HIV+"]["HM"] + deaths["HIV+"]["HF"] + deaths["HIV+"]["MSM"]
    )
    deaths_HIV_HM = deaths["HIV+"]["HM"]
    deaths_HIV_MSM = deaths["HIV+"]["MSM"]
    deaths_HIV_HF = deaths["HIV+"]["HF"]

    W_rslts = rsltdic["WHITE"]
    B_rslts = rsltdic["BLACK"]

    # Sum 'ALL' categories for race/SO bins
    for race in rsltdic:
        for param in rc_template:
            rsltdic[race]["ALL"][param] = (
                rsltdic[race]["MSM"][param]
                + rsltdic[race]["HM"][param]
                + rsltdic[race]["HF"][param]
            )
    for race in rsltdic:
        for param in rc_template:
            tot_rsltdic["ALL"]["ALL"][param] += rsltdic[race]["ALL"][param]
            tot_rsltdic["ALL"]["HM"][param] += rsltdic[race]["HM"][param]
            tot_rsltdic["ALL"]["HF"][param] += rsltdic[race]["HF"][param]
            tot_rsltdic["ALL"]["MSM"][param] += rsltdic[race]["MSM"][param]
            tot_rsltdic["ALL"]["IDU"][param] += rsltdic[race]["IDU"][param]
            rsltdic["ALL"]["ALL"][param] += rsltdic[race]["ALL"][param]
            rsltdic["ALL"]["HM"][param] += rsltdic[race]["HM"][param]
            rsltdic["ALL"]["HF"][param] += rsltdic[race]["HF"][param]
            rsltdic["ALL"]["MSM"][param] += rsltdic[race]["MSM"][param]
            rsltdic["ALL"]["IDU"][param] += rsltdic[race]["IDU"][param]
    for agentRace in ["WHITE", "BLACK", "ALL"]:
        for agentTypes in params.agentPopulations:
            name = "basicReport_" + agentTypes + "_" + agentRace
            tmpReport = open("results/" + name + ".txt", "a")
            tmpReport.write(
                (
                    "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                        self.runseed,
                        self.popseed,
                        self.netseed,
                        t,
                        rsltdic[agentRace][agentTypes]["numAgents"],
                        rsltdic[agentRace][agentTypes]["numHIV"],
                        rsltdic[agentRace][agentTypes]["numAIDS"],
                        rsltdic[agentRace][agentTypes]["numTested"],
                        rsltdic[agentRace][agentTypes]["numART"],
                        rsltdic[agentRace][agentTypes]["numHR"],
                        rsltdic[agentRace][agentTypes]["inf_newInf"],
                        rsltdic[agentRace][agentTypes]["inf_HR6m"],
                        rsltdic[agentRace][agentTypes]["inf_HRever"],
                        rsltdic[agentRace][agentTypes]["newlyTested"],
                        rsltdic[agentRace][agentTypes]["deaths"],
                        rsltdic[agentRace][agentTypes]["numPrEP"],
                        rsltdic[agentRace][agentTypes]["iduPartPrep"],
                        rsltdic[agentRace][agentTypes]["msmwPartPrep"],
                        rsltdic[agentRace][agentTypes]["testedPartPrep"],
                    )
                )
            )
            tmpReport.close()

    for demographicTypes in list(params.DemographicParams.keys()):
        name = "basicReport_" + demographicTypes
        tmpReport = open("results/" + name + ".txt", "a")
        tmpReport.write(
            (
                "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
                % (
                    self.runseed,
                    self.popseed,
                    self.netseed,
                    t,
                    totalAgents._subset["Race"]._subset[demographicTypes].num_members(),
                    rsltdic[demographicTypes]["ALL"]["numHIV"],
                    rsltdic[demographicTypes]["ALL"]["numAIDS"],
                    rsltdic[demographicTypes]["ALL"]["numTested"],
                    rsltdic[demographicTypes]["ALL"]["numART"],
                    rsltdic[demographicTypes]["ALL"]["numHR"],
                    rsltdic[demographicTypes]["ALL"]["inf_newInf"],
                    rsltdic[demographicTypes]["ALL"]["inf_HR6m"],
                    rsltdic[demographicTypes]["ALL"]["inf_HRever"],
                    rsltdic[demographicTypes]["ALL"]["newlyTested"],
                    rsltdic[demographicTypes]["ALL"]["deaths"],
                    rsltdic[demographicTypes]["ALL"]["numPrEP"],
                )
            )
        )
        tmpReport.close()

    incidenceReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            NewInfections.num_members(),
            rsltdic["WHITE"]["HM"]["inf_newInf"],
            rsltdic["BLACK"]["HM"]["inf_newInf"],
            tot_rsltdic["ALL"]["HM"]["inf_newInf"],
            rsltdic["WHITE"]["HF"]["inf_newInf"],
            rsltdic["BLACK"]["HF"]["inf_newInf"],
            tot_rsltdic["ALL"]["HF"]["inf_newInf"],
            rsltdic["WHITE"]["MSM"]["inf_newInf"],
            rsltdic["BLACK"]["MSM"]["inf_newInf"],
            tot_rsltdic["ALL"]["MSM"]["inf_newInf"],
        )
    )

    prevalenceReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            totalAgents.num_members(),
            totalAgents._subset["SO"]._subset["HM"].num_members(),
            totalAgents._subset["SO"]._subset["HF"].num_members(),
            HIVAgents.num_members(),
            tot_rsltdic["ALL"]["HM"]["numHIV"],
            tot_rsltdic["ALL"]["HF"]["numHIV"],
        )
    )

    deathReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            deaths_total,
            deaths_HM,
            deaths_MSM,
            deaths_HF,
            deaths_HIV_total,
            deaths_HIV_HM,
            deaths_HIV_MSM,
            deaths_HIV_HF,
        )
    )

    highriskReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            tot_rsltdic["ALL"]["ALL"]["inf_HRever"],
            tot_rsltdic["ALL"]["HM"]["inf_HRever"],
            tot_rsltdic["ALL"]["HF"]["inf_HRever"],
            tot_rsltdic["ALL"]["ALL"]["inf_HR6m"],
            tot_rsltdic["ALL"]["HM"]["inf_HR6m"],
            tot_rsltdic["ALL"]["HF"]["inf_HR6m"],
        )
    )

    newlyhighriskReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            newHR_HM,
            newHR_HIV_HM,
            newHR_AIDS_HM,
            newHR_Tested_HM,
            newHR_ART_HM,
            newHR_HF,
            newHR_HIV_HF,
            newHR_AIDS_HF,
            newHR_Tested_HF,
            newHR_ART_HF,
        )
    )

    femaleReport.write(
        (
            "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
            % (
                rseed,
                t,
                totalAgents._subset["SO"]._subset["HF"].num_members(),
                tot_rsltdic["ALL"]["HF"]["numHIV"],
                tot_rsltdic["ALL"]["HF"]["numAIDS"],
                tot_rsltdic["ALL"]["HF"]["numTested"],
                tot_rsltdic["ALL"]["HF"]["numART"],
                tot_rsltdic["ALL"]["HF"]["inf_newInf"],
                tot_rsltdic["ALL"]["HF"]["inf_HR6m"],
                tot_rsltdic["ALL"]["HF"]["inf_HRever"],
                tot_rsltdic["ALL"]["HF"]["newlyTested"],
                tot_rsltdic["ALL"]["HF"]["deaths"],
                tot_rsltdic["ALL"]["HF"]["numPrEP"],
            )
        )
    )

    maleReport.write(
        (
            "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
            % (
                rseed,
                t,
                totalAgents._subset["SO"]._subset["HM"].num_members(),
                tot_rsltdic["ALL"]["HM"]["numHIV"],
                tot_rsltdic["ALL"]["HM"]["numAIDS"],
                tot_rsltdic["ALL"]["HM"]["numTested"],
                tot_rsltdic["ALL"]["HM"]["numART"],
                tot_rsltdic["ALL"]["HM"]["inf_newInf"],
                tot_rsltdic["ALL"]["HM"]["inf_HR6m"],
                tot_rsltdic["ALL"]["HM"]["inf_HRever"],
                tot_rsltdic["ALL"]["HM"]["newlyTested"],
                tot_rsltdic["ALL"]["HM"]["deaths"],
                tot_rsltdic["ALL"]["HM"]["numPrEP"],
            )
        )
    )

    msmReport.write(
        (
            "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
            % (
                rseed,
                t,
                totalAgents._subset["SO"]._subset["HM"].num_members(),
                tot_rsltdic["ALL"]["MSM"]["numHIV"],
                tot_rsltdic["ALL"]["MSM"]["numAIDS"],
                tot_rsltdic["ALL"]["MSM"]["numTested"],
                tot_rsltdic["ALL"]["MSM"]["numART"],
                tot_rsltdic["ALL"]["MSM"]["inf_newInf"],
                tot_rsltdic["ALL"]["MSM"]["inf_HR6m"],
                tot_rsltdic["ALL"]["MSM"]["inf_HRever"],
                tot_rsltdic["ALL"]["MSM"]["newlyTested"],
                tot_rsltdic["ALL"]["MSM"]["deaths"],
                tot_rsltdic["ALL"]["MSM"]["numPrEP"],
            )
        )
    )

    nalReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            self.runseed,
            self.popseed,
            self.netseed,
            t,
            Naltrex_NIDU_M,
            Naltrex_NIDU_F,
            Naltrex_IDU_M,
            Naltrex_IDU_F,
            DOC_Naltrex_M,
            DOC_Naltrex_F,
            Prior_Year_Naltrex,
        )
    )
    oatReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            self.runseed,
            self.popseed,
            self.netseed,
            t,
            OAT_NIDU_M,
            OAT_NIDU_F,
            OAT_IDU_M,
            OAT_IDU_F,
            DOC_OAT_M,
            DOC_OAT_F,
            Prior_Year_OAT,
        )
    )

    incarReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            IncarAgents.num_members(),
            rsltdic["WHITE"]["HM"]["incar"],
            rsltdic["BLACK"]["HM"]["incar"],
            rsltdic["WHITE"]["HF"]["incar"],
            rsltdic["BLACK"]["HF"]["incar"],
            rsltdic["WHITE"]["MSM"]["incar"],
            rsltdic["BLACK"]["MSM"]["incar"],
            rsltdic["WHITE"]["ALL"]["incarHIV"],
            rsltdic["BLACK"]["ALL"]["incarHIV"],
            rsltdic["WHITE"]["HM"]["newRelease"],
            rsltdic["WHITE"]["HF"]["newRelease"],
            rsltdic["BLACK"]["HM"]["newRelease"],
            rsltdic["BLACK"]["HF"]["newRelease"],
            rsltdic["WHITE"]["ALL"]["newReleaseHIV"],
            rsltdic["BLACK"]["ALL"]["newReleaseHIV"],
        )
    )

    iduReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            totalAgents._subset["DU"]._subset["IDU"].num_members(),
            tot_rsltdic["ALL"]["IDU"]["numHIV"],
            tot_rsltdic["ALL"]["IDU"]["numAIDS"],
            tot_rsltdic["ALL"]["IDU"]["numART"],
            tot_rsltdic["ALL"]["IDU"]["numTested"],
        )
    )

    whiteReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            rsltdic["WHITE"]["ALL"]["numHIV"],
            rsltdic["WHITE"]["MSM"]["numHIV"],
            rsltdic["WHITE"]["ALL"]["numTested"],
            rsltdic["WHITE"]["ALL"]["numART"],
        )
    )

    blackReport.write(
        "%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            rseed,
            t,
            rsltdic["BLACK"]["ALL"]["numHIV"],
            rsltdic["BLACK"]["MSM"]["numHIV"],
            rsltdic["BLACK"]["ALL"]["numTested"],
            rsltdic["BLACK"]["ALL"]["numART"],
        )
    )

    if t == 0:
        cumulativeI = len(NewInfections._members)
        ResultDict["Inc_c_HM"].update({t: rsltdic["WHITE"]["HM"]["inf_newInf"]})
        ResultDict["Inc_c_HF"].update({t: rsltdic["WHITE"]["HF"]["inf_newInf"]})
    else:
        cumulativeI = 0
        pass

    ResultDict["Inc_c_Tot"].update({t: cumulativeI})

    try:
        ResultDict["Prv_HIV"].update(
            {t: (1.0 * tot_rsltdic["ALL"]["ALL"]["numHIV"] / totalAgents.num_members())}
        )
    except:
        ResultDict["Prv_HIV"].update({t: (0.0)})
    try:
        ResultDict["Prv_AIDS"].update(
            {
                t: (
                    1.0
                    * tot_rsltdic["ALL"]["ALL"]["numAIDS"]
                    / tot_rsltdic["ALL"]["ALL"]["numHIV"]
                )
            }
        )
    except:
        ResultDict["Prv_AIDS"].update({t: (0.0)})
    try:
        ResultDict["Prv_Test"].update(
            {
                t: (
                    1.0
                    * tot_rsltdic["ALL"]["ALL"]["numTested"]
                    / max(tot_rsltdic["ALL"]["ALL"]["numHIV"], 1)
                )
            }
        )
    except:
        ResultDict["Prv_Test"].update({t: (0.0)})
    try:
        ResultDict["Prv_ART"].update(
            {
                t: (
                    1.0
                    * tot_rsltdic["ALL"]["ALL"]["numART"]
                    / tot_rsltdic["ALL"]["ALL"]["numTested"]
                )
            }
        )
    except:
        ResultDict["Prv_ART"].update({t: (0.0)})

    ResultDict["n_Relations"].update({t: Relationships.num_members()})
    ResultDict["Inc_t_HM"].update({t: rsltdic["WHITE"]["HM"]["inf_newInf"]})
    ResultDict["Inc_t_HF"].update({t: rsltdic["WHITE"]["HF"]["inf_newInf"]})

    num_partners = []
    num_partners_hr = []
    ann_num_partners = []
    turnover = []
    ann_turnover = []
