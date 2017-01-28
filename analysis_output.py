__author__ = 'MaximilianKing'
#!/usr/bin/env python

# HIVABM assessment library.
# Needed for Evolution class.
#
# Lars Seemann
# lseemann [at] uh.edu

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

try:
    from agent import *
except ImportError:
    raise ImportError("Can't import Agent class")

def initiate_ResultDict():
    # nested dictionary for results (inner dictionary has the form: time:result)
    #incidenceMatrix = [[[0 for x in range(2)] for x in range(2)] for x in range(2)]
    ResultDict = {  'ND XMS UnDX':{},
                    'ND XMS DX':{},
                    'ND XMS HAART 5':{},
                    'ND XMS HAART Other':{},
                    'PWID XMS UnDX':{},
                    'PWID XMS DX':{},
                    'PWID XMS HAART 5':{},
                    'PWID XMS HAART Other':{}}


    """
    whiteReport = open('Results/W_pop_report.txt', 'w')
    whiteReport.write("t\tTotal-HIV\tMSM\tTested+\tHAART\n")
    whiteReport.write("0\t0\t0\t0\t0\n")
    whiteReport.close()

    blackReport = open('Results/B_pop_report.txt', 'a')
    blackReport.write("t\tTotal-HIV\tMSM\tTested+\tHAART\n")
    blackReport.write("0\t0\t0\t0\t0\n")
    blackReport.close()

    incidenceReport = open('Results/IncidenceReport.txt', 'w')
    incidenceReport.write("t\ttotalIncidence\tIDU-incidence\tAcute-IDU\n")
    incidenceReport.write("0\t0\t0\t0\n")
    incidenceReport.close()


    prevalenceReport = open('Results/PrevalenceReport.txt', 'w')
    prevalenceReport.write("t\tND\tIDU\tHM\tHF\tMSM\n")
    prevalenceReport.write("0\t0\t0\t0\t0\t0\n")
    prevalenceReport.close()

    iduReport = open('Results/iduReport.txt', 'w')
    iduReport.write("t\tTotal-IDU\tIDU-HIV\tIDU-AIDS\tIDU-HAART\tIDU-tested\n")
    iduReport.write("0\t0\t0\t0\t0\t0\n")
    iduReport.close()
    """


    return ResultDict

def print_stats(t, totalAgents, HIVAgents, IncarAgents,PrEPAgents, NewInfections, deaths):

    incidenceReport = open('Results/IncidenceReport.txt', 'a')
    prevalenceReport = open('Results/PrevalenceReport.txt', 'a')
    deathReport = open('Results/DeathReport.txt', 'a')
    incarReport = open('Results/IncarReport.txt', 'a')
    PrEPReport = open('Results/PrEPReport.txt', 'a')
    iduReport = open('Results/iduReport.txt', 'a')
    highriskReport = open('Results/highriskReport.txt', 'a')


    whiteReport = open('Results/W_pop_report.txt', 'a')
    blackReport = open('Results/B_pop_report.txt', 'a')
    num_SEP = 0

    numToHM = 0
    numToHF = 0
    numHIV_HM = 0
    numHIV_HF = 0

    infHM_HRever = 0
    infHM_HR6m = 0
    infHF_HRever = 0
    infHF_HR6m = 0

    for tmpA in NewInfections.iter_agents():
        if tmpA._SO == "HM":
            numToHM +=1
            if tmpA._everhighrisk_bool:
                infHM_HRever += 1
                if tmpA._highrisk_bool:
                    infHM_HR6m += 1

        elif tmpA._SO == "HF":
            numToHF += 1
            if tmpA._everhighrisk_bool:
                infHF_HRever += 1
                if tmpA._highrisk_bool:
                    infHF_HR6m += 1

    for tmpA in HIVAgents.iter_agents():
        if tmpA._SO == "HM":
            numHIV_HM += 1
        elif tmpA._SO == "HF":
            numHIV_HF += 1

    deaths_total = deaths["Total"]["HM"]+deaths["Total"]["HF"]+deaths["Total"]["MSM"]
    deaths_HM = deaths["Total"]["HM"]
    deaths_HF = deaths["Total"]["HF"]
    deaths_HIV_total = deaths["HIV+"]["HM"]+deaths["HIV+"]["HF"]+deaths["HIV+"]["MSM"]
    deaths_HIV_HM = deaths["HIV+"]["HM"]
    deaths_HIV_HF = deaths["HIV+"]["HF"]

    incidenceReport.write("%d\t%d\t%d\t%d\n" % (t, NewInfections.num_members(),numToHM,numToHF))
    prevalenceReport.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (t, totalAgents.num_members(), totalAgents._subset["HM"].num_members(), totalAgents._subset["HF"].num_members(), HIVAgents.num_members(), numHIV_HM, numHIV_HF))
    deathReport.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (t, deaths_total, deaths_HM, deaths_HF, deaths_HIV_total, deaths_HIV_HM,deaths_HIV_HF))
    highriskReport.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (t, infHM_HRever+infHF_HRever, infHM_HRever, infHF_HRever, infHM_HR6m+infHF_HR6m, infHM_HR6m, infHF_HR6m))
    #iduReport.write("%d\t%d\t%d\t%d\t%d\t%d\n" % (t, numIDU, numHIV_IDU_p, numAIDS_IDU_p, numHIDU, numTested_IDU))

    #whiteReport.write("%d\t%d\t%d\t%d\t%d\n" % (t, numHIV_White, numHIV_MSM_White, numTested_White, numHAART_White))
    #blackReport.write("%d\t%d\t%d\t%d\t%d\n" % (t, numHIV_Black, numHIV_MSM_Black, numTested_Black, numHAART_Black))

    num_partners = []
    num_partners_hr = []
    ann_num_partners = []
    turnover = []
    ann_turnover = []


def assess_before_update(t,
                         ResultDict,
                         Agents,
                         HIV_agents, tmp_HIV_agents,
                         AIDS_agents,tmp_AIDS_agents,
                         SEPAgents,
                         HAART_agents, tmp_HAART_agents,
                         AdherenceAgents,
                         num_Deaths,
                         AdjMats_by_time,
                         Acute_agents,
                         Transmit_from_agents,
                         Transmit_to_agents,
                         Transmission_tracker,
                         high_risk_agents,
                         Incarcerated_agents,
                         HIVIdentified_agents):
    """
    Assess population before update.
    Called in HIVABM Evolution class.
    Usage: assess_before_update(time,
                                self.Agents,
                                self.HIV_agents,
                                self.tmp_HIV_agents,
                                self.AIDS_agents,
                                self.tmp_AIDS_agents,
                                self.SEPAgents,
                                self.HAART_agents,
                                self.tmp_HAART_agents,
                                self.AdherenceAgents,
                                self.num_Deaths,
                                self.adjMat,
                                self.Acute_agents,
                                self.Transmit_from_agents)
    Output: ResultDict
    """

    incidenceReport = open('Results/IncidenceReport.txt', 'a')
    prevalenceReport = open('Results/PrevalenceReport.txt', 'a')
    deathReport = open('Results/DeathReport.txt', 'a')
    iduReport = open('Results/iduReport.txt', 'a')

    whiteReport = open('Results/W_pop_report.txt', 'a')
    blackReport = open('Results/B_pop_report.txt', 'a')
    num_SEP = 0


    SexTrans = {'0':{}, '1':{}}
    SexTrans['0'] = {'0':0.005, '1':0.005, '2':0.004, '3':0.002, '4':0.001, '5':0.0001}
    SexTrans['1'] = {'0':0.001, '1':0.001, '2':0.0008, '3':0.0004, '4':0.0002, '5':0.0001}

    NeedleTrans = {'0':0.007, '1':0.007, '2':0.0056, '3':0.0028, '4':0.0014, '5':0.0002}


    transmissionTotalCounter = {'0':SexTrans, '1':NeedleTrans}
    transmissionAcuteCounter = {'SEX':SexTrans, 'NEEDLE':NeedleTrans}

    demographicMatrix = [[[0 for x in range(2)] for x in range(2)] for x in range(2)]
    incidenceMatrix = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    transmittedToMatrix = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    transmittedIDUArray = 0
    transmittedAcutedArray = 0
    numHIV_MSM = 0
    numPWID = 0
    numIDU = 0
    numHIV_IDU_i = 0
    numHIV_IDU_p = 0
    numHIV_White = 0
    numHIV_Black = 0
    numAIDS_IDU_i = 0
    numAIDS_IDU_p = 0
    numHIDU = 0
    numHIDUa = 0
    numDIDU = 0

    numToND = 0
    numToIDU = 0
    numToHM = 0
    numToHF = 0
    numToMSM = 0

    num_White = 0
    numMSM_White = 0
    numHM_White = 0
    numHF_White = 0
    numPWID_White = 0
    numHIV_MSM_White = 0

    num_Black = 0
    numMSM_Black = 0
    numHM_Black = 0
    numHF_Black = 0
    numPWID_Black = 0
    numHIV_MSM_Black = 0

    numMSM = 0
    numWSW = 0
    numHF = 0
    numHM = 0

    numMSM_Black = 0

    numACUTE_IDU = 0
    numACUTE_NIDU = 0
    numACUTE_ND = 0
    numACUTE_HM = 0
    numACUTE_HF = 0
    numACUTE_MSM = 0
    numACUTE_WSW = 0

    numFROM_IDU = 0
    numFROM_NIDU = 0
    numFROM_ND = 0
    numFROM_HM = 0
    numFROM_HF = 0
    numFROM_MSM = 0
    numFROM_WSW = 0

    num_Incar_White = 0
    numMSM_Incar_White = 0
    numHM_Incar_White = 0
    numHF_Incar_White = 0
    numPWID_Incar_White = 0

    num_Incar_Black = 0
    numMSM_Incar_Black = 0
    numHM_Incar_Black = 0
    numHF_Incar_Black = 0
    numPWID_Incar_Black = 0

    numTested_White = 0
    numHAART_White = 0

    numTested_Black = 0
    numHAART_Black = 0


    numTested_IDU = 0

    num_SEP = len(list(set(SEPAgents.keys())))   #number of IDU in SEP
    #print Acute_agents
    for agent in Agents:
        agent_sex_type = Agents[agent]['Sex Type']
        agent_drug_type = Agents[agent]['Drug Type']
        if agent_drug_type == 'IDU':
            agent_drug_bool = 1
        else:
            agent_drug_bool = 0

        agent_race_type = Agents[agent]['Race']
        agent_HAART_adh = AdherenceAgents[agent]
        # Booleans:
        agent_HIV_bool = Agents[agent]['HIV']#agent in HIV_agents
        agent_HIV_tmp_bool = 0#agent in tmp_HIV_agents
        agent_AIDS_bool = Agents[agent]['AIDS']#agent in AIDS_agents
        agent_AIDS_tmp_bool = 0#agent in tmp_AIDS_agents
        agent_HAART_bool = Agents[agent]['HAARTa']#agent in HAART_agents
        agent_HAART_tmp_bool = 0#agent in tmp_HAART_agents
        agent_Incar_bool = 0#agent in Incarcerated_agents
        agent_Test_bool = Agents[agent]['Tested']#agent in HIVIdentified_agents

        #print "DT: %d\tST: %d\t Test: %d"%(agent_drug_bool,agent_Test_bool, agent_HAART)


        #transmissionTotalCounter[agent_drug_bool][agent_sex_type][agent_Test_bool][agent_HAART_bool]
        #transmittedTotalArray[agent_Test_bool][agent_HIV_bool][agent_HAART_bool] += 1
        #print agent_HAART_adh
        if agent in Acute_agents:
            if agent_drug_type == 'ND': numACUTE_ND += 1
            elif agent_drug_type == 'IDU': numACUTE_IDU += 1
            elif agent_drug_type == 'NIDU': numACUTE_NIDU += 1
            if agent_sex_type=='HM': numACUTE_HM += 1
            elif agent_sex_type=='HF': numACUTE_HF += 1
            elif agent_sex_type=='MSM': numACUTE_MSM += 1
            elif agent_sex_type=='WSW': numACUTE_WSW += 1

        if agent in Transmit_from_agents:
            #Determine adh class. 0: nohaart, 1: haart w/ adh 1-4, 2: HAART w/ adh 5
            if agent_HAART_adh == 5:
                agent_HAART_class = 2
            elif agent_HAART_adh > 0:
                agent_HAART_class = 1
            elif agent_HAART_adh == 0:
                agent_HAART_class = 0
            else:
                print "BROKE"

            incidenceMatrix[agent_drug_bool][agent_Test_bool][agent_HAART_class] += 1
            if agent_drug_type == 'ND': numFROM_ND += 1
            elif agent_drug_type == 'IDU': numFROM_IDU += 1
            elif agent_drug_type == 'NIDU': numFROM_NIDU += 1

            if agent_sex_type=='HM': numFROM_HM += 1
            elif agent_sex_type=='HM': numFROM_HF += 1
            elif agent_sex_type=='MSM': numFROM_MSM += 1
            elif agent_sex_type=='WSW': numFROM_WSW += 1


        elif agent in Transmit_to_agents:
            if agent_drug_type == 'ND': numToND += 1
            elif agent_drug_type == 'IDU': numToIDU += 1

            if agent_sex_type=='HM': numToHM += 1
            elif agent_sex_type=='HF': numToHF += 1
            elif agent_sex_type=='MSM': numToMSM += 1


        #if agent_sex_type=='MSM': numMSM += 1
        #elif agent_sex_type=='WSW': numWSW += 1
        #elif agent_sex_type=='HF': numHF += 1
        #elif agent_sex_type=='HM': numHM += 1
        if agent_drug_type == 'IDU':
            numIDU += 1

            if agent_HIV_bool:
                numHIV_IDU_p += 1

                if agent_AIDS_bool:
                    numAIDS_IDU_p +=1

                if agent_HAART_bool:
                    numHIDU += 1

                if agent_Test_bool:
                    numTested_IDU +=1



        if agent_race_type=='WHITE':

            num_White +=1
            if agent_sex_type=='MSM':
                numMSM_White +=1
                numMSM +=1


                if agent_Incar_bool:
                    numMSM_Incar_White +=1
                    num_Incar_White +=1

            elif agent_sex_type=='HM':
                numHM_White +=1
                numHM +=1

                if agent_Incar_bool:
                    numHM_Incar_White +=1
                    num_Incar_White +=1

            elif agent_sex_type=='HF':
                numHF_White +=1
                numHF +=1

                if agent_Incar_bool:
                    numHF_Incar_White +=1
                    num_Incar_White +=1

            if agent_sex_type=='PWID':
                numPWID_White +=1
                numPWID +=1

                if agent_Incar_bool:
                    numPWID_Incar_White +=1

            if agent_HIV_bool:
                numHIV_White +=1
                if agent_sex_type == 'MSM':
                    numHIV_MSM_White +=1

            if agent_HAART_bool:
                numHAART_White +=1
            if agent_Test_bool:
                numTested_White +=1

        elif agent_race_type=='BLACK':

            num_Black +=1
            if agent_sex_type=='MSM':
                numMSM_Black +=1
                numMSM +=1


                if agent_Incar_bool:
                    numMSM_Incar_Black +=1
                    num_Incar_Black +=1

            elif agent_sex_type=='HM':
                numHM_Black +=1
                numHM +=1

                if agent_Incar_bool:
                    numHM_Incar_Black +=1
                    num_Incar_Black +=1

            elif agent_sex_type=='HF':
                numHF_Black +=1
                numHF +=1

                if agent_Incar_bool:
                    numHF_Incar_Black +=1
                    num_Incar_Black +=1

            if agent_sex_type=='PWID':
                numPWID_Black +=1
                numPWID +=1

                if agent_Incar_bool:
                    numPWID_Incar_Black +=1

            if agent_HIV_bool:
                numHIV_Black +=1
                if agent_sex_type == 'MSM':
                    numHIV_MSM_Black +=1

            if agent_HAART_bool:
                numHAART_Black +=1
            if agent_Test_bool:
                numTested_Black +=1

    """
    #PLOTTING FOR RUN
    #plt.ion()
    plt.subplot(2,2,1)
    plt.plot(t,numHIV_White,'rs', color='cornflowerblue', linewidth=3)
    plt.title('Vampire Infection')
    plt.ylabel('Vampire Population')
    plt.subplot(2,1,2)
    plt.plot(t,numHIV_MSM,'bs')
    plt.xlabel('time (hours)')
    plt.ylabel('Human Population')
    plt.draw()

    #plt.show()
    plt.pause(.1)
    """
    numNEEDLE = Transmission_tracker['NEEDLE'][t]

    incidenceReport.write("%d\t%d\t%d\t%d\n" % (t, numFROM_ND + numFROM_IDU + numFROM_NIDU, numFROM_IDU, numACUTE_IDU)) #3rd was numHIV_IDU_MSM_i + numHIV_ND_MSM_i
    prevalenceReport.write("%d\t%d\t%d\t%d\t%d\t%d\n" % (t, numToND, numToIDU, numToHM, numToHF, numToMSM))
    iduReport.write("%d\t%d\t%d\t%d\t%d\t%d\n" % (t,numIDU, numHIV_IDU_p, numAIDS_IDU_p, numHIDU, numTested_IDU))

    whiteReport.write("%d\t%d\t%d\t%d\t%d\n" % (t, numHIV_White, numHIV_MSM_White, numTested_White, numHAART_White))
    blackReport.write("%d\t%d\t%d\t%d\t%d\n" % (t, numHIV_Black, numHIV_MSM_Black, numTested_Black, numHAART_Black))



    num_partners = []
    num_partners_hr = []
    ann_num_partners = []
    turnover = []
    ann_turnover = []

    """
    for agent in HIV_agents:#Agents:
        current_partners = list(AdjMats_by_time[t].rows[agent])
        last_month_partners = list(AdjMats_by_time[t-1].rows[agent])
        if t > 12:
            last_year_partners = list(AdjMats_by_time[t-12].rows[agent])
            whole_year_partners = []
            for tt in range(t-11,t+1):
                whole_year_partners += list(AdjMats_by_time[tt].rows[agent])
            whole_year_partners = list(set(whole_year_partners))
        else:
            last_year_partners = []
            whole_year_partners = []

        new_partners = [p for p in current_partners if not p in last_month_partners]
        ann_new_partners = [p for p in current_partners if not p in last_year_partners]

        num_partners += [len(current_partners)]
        if agent in high_risk_agents: num_partners_hr += [len(current_partners)]
        ann_num_partners += [len(whole_year_partners)]

        if len(current_partners) > 0:
            turnover += [float(len(new_partners)) / len(current_partners)]
            ann_turnover += [float(len(ann_new_partners)) / len(current_partners)]


    if t == 1 or t == 12:
    	bc = np.bincount(num_partners)
        print np.mean(num_partners)
    	plt.bar(np.arange(len(bc)), bc, width=.5)
    	plt.show()
    """

    #ResultDict['Number Agents'].update({t:num_White})
    #ResultDict['Total Incidence'].update({t:numMSM_White})
    #ResultDict['Total IDU Incidence'].update({t:incidenceMatrix[1][0][1]})

    ResultDict['ND XMS UnDX'].update({t:incidenceMatrix[0][0][0]})
    ResultDict['ND XMS DX'].update({t:incidenceMatrix[0][1][0]})
    ResultDict['ND XMS HAART 5'].update({t:incidenceMatrix[0][1][2]})
    ResultDict['ND XMS HAART Other'].update({t:incidenceMatrix[0][1][1]})
    ResultDict['PWID XMS UnDX'].update({t:incidenceMatrix[1][0][0]})
    ResultDict['PWID XMS DX'].update({t:incidenceMatrix[1][1][0]})
    ResultDict['PWID XMS HAART 5'].update({t:incidenceMatrix[1][1][2]})
    ResultDict['PWID XMS HAART Other'].update({t:incidenceMatrix[1][1][1]})
    #ResultDict['HAART Proportion'].update({t:np.mean(ann_turnover)})
    #ResultDict['DX Proportion'].update({t:Transmission_tracker['SEX_MSM'][t]})
    #ResultDict['SEX_NMSM_transmissions'].update({t:Transmission_tracker['SEX_NMSM'][t]})
    #ResultDict['NEEDLE_transmissions'].update({t:Transmission_tracker['NEEDLE'][t]})
    return ResultDict

def assess_interaction_distribution(Agents,
                                    agent, agent_drug_type, agent_sex_type,
                                    partners,
                                    tmp_ND_NumPartners_Count,
                                    tmp_NIDU_NumPartners_Count,
                                    tmp_IDU_NumPartners_Count,
                                    tmp_MSM_NumPartners_Count,
                                    tmp_WSW_NumPartners_Count):
    """
    :Purpose:
        Assess distribution of connection type.
    :Usage:
        assess_interaction_distribution(self.Agents,
                                        agent, agent_drug_type, agent_sex_type,
                                        partners,
                                        self.tmp_ND_NumPartners_Count,
                                        self.tmp_NIDU_NumPartners_Count,
                                        self.tmp_IDU_NumPartners_Count,
                                        self.tmp_MSM_NumPartners_Count,
                                        self.tmp_WSW_NumPartners_Count
                                        )

    """
    if agent_drug_type == 'ND':
        for p in partners:
            p_drug_type=Agents[p]['Drug Type']
            tmp_ND_NumPartners_Count.setdefault(agent,
                {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[p_drug_type] += 1
            if Agents[p]['Sex Type']=='MSM':
                tmp_ND_NumPartners_Count[agent]['MSM'] += 1
    elif agent_drug_type == 'NIDU':
        for p in partners:
            p_drug_type=Agents[p]['Drug Type']
            tmp_NIDU_NumPartners_Count.setdefault(agent,
                                                       {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[p_drug_type] += 1
            if Agents[p]['Sex Type']=='MSM':
                tmp_NIDU_NumPartners_Count[agent]['MSM'] += 1
    elif agent_drug_type == 'IDU':
        for p in partners:
            p_drug_type=Agents[p]['Drug Type']
            tmp_IDU_NumPartners_Count.setdefault(agent,
                                                      {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[p_drug_type] += 1
            if Agents[p]['Sex Type']=='MSM':
                tmp_IDU_NumPartners_Count[agent]['MSM'] += 1
    else:
        raise ValueError("Invalid drug type! %s"%agent_drug_type)

    if agent_sex_type == 'MSM':
        for p in partners:
            p_drug_type=Agents[p]['Drug Type']
            tmp_MSM_NumPartners_Count.setdefault(agent,
                                                      {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[p_drug_type] += 1
    elif agent_sex_type == 'WSW':
        for p in partners:
            p_drug_type=Agents[p]['Drug Type']
            tmp_WSW_NumPartners_Count.setdefault(agent,
                                                      {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[p_drug_type] += 1
    elif agent_sex_type not in ['HM','HF']:
        raise ValueError("Invalid sex type! %s"%agent_sex_type)

    for p in partners:
        partner_drug_type = Agents[p]['Drug Type']
        if partner_drug_type == 'ND':
            tmp_ND_NumPartners_Count.setdefault(p,
                                                     {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[agent_drug_type] += 1
            if agent_sex_type=='MSM':
                tmp_ND_NumPartners_Count[p]['MSM'] += 1
        elif partner_drug_type == 'NIDU':
            tmp_NIDU_NumPartners_Count.setdefault(p,
                                                       {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[agent_drug_type] += 1
            if agent_sex_type=='MSM':
                tmp_NIDU_NumPartners_Count[p]['MSM'] += 1
        elif partner_drug_type == 'IDU':
            tmp_IDU_NumPartners_Count.setdefault(p,
                                                      {'ND':0, 'NIDU':0, 'IDU':0, 'MSM':0})[agent_drug_type] += 1
            if agent_sex_type=='MSM':
                tmp_IDU_NumPartners_Count[p]['MSM'] += 1
        else:
            raise ValueError("Invalid partner drug type! %s"%
                             partner_drug_type)

    return [tmp_ND_NumPartners_Count,
            tmp_NIDU_NumPartners_Count,
            tmp_IDU_NumPartners_Count,
            tmp_MSM_NumPartners_Count,
            tmp_WSW_NumPartners_Count]
