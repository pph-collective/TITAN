import pytest

from src.ABM_core import *

# #REVIEW not used anywhere, but may be useful in writing tests?
# def _check_population(self):
#     """
#     :Purpose:
#         Check consistency of population.
#         Only called in unittest.
#
#     """
#
#     # Check consistency of last partners
#     if not (
#         np.all(self.AdjMat.sum(0) == self.AdjMat.conj().sum(0))
#         and np.all(self.AdjMat.sum(1) == self.AdjMat.conj().sum(1))
#     ):
#         raise ValueError("Adjacency matrix not symmetric!")
#
#     # Check consistency of real population
#     count_HF = 0
#     count_HM = 0
#     count_MSM = 0
#     count_WSW = 0
#     count_ND = 0
#     count_NIDU = 0
#     count_IDU = 0
#     count_HIV = 0
#     count_AIDS = 0
#     for (agent, d) in list(self.Agents.items()):
#         agent_dict = d
#         # Sex type
#         sex_type = agent_dict["Sex Type"]
#         if sex_type == "HF":
#             if agent not in self.HF_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agents HF Sex type %d" % agent)
#             else:
#                 count_HF += 1
#         elif sex_type == "HM":
#             if agent not in self.HM_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agents HM Sex type %d" % agent)
#             else:
#                 count_HM += 1
#         elif sex_type == "MSM":
#             if agent not in self.MSM_agents:
#                 raise ValueError("Check agents MSM Sex type %d" % agent)
#             else:
#                 count_MSM += 1
#         elif sex_type == "WSW":
#             if agent not in self.WSW_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agents WSW Sex type %d" % agent)
#             else:
#                 count_WSW += 1
#         else:
#             raise ValueError("Invalid sex type %s" % str(sex_type))
#
#         # Drug type
#         drug_type = agent_dict["Drug Type"]
#         if drug_type == "ND":
#             if agent not in self.ND_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agents ND Drug type %d" % agent)
#             else:
#                 count_ND += 1
#         elif drug_type == "NIDU":
#             if agent not in self.NIDU_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agents NIDU Drug type %d" % agent)
#             else:
#                 count_NIDU += 1
#         elif drug_type == "IDU":
#             if agent not in self.IDU_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agents IDU Drug type %d" % agent)
#             else:
#                 count_IDU += 1
#         else:
#             raise ValueError("Invalid drug type %s" % str(drug_type))
#
#         # HIV
#         HIVstatus = agent_dict["HIV"]
#         if HIVstatus != 0:
#             if agent not in self.HIV_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agent HIV %d" % agent)
#             else:
#                 count_HIV += 1
#         # AIDS
#         AIDSstatus = agent_dict["AIDS"]
#         if AIDSstatus != 0:
#             if agent not in self.AIDS_agents:
#                 print((self.Agents[agent]))
#                 raise ValueError("Check agent AIDS %d" % agent)
#             else:
#                 count_AIDS += 1
#
#     if len(self.HF_agents) != count_HF:
#         raise ValueError("self.HF agents contains too many agents!")
#     if len(self.HM_agents) != count_HM:
#         print(("len(self.HM_agents)=%d" % len(self.HM_agents)))
#         print(("count_HM=%d" % count_HM))
#         raise ValueError("self.HM agents contains too many agents!")
#     if len(self.MSM_agents) != count_MSM:
#         raise ValueError("self.MSM agents contains too many agents!")
#     if len(self.WSW_agents) != count_WSW:
#         raise ValueError("self.WSW agents contains too many agents!")
#
#     if len(self.NIDU_agents) != count_NIDU:
#         raise ValueError("self.NIDU_agents contains too many agents!")
#     if len(self.ND_agents) != count_ND:
#         raise ValueError("self.ND agents contains too many agents!")
#     if len(self.IDU_agents) != count_IDU:
#         mssg = "self.IDU agents contains too many agents!\
#                 \nlen(self.IDU_agents)=%d\ncount_IDU=%d\n"
#         raise ValueError(mssg % (len(self.IDU_agents), count_IDU))
#
#     if len(self.HIV_agents) != count_HIV:
#         raise ValueError(
#             "self.HIV_agents contains too many agents!\
#             \nlen(self.HIV_agents) = %d\ncount_HIV = %d\n"
#             % (len(self.HIV_agents), count_HIV)
#         )
#     if len(self.AIDS_agents) != count_AIDS:
#         raise ValueError("self.AIDS agents contains too many agents!")
#
#     # Check consistency of tmp population
#     count_HF = 0
#     count_HM = 0
#     count_MSM = 0
#     count_WSW = 0
#     count_ND = 0
#     count_NIDU = 0
#     count_IDU = 0
#     count_HIV = 0
#     count_AIDS = 0
#     for (agent, d) in list(self.tmp_Agents.items()):
#         agent_dict = d
#         # Sex type
#         sex_type = agent_dict["Sex Type"]
#         if sex_type == "HF":
#             if agent not in self.tmp_HF_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Sex type %d" % agent)
#             else:
#                 count_HF += 1
#         elif sex_type == "HM":
#             if agent not in self.tmp_HM_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Sex type %d" % agent)
#             else:
#                 count_HM += 1
#         elif sex_type == "MSM":
#             if agent not in self.tmp_MSM_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Sex type %d" % agent)
#             else:
#                 count_MSM += 1
#         elif sex_type == "WSW":
#             if agent not in self.tmp_WSW_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Sex type %d" % agent)
#             else:
#                 count_WSW += 1
#         else:
#             raise ValueError("Invalid sex type %s" % str(sex_type))
#
#         # Drug type
#         drug_type = agent_dict["Drug Type"]
#         if drug_type == "ND":
#             if agent not in self.tmp_ND_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Drug type %d" % agent)
#             else:
#                 count_ND += 1
#         elif drug_type == "NIDU":
#             if agent not in self.tmp_NIDU_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Drug type %d" % agent)
#             else:
#                 count_NIDU += 1
#         elif drug_type == "IDU":
#             if agent not in self.tmp_IDU_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agents Drug type %d" % agent)
#             else:
#                 count_IDU += 1
#         else:
#             raise ValueError("Invalid drug type %s" % str(drug_type))
#
#         # HIV
#         HIVstatus = agent_dict["HIV"]
#         if HIVstatus != 0:
#             if agent not in self.tmp_HIV_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check tmp_agent HIV %d" % agent)
#             else:
#                 count_HIV += 1
#         # AIDS
#         AIDSstatus = agent_dict["AIDS"]
#         if AIDSstatus != 0:
#             if agent not in self.tmp_AIDS_agents:
#                 print((self.tmp_Agents[agent]))
#                 raise ValueError("Check agent AIDS %d" % agent)
#             else:
#                 count_AIDS += 1
#
#     if len(self.tmp_HF_agents) != count_HF:
#         raise ValueError("self.tmp_HF agents contains too many agents!")
#     if len(self.tmp_HM_agents) != count_HM:
#         raise ValueError("self.tmp_HM agents contains too many agents!")
#     if len(self.tmp_MSM_agents) != count_MSM:
#         raise ValueError("self.tmp_MSM agents contains too many agents!")
#     if len(self.tmp_WSW_agents) != count_WSW:
#         raise ValueError("self.tmp_WSW agents contains too many agents!")
#
#     if len(self.tmp_NIDU_agents) != count_NIDU:
#         raise ValueError("self.tmp_NIDU_agents contains too many agents!")
#     if len(self.tmp_ND_agents) != count_ND:
#         raise ValueError("self.tmp_ND agents contains too many agents!")
#     if len(self.tmp_IDU_agents) != count_IDU:
#         mssg = "self.tmp_IDU agents contains too many agents!\
#                 \nlen(self.tmp_IDU_agents)=%d\ncount_IDU=%d\n"
#         raise ValueError(mssg % (len(self.IDU_agents), count_IDU))
#
#     if len(self.tmp_HIV_agents) != count_HIV:
#         raise ValueError("self.tmp_HIV_agents contains too many agents!")
#     if len(self.tmp_AIDS_agents) != count_AIDS:
#         print(("len(self.tmp_AIDS_agents)=%d" % len(self.tmp_AIDS_agents)))
#         print(("count_AIDS=%d" % count_AIDS))
#         raise ValueError("self.tmp_AIDS agents contains too many agents!")
