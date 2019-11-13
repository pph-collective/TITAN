# TODO convert to pytest, update

# class TestClassMethods(unittest.TestCase):
#     """
#     :Purpose:
#         unittest
#     """
#
#     def setUp(self):
#         """
#         :Purpose:
#             Tests that all models from setup pass inspection. ``setUp`` is perfomed before each method.
#         """
#         self.N_pop = 10000
#
#     def test_SexType(self):
#         """ Test: Testing consistency of Sex type agents"""
#         print("\t__Testing HM agents list")
#         myPopulation = PopulationClass(n=self.N_pop)
#         for a in list(myPopulation.Agents.keys()):
#             SexType = myPopulation.get_agent_characteristic(a, "Sex Type")
#             if SexType == "HM":
#                 self.assertTrue(a in myPopulation.HM_agents)
#             elif SexType == "HF":
#                 self.assertTrue(a in myPopulation.HF_agents)
#             elif SexType == "MSM":
#                 self.assertTrue(a in myPopulation.MSM_agents)
#             elif SexType == "WSW":
#                 self.assertTrue(a in myPopulation.WSW_agents)
#
#         for agent in myPopulation.HM_agents:
#             self.assertTrue(myPopulation.get_agent_characteristic(agent, "Sex Type") == "HM")
#         for agent in myPopulation.HF_agents:
#             self.assertTrue(myPopulation.get_agent_characteristic(agent, "Sex Type") == "HF")
#         for agent in myPopulation.MSM_agents:
#             self.assertTrue(myPopulation.get_agent_characteristic(agent, "Sex Type") == "MSM")
#         for agent in myPopulation.WSW_agents:
#             self.assertTrue(myPopulation.get_agent_characteristic(agent, "Sex Type") == "WSW")
#
#     def test_HIV(self):
#         """ Test: Testing HIV agent array"""
#         print("\t__Testing the HIV agent list")
#         tmpCount = 0
#         myPopulation = PopulationClass(n=self.N_pop)
#         for a in list(myPopulation.Agents.keys()):
#             HIVstatus = myPopulation.get_agent_characteristic(a, "HIV")
#             if HIVstatus != 0:
#                 tmpCount += 1
#                 self.assertTrue(a in myPopulation.HIV_agents)
#         self.assertTrue(len(myPopulation.HIV_agents) == tmpCount)
#
#     def test_AIDS(self):
#         """ Test: Testing AIDS agent array"""
#         print("\t__Testing the AIDS agent list")
#         tmpCount = 0
#         myPopulation = PopulationClass(n=self.N_pop)
#         for a in list(myPopulation.Agents.keys()):
#             AIDSstatus = myPopulation.get_agent_characteristic(a, "AIDS")
#             if AIDSstatus != 0:
#                 tmpCount += 1
#                 self.assertTrue(a in myPopulation.AIDS_agents)
#                 self.assertTrue(a in myPopulation.HIV_agents)
#         self.assertTrue(len(myPopulation.AIDS_agents) == tmpCount)
#
#     def test_consistency(self):
#         """ Test: Testing consistency"""
#         print("\t__Testing consistency of agent lists")
#         myPop = PopulationClass(n=self.N_pop)
#         NormalAgents = list(
#             set(range(myPop.PopulationSize)).difference(
#                 set(myPop.IDU_agents).union(set(myPop.MSM_agents))
#             )
#         )
#         MSMAgents = list(set(myPop.MSM_agents).difference(set(myPop.IDU_agents)))
#         IDUagents = myPop.IDU_agents
#         NumAgents = len(NormalAgents) + len(MSMAgents) + len(IDUagents)
#         self.assertTrue(
#             NumAgents == self.N_pop,
#             "NumAgents=%d, \
# 				PopulationSize = %d"
#             % (NumAgents, self.N_pop),
#         )
#
#     def test_MSM(self):
#         """ Test: Testing MSM agent array"""
#         print("\t__Testing the MSM agent list")
#         tmpCount = 0
#         myPopulation = PopulationClass(n=self.N_pop)
#         for a in list(myPopulation.Agents.keys()):
#             SEXstatus = myPopulation.get_agent_characteristic(a, "Sex Type")
#             if SEXstatus == "MSM":
#                 tmpCount += 1
#                 self.assertTrue(a in myPopulation.MSM_agents)
#         self.assertTrue(len(myPopulation.MSM_agents) == tmpCount)
#
#         for agent in myPopulation.MSM_agents:
#             agent_sex_type = myPopulation.get_agent_characteristic(agent, "Sex Type")
#             self.assertTrue(agent_sex_type == "MSM")
#
#     def test_IDU(self):
#         """ Test: Testing IDU agent array"""
#         print("\t__Testing the IDU agent list")
#         tmpCount = 0
#         myPopulation = PopulationClass(n=self.N_pop)
#         for agent in list(myPopulation.Agents.keys()):
#             agent_drug_status = myPopulation.get_agent_characteristic(agent, "Drug Type")
#             if agent_drug_status == "IDU":
#                 tmpCount += 1
#                 self.assertTrue(agent in myPopulation.IDU_agents)
#         self.assertTrue(len(myPopulation.IDU_agents) == tmpCount)
#
#         for agent in myPopulation.IDU_agents:
#             agent_drug_type = myPopulation.get_agent_characteristic(agent, "Drug Type")
#             self.assertTrue(agent_drug_type == "IDU")
#
#     def test_NIDU(self):
#         """ Test: Testing INDU agent array"""
#         print("\t__Testing the NIDU agent list")
#         tmpCount = 0
#         myPopulation = PopulationClass(n=self.N_pop)
#         for agent in list(myPopulation.Agents.keys()):
#             agent_drug_status = myPopulation.get_agent_characteristic(agent, "Drug Type")
#             if agent_drug_status == "NIDU":
#                 tmpCount += 1
#                 self.assertTrue(agent in myPopulation.NIDU_agents)
#         self.assertTrue(len(myPopulation.NIDU_agents) == tmpCount)
#
#         for agent in myPopulation.NIDU_agents:
#             agent_drug_type = myPopulation.get_agent_characteristic(agent, "Drug Type")
#             self.assertTrue(agent_drug_type == "NIDU")
#
#     def test_Population(self):
#         """ Test: Testing the population"""
#         print("\t__Testing the population")
#         myPopulation = PopulationClass(n=self.N_pop)
#         for agent in list(myPopulation.Agents.keys()):
#             tmp_DrugType = myPopulation.get_agent_characteristic(agent, "Drug Type")
#             self.assertTrue(tmp_DrugType in ["NIDU", "IDU", "ND"])
#             if tmp_DrugType == "NIDU":
#                 self.assertTrue(agent in myPopulation.NIDU_agents)
#             elif tmp_DrugType == "IDU":
#                 self.assertTrue(agent in myPopulation.IDU_agents)
#             else:
#                 self.assertTrue(agent in myPopulation.ND_agents)
#         self.assertTrue(len(myPopulation.NIDU_agents) == myPopulation.numNIDU)
#         self.assertTrue(len(myPopulation.IDU_agents) == myPopulation.numIDU)
#         self.assertTrue(len(myPopulation.ND_agents) == myPopulation.numND)
