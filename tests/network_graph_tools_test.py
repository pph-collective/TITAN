import pytest

from src.network_graph_tools import *

n_pop = 10

def test_NormalAgents():
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    myNetworkObj = NetworkClass(N=n_pop, m_0=3)
    assert n_pop == myNetworkObj.All_agentSet.num_members()
    for agent in myNetworkObj.All_agentSet.get_agents():
        agent_sex_type = agent._SO
        agent_drug_type = agent._DU
        assert agent_drug_type in ["IDU", "NIDU", "NDU"]
        assert agent_sex_type in ["MSM", "MSW", "HM", "HF"]


def test_PartialNetwork():
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    myNetworkObj = NetworkClass(N=n_pop, m_0=3)

    for agent in myNetworkObj.All_agentSet.get_agents():
        assert agent in myNetworkObj.G.nodes()

    for agent in myNetworkObj.G:
        agent_sex_type = agent._SO
        agent_drug_type = agent._DU
        assert agent_drug_type in ["IDU", "NIDU", "NDU"]
        assert agent_sex_type in ["MSM", "MSW", "HM", "HF"]


def test_PopulationConsistency():
    """Test if Drug users add up"""
    myNetworkObj = NetworkClass(N=n_pop, m_0=3)
    CheckSumDrug = (
        myNetworkObj.DU_IDU_agentSet.num_members()
        + myNetworkObj.DU_NIDU_agentSet.num_members()
        + myNetworkObj.DU_NDU_agentSet.num_members()
    )
    assert myNetworkObj.drugUse_agentSet.num_members() == CheckSumDrug
    assert myNetworkObj.PopulationSize == CheckSumDrug


def test_HIVConsistency():
    """Test HIV consistency"""
    myNetworkObj = NetworkClass(N=n_pop, m_0=3)
    for agent in myNetworkObj.All_agentSet.get_agents():
        HIV_status = agent._HIV_bool
        if HIV_status:
            assert agent in myNetworkObj.HIV_agentSet.get_agents()

    for agent in myNetworkObj.HIV_agentSet.get_agents():
        HIV_status = agent._HIV_bool
        assert HIV_status
