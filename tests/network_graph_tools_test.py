import pytest

from titan.network_graph_tools import *
from titan import params

n_pop = 10


def test_network_init():
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    myNetworkObj = NetworkClass(N=n_pop)
    assert n_pop == myNetworkObj.All_agentSet.num_members()

    for agent in myNetworkObj.All_agentSet.get_agents():
        assert agent in myNetworkObj.G.nodes()

    for agent in myNetworkObj.All_agentSet.get_agents():
        assert agent._DU in ["IDU", "NIDU", "NDU"]
        assert agent._SO in params.agentSexTypes


def test_population_consistency_DU():
    """Test if Drug users add up"""
    myNetworkObj = NetworkClass(N=n_pop)
    CheckSumDrug = (
        myNetworkObj.DU_IDU_agentSet.num_members()
        + myNetworkObj.DU_NIDU_agentSet.num_members()
        + myNetworkObj.DU_NDU_agentSet.num_members()
    )
    assert myNetworkObj.drugUse_agentSet.num_members() == CheckSumDrug
    assert myNetworkObj.PopulationSize == CheckSumDrug


def test_population_consistency_HIV():
    """Test HIV consistency"""
    myNetworkObj = NetworkClass(N=n_pop,)
    for agent in myNetworkObj.All_agentSet.get_agents():
        HIV_status = agent._HIV_bool
        if HIV_status:
            assert agent in myNetworkObj.HIV_agentSet.get_agents()

    for agent in myNetworkObj.HIV_agentSet.get_agents():
        HIV_status = agent._HIV_bool
        assert HIV_status
