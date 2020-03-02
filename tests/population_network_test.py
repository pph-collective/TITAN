import pytest
from titan.population_network import PopulationClass
from titan.agent import Agent
from titan import params

n_pop = 100


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="WHITE", DU="NDU"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_population():
    def _make_population(n=n_pop):
        return PopulationClass(n, 0, None, enable_nx_graph=True)

    return _make_population


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float):
        self.num = num

    def random(self):
        return self.num

    def randrange(self, start, stop, step=1):
        return start

    def randint(self, start, stop):
        return start


def test_pop_init(make_population):
    n_pop = 100
    pop = make_population(n=n_pop)

    assert pop.numWhite + pop.numBlack == n_pop
    assert pop.All_agentSet.num_members() == n_pop

    # test umbrella sets are all consistent
    parent_sets = ["drugUse_agentSet", "SO_agentSet", "racial_agentSet"]
    for set_name in parent_sets:
        set = getattr(pop, set_name)

        assert set.num_members() == n_pop

        child_pops = 0
        for child_set in set.iter_subset():
            child_pops += child_set.num_members()

        assert child_pops == n_pop


def test_create_agent(make_population):
    pop = make_population()

    a1 = pop.create_agent("WHITE")
    assert a1._race == "WHITE"
    assert a1.opinion in range(
        5
    ), f"Agents opinion of injectible PrEP is out of bounds {a1.opinion}"

    a2 = pop.create_agent("BLACK")
    assert a2._race == "BLACK"

    a3 = pop.create_agent("WHITE", "HM")
    assert a3._SO == "HM"
    assert a3._race == "WHITE"

    # check IDU and HIV and high risk
    pop.popRandom = FakeRandom(-0.1)
    a4 = pop.create_agent("WHITE")
    assert a4._DU == "IDU"
    assert a4._HIV_bool
    assert a4._AIDS_bool
    assert a4._tested
    assert a4._HAART_bool
    assert a4._HAART_adh == 5
    assert a4._HAART_time == 0
    assert a4._treatment_bool
    assert a4._highrisk_bool
    assert a4._everhighrisk_bool

    # check not IDU and HIV
    pop.popRandom = FakeRandom(0.999)
    a4 = pop.create_agent("WHITE")
    assert a4._DU == "NDU"
    assert a4._HIV_bool is False
    assert a4._PrEP_bool is False
    assert a4._treatment_bool is False


def test_create_agent_proportions(make_population):
    pop = make_population()

    n = 1000
    race = "WHITE"
    # check proportions
    pop.pop_weights[race] = {"values": ["HM", "HF"], "weights": [0.1, 0.9]}
    prop_idu = round(params.DemographicParams[race]["IDU"]["POP"] * n)
    num_HM = 0
    num_HF = 0
    num_IDU = 0
    for i in range(n):
        a = pop.create_agent(race)
        if a._DU == "IDU":
            num_IDU += 1

        if a._SO == "HF":
            num_HF += 1
        elif a._SO == "HM":
            num_HM += 1
        else:
            assert False

    assert num_HM > 80 and num_HM < 120
    assert num_HF > 880 and num_HF < 920
    assert num_IDU > prop_idu - 20 and num_IDU < prop_idu + 20


def test_add_agent_to_pop(make_population):
    pop = make_population()
    agent = pop.create_agent("WHITE", "HM")
    agent._DU = "IDU"
    agent._HIV_bool = True
    agent._AIDS_bool = True
    agent._treatment_bool = True
    agent._HAART_bool = True
    agent._PrEP_bool = True
    agent._tested = True
    agent._incar_bool = True
    agent._highrisk_bool = True

    pop.add_agent_to_pop(agent)

    assert agent in pop.All_agentSet._members
    assert agent in pop.racial_agentSet._members
    assert agent in pop.Race_WHITE_agentSet._members
    assert agent in pop.SO_agentSet._members
    assert agent in pop.SO_HM_agentSet._members
    assert agent in pop.drugUse_agentSet._members
    assert agent in pop.DU_IDU_agentSet._members
    assert agent in pop.HIV_agentSet._members
    assert agent in pop.HIV_AIDS_agentSet._members
    assert agent in pop.treatment_agentSet._members
    assert agent in pop.Trt_ART_agentSet._members
    assert agent in pop.Trt_ART_agentSet._members
    assert agent in pop.Trt_PrEP_agentSet._members
    assert agent in pop.Trt_Tstd_agentSet._members
    assert agent in pop.incarcerated_agentSet._members
    assert agent in pop.highrisk_agentsSet._members

    # check not in all agent sets
    assert agent not in pop.Race_BLACK_agentSet._members


def test_get_age(make_population):
    pop = make_population()

    race = "BLACK"

    expected_ages = [15, 25, 35, 45]
    for i in range(1, 5):
        # make sure rand is less tan the setting
        pop.popRandom = FakeRandom(params.ageMatrix[race]["Prop"][i] - 0.001)
        age, ageBin = pop.get_age(race)
        assert age == expected_ages[i - 1]
        assert ageBin == i

    # test else case
    pop.popRandom = FakeRandom(1.1)
    age, ageBin = pop.get_age(race)
    assert age == 55
    assert ageBin == 5


def test_update_agent_partners_no_match(make_population):
    pop = PopulationClass(n=1, enable_nx_graph=True)
    G = pop.get_Graph()

    agent = pop.All_agentSet._members[0]  # the only agent in the pop

    assert pop.update_agent_partners(agent)  # noMatch == True
    assert agent in G.nodes()


def test_update_agent_partners_match(make_population):
    # pop = make_population(n=0)
    pop = PopulationClass(n=0, enable_nx_graph=True)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    pop.add_agent_to_pop(a, update_nx=True)
    pop.add_agent_to_pop(p, update_nx=True)

    assert pop.update_agent_partners(a, update_nx=True) is False  # noMatch == False
    assert a in pop.nx_graph.nodes()
    assert p in pop.nx_graph.nodes()
    assert len(pop.nx_graph.edges()) == 1


def test_update_partner_assignments_match(make_population):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    pop.add_agent_to_pop(a, update_nx=True)
    pop.add_agent_to_pop(p, update_nx=True)
    pop.update_partner_assignments(100.0, update_nx=True) is True  # noMatch == False
    assert a in pop.nx_graph.nodes()
    assert p in pop.nx_graph.nodes()
    assert len(pop.nx_graph.edges()) == 1


def test_update_partner_assignments_no_match(make_population):
    pop = PopulationClass(n=0, enable_nx_graph=True)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    pop.add_agent_to_pop(a, update_nx=True)
    pop.add_agent_to_pop(p, update_nx=True)

    pop.update_partner_assignments(0.0, pop.nx_graph) is False  # noMatch == False
    assert a in pop.nx_graph.nodes()
    assert p in pop.nx_graph.nodes()
    assert len(pop.nx_graph.edges()) == 0


def test_network_init_scale_free(make_population):
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    net = make_population()
    assert n_pop == net.All_agentSet.num_members()

    for agent in net.All_agentSet.get_agents():
        assert agent in net.nx_graph.nodes()

    for agent in net.All_agentSet.get_agents():
        assert agent._DU in ["IDU", "NIDU", "NDU"]
        assert agent._SO in params.agentSexTypes

    assert net.get_Graph() == net.nx_graph


def test_network_init_max_k(make_population):
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    net = make_population()
    net.create_network(network_type="max_k_comp_size")
    assert n_pop == net.All_agentSet.num_members()

    for agent in net.All_agentSet.get_agents():
        assert agent in net.nx_graph.nodes()

    for agent in net.All_agentSet.get_agents():
        assert agent._DU in ["IDU", "NIDU", "NDU"]
        assert agent._SO in params.agentSexTypes


def test_population_consistency_DU(make_population):
    """Test if Drug users add up"""
    net = make_population()
    check_sum_DU = (
        net.DU_IDU_agentSet.num_members()
        + net.DU_NIDU_agentSet.num_members()
        + net.DU_NDU_agentSet.num_members()
    )

    assert net.drugUse_agentSet.num_members() == check_sum_DU
    assert net.PopulationSize == check_sum_DU


def test_population_consistency_HIV(make_population):
    """Test HIV consistency"""
    net = make_population()
    for agent in net.All_agentSet.get_agents():
        if agent._HIV_bool:
            assert agent in net.HIV_agentSet.get_agents()

    for agent in net.HIV_agentSet.get_agents():
        assert agent._HIV_bool


def test_create_graph_from_agents(make_agent):
    net = PopulationClass(n=n_pop)

    assert net.nx_graph.number_of_nodes() == 0

    net.create_graph_from_agents()

    assert net.nx_graph.number_of_nodes() == n_pop
