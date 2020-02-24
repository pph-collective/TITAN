import pytest

from copy import deepcopy
import os

from titan.ABM_core import *
from titan.agent import Agent, Relationship
from titan.params_parse import create_params


@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )
    return create_params(None, param_file, tmpdir)


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="None"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_model(params):
    def _make_model():
        return HIVModel(params)

    return _make_model


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float):
        self.num = num

    def random(self):
        return self.num

    def randrange(self, start, stop, step):
        return start

    def sample(self, seq, rate):
        return seq

    def choice(self, seq):
        return seq[0]

    def randint(self, start, stop):
        return start


# ================================ MODEL TESTS =================================


def test_model_init_error(params):
    params.model.seed.run = 0.5
    with pytest.raises(ValueError):
        HIVModel(params)


def test_model_init(params):
    model = HIVModel(params)

    assert model.runseed > 0
    assert model.popseed > 0
    assert model.netseed > 0

    assert model.NewInfections.num_members() == 0
    assert model.NewDiagnosis.num_members() == 0
    assert model.NewIncarRelease.num_members() == 0
    assert model.NewHRrolls.num_members() == 0

    assert model.totalDiagnosis == 0
    assert model.needle_exchange == False


@pytest.mark.skip("too parameter dependent to test at this point")
def test_update_AllAgents():
    pass


def test_agents_interact(make_model, make_agent):
    model = make_model()
    a = make_agent(race="WHITE", SO="HM")
    p = make_agent(race="WHITE", SO="HF")
    rel = Relationship(a, p, 10, bond_type="sexOnly")

    model.runRandom = FakeRandom(0.6)

    a.incar = True
    assert model._agents_interact(0, rel) is False

    a.incar = False
    assert model._agents_interact(0, rel) is False  # neither HIV

    a.hiv = True
    p.hiv = True
    assert model._agents_interact(0, rel) is False  # both HIV

    p.hiv = False

    assert model._agents_interact(0, rel)  # sex transmssion
    assert p.hiv is False  # but nothing happened (see skipped test)

    a.drug_use = "Inj"
    p.drug_use = "Inj"

    model.runRandom = FakeRandom(-0.1)

    assert model._agents_interact(0, rel)  # needle transmission
    assert p.hiv

    p.hiv = False
    model.runRandom = FakeRandom(1.1)

    assert model._agents_interact(0, rel)  # needle and sex
    assert p.hiv is False  # but nothing happened


def test_needle_transmission(make_model, make_agent):
    model = make_model()
    a = make_agent(race="WHITE", DU="Inj", SO="HM")
    p = make_agent(race="WHITE", DU="Inj", SO="HF")

    with pytest.raises(AssertionError):
        model._needle_transmission(a, p, time=0)

    a.hiv = True
    a.hiv_time = 1  # acute

    model.runRandom = FakeRandom(-0.1)

    model._needle_transmission(a, p, time=0)

    assert p.hiv


def test_sex_transmission(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    rel = Relationship(a, p, 10, bond_type="sexOnly")

    a.hiv = True
    a.hiv_time = 1  # acute

    rel.total_sex_acts = 0

    model.runRandom = FakeRandom(0.6)

    # test partner becomes
    model._sex_transmission(rel, 0)

    assert p.hiv


def test_sex_transmission_do_nothing(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    rel = Relationship(a, p, 10, bond_type="sexOnly")

    with pytest.raises(ValueError):
        model._sex_transmission(rel, 0)

    a.hiv = True
    p.hiv = True

    # test nothing happens
    model._sex_transmission(rel, 0)


def test_pca_interaction(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.opinion = 4  # TO_REVIEW opinino and awareness are both prep things right? should the be prepended with prep_?
    p.opinion = 2
    a.awareness = True

    model.runRandom = FakeRandom(1.0)

    model.G.add_edge(a, p)
    model.G.add_edge(a, "edge")

    rel = Relationship(a, p, 10, bond_type="multiplex")
    model._pca_interaction(rel, 5, force=True)

    assert p.awareness

    model._pca_interaction(rel, 6, force=True)

    assert p.opinion == 3


def test_become_HIV(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.prep = True

    model.runRandom = FakeRandom(-0.1)

    model._become_HIV(a)

    assert a.hiv
    assert a.hiv_time == 1
    assert a in model.NewInfections.members
    assert a in model.HIV_agentSet.members
    assert a.prep is False


def testenroll_needle_exchange(make_model):
    model = make_model()
    model.runRandom = FakeRandom(-0.1)  # all "Inj" agents will be _SNE_bool

    # make at least one agent PWID
    model.All_agentSet.members[0].drug_use = "Inj"

    assert model.needle_exchange is False

    model.enroll_needle_exchange()

    assert model.needle_exchange is True

    for a in model.All_agentSet.members:
        if a.drug_use == "Inj":
            assert a.sne


def test_becomeHighRisk(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model._become_high_risk(a, 10)

    assert a in model.highrisk_agentsSet.members
    assert a in model.NewHRrolls.members
    assert a.high_risk
    assert a.high_risk_ever
    assert a.high_risk_time == 10


def test_incarcerate_tested(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="WHITE")  # incarceration only for HM and HF?
    a.hiv = True
    a.hiv_dx = True

    model.runRandom = FakeRandom(-0.1)  # always less than params

    model._incarcerate(a, 10)

    assert a.incar
    assert a.incar_time == 1
    assert a in model.incarcerated_agentSet.members
    assert a.haart
    assert a.haart_adherence == 5
    assert a.haart_time == 10
    assert a in model.Trt_ART_agentSet.members


def test_incarcerate_not_tested(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="WHITE")  # incarceration only for HM and HF?
    a.hiv = True

    p = make_agent(SO="HF")
    rel = Relationship(a, p, 10, bond_type="sexOnly")

    model.runRandom = FakeRandom(-0.1)  # always less than params

    model._incarcerate(a, 0)

    assert a.incar
    assert a.incar_time == 1
    assert a in model.incarcerated_agentSet.members
    assert a.hiv_dx

    assert p in model.highrisk_agentsSet.members
    assert p in model.NewHRrolls.members
    assert p.high_risk
    assert p.high_risk_ever
    assert p.high_risk_time > 0


def test_incarcerate_unincarcerate(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.incar = True
    a.incar_time = 2
    model.incarcerated_agentSet.add_agent(a)

    model._incarcerate(a, 0)

    assert a.incar
    assert a.incar_time == 1
    assert a in model.incarcerated_agentSet.members

    model._incarcerate(a, 0)

    assert a.incar is False
    assert a.incar_time == 0
    assert a not in model.incarcerated_agentSet.members
    assert a in model.NewIncarRelease.members
    assert a.incar_ever


def test_HIVtest(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.runRandom = FakeRandom(1.1)  # always greater than param
    model._HIVtest(a, 0)

    assert a.hiv_dx is False
    assert a not in model.NewDiagnosis.members
    assert a not in model.Trt_Tstd_agentSet.members

    model.runRandom = FakeRandom(-0.1)  # always less than param
    model._HIVtest(a, 0)

    assert a.hiv_dx
    assert a in model.NewDiagnosis.members
    assert a in model.Trt_Tstd_agentSet.members


def test_HIVtest_already_tested(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.hiv_dx = True

    model.runRandom = FakeRandom(-0.1)  # always less than param
    model._HIVtest(a, 0)

    assert a.hiv_dx
    assert a not in model.NewDiagnosis.members
    assert a not in model.Trt_Tstd_agentSet.members


def test_update_HAART_t1(make_model, make_agent):
    model = make_model()
    a = make_agent(race="WHITE")

    a.hiv = True

    # nothing happens, not tested
    model._update_HAART(a, 1)
    assert a.haart_adherence == 0
    assert a.haart is False
    assert a not in model.Trt_ART_agentSet.members

    # t0 agent initialized HAART
    a.hiv_dx = True

    # go on haart
    model.runRandom = FakeRandom(
        -0.1
    )  # means this will always be less than params even though not physically possible in reality
    model._update_HAART(a, 1)

    assert a.haart_adherence == 5
    assert a.haart_time == 1
    assert a.haart
    assert a in model.Trt_ART_agentSet.members

    # go off haart
    model._update_HAART(a, 1)

    assert a.haart_adherence == 0
    assert a.haart_time == 0
    assert a.haart is False
    assert a not in model.Trt_ART_agentSet.members


def test_discont_PrEP_force(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_reason = ["blah"]
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a, True)

    assert a.prep is False
    assert a.prep_reason == []
    assert a not in model.Trt_PrEP_agentSet.members


def test_discont_PrEP_decrement_time(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_reason = ["blah"]
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a)

    assert a.prep
    assert a.prep_reason == ["blah"]


def test_discont_PrEP_decrement_end(make_model, make_agent):
    model = make_model()
    a = make_agent(race="WHITE")

    model.runRandom = FakeRandom(-0.1)

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_reason = ["blah"]
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a)

    assert a.prep is False
    assert a.prep_reason == []
    assert a not in model.Trt_PrEP_agentSet.members


def test_discont_PrEP_decrement_not_end(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.runRandom = FakeRandom(1.1)

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_reason = ["blah"]
    a.prep_last_dose = 3
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a)

    assert a.prep
    assert a.prep_reason == ["blah"]
    assert a.prep_last_dose == -1  # 3 -> -1 -> +1 == 0 # Inj no longer in PrEP types
    assert a in model.Trt_PrEP_agentSet.members
    # assert a.prep_load > 0 # Inj no longer in PrEP types


def test_initiate_PrEP_assertions(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # no PreP if already PreP
    a.prep = True
    assert model._initiate_PrEP(a, 0) is None

    # no PrEP if already HIV
    a.prep = False
    a.hiv = True
    assert model._initiate_PrEP(a, 0) is None


def test_initiate_PrEP_force_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # forcing, adherant, inj
    model.runRandom = FakeRandom(-0.1)
    model._initiate_PrEP(a, 0, True)
    assert a.prep
    assert a in model.Trt_PrEP_agentSet.members
    assert a in model.newPrEPagents.members
    assert a.prep_adherence == 1
    # assert a.prep_load > 0.0 # Inj no longer in PrEP types
    assert a.prep_last_dose == 0


def test_initiate_PrEP_force_non_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()
    # forcing, non-adherant, inj
    model.runRandom = FakeRandom(1.0)
    model._initiate_PrEP(a, 0, True)
    assert a.prep
    assert a in model.Trt_PrEP_agentSet.members
    assert a in model.newPrEPagents.members
    assert a.prep_adherence == 0
    # assert a.prep_load > 0.0 # Inj no longer in PrEP types
    assert a.prep_last_dose == 0


def test_initiate_PrEP_eligible(make_model, make_agent):
    model = make_model()

    # make sure there's room to add more prep agents
    model.Trt_PrEP_agentSet.members = []
    a = make_agent(SO="HF")  # model is "CDCwomen"
    p = make_agent(DU="Inj")
    p.hiv_dx = True
    p.msmw = True
    rel = Relationship(a, p, 10, bond_type="sexOnly")
    # non-forcing, adherant, inj
    model.runRandom = FakeRandom(-0.1)
    model._initiate_PrEP(a, 0)
    assert a.prep
    assert a in model.Trt_PrEP_agentSet.members
    assert a in model.newPrEPagents.members
    assert a.prep_adherence == 1
    # assert a.prep_load > 0.0 # Inj not in params prep_type anymore
    assert a.prep_last_dose == 0
    assert "PWID" in a.prep_reason
    assert "HIV test" in a.prep_reason
    assert "MSMW" in a.prep_reason


def test_progress_to_AIDS_error(make_agent, make_model):
    a = make_agent()
    model = make_model()
    num_aids = model.HIV_AIDS_agentSet.num_members()  # get baseline

    # test error case, agent must be HIV+
    with pytest.raises(ValueError):
        model._progress_to_AIDS(a)

    assert model.HIV_AIDS_agentSet.num_members() == num_aids


def test_progress_to_AIDS_nothing(make_agent, make_model):
    a = make_agent()
    model = make_model()
    num_aids = model.HIV_AIDS_agentSet.num_members()  # get baseline

    # test nothing case
    a.hiv = True
    a.haart_adherence = 1  # .0051 prob

    model.runRandom = FakeRandom(0.9)  # no AIDS

    assert model._progress_to_AIDS(a) is None
    assert model.HIV_AIDS_agentSet.num_members() == num_aids
    assert a.aids is False


def test_progress_to_AIDS_progress(make_agent, make_model):
    a = make_agent()
    model = make_model()
    num_aids = model.HIV_AIDS_agentSet.num_members()  # get baseline

    a.hiv = True
    a.haart_adherence = 1  # .0051 prob

    # test progress case
    model.runRandom = FakeRandom(0.001)  # AIDS

    assert model._progress_to_AIDS(a) is None
    assert model.HIV_AIDS_agentSet.num_members() == num_aids + 1
    assert a in model.HIV_AIDS_agentSet.members
    assert a.aids is True


def test_die_and_replace_none(make_model):
    model = make_model()
    model.runRandom = FakeRandom(0.999)  # always greater than death rate
    baseline_pop = deepcopy(model.All_agentSet.members)

    model._die_and_replace()

    ids = [a.id for a in baseline_pop]
    for agent in model.All_agentSet.members:
        assert agent.id in ids


def test_die_and_replace_all(make_model):
    model = make_model()
    model.runRandom = FakeRandom(0.0000001)  # always lower than death rate

    # un-incarcerate everyone
    for agent in model.All_agentSet.members:
        agent.incar = False

    baseline_pop = deepcopy(model.All_agentSet.members)
    old_ids = [a.id for a in baseline_pop]

    num_HM = len([x for x in baseline_pop if x.so == "HM"])
    num_WHITE = len([x for x in baseline_pop if x.race == "WHITE"])

    model._die_and_replace()

    assert num_HM == len([x for x in model.All_agentSet.members if x.so == "HM"])
    assert num_WHITE == len(
        [x for x in model.All_agentSet.members if x.race == "WHITE"]
    )

    new_ids = [a.id for a in model.All_agentSet.members]
    death_ids = [a.id for a in model.deathSet]

    for agent in model.All_agentSet.members:
        assert agent.id not in old_ids
        assert agent in model.get_Graph().nodes()

    for agent in baseline_pop:
        assert agent.id not in new_ids
        assert agent not in model.get_Graph().nodes()
        assert agent.id in death_ids


def test_die_and_replace_incar(make_model):
    model = make_model()
    model.runRandom = FakeRandom(0.0000001)  # always lower than death rate
    baseline_pop = deepcopy(model.All_agentSet.members)
    old_ids = [a.id for a in baseline_pop]

    model.All_agentSet.members[0].incar = True
    agent_id = model.All_agentSet.members[0].id

    model._die_and_replace()

    new_ids = [a.id for a in model.All_agentSet.members]
    death_ids = [a.id for a in model.deathSet]

    assert agent_id in old_ids
    assert agent_id not in death_ids
    assert agent_id in new_ids
