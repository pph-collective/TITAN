import pytest

from copy import deepcopy

from titan.ABM_core import *
from titan.agent import Agent, Relationship
from titan import params


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="NDU"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_model():
    def _make_model():
        return HIVModel(100, 10, 0, 0, 0, "scale_free")

    return _make_model


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float):
        assert num >= 0 and num <= 1
        self.num = num

    def random(self):
        return self.num

    def randrange(self, start, stop, step):
        return start

    def sample(self, seq, rate):
        return seq

    def choice(self, seq):
        return seq[0]


# ================================ MODEL TESTS =================================


def test_model_init_error():
    with pytest.raises(ValueError):
        HIVModel(100, 0.5, 0, 0, 0, "scale_free")


def test_model_init():
    model = HIVModel(100, 10, 0, 0, 0, "scale_free")

    assert model.runseed > 0
    assert model.popseed > 0
    assert model.netseed > 0

    assert model.NewInfections.num_members() == 0
    assert model.NewDiagnosis.num_members() == 0
    assert model.NewIncarRelease.num_members() == 0
    assert model.NewHRrolls.num_members() == 0

    assert len(model.Acute_agents) == 0
    assert len(model.Transmit_from_agents) == 0
    assert len(model.Transmit_to_agents) == 0

    assert model.totalDiagnosis == 0
    assert model.totalIncarcerated == 0
    assert model.treatmentEnrolled == False


@pytest.mark.skip("to do")
def test_update_AllAgents():
    pass


@pytest.mark.skip("to do")
def test_agents_interact():
    pass


@pytest.mark.skip("to do")
def test_needle_transmission():
    pass


@pytest.mark.skip("to do")
def test_sex_transmission():
    pass


@pytest.mark.skip("to do")
def test_become_HIV():
    pass


@pytest.mark.skip("to do")
def test_enroll_treatment():
    pass


@pytest.mark.skip("to do")
def test_becomeHighRisk():
    pass


@pytest.mark.skip("to do")
def test_incarcerate():
    pass


@pytest.mark.skip("to do")
def test_HIVtest():
    pass


@pytest.mark.skip("to do")
def test_update_HAART():
    pass


def test_discont_PrEP_force(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # set up so the agent appears to be on PrEP
    a._PrEP_bool = True
    a._PrEP_reason = ["blah"]
    a._PrEP_time = 1000
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a, 0, True)

    assert a._PrEP_bool is False
    assert a._PrEP_reason == []
    assert a._PrEP_time == 0
    assert a not in model.Trt_PrEP_agentSet._members


def test_discont_PrEP_decrement_time(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # set up so the agent appears to be on PrEP
    a._PrEP_bool = True
    a._PrEP_reason = ["blah"]
    a._PrEP_time = 1000
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a, 0)

    assert a._PrEP_bool
    assert a._PrEP_reason == ["blah"]
    assert a._PrEP_time == 999


def test_discont_PrEP_decrement_end(make_model, make_agent):
    model = make_model()
    a = make_agent(race="WHITE")

    model.runRandom = FakeRandom(0.00000001)

    # set up so the agent appears to be on PrEP
    a._PrEP_bool = True
    a._PrEP_reason = ["blah"]
    a._PrEP_time = 0
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a, 0)

    assert a._PrEP_bool
    assert a._PrEP_reason == ["blah"]
    assert a._PrEP_time == params.PrEP_falloutT
    assert a not in model.Trt_PrEP_agentSet._members


def test_discont_PrEP_decrement_not_end(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.runRandom = FakeRandom(0.9999)

    # set up so the agent appears to be on PrEP
    a._PrEP_bool = True
    a._PrEP_reason = ["blah"]
    a._PrEP_time = 0
    a._PrEP_lastDose = 3
    model.Trt_PrEP_agentSet.add_agent(a)

    model._discont_PrEP(a, 0)

    assert a._PrEP_bool
    assert a._PrEP_reason == ["blah"]
    assert a._PrEP_time == 0
    assert a._PrEP_lastDose == 0  # 3 -> -1 -> +1 == 0
    assert a in model.Trt_PrEP_agentSet._members
    assert a._PrEP_load > 0


def test_initiate_PrEP_assertions(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # no PreP if already PreP
    a._PrEP_bool = True
    assert model._initiate_PrEP(a, 0) is None

    # no PrEP if already HIV
    a._PrEP_bool = False
    a._HIV_bool = True
    assert model._initiate_PrEP(a, 0) is None


def test_initiate_PrEP_force_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # forcing, adherant, inj
    model.runRandom = FakeRandom(0.00001)
    model._initiate_PrEP(a, 0, True)
    assert a._PrEP_bool
    assert a._PrEP_time == 0
    assert a in model.Trt_PrEP_agentSet._members
    assert a in model.newPrEPagents._members
    assert a._PrEP_adh == 1
    assert a._PrEP_load > 0.0
    assert a._PrEP_lastDose == 0


def test_initiate_PrEP_force_non_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()
    # forcing, non-adherant, inj
    model.runRandom = FakeRandom(1.0)
    model._initiate_PrEP(a, 0, True)
    assert a._PrEP_bool
    assert a._PrEP_time == 0
    assert a in model.Trt_PrEP_agentSet._members
    assert a in model.newPrEPagents._members
    assert a._PrEP_adh == 0
    assert a._PrEP_load > 0.0
    assert a._PrEP_lastDose == 0


def test_initiate_PrEP_eligible(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HF")  # model is "CDCwomen"
    p = make_agent(DU="IDU")
    p._tested = True
    p._MSMW = True
    rel = Relationship(a, p, 10)
    # non-forcing, adherant, inj
    model.runRandom = FakeRandom(0.00001)
    model._initiate_PrEP(a, 0)
    assert a._PrEP_bool
    assert a._PrEP_time == 0
    assert a in model.Trt_PrEP_agentSet._members
    assert a in model.newPrEPagents._members
    assert a._PrEP_adh == 1
    assert a._PrEP_load > 0.0
    assert a._PrEP_lastDose == 0
    assert "IDU" in a._PrEP_reason
    assert "HIV test" in a._PrEP_reason
    assert "MSMW" in a._PrEP_reason


def test_get_clinic_agent_none(make_model, make_agent):
    model = make_model()
    clinic_cat = "Mid"

    a = make_agent()
    a._mean_num_partners = 2
    pool = [a]

    # probability of matching is on bin 1 0.054
    # min is 0 max is 1
    model.runRandom = FakeRandom(0.0001)

    assert model._get_clinic_agent(clinic_cat, pool) is None


def test_get_clinic_agent_match(make_model, make_agent):
    model = make_model()
    clinic_cat = "Mid"

    a = make_agent()
    a._mean_num_partners = 1
    pool = [a]

    # probability of matching is on bin 1 0.054
    # min is 0 max is 1
    model.runRandom = FakeRandom(0.0001)

    assert model._get_clinic_agent(clinic_cat, pool) == a


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
    a._HIV_bool = True
    a._HAART_adh = 1  # .0051 prob

    model.runRandom = FakeRandom(0.9)  # no AIDS

    assert model._progress_to_AIDS(a) is None
    assert model.HIV_AIDS_agentSet.num_members() == num_aids
    assert a._AIDS_bool is False


def test_progress_to_AIDS_progress(make_agent, make_model):
    a = make_agent()
    model = make_model()
    num_aids = model.HIV_AIDS_agentSet.num_members()  # get baseline

    a._HIV_bool = True
    a._HAART_adh = 1  # .0051 prob

    # test progress case
    model.runRandom = FakeRandom(0.001)  # AIDS

    assert model._progress_to_AIDS(a) is None
    assert model.HIV_AIDS_agentSet.num_members() == num_aids + 1
    assert a in model.HIV_AIDS_agentSet._members
    assert a._AIDS_bool is True


def test_die_and_replace_none(make_model):
    model = make_model()
    model.runRandom = FakeRandom(0.999)  # always greater than death rate
    baseline_pop = deepcopy(model.All_agentSet._members)

    model._die_and_replace(0)  # add arbitrarty time

    ids = [a._ID for a in baseline_pop]
    for agent in model.All_agentSet._members:
        assert agent._ID in ids


def test_die_and_replace_all(make_model):
    model = make_model()
    model.runRandom = FakeRandom(0.0000001)  # always lower than death rate

    # un-incarcerate everyone
    for agent in model.All_agentSet._members:
        agent._incar_bool = False

    baseline_pop = deepcopy(model.All_agentSet._members)
    old_ids = [a._ID for a in baseline_pop]

    num_HM = len([x for x in baseline_pop if x._SO == "HM"])
    num_WHITE = len([x for x in baseline_pop if x._race == "WHITE"])

    model._die_and_replace(0)  # add arbitrarty time

    assert num_HM == len([x for x in model.All_agentSet._members if x._SO == "HM"])
    assert num_WHITE == len(
        [x for x in model.All_agentSet._members if x._race == "WHITE"]
    )

    new_ids = [a._ID for a in model.All_agentSet._members]
    death_ids = [a._ID for a in model.deathSet]

    for agent in model.All_agentSet._members:
        assert agent._ID not in old_ids
        assert agent in model.get_Graph().nodes()

    for agent in baseline_pop:
        assert agent._ID not in new_ids
        assert agent not in model.get_Graph().nodes()
        assert agent._ID in death_ids


def test_die_and_replace_incar(make_model):
    model = make_model()
    model.runRandom = FakeRandom(0.0000001)  # always lower than death rate
    baseline_pop = deepcopy(model.All_agentSet._members)
    old_ids = [a._ID for a in baseline_pop]

    model.All_agentSet._members[0]._incar_bool = True
    agent_id = model.All_agentSet._members[0]._ID

    model._die_and_replace(0)  # add arbitrarty time

    new_ids = [a._ID for a in model.All_agentSet._members]
    death_ids = [a._ID for a in model.deathSet]

    assert agent_id in old_ids
    assert agent_id not in death_ids
    assert agent_id in new_ids
