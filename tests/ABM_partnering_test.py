import pytest

from src.ABM_partnering import *
from src.agent import Agent


@pytest.fixture
def make_agent():
    def _make_agent(id, SO="MSM", age=30, race="BLACK", DU="NDU"):
        return Agent(id, SO, age, race, DU)

    return _make_agent


@pytest.mark.skip(reason="#REVIEW - EligSE_PartnerType")
def test_get_IDU_partner_no_IDU(make_agent):
    idu_agent = make_agent(1, DU="IDU")
    nd_agent = make_agent(2)
    avail_partners = [nd_agent]
    assert get_partner(idu_agent, avail_partners) == nd_agent


def test_sex_possible():
    # agent sex types are ["HM", "MSM", "WSW", "HF", "MTF"]
    assert sex_possible("HM", "HM") == False
    assert sex_possible("HM", "MSM") == False
    assert sex_possible("HM", "HF") == True
    assert sex_possible("HM", "WSW") == True
    assert sex_possible("HM", "MTF") == True

    assert sex_possible("MSM", "HM") == False
    assert sex_possible("MSM", "MSM") == True
    assert sex_possible("MSM", "HF") == True
    assert sex_possible("MSM", "WSW") == True
    assert sex_possible("MSM", "MTF") == True

    assert sex_possible("WSW", "HM") == True
    assert sex_possible("WSW", "MSM") == True
    assert sex_possible("WSW", "HF") == False
    assert sex_possible("WSW", "WSW") == True
    assert sex_possible("WSW", "MTF") == False

    assert sex_possible("HF", "HM") == True
    assert sex_possible("HF", "MSM") == True
    assert sex_possible("HF", "HF") == False
    assert sex_possible("HF", "WSW") == False
    assert sex_possible("HF", "MTF") == False

    assert sex_possible("MTF", "HM") == True
    assert sex_possible("MTF", "MSM") == True
    assert sex_possible("MTF", "HF") == False
    assert sex_possible("MTF", "WSW") == False
    assert sex_possible("MTF", "MTF") == False

    with pytest.raises(ValueError, match=r"Invalid .*_sex_type.*"):
        sex_possible("HM", "XYZ")
        sex_possible("XYZ", "HM")
