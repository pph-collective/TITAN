import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_update_haart_t1(make_model, make_agent):
    model = make_model()
    model.time = 1
    a = make_agent(race="white")

    a.hiv = True

    # nothing happens, not tested
    a.haart.update_agent(model)
    assert a.haart.adherence == 0
    assert a.haart.active is False

    # t0 agent initialized HAART
    a.hiv_dx = True

    # go on haart
    model.run_random = FakeRandom(
        -0.1
    )  # means this will always be less than params even though not physically possible in reality
    a.haart.update_agent(model)

    assert a.haart.adherence == 5
    assert a.haart.time == 1
    assert a.haart.active

    # go off haart
    a.haart.update_agent(model)

    assert a.haart.adherence == 0
    assert a.haart.time == 0
    assert a.haart.active is False
