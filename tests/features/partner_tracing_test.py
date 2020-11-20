import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_partner_tracing(make_model, make_agent):
    model = make_model()
    model.time = model.params.partner_tracing.start_time + 1
    a = make_agent()
    p = make_agent()
    a.hiv.active = True
    p.hiv.active = True
    p.location.params.partner_tracing.trace_duration = 2
    a.partners["Sex"].add(p)

    model.run_random = FakeRandom(-0.1)  # always less than param
    a.hiv.diagnose(model)

    assert a.hiv.dx
    assert a.hiv.dx_time == model.time

    model.time += 1

    a.partner_tracing.update_agent(model)

    assert p in a.get_partners()
    assert p.partner_tracing.active
    assert p.partner_tracing.time == model.time

    assert p.hiv.dx is False
    model.params.demographics[p.race][p.sex_type].hiv.dx.prob = 0

    model.time += 1
    p.partner_tracing.update_agent(model)
    assert p.hiv.dx
    assert p.partner_tracing.active

    model.time += 1
    p.partner_tracing.update_agent(model)
    assert p.partner_tracing.active is False
    assert p.partner_tracing.time is None
