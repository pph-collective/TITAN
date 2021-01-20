import pytest

from titan.run_titan import *


@pytest.mark.unit
def test_sweep_range():
    assert sweep_range("param:1:2") == {
        "start": 1,
        "stop": 2,
        "step": 1,
        "param": "param",
    }
    assert sweep_range("param:1.5:2:0.5") == {
        "start": 1.5,
        "stop": 2.0,
        "step": 0.5,
        "param": "param",
    }

    with pytest.raises(ValueError):
        sweep_range("param:a:b:c")


@pytest.mark.unit
def test_drange():
    start = 1
    step = 0.1
    stop = 2
    i = 0
    for val in drange(start, stop, step):
        i += 1
        assert val >= start and val < stop

    assert i == 10


@pytest.mark.unit
def test_setup_sweeps():
    sweeps = [
        {"start": 1.5, "stop": 2.3, "step": 0.5, "param": "param1"},
        {"start": 1, "stop": 4, "step": 1, "param": "param2"},
    ]

    defs = setup_sweeps(sweeps)

    assert len(defs) == 6


@pytest.mark.unit
def test_setup_seepfile():
    sweepfile = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "sweepfile.csv"
    )
    rows = None
    defs = setup_sweeps_file(sweepfile, rows)

    assert len(defs) == 7
    assert defs[0] == {"model.seed.run": 0, "model.seed.ppl": 0}

    rows = "2:3"
    defs = setup_sweeps_file(sweepfile, rows)

    assert len(defs) == 2
    assert defs[0] == {"model.seed.run": 1, "model.seed.ppl": 1}
