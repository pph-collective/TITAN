import pytest

import os

def test_model_runs():
    g = globals().copy()
    g["__name__"] = "__main__"

    exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")).read(), g)
    assert True
