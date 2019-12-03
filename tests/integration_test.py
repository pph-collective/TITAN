import pytest

import os


def test_model_runs():
    g = globals().copy()
    g["__name__"] = "__main__"

    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    exec(open(f).read(), g)
    assert True
