import pytest

import os
import subprocess

@pytest.mark.integration
def test_model_runs():
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml")

    subprocess.call([f, f"-p {param_file}"])
    assert True
