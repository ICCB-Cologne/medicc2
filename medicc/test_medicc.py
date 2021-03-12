import pathlib
import subprocess
import time

import numpy as np
import pytest

print(pathlib.Path(__file__).parent.parent.absolute())

def test_medicc_help_box():
    "Just testing that medicc can be started"
    process = subprocess.Popen(['python', "medicc2.py", "--help"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet, let's wait some
        time.sleep(0.5)

    assert process.returncode == 0


def test_medicc_with_example():
    "Just testing that medicc can be started"
    process = subprocess.Popen(['python', "medicc2.py", "examples/example1/example1.tsv", "examples/test_output"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet, let's wait some
        time.sleep(0.5)

    # TODO: Check which files are created
    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0
