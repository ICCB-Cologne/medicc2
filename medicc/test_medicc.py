import os
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
    "Testing small example"
    process = subprocess.Popen(['python', "medicc2.py", "examples/simple_example/simple_example.tsv", "examples/test_output"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet, let's wait some
        time.sleep(0.5)

    # TODO: Check which files are created
    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0


def test_medicc_with_bootstrap():
    "Testing bootstrap workflow"
    process = subprocess.Popen(['python', "medicc2.py",
                                "examples/simple_example/simple_example.tsv", "examples/test_output",
                                "--bootstrap-nr", "5"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet, let's wait some
        time.sleep(0.5)

    support_tree_exists = os.path.isfile('examples/test_output/simple_example_support_tree.new')

    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert support_tree_exists, "Support tree file was not created"
