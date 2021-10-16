import os
import pathlib
import subprocess
import time

import numpy as np


def test_medicc_help_box():
    "Just testing that medicc can be started"
    process = subprocess.Popen(['python', "medicc2", "--help"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    assert process.returncode == 0


def test_medicc_with_example():
    "Testing small example"
    process = subprocess.Popen(['python', "medicc2", "examples/simple_example/simple_example.tsv", 
                                "examples/test_output"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['simple_example_cn_profiles.pdf', 'simple_example_final_cn_profiles.tsv',
                      'simple_example_final_tree.new', 'simple_example_final_tree.png',
                      'simple_example_final_tree.txt', 'simple_example_final_tree.xml',
                      'simple_example_nj_tree.new', 'simple_example_nj_tree.png',
                      'simple_example_nj_tree.txt', 'simple_example_nj_tree.xml',
                      'simple_example_pdm_total.tsv', 'simple_example_summary.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output/', f)) for f in expected_files]
    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! \nMissing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])


def test_medicc_with_example_total_copy_numbers():
    "Testing small example"
    process = subprocess.Popen(['python', "medicc2", "examples/simple_example/simple_example.tsv", 
                                "examples/test_output", "--total-copy-numbers", 
                                "--input-allele-columns", "cn_a"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['simple_example_cn_profiles.pdf', 'simple_example_final_cn_profiles.tsv',
                      'simple_example_final_tree.new', 'simple_example_final_tree.png',
                      'simple_example_final_tree.txt', 'simple_example_final_tree.xml',
                      'simple_example_nj_tree.new', 'simple_example_nj_tree.png',
                      'simple_example_nj_tree.txt', 'simple_example_nj_tree.xml',
                      'simple_example_pdm_total.tsv', 'simple_example_summary.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output/', f))
                       for f in expected_files]
    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! \nMissing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])



def test_medicc_with_bootstrap():
    "Testing bootstrap workflow"
    process = subprocess.Popen(['python', "medicc2",
                                "examples/simple_example/simple_example.tsv", "examples/test_output",
                                "--bootstrap-nr", "5"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    support_tree_exists = os.path.isfile('examples/test_output/simple_example_support_tree.new')
    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert support_tree_exists, "Support tree file was not created"
