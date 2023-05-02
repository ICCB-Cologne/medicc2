import os
import pathlib
import subprocess
import time

import numpy as np
import pandas as pd
import pytest


def test_medicc_help_box():
    "Just testing that medicc can be started"
    process = subprocess.Popen(['python', "medicc2", "--help"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    assert process.returncode == 0


def test_medicc_with_simple_example():
    "Testing small example"
    output_dir = 'examples/test_output'
    process = subprocess.Popen(['python', "medicc2", "examples/simple_example/simple_example.tsv", 
                                output_dir, "--plot", "both", "--events", "--chromosomes-bed",
                                "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['simple_example_cn_profiles.pdf', 'simple_example_final_cn_profiles.tsv',
                      'simple_example_final_tree.new', 'simple_example_final_tree.png',
                      'simple_example_final_tree.xml', 'simple_example_pairwise_distances.tsv',
                      'simple_example_summary.tsv', 'simple_example_copynumber_events_df.tsv',
                      'simple_example_events_overlap.tsv', 'simple_example_branch_lengths.tsv',
                      'simple_example_cn_profiles_heatmap.pdf']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output/', f)) for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'simple_example')
    output_df = pd.read_csv(os.path.join(output_dir, "simple_example_final_cn_profiles.tsv"), sep='\t')
    events_df = pd.read_csv(os.path.join(output_dir, "simple_example_copynumber_events_df.tsv"), sep='\t')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"

    assert output_df['is_gain'].sum() == 7, f"Number of gained segments in _final_cn_profiles.tsv is not 7 but {output_df['is_gain'].sum()}"
    assert output_df['is_loss'].sum() == 5, f"Number of lost segments in _final_cn_profiles.tsv is not 5 but {output_df['is_loss'].sum()}"

    assert (events_df['type'] == 'gain').sum() == 4, f"Number of gains in events_df is not 4 but {(events_df['type'] == 'gain').sum()}"
    assert (events_df['type'] == 'loss').sum() == 3, f"Number of losses in events_df is not 3 but {(events_df['type'] == 'loss').sum()}"

def test_medicc_with_testing_example():
    "Testing testing example"
    output_dir = 'examples/test_output'
    process = subprocess.Popen(['python', "medicc2", "examples/testing_example/testing_example.tsv", 
                                output_dir, "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['testing_example_cn_profiles.pdf', 'testing_example_final_cn_profiles.tsv',
                      'testing_example_final_tree.new', 'testing_example_final_tree.png',
                      'testing_example_final_tree.xml', 'testing_example_pairwise_distances.tsv',
                      'testing_example_summary.tsv', 'testing_example_copynumber_events_df.tsv',
                      'testing_example_events_overlap.tsv', 'testing_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output/', f)) for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'testing_example')
    output_df = pd.read_csv(os.path.join(output_dir, "testing_example_final_cn_profiles.tsv"), sep='\t')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"

    assert output_df['is_gain'].sum() == 187, f"Number of gains in _final_cn_profiles.tsv is not 187 but {output_df['is_gain'].sum()}"
    assert output_df['is_loss'].sum() == 170, f"Number of losses in _final_cn_profiles.tsv is not 170 but {output_df['is_loss'].sum()}"


def test_medicc_with_OV03_04():
    "Testing testing example"
    output_dir = 'examples/test_output'
    process = subprocess.Popen(['python', "medicc2", "examples/OV03-04/OV03-04_descr.txt", 
                                output_dir, "-i", "fasta", "--normal-name", "OV03-04_diploid",
                                "--plot", "both", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['OV03-04_descr_cn_profiles.pdf', 'OV03-04_descr_final_cn_profiles.tsv',
                      'OV03-04_descr_final_tree.new', 'OV03-04_descr_final_tree.png',
                      'OV03-04_descr_final_tree.xml', 'OV03-04_descr_pairwise_distances.tsv',
                      'OV03-04_descr_summary.tsv', 'OV03-04_descr_copynumber_events_df.tsv',
                      'OV03-04_descr_events_overlap.tsv', 'OV03-04_descr_branch_lengths.tsv',
                      'OV03-04_descr_cn_profiles_heatmap.pdf']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output/', f)) for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'OV03-04_descr')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_testing_example_total_copy_numbers():
    "Testing small example"
    output_dir = 'examples/test_output_total_cn'
    process = subprocess.Popen(['python', "medicc2", "examples/testing_example/testing_example.tsv", 
                                output_dir, "--total-copy-numbers", 
                                "--input-allele-columns", "cn_a", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['testing_example_cn_profiles.pdf', 'testing_example_final_cn_profiles.tsv',
                      'testing_example_final_tree.new', 'testing_example_final_tree.png',
                      'testing_example_final_tree.xml', 'testing_example_pairwise_distances.tsv',
                      'testing_example_summary.tsv', 'testing_example_copynumber_events_df.tsv',
                      'testing_example_events_overlap.tsv', 'testing_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join(output_dir, f))
                       for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'testing_example')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_testing_example_parallelization():
    "Testing small example"
    output_dir = 'examples/test_output_parallelization'
    process = subprocess.Popen(['python', "medicc2", "examples/testing_example/testing_example.tsv", 
                                output_dir, "--n-cores", "4", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['testing_example_cn_profiles.pdf', 'testing_example_final_cn_profiles.tsv',
                      'testing_example_final_tree.new', 'testing_example_final_tree.png',
                      'testing_example_final_tree.xml', 'testing_example_pairwise_distances.tsv',
                      'testing_example_summary.tsv', 'testing_example_copynumber_events_df.tsv',
                      'testing_example_events_overlap.tsv', 'testing_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join(output_dir, f))
                       for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'testing_example')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_testing_example_parallelization():
    "Testing small example"
    output_dir = 'examples/test_output_parallelization'
    process = subprocess.Popen(['python', "medicc2", "examples/testing_example/testing_example.tsv", 
                                output_dir, "--n-cores", "4", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['testing_example_cn_profiles.pdf', 'testing_example_final_cn_profiles.tsv',
                      'testing_example_final_tree.new', 'testing_example_final_tree.png',
                      'testing_example_final_tree.xml', 'testing_example_pairwise_distances.tsv',
                      'testing_example_summary.tsv', 'testing_example_copynumber_events_df.tsv',
                      'testing_example_events_overlap.tsv', 'testing_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join(output_dir, f))
                       for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'testing_example')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_testing_example_nowgd():
    "Testing small example"
    output_dir = 'examples/test_output_nowgd'
    process = subprocess.Popen(['python', "medicc2", "examples/testing_example/testing_example.tsv", 
                                output_dir, "--no-wgd", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['testing_example_cn_profiles.pdf', 'testing_example_final_cn_profiles.tsv',
                      'testing_example_final_tree.new', 'testing_example_final_tree.png',
                      'testing_example_final_tree.xml', 'testing_example_pairwise_distances.tsv',
                      'testing_example_summary.tsv', 'testing_example_copynumber_events_df.tsv',
                      'testing_example_events_overlap.tsv', 'testing_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join(output_dir, f))
                        for f in expected_files]

    nr_events, tree_size = get_number_of_events(output_dir, 'testing_example')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_testing_example_WGD_x2():
    "Testing small example"
    output_dir = 'examples/test_output_wgd_x2'
    process = subprocess.Popen(['python', "medicc2", "examples/testing_example/testing_example.tsv", 
                                output_dir, "--wgd-x2", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['testing_example_cn_profiles.pdf', 'testing_example_final_cn_profiles.tsv',
                      'testing_example_final_tree.new', 'testing_example_final_tree.png',
                      'testing_example_final_tree.xml', 'testing_example_pairwise_distances.tsv',
                      'testing_example_summary.tsv', 'testing_example_copynumber_events_df.tsv',
                      'testing_example_events_overlap.tsv', 'testing_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join(output_dir, f))
                       for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'testing_example')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_multiple_cores():
    "Testing small example"
    output_dir = 'examples/test_output_multiple_cores'
    process = subprocess.Popen(['python', "medicc2", "examples/simple_example/simple_example.tsv", 
                                output_dir, "--n-cores", "4", "--events", "--chromosomes-bed", "default", "--regions-bed", "default"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['simple_example_cn_profiles.pdf', 'simple_example_final_cn_profiles.tsv',
                      'simple_example_final_tree.new', 'simple_example_final_tree.png',
                      'simple_example_final_tree.xml', 'simple_example_pairwise_distances.tsv',
                      'simple_example_summary.tsv', 'simple_example_copynumber_events_df.tsv',
                      'simple_example_events_overlap.tsv', 'simple_example_branch_lengths.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output_multiple_cores/', f))
                       for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, 'simple_example')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


def test_medicc_with_bootstrap():
    "Testing bootstrap workflow"
    output_dir = 'examples/test_output_bootstrap'
    process = subprocess.Popen(['python', "medicc2",
                                "examples/simple_example/simple_example.tsv",
                                output_dir,
                                "--bootstrap-nr", "5"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    support_tree_exists = os.path.isfile('examples/test_output_bootstrap/simple_example_support_tree.new')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert support_tree_exists, "Support tree file was not created"


gundem_et_al_2015_patients = ['PTX004', 'PTX005', 'PTX006', 'PTX007', 'PTX008', 
                              'PTX009', 'PTX010', 'PTX011', 'PTX012', 'PTX013']
extra_condition = ['normal', 'no_wgd', 'total_cn', 'wgd_x2']
@pytest.mark.parametrize("patient", gundem_et_al_2015_patients)
@pytest.mark.parametrize("extra_condition", extra_condition)
def test_gundem_et_al_2015(patient, extra_condition):
    "Testing if running of all Gundem data works"

    output_dir = f"examples/test_output_{patient}"
    command = ['python', "medicc2", f"examples/gundem_et_al_2015/{patient}_input_df.tsv", output_dir,
               "--events", "--chromosomes-bed", "default", "--regions-bed", "default"]
    if extra_condition == 'normal':
        pass
    elif extra_condition == 'no_wgd':
        command.append('--no-wgd')
    elif extra_condition == 'total_cn':
        command += ['--total-copy-numbers', '--input-allele-columns', 'cn_a']
    elif extra_condition == 'wgd_x2':
        command.append('--wgd-x2')

    command += ["--events", "--chromosomes-bed", "default", "--regions-bed", "default"]
        
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = [f'{patient}_input_df_cn_profiles.pdf', f'{patient}_input_df_final_cn_profiles.tsv',
                      f'{patient}_input_df_final_tree.new', f'{patient}_input_df_final_tree.png',
                      f'{patient}_input_df_final_tree.xml', f'{patient}_input_df_pairwise_distances.tsv',
                      f'{patient}_input_df_summary.tsv', f'{patient}_input_df_copynumber_events_df.tsv',
                      f'{patient}_input_df_events_overlap.tsv', f'{patient}_input_df_branch_lengths.tsv']

    all_files_exist = [os.path.isfile(os.path.join(output_dir, f)) for f in expected_files]
    nr_events, tree_size = get_number_of_events(output_dir, f'{patient}_input_df')
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, f'Error while running MEDICC for Gundem et al patient {patient}'
    assert np.all(all_files_exist), "Some files were not created! Missing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])
    assert nr_events == tree_size, f"Number of events is {nr_events}, but tree size is {tree_size}"


all_ipynb_notebooks = [x for x in os.listdir('notebooks') if '.ipynb' in x]
@pytest.mark.parametrize("notebook", all_ipynb_notebooks)
def test_all_ipynb_notebooks(notebook):
    "Testing if all notebooks (with ending .ipynb) work"

    process = subprocess.Popen([f'ipython -c "%run {notebook}"'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               shell=True,
                               cwd=os.path.join(pathlib.Path(__file__).parent.parent.absolute(), 'notebooks'))

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    assert process.returncode == 0, f'Error while running notebook {notebook}: {process.stderr.read()}'


def get_number_of_events(output_dir, file_prefix):
    with open(os.path.join(output_dir, f"{file_prefix}_copynumber_events_df.tsv"), 'r') as f:
        events = f.readlines()
    nr_events = len(events) - 1

    with open(os.path.join(output_dir, f"{file_prefix}_summary.tsv"), 'r') as f:
        summary = f.readlines()
    tree_size = float(summary[2].split('\t')[1].rstrip())

    return nr_events, tree_size
