import os
import pathlib
import subprocess
import time

import numpy as np
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
                      'simple_example_final_tree.xml', 'simple_example_pairwise_distances.tsv',
                      'simple_example_summary.tsv', 'simple_example_copynumber_events_df.tsv',
                      'simple_example_events_overlap.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output/', f)) for f in expected_files]
    subprocess.Popen(["rm", "examples/test_output", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! \nMissing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])


def test_medicc_with_example_total_copy_numbers():
    "Testing small example"
    process = subprocess.Popen(['python', "medicc2", "examples/simple_example/simple_example.tsv", 
                                "examples/test_output_total_cn", "--total-copy-numbers", 
                                "--input-allele-columns", "cn_a"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['simple_example_cn_profiles.pdf', 'simple_example_final_cn_profiles.tsv',
                      'simple_example_final_tree.new', 'simple_example_final_tree.png',
                      'simple_example_final_tree.xml', 'simple_example_pairwise_distances.tsv',
                      'simple_example_summary.tsv', 'simple_example_copynumber_events_df.tsv',
                      'simple_example_events_overlap.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output_total_cn/', f))
                       for f in expected_files]
    subprocess.Popen(["rm", "examples/test_output_total_cn", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! \nMissing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])


def test_medicc_with_multiple_cores():
    "Testing small example"
    process = subprocess.Popen(['python', "medicc2", "examples/simple_example/simple_example.tsv", 
                                "examples/test_output_multiple_cores", "-j 4"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = ['simple_example_cn_profiles.pdf', 'simple_example_final_cn_profiles.tsv',
                      'simple_example_final_tree.new', 'simple_example_final_tree.png',
                      'simple_example_final_tree.xml', 'simple_example_pairwise_distances.tsv',
                      'simple_example_summary.tsv', 'simple_example_copynumber_events_df.tsv',
                      'simple_example_events_overlap.tsv']
    all_files_exist = [os.path.isfile(os.path.join('examples/test_output_multiple_cores/', f))
                       for f in expected_files]
    subprocess.Popen(["rm", "examples/test_output_multiple_cores", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert np.all(all_files_exist), "Some files were not created! \nMissing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])



def test_medicc_with_bootstrap():
    "Testing bootstrap workflow"
    process = subprocess.Popen(['python', "medicc2",
                                "examples/simple_example/simple_example.tsv",
                                "examples/test_output_bootstrap",
                                "--bootstrap-nr", "5"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    support_tree_exists = os.path.isfile('examples/test_output_bootstrap/simple_example_support_tree.new')
    subprocess.Popen(["rm", "examples/test_output_bootstrap", "-rf"])

    assert process.returncode == 0, 'Error while running MEDICC'
    assert support_tree_exists, "Support tree file was not created"


gundem_et_al_2015_patients = ['PTX004', 'PTX005', 'PTX006', 'PTX007', 'PTX008', 
                              'PTX009', 'PTX010', 'PTX011', 'PTX012', 'PTX013']
@pytest.mark.parametrize("patient", gundem_et_al_2015_patients)
def test_gundem_et_al_2015(patient):
    "Testing if running of all Gundem data works"

    output_dir = f"examples/test_output_{patient}"

    process = subprocess.Popen(['python', "medicc2",
                                f"examples/gundem_et_al_2015/{patient}_input_df.tsv",
                                output_dir],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    expected_files = [f'{patient}_input_df_cn_profiles.pdf', f'{patient}_input_df_final_cn_profiles.tsv',
                      f'{patient}_input_df_final_tree.new', f'{patient}_input_df_final_tree.png',
                      f'{patient}_input_df_final_tree.xml', f'{patient}_input_df_pairwise_distances.tsv',
                      f'{patient}_input_df_summary.tsv', f'{patient}_input_df_copynumber_events_df.tsv',
                      f'{patient}_input_df_events_overlap.tsv']

    all_files_exist = [os.path.isfile(os.path.join(output_dir, f))
                       for f in expected_files]
    subprocess.Popen(["rm", output_dir, "-rf"])

    assert process.returncode == 0, f'Error while running MEDICC for Gundem et al patient {patient}'
    assert np.all(all_files_exist), "Some files were not created! \nMissing files are: {}".format(
        np.array(expected_files)[~np.array(all_files_exist)])


all_notebooks = [x for x in os.listdir('notebooks') if '.py' in x]
@pytest.mark.parametrize("notebook", all_notebooks)
def test_all_notebooks(notebook):
    "Testing if all notebooks work"

    process = subprocess.Popen(['python', 
                                f"notebooks/{notebook}"],
                               stdout=subprocess.PIPE,
                               cwd=pathlib.Path(__file__).parent.parent.absolute())

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    assert process.returncode == 0, f'Error while running notebook {notebook}'


all_files = [x for x in os.listdir('Figures_Kaufmann_et_al_2021') if '.py' in x]
@pytest.mark.parametrize("cur_file", all_files)
def test_figures_for_publication(cur_file):
    "Testing if all Figure scripts work"
    print(os.path.join(pathlib.Path(__file__).parent.parent.absolute(),
                       "Figures_Kaufmann_et_al_2021"))

    process = subprocess.Popen(['python',
                                f"{cur_file}"],
                                stdout=subprocess.PIPE,
                                cwd=os.path.join(pathlib.Path(__file__).parent.parent.absolute(),
                                                 "Figures_Kaufmann_et_al_2021"))

    while process.poll() is None:
        # Process hasn't exited yet
        time.sleep(0.5)

    assert process.returncode == 0, f'Error while running {cur_file}'
