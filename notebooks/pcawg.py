#%% imports
import os
import glob
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import imp
import numpy as np
import pandas as pd
import scipy as sp
import sklearn as skl
import sklearn.metrics
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
import plotnine as pn
import fstlib
import medicc


## REQUIRED FILES TO RUN THIS:
## https://dcc.icgc.org/releases/PCAWG/consensus_cnv/consensus.20170119.somatic.cna.annotated.tar.gz
## https://dcc.icgc.org/releases/PCAWG/consensus_cnv/consensus.20170217.purity.ploidy.txt.gz
## https://dcc.icgc.org/releases/PCAWG/data_releases/latest/pcawg_sample_shee.tsv
## Meta data from: https://dcc.icgc.org/releases/PCAWG/clinical_and_histology

#%%
LOAD = True

#%% FSTs
T_wgd_asymm = fstlib.Fst.read('../objects/wgd_asymm.fst')
T_no_wgd_asymm = fstlib.Fst.read('../objects/no_wgd_asymm.fst')

#%% folders and HDF5 store
pcawg_folder = '/mnt/d/data/Documents/Projects/project-pcawg'
data_folder = os.path.join(pcawg_folder, 'consensus.20170119.somatic.cna.annotated')
metadata_folder = os.path.join(pcawg_folder, 'metadata')

#%% open HDF5 store
hdfstore = pd.HDFStore(os.path.join(pcawg_folder, 'pcawg_scna.hdf5'))

#%% metadata
if LOAD:
    meta = hdfstore['metadata/simple']
    tumour_types = hdfstore['metadata/tumour_types']
else:
    ## sample master table
    pmd_samples = pd.read_csv(os.path.join(metadata_folder, 'pcawg_sample_sheet.tsv'), sep='\t')
    pmd_samples.set_index('aliquot_id', inplace=True)
    pmd_samples.sort_index(inplace=True)
    specimen_details = pmd_samples.dcc_specimen_type.str.split("-", expand=True)
    specimen_details.columns=["dcc_specimen_type1", "dcc_specimen_type2"]
    specimen_details['dcc_specimen_type1'] = specimen_details['dcc_specimen_type1'].str.strip()
    specimen_details['dcc_specimen_type2'] = specimen_details['dcc_specimen_type2'].str.strip()
    pmd_samples = pmd_samples.join(specimen_details)

    ## histology table
    pmd_histology = pd.read_excel(os.path.join(metadata_folder, 'pcawg_specimen_histology_August2016_v9.xlsx'), dtype='str')
    pmd_histology['percentage_cellularity'] = pmd_histology['percentage_cellularity'].astype('float')
    pmd_histology.set_index('icgc_sample_id', inplace=True)
    pmd_histology.sort_index(inplace=True)
    pmd_histology = pmd_histology.drop(pmd_histology.columns[pmd_histology.columns.str.startswith('Unnamed')], axis=1)

    ## clinical table
    pmd_clinical = pd.read_excel(os.path.join(metadata_folder, 'pcawg_donor_clinical_August2016_v9.xlsx'), dtype='str')
    pmd_clinical = pmd_clinical.drop(pmd_clinical.columns[pmd_clinical.columns.str.startswith('Unnamed')], axis=1)
    pmd_clinical.columns = pmd_clinical.columns.str.replace('#', '').str.strip()
    pmd_clinical.set_index('donor_unique_id', inplace=True)
    pmd_clinical.sort_index(inplace=True)

    ## purity / ploidy table from consensus SCNA callset
    pmd_purity = pd.read_csv(os.path.join(metadata_folder, 'consensus.20170217.purity.ploidy.txt'), sep='\t')
    pmd_purity.rename({'samplename':'aliquot_id'}, inplace=True, axis=1)
    pmd_purity.set_index('aliquot_id', inplace=True)

    ## glossary (needs fixing...)
    pmd_glossary = pd.read_excel(os.path.join(metadata_folder, 'pcawg-glossary-colour-references.xlsx'), dtype='str', sheet_name='Glossary')    
    pmd_glossary = pmd_glossary.drop(pmd_glossary.columns[pmd_glossary.columns.str.startswith('Unnamed')], axis=1)
    pmd_glossary['Term'] = pmd_glossary['Term'].str.replace('Aeno', 'Adeno')
    pmd_glossary['Term'] = pmd_glossary['Term'].str.replace('Bladder', 'Bladder-TCC')
    pmd_glossary['Term'] = pmd_glossary['Term'].str.replace('Bone-Leiomyo', 'SoftTissue-Leiomyo')
    pmd_glossary['Term'] = pmd_glossary['Term'].str.replace('Bone-Cart', 'Bone-Benign')
    pmd_glossary.loc[pmd_glossary.shape[0], :] = ['Tumour Subtype', 'SoftTissue-Liposarc', '-', 'SoftTissue-Liposarc', '?']

    ## colors
    pmd_colors = pd.read_excel(os.path.join(metadata_folder, 'pcawg-glossary-colour-references.xlsx'), dtype='str', sheet_name='ColourReferenceTable_2017-08-23')
    pmd_colors = pmd_colors.drop(pmd_colors.columns[pmd_colors.columns.str.startswith('Unnamed')], axis=1)
    pmd_colors['Scheme'] = pmd_colors['Scheme'].fillna(method='ffill')
    pmd_colors = pmd_colors.loc[~pmd_colors.iloc[:,1:].isnull().all(axis=1),:]
    pmd_colors['Term'] = pmd_colors['Term'].str.replace('Bone-Cart', 'Bone-Benign')

    ## extract tumour types from glossary and colors
    tumour_types = pmd_glossary.query("Context=='Tumour Subtype'")[['Term', 'Definition']]
    tumour_types.set_index('Term', inplace=True)
    tumour_types = tumour_types.join(pmd_colors.query("Scheme=='Tumour.Subtype' or Scheme=='Low Sample Size Cohorts'").drop('Scheme', axis=1).set_index('Term'))
    tumour_types.index.name = 'histology_abbreviation'
    tumour_types.sort_index(inplace=True)
    tumour_types['histology_useable'] = True
    tumour_types['histology_useable'] = tumour_types['histology_useable'].astype('bool')
    tumour_types.loc['Bone-Benign', 'histology_useable'] = False
    tumour_types.loc['Breast-LobularCA', 'histology_useable'] = False
    tumour_types.loc['Breast-DCIS', 'histology_useable'] = False
    tumour_types.loc['Lymph-NOS', 'histology_useable'] = False
    tumour_types.loc['Myeloid-MDS', 'histology_useable'] = False
    tumour_types.loc['Cervix-AdenoCA', 'histology_useable'] = False

    meta = pmd_samples[['icgc_sample_id', 'donor_unique_id', 'library_strategy', 'dcc_specimen_type1', 'dcc_specimen_type2']]
    meta = meta.merge(pmd_histology[['histology_abbreviation', 'histology_tier1', 'histology_tier2', 'histology_tier3', 'histology_tier4', 'percentage_cellularity']], left_on='icgc_sample_id', right_index=True, how='inner')
    meta = meta.loc[~meta.histology_abbreviation.isnull(),:]
    meta = meta.merge(tumour_types, left_on='histology_abbreviation', right_index=True, how='left')
    meta.rename({'Definition':'histology_full_name'}, inplace=True, axis=1)
    meta = meta.merge(pmd_clinical, how='left', left_on='donor_unique_id', right_index=True)
    meta = meta.join(pmd_purity)
    meta['wgd_uncertain'] = meta['wgd_uncertain'].astype('bool')

    hdfstore.put('metadata/samples', pmd_samples, format='table')
    hdfstore.put('metadata/histology', pmd_histology, format='table')
    hdfstore.put('metadata/clinical', pmd_clinical, format='table')
    hdfstore.put('metadata/purity', pmd_purity, format='table')
    hdfstore.put('metadata/colors', pmd_colors, format='table')
    hdfstore.put('metadata/tumour_types', tumour_types, format='table')
    hdfstore.put('metadata/simple', meta, format='table')

#%% WGD (Haase)
if LOAD:
    wgd = hdfstore['wgd']
else:
    wgd = pd.read_csv(os.path.join(pcawg_folder, 'wgd.status.txt'), sep='\t')
    wgd.rename({'samplename':'sample_id', 'ploidy':'ploidy_haase'}, inplace=True, axis=1)
    wgd.set_index('sample_id', inplace=True)

    hdfstore.put('wgd', wgd, format='table')

# %% CNPs
if LOAD:
    dat = hdfstore['scna']
else:
    ## load CNPs
    files = glob.glob(data_folder + '/*.txt')
    index = []
    data = []
    for f in files:
        pid = os.path.basename(f).split('.')[0]
        try:
            df = pd.read_csv(f, sep='\t', index_col=['chromosome', 'start', 'end'])
            index.append(pid)
            data.append(df)
        except ValueError:
            print(f)

    dat = pd.concat(data, keys=index, names=['sample_id', 'chrom', 'start', 'end'])
    dat = dat[['major_cn', 'minor_cn']]
    dat = dat.loc[~dat.isnull().any(axis=1), :]
    dat['major_cn'] = np.fmin(dat['major_cn'], 8)
    dat['minor_cn'] = np.fmin(dat['minor_cn'], 8)
    dat['major_cn'] = dat['major_cn'].astype('int').astype('str')
    dat['minor_cn'] = dat['minor_cn'].astype('int').astype('str')

    ## phase CNPs
    diploid = fstlib.Fst()
    diploid.set_input_symbols(T_wgd_asymm.input_symbols())
    diploid.set_output_symbols(T_wgd_asymm.output_symbols())
    diploid.add_state()
    diploid.set_start(0)
    diploid.set_final(0, 0)
    diploid.add_arc(0, ('1','1',0,0))
    diploid.add_arc(0, ('X','X',0,0))

    dat = dat.query("chrom != 'X' and chrom != 'Y'") ## filter X/Y chromosomes

    phasing_scores = {}
    phased_dfs = []
    for idx, df in dat.groupby('sample_id'):
        phasing_dict = medicc.create_phasing_fsa_dict_from_df(df, T_wgd_asymm.input_symbols(), 'X')
        fsa_dict_a, fsa_dict_b, scores = medicc.phase_dict(phasing_dict, T_wgd_asymm, diploid)
        phased = medicc.create_df_from_fsa_dicts(df, [fsa_dict_a, fsa_dict_b], 'X')
        phased.columns = ['cn_a', 'cn_b']
        phased_dfs.append(phased)
        phasing_scores.update(scores)

    phasing_scores = pd.Series(phasing_scores, name='dist_to_diploid')
    phased_df = pd.concat(phased_dfs)

    dat = dat.join(phased_df, how='left')

    ## store CNPs
    hdfstore.put('scna', dat, format='table')

#%% determine no_wgd distances
if LOAD:
    distances = hdfstore['distances']
else:
    dist_no_wgd_a = pd.Series({aliquot_id:float(fstlib.score(T_no_wgd_asymm, diploid, fsa)) for aliquot_id, fsa in fsa_dict_a.items()}, name='dist_no_wgd_a')
    dist_no_wgd_b = pd.Series({aliquot_id:float(fstlib.score(T_no_wgd_asymm, diploid, fsa)) for aliquot_id, fsa in fsa_dict_b.items()}, name='dist_no_wgd_b')
    dist_no_wgd_total = dist_no_wgd_a + dist_no_wgd_b
    dist_no_wgd_total.name='dist_no_wgd_total'

    ## determine wgd distances (are in sum equivalent to the phasing scores)
    dist_wgd_a = pd.Series({aliquot_id:float(fstlib.score(T_wgd_asymm, diploid, fsa)) for aliquot_id, fsa in fsa_dict_a.items()}, name='dist_wgd_a')
    dist_wgd_b = pd.Series({aliquot_id:float(fstlib.score(T_wgd_asymm, diploid, fsa)) for aliquot_id, fsa in fsa_dict_b.items()}, name='dist_wgd_b')
    dist_wgd_total = dist_wgd_a + dist_wgd_b
    dist_wgd_total.name='dist_wgd_total'

    distances = pd.concat([dist_no_wgd_a, dist_no_wgd_b, dist_no_wgd_total, dist_wgd_a, dist_wgd_b, dist_wgd_total], axis=1)
    distances['dist_a_diff'] = distances['dist_no_wgd_a'] - distances['dist_wgd_a']
    distances['dist_b_diff'] = distances['dist_no_wgd_b'] - distances['dist_wgd_b']
    distances['dist_total_diff'] = distances['dist_no_wgd_total'] - distances['dist_wgd_total']

    hdfstore.put('distances', distances, format='table')

# %%
hdfstore.close()
result = meta.join(distances, how='inner').join(wgd[['ploidy_haase','hom']])
result['wgd_status'] = result['wgd_status'].map({'wgd':'WGD', 'no_wgd':'No WGD'})
result['wgd_status_medicc'] = (result['dist_total_diff']>=2).map({True:'WGD', False:'No WGD'})
linex = np.linspace(0, result.hom.max())
linea = -2
lineb = -1
linec = 2.9
liney = linea * linex + linec
## distance from line:
#linedist = np.abs(linea * result.hom + lineb * result.ploidy_haase + linec) / np.sqrt(linea**2 + lineb**2)
linedist = -(linea * result.hom + lineb * result.ploidy_haase + linec) / np.sqrt(linea**2 + lineb**2)
result['linedist'] = linedist
scale = 0.8

#%% test for association with uncertainty
confusion = pd.crosstab(result['wgd_status'], result['wgd_status_medicc'])
sp.stats.chi2_contingency(pd.crosstab(result.eval("wgd_status != wgd_status_medicc"), result['wgd_uncertain']))

# %% plot PCAWG MED distributions
palette = dict(zip(tumour_types.index, tumour_types.Hexadecimal))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,10))
plotdat = result.query("histology_useable").sort_values('histology_abbreviation')
sns.boxplot(x='dist_wgd_total', y='histology_abbreviation', palette=palette, data=plotdat, ax=ax)
#sns.violinplot(x='dist_wgd_total', y='histology_abbreviation', palette=palette, data=plotdat, ax=ax, scale='width')
ax.set_xlabel('MED (WGD) from diploid normal')
ax.set_ylabel('PCAWG histology')
fig.show()
fig.savefig('figures/pcawg_supp_MED_per_cancer_type.pdf', bbox_inches='tight')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,10))
plotdat = result.query("histology_useable").sort_values('histology_abbreviation')
sns.boxplot(x='dist_total_diff', y='histology_abbreviation', palette=palette, data=plotdat, ax=ax)
#sns.violinplot(x='dist_wgd_total', y='histology_abbreviation', palette=palette, data=plotdat, ax=ax, scale='width')
ax.set_xlabel('MEDICC2 WGD score')
ax.set_ylabel('PCAWG histology')
fig.show()
fig.savefig('figures/pcawg_supp_MEDICC2_WGD_score_per_cancer_type.pdf', bbox_inches='tight')


# %% PCAWG Figure
fig, ax = plt.subplots(figsize=(8 * scale,6 * scale))
sns.scatterplot('hom', 'ploidy_haase', data=result, hue='wgd_status', ax=ax)
ax.get_legend().set_title('WGD status')
ax.set_xlabel('Fraction of genome with LOH')
ax.set_ylabel('Ploidy')
ax.plot(linex, liney, '--', color='grey')
fig.show()
fig.savefig('figures/pcawg_supp_ploidy_vs_loh_wgd.pdf', bbox_inches='tight')

# %% PCAWG Figure with our score
fig, ax = plt.subplots(figsize=(8 * scale, 6 * scale))
sns.scatterplot('hom', 'ploidy_haase', data=result, hue='dist_total_diff', palette="viridis", ax=ax)
ax.get_legend().set_title('MEDICC2\nWGD score')
ax.set_xlabel('Fraction of genome with LOH')
ax.set_ylabel('Ploidy')
ax.plot(linex, liney, '--', color='grey')
fig.show()
fig.savefig('figures/pcawg_ploidy_vs_loh_our_score.pdf', bbox_inches='tight')

#%% distance from line:
plotdat = result
fig = plt.figure(figsize=(8 * scale, 8 * scale))
gs = fig.add_gridspec(2, 1, height_ratios=[0.2, 0.8])
gs.update(wspace=0, hspace=0.03) ## doesn't work with constrained layout
axhist = fig.add_subplot(gs[0])
#axhist.axes.get_yaxis().set_visible(False)
#axhist.spines["left"].set_visible(False)
axhist.spines["top"].set_visible(False)
axhist.spines["right"].set_visible(False)
sns.histplot(x='dist_total_diff', hue='wgd_status', legend=True, data=plotdat, ax=axhist)
axhist.get_legend().set_title('WGD status')

ax = fig.add_subplot(gs[1])
ax.spines['top'].set_visible(False)
sns.scatterplot('dist_total_diff', 'linedist', data=plotdat, hue='wgd_status', ax=ax, legend=False)
ax.set_xlabel('MEDICC2 WGD score')
ax.set_ylabel('PCAWG WGD score')
r = sp.stats.pearsonr(plotdat['dist_total_diff'], plotdat['linedist'])[0]**2
ax.text(0.98, 0.02, r"$r^2 = %.2f$" % r, transform=ax.transAxes, ha='right', va='bottom')
fig.show()
fig.savefig('figures/pcawg_our_score_vs_linedist.pdf', bbox_inches='tight')

