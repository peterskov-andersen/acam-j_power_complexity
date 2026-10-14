# ACAM-J EEG analysis code
Data are not included.

## Layout

```
preprocessing/          cleaning, epoching, Laplacian (MATLAB)
eeg_measures/           complexity and power measures (MATLAB)
statistics/             cluster tests, mixed models, correlation matrices (MATLAB, R)
epoch_stability/        epoch-stability / ICC analysis (Python)
neurophenomenology/
  code/                 PLSC and factor correlations (Python)
  inputs/               not included
  outputs/              created when the analysis runs
```

## 1. Preprocessing (`preprocessing/`)

| Script | Step |
|---|---|
| `RELAX_SET_PARAMETERS_SLURM.m` | Cleaning of continuous data (RELAX 2.0.1) |
| `RELAX_EPOCH_CLEAN_DATA_FOR_ANALYSIS1000_2s.m` | Epoching and epoch rejection (RELAX) |
| `laplacian_35hz.m` | Low-pass filter and surface Laplacian on the epoched data |

**Cleaning:**

- 1–80 Hz band-pass (EEGLAB `pop_eegfiltnew`, acausal) and 60 Hz Butterworth notch. Data stay at 1000 Hz.
- Bad electrodes are removed with PREP and RELAX criteria and interpolated after cleaning.
- No MWF cleaning. ICA uses PICARD; ICLabel classifies the components, and artefact components are cleaned with wavelet-enhanced ICA.

**Epoching:**

- An epoch is rejected if any of these hold:
  - absolute amplitude > 120 µV
  - improbable data > 5 MAD (single channel) or > 3 SD (all channels)
  - kurtosis > 5 SD (single channel) or > 3 SD (all channels)
- Epochs with muscle activity are not removed.

**Laplacian:** spline surface Laplacian (`pop_currentdensity`, λ = 1e-5), applied after the low-pass filter.

The epoching and Laplacian scripts were run once for each data set below, with the epoch length, input folder and low-pass changed. The files contain the settings for the 2 s data set.

| Epoch length | Low-pass before Laplacian | Used for |
|---|---|---|
| 2.048 s | 35 Hz | HFD, LZ, PE, spectral entropy |
| 5 s | none (80 Hz from cleaning) | IRASA power, aperiodic parameters |
| 10 s | 35 Hz | DFA |

## 2. EEG measures (`eeg_measures/`)

Input files are named `sub-<ID>_ses-<X>_task-<task>_run-<N>_eeg*.set`. The state is taken from the task label: `jhana1`–`jhana8`, `counting` or `memory`. Files labelled `jhana` with no number are skipped.

| Script | Measure | Output |
|---|---|---|
| `dfa.m` (calls `dfa_custom_box_size.m`) | DFA exponent | `dfaData.mat` |
| `hdf_struct.m` | Higuchi fractal dimension | `hfd_singlechan_kmax200_higuchi.mat` |
| `lzc_lz76.m` | Lempel-Ziv complexity | `lz76_singlechan2norm.mat` |
| `PeEn.m` | Permutation entropy | `pe_singlechan_m3_tau4s2_tau8s4_entropyhub.mat` |
| `spectral_entropy.m` | Spectral entropy | `specen_singlechan_broadband_spectralentropy.mat` |
| `IRASA_pwr_struct.m` | Band power (mixed, fractal, oscillatory) | `powerData_broadband_irasa.mat`, `powerData_full_spectrum_irasa.mat` |
| `slope_offset.m` | Aperiodic exponent and offset | `aperiodic_params_irasa.mat` |

Run `IRASA_pwr_struct.m` before `slope_offset.m`. The other scripts are independent.

All measures are computed per epoch and channel.

**DFA.** Each signal is mean-centred and integrated, then passed to `dfa_custom_box_size.m`.

- Detrending is quadratic (DFA2).
- Box sizes: 50 log-spaced values from 100 to 1000 samples (0.1–1 s), non-overlapping. They are set inside `dfa_custom_box_size.m`; its `box_sizes` argument is not used.
- The exponent is the slope of log10 F(n) against log10 n.

**HFD.** `Higuchi_FD` with Kmax = 200.

**LZ.** The signal is binarised at its median and LZ76 is computed with EntRate. Six features are stored:

| Feature | Definition |
|---|---|
| `C0` | raw LZ76 count |
| `c0` | C0 · log2(T) / T |
| `Cshuf_mean` | mean LZ76 of 50 time-shuffled copies |
| `LZsN` | C0 / Cshuf_mean |
| `Cphase_mean` | mean LZ76 of 100 phase-randomised surrogates |
| `LZphaseN` | C0 / Cphase_mean |

The paper uses `c0` (reported as LZsum) and `LZphaseN` (reported as LZsumN).

**PE.** EntropyHub `PermEn` with m = 3, normalised, log base 2. Two settings:

- τ = 4 samples, implemented as downsampling by 2 with τ' = 2. This is the setting used in the paper.
- τ = 8 samples, implemented as downsampling by 4 with τ' = 2.

**Spectral entropy.** MATLAB `spectralEntropy`, with these settings:

- the power spectrum of the whole epoch
- a Hamming window the length of the epoch, no overlap
- 1–35 Hz

**IRASA power.** FieldTrip `ft_freqanalysis` (`method = 'irasa'`), 1–70 Hz, `nextpow2` padding.

- Oscillatory power = mixed − fractal, with negative values set to 0. Exact zeros occurred in fewer than 0.005 % of epoch × channel band-power values.
- Band power is the sum over frequency bins: delta 1–4, theta 4–8, alpha 8–13, low beta 13–20, high beta 20–35 Hz.
- Epoch averages are stored as mean of dB (`powspctrm`) and as linear mean (`powspctrm_linear_mean`).

**Aperiodic parameters.** Linear fit of log10 power against log10 frequency. The fit uses the participant × state linear-mean fractal spectrum. Stored features:

- exponent (−slope) for 3–35, 3–13 and 13–30 Hz
- offset (intercept) for 3–35 Hz
- R² for each fit

**Output format.** FieldTrip-style structs with dimensions `rpt × chan × feature`:

- `d_<measure>.(subject).(state)`: epoch-averaged
- `d_<measure>_trials.(subject).(state)`: epoch-level

Composite states:

| Composite | Contents |
|---|---|
| `form_jhana` | jhana 1–4 |
| `formless_jhana` | jhana 5–8 |
| `all_jhana` | jhana 1–8 |
| `ctrl_combined` | counting and memory |

`form_jhana`, `formless_jhana` and `all_jhana` pool epochs across their states.

For `ctrl_combined`:

- averaged structs: the equal-weight mean of the counting and memory averages
- trial-level structs: the pooled epochs

**Epoch counts.** Each script writes a table of epoch counts per participant × state and flags counts below a cutoff. The flags are for inspection only; these scripts exclude no data.

## 3. Statistics (`statistics/`)

The R scripts read a master table with one row per participant × state and one column per measure, averaged across channels. The states are `ACAM1`–`ACAM8`, `c` (counting), `m` (memory) and `cm` (combined control). The table is not included.

**Cluster-based permutation tests (`CBPT_example.m`).** The script shows the test for spectral entropy. The same settings were used for every complexity and power measure; only the input file and feature changed.

- Contrasts: `ctrl_combined` against each of jhana 1–8, using participant-level channel values.
- Test settings:

  | Setting | Value |
  |---|---|
  | Statistic | paired t (`depsamplesT`) |
  | Cluster statistic | `maxsum` |
  | Cluster-forming threshold | α = 0.05, two-sided |
  | Cluster α | 0.025 per tail |
  | Neighbours | triangulation, at least 2 neighbouring channels |
  | Permutations | all possible |

- For each contrast the smallest cluster p is kept. Holm and Benjamini-Hochberg corrections are applied across the 8 contrasts.

**Mixed models (`power_entropy_lmm.R`).**

- Model per measure: `measure ~ state + (1 | subject)`, fitted with lme4/lmerTest (REML, bobyqa optimiser).
- Contrasts (emmeans, Satterthwaite df):
  - ACAM-J1–8 vs `cm`, with BH-FDR across the 8 contrasts
  - `c − m`, uncorrected
- Effect sizes: emmeans `eff_size`, using the residual SD and residual df.
- Also reported: likelihood-ratio test against an intercept-only model.
- The model includes the `c`, `m` and `cm` rows; `cm` is the mean of `c` and `m`.
- Output: raincloud figures and supplementary tables.

**Correlation matrices (`corr_matrix_script.R`).**

- Spearman correlations between measures, computed on ACAM-J − control differences.
- Run for each ACAM-J state separately and for the mean of all ACAM-J states.
- BH-FDR across the unique pairs in each matrix.

## 4. Epoch stability (`epoch_stability/`)

`epoch_stability_icc.py` reads the epoch-level structs from section 2. For each measure it uses the global value per epoch (mean across 128 channels).

- Reliability: ICC(1,1), projected to a K-epoch average with the Spearman-Brown formula.
- K is the minimum number of epochs per participant × state: 6 for power, 3 for DFA, 15 for the other complexity measures.

`icc2_1_analysis.py` repeats the analysis with ICC(2,1) (two-way random effects, absolute agreement). It uses the first k epochs of each participant (k = the smallest epoch count in the condition) and bootstrap CIs (300 resamples). It imports `epoch_stability_icc.py`.

## 5. Neurophenomenology (`neurophenomenology/code/`)

All settings are in `config.py`. Paths are relative to `neurophenomenology/`. Run from `neurophenomenology/code/`:

```
python 01_parse_phenomenology.py
python 02_build_analysis_table.py
python 03_plsc_pyls.py
python 04_univariate_factors.py
python 05_figures.py
```

To run a single PLSC, pass its name, e.g. `python 03_plsc_pyls.py within`.

Inputs (in `neurophenomenology/inputs/`):

- `ratings_raw/*Jhana*_ratings.txt`: ACAM-J rating files, one per session.
- `eeg_measures.csv`: one row per participant × state, with `subject`, `state` (`ACAM1`–`ACAM8`) and the EEG columns listed in `config.py`.

**Parsing (01).**

- Scores above 10 are missing-value codes and are dropped before averaging across runs.
- `Adept-035-Jhana-1_ratings.txt` is excluded, because only that participant's second session has retained EEG.

**Merging (02).** Early sensation ratings are assigned to ACAM-J 1–4 and late ratings to ACAM-J 5–8.

**PLSC (03).** pyls behavioural PLS.

- Variables:
  - X: 16 EEG measures.
  - Y: 7 state-level items.
- Data handling:
  - Missing Y values are imputed with scikit-learn `IterativeImputer` (20 iterations, seed 0).
  - Rows with missing EEG are removed.
  - EEG values are rounded to 4 decimals, as in the original analysis table.
- Permutations and bootstraps: 10,000 each, seed 999.
- Three analyses:

  | Analysis | Data | Permutation | Bootstrap |
  |---|---|---|---|
  | Between-participant | z-scored observations | unrestricted | over observations |
  | Within-participant | X and Y centred within participant | within participant | over participants (resamples truncated or repeated to the original row count, as pyls needs a fixed sample size) |
  | Run-level quality | EEG averaged across ACAM-J states, Y = 4 run-level items, one row per participant | unrestricted | over observations |

- An LV is significant if permutation p < .05 and its variance explained exceeds 100 / number of LVs.
- pyls 0.0.2 needs `permindices=True`; without it, the permutation indices are passed in place of Y.

**Factors (04).**

- Spearman correlations between each EEG measure and each ACAM-J-specific factor.
- Both are averaged over the states where the factor is rated, giving one value per participant.
- Benjamini-Hochberg FDR is applied within each factor across the 16 measures.

**Figures (05).** LV signs are arbitrary. The main LVs are flipped so that the positive pole corresponds to lower spectral power.

## Notes

- **Paths.** They are hard-coded in the MATLAB and R scripts and in `icc2_1_analysis.py`. `epoch_stability_icc.py` reads `ICC_DATA_DIR` and `ICC_OUT_DIR` from the environment.
- **Checkpoints.** HFD, LZ, PE and spectral entropy save a checkpoint per run and resume from existing checkpoints.
- **LZ random seed.** `rng(1)` seeds only the main MATLAB session, not `parfor` workers, so the surrogate-based LZ features differ slightly between runs. `C0` and `c0` are deterministic.
- **Comments.** Comments were revised after the analysis to match the code; the code itself is unchanged.

## Dependencies

MATLAB (R2023b or later for `spectralEntropy`):

- EEGLAB 2025.1.0, with the plugins:
  - RELAX 2.0.1
  - PREP 0.57.0
  - ICLabel
  - PICARD 2.0
  - currentdensity
  - FieldTrip-lite 250523
- Signal Processing, Statistics and Machine Learning, and Wavelet Toolboxes.
- Parallel Computing Toolbox, optional.
- [EntRate](https://github.com/pmediano/EntRate), for LZ
- [EntropyHub](https://github.com/MattWillFlood/EntropyHub) v2.0.0, for PE
- [Higuchi and Katz fractal dimension measures](https://www.mathworks.com/matlabcentral/fileexchange/50290) (File Exchange 50290), for HFD

R:

- readr, dplyr, tidyr, tibble, purrr, stringr
- lme4, lmerTest, emmeans, Hmisc
- ggplot2, ggdist, patchwork

Python 3:

- numpy, pandas, scipy, scikit-learn, matplotlib, h5py, mne
- [pyls](https://github.com/netneurolab/pypyls) v0.0.2. Install with `pip install git+https://github.com/netneurolab/pypyls.git`; `pip install pyls` installs an unrelated package.

## Licence

Code: MIT (see `LICENSE`), with these exceptions:

- RELAX scripts in `preprocessing/`: GPL-3.0 (Neil Bailey). Their licence headers are kept.
- `dfa_custom_box_size.m`: modified from DFluC by Yishul Wei. It remains under the original terms in its header.

Cite when using:

- RELAX: Bailey NW et al. Introducing RELAX: an automated pre-processing pipeline for cleaning EEG data, parts 1 and 2. *Clin Neurophysiol* 149 (2023).
- DFluC: Colombo MA, Wei Y, Ramautar JR, Linkenkaer-Hansen K, Tagliazucchi E, Van Someren EJW. *Front Physiol* 7:576 (2016). doi:10.3389/fphys.2016.00576


