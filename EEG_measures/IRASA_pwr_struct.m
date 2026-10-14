% IRASA power per epoch and channel (5 s epochs).
% Pools epochs across runs (manual concatenation along trial dim),
% computes band & full-spectrum power using IRASA, stores:
%   - original/mixed      (conventional PSD)
%   - fractal/aperiodic   (1/f)
%   - oscillatory         (original - fractal, rectified at 0)
%
% For each component:
%   - mean-of-dB and linear means
%   - epoch counts, weighted composites, canonicalized labels
%   - QC table
%
% Backwards-compat aliases:
%   d_band        -> original/mixed band power
%   d_full        -> original/mixed full-spectrum power
%   d_band_trials -> original/mixed trial-level band power
%   d_full_trials -> original/mixed trial-level full-spectrum power

clear; clc;

VERBOSE    = true;   % set false to silence per-run summaries
USE_PARFOR = true;   % toggle parallel processing across files

%% 0) Initialize FieldTrip
if exist('ft_defaults','file'), ft_defaults; end
% EEGLAB
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab
%% 1) User parameters
dataDir  = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_5s/80Hz_laplacian';
% 5 s epochs, surface Laplacian, no low-pass beyond the 80 Hz applied during cleaning
saveBand = 'powerData_broadband_irasa.mat';
saveFull = 'powerData_full_spectrum_irasa.mat';
epsP     = 1e-20;   % tiny guard to avoid log(0) in dB

%% 2) Frequency bands (Hz)
freq_bands = { ...
  'delta',     [1 4];
  'theta',     [4 8];
  'alpha',     [8 13];
  'low_beta',  [13 20];
  'high_beta', [20 35];
};
bandNames = freq_bands(:,1);
nBands    = numel(bandNames);

%% 3) Valid state labels (canonical targets)  (for reference/order)
validStates = { ...
  'counting','memory', ...
  '1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
  '5th_jhana','6th_jhana','7th_jhana','8th_jhana'
}; %#ok<NASGU> % not used; listed for reference

%% 4) Prepare per-run maps (separate for orig / fractal / oscillatory)
runsMapBand_orig = struct();   % broadband per-run, trial-level (linear)
runsMapBand_frac = struct();
runsMapBand_osc  = struct();
runsMapFull_orig = struct();   % full-spectrum per-run, trial-level (linear)
runsMapFull_frac = struct();
runsMapFull_osc  = struct();

fileList = dir(fullfile(dataDir,'*.set'));
if isempty(fileList)
    error('No .set files found under: %s', dataDir);
end

nFiles = numel(fileList);
[~, ord] = sort({fileList.name});
fileList = fileList(ord);  % stable order

%% 5) Loop files (possibly in parallel) -> per-run IRASA results in cell array
fileResults = cell(nFiles,1);

usePar = USE_PARFOR && ~isempty(ver('parallel')) && nFiles >= 2;
if usePar
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool;
        catch
            warning('Could not start parpool; falling back to serial loop.');
            usePar = false;
        end
    end
end

if usePar
    if VERBOSE
        fprintf('Using PARFOR over %d files\n', nFiles);
    end
    parfor iF = 1:nFiles
        fileInfo = fileList(iF);
        fileResults{iF} = process_irasa_one_file(fileInfo, freq_bands, bandNames, VERBOSE);
    end
else
    if VERBOSE
        fprintf('Using serial FOR loop over %d files\n', nFiles);
    end
    for iF = 1:nFiles
        fileInfo = fileList(iF);
        fileResults{iF} = process_irasa_one_file(fileInfo, freq_bands, bandNames, VERBOSE);
    end
end

%% 5b) Merge per-file results into runsMap* structs (serial, no parfor)
for iF = 1:nFiles
    R = fileResults{iF};
    if isempty(R) || ~isfield(R,'hasData') || ~R.hasData
        continue;
    end

    subjFld = R.subjFld;
    stateFld = R.stateFld;

    % ensure subject fields exist
    if ~isfield(runsMapBand_orig,subjFld)
        runsMapBand_orig.(subjFld) = struct();
        runsMapBand_frac.(subjFld) = struct();
        runsMapBand_osc.(subjFld)  = struct();
        runsMapFull_orig.(subjFld) = struct();
        runsMapFull_frac.(subjFld) = struct();
        runsMapFull_osc.(subjFld)  = struct();
    end

    % stash per-run structs as cell array entries
    if ~isfield(runsMapBand_orig.(subjFld),stateFld)
        runsMapBand_orig.(subjFld).(stateFld) = {R.frB_orig};
        runsMapBand_frac.(subjFld).(stateFld) = {R.frB_frac};
        runsMapBand_osc.(subjFld).(stateFld)  = {R.frB_osc};
        runsMapFull_orig.(subjFld).(stateFld) = {R.frF_orig};
        runsMapFull_frac.(subjFld).(stateFld) = {R.frF_frac};
        runsMapFull_osc.(subjFld).(stateFld)  = {R.frF_osc};
    else
        runsMapBand_orig.(subjFld).(stateFld){end+1} = R.frB_orig;
        runsMapBand_frac.(subjFld).(stateFld){end+1} = R.frB_frac;
        runsMapBand_osc.(subjFld).(stateFld){end+1}  = R.frB_osc;
        runsMapFull_orig.(subjFld).(stateFld){end+1} = R.frF_orig;
        runsMapFull_frac.(subjFld).(stateFld){end+1} = R.frF_frac;
        runsMapFull_osc.(subjFld).(stateFld){end+1}  = R.frF_osc;
    end
end

%% 6) Aggregate per subject×state for EACH IRASA component
% Composite definitions reused for all components
composites = {
  'form_jhana',     {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana'};
  'formless_jhana', {'5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
  'all_jhana',      {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
                     '5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
  'ctrl_combined',  {'counting','memory'}          
};

[d_band_orig, d_full_orig, d_band_trials_orig, d_full_trials_orig] = ...
    aggregate_irasa_component(runsMapBand_orig, runsMapFull_orig, epsP, VERBOSE, composites);

[d_band_frac, d_full_frac, d_band_trials_frac, d_full_trials_frac] = ...
    aggregate_irasa_component(runsMapBand_frac, runsMapFull_frac, epsP, VERBOSE, composites);

[d_band_osc,  d_full_osc,  d_band_trials_osc,  d_full_trials_osc ] = ...
    aggregate_irasa_component(runsMapBand_osc,  runsMapFull_osc,  epsP, VERBOSE, composites);

%% 7) ctrl_combined: equal-weight mean of counting and memory (all components;
%     replaces the pooled-epoch version computed in step 6)
[d_band_orig, d_full_orig] = add_ctrl_combined(d_band_orig, d_full_orig);
[d_band_frac, d_full_frac] = add_ctrl_combined(d_band_frac, d_full_frac);
[d_band_osc,  d_full_osc ] = add_ctrl_combined(d_band_osc,  d_full_osc);

%% 8) Backwards-compatible aliases (original/mixed only)
d_band        = d_band_orig;
d_full        = d_full_orig;
d_band_trials = d_band_trials_orig;
d_full_trials = d_full_trials_orig;

%% 9) Save outputs
save(fullfile(dataDir,saveBand), ...
     'd_band','d_band_frac','d_band_osc', ...
     'd_band_trials','d_band_trials_frac','d_band_trials_osc', '-v7.3');
fprintf('Saved broadband IRASA data to %s\n', fullfile(dataDir,saveBand));

save(fullfile(dataDir,saveFull), ...
     'd_full','d_full_frac','d_full_osc', ...
     'd_full_trials','d_full_trials_frac','d_full_trials_osc', '-v7.3');
fprintf('Saved full-spectrum IRASA data to %s\n', fullfile(dataDir,saveFull));

%% 10) QC: subject×state epoch counts (from original/mixed d_band)
cutoff = 25;   % flag subject×state cells below this epoch count (QC table only)
qcTable = qc_ntrials_table(d_band, cutoff);
disp(qcTable)

bad = qcTable.flag_below_cutoff;
if any(bad)
    fprintf('\nStages below %d epochs:\n', cutoff);
    disp(qcTable(bad, :));
else
    fprintf('\nNo stages below %d epochs.\n', cutoff);
end
writetable(qcTable, fullfile(dataDir, sprintf('ntrials_qc_cutoff%d.csv', cutoff)));

%% ----------------- per-file IRASA helper (used by parfor/for) -----------------
function R = process_irasa_one_file(fileInfo, freq_bands, bandNames, VERBOSE)
% Process a single .set file with IRASA and return a struct with:
%  - hasData (logical)
%  - subjFld, stateFld (sanitized field names)
%  - frB_orig, frB_frac, frB_osc (broadband trial-level)
%  - frF_orig, frF_frac, frF_osc (full-spectrum trial-level)

    R = struct();
    R.hasData = false;  % default

    fname = fileInfo.name;
    fpath = fileInfo.folder;

    % ---- Parse subject, task and run from the filename
    tok = regexp(fname, ...
        'sub-(\d+)_ses-[^_]+_task-([^_]+)_run-(\d+)_eeg.*\.set$', 'tokens','once');
    if isempty(tok)
        warning('Skipping %s (name parse failed)', fname);
        return;
    end
    rawSub  = ['sub-' tok{1}];
    rawTask = tok{2};
    rawRun  = tok{3};
    subjFld = matlab.lang.makeValidName(rawSub);

    % ---- Infer a single label for ALL epochs from filename task
    inferredState = task_to_state(rawTask);  % -> "1st_jhana", "memory", etc., or ""
    if inferredState == ""
        if VERBOSE
            fprintf('[SKIP ] %s | task=%s (no usable state label)\n', rawSub, rawTask);
        end
        return;
    end

    % ---- Load EEG & convert to FieldTrip
    EEG     = pop_loadset(fname, fpath);
    ft_data = eeglab2fieldtrip(EEG,'preprocessing','none');

    % ---- IRASA: original & fractal spectra, then oscillatory = original - fractal
    cfg_ir            = [];
    cfg_ir.method     = 'irasa';
    cfg_ir.foilim     = [1 70];        % output range; band power uses 1-35 Hz
    cfg_ir.pad        = 'nextpow2';
    cfg_ir.keeptrials = 'yes';         % per-epoch spectra

    cfg_ir.output     = 'original';
    freq_orig         = ft_freqanalysis(cfg_ir, ft_data);

    cfg_ir.output     = 'fractal';
    freq_frac         = ft_freqanalysis(cfg_ir, ft_data);

    % oscillatory = original - fractal
    cfg_math             = [];
    cfg_math.parameter   = 'powspctrm';
    cfg_math.operation   = 'x2-x1';    % original - fractal (note argument order)
    freq_osc             = ft_math(cfg_math, freq_frac, freq_orig);

    % Set negative oscillatory power to 0
    freq_osc.powspctrm(freq_osc.powspctrm < 0) = 0;

    freqsFull = freq_orig.freq;   % shared freq axis

    % ---- Per-epoch bandpower (linear) for each component
    nE = size(freq_orig.powspctrm,1);
    nC = size(freq_orig.powspctrm,2);
    nBands = numel(bandNames);

    absPow_orig = zeros(nE,nC,nBands);
    absPow_frac = zeros(nE,nC,nBands);
    absPow_osc  = zeros(nE,nC,nBands);

    for b = 1:nBands
        fRng = freq_bands{b,2};
        idxF = freqsFull>=fRng(1) & freqsFull<=fRng(2);

        absPow_orig(:,:,b) = sum(freq_orig.powspctrm(:,:,idxF),3);
        absPow_frac(:,:,b) = sum(freq_frac.powspctrm(:,:,idxF),3);
        absPow_osc(:,:,b)  = sum(freq_osc.powspctrm(:,:,idxF),3);
    end

    labels = repmat(inferredState, nE, 1);  % one label per epoch in this file

    powFull_orig = freq_orig.powspctrm;
    powFull_frac = freq_frac.powspctrm;
    powFull_osc  = freq_osc.powspctrm;

    if VERBOSE
        fprintf('[IRASA] %s | %s run-%s : %d epochs -> %s=%d\n', ...
            rawSub, rawTask, rawRun, numel(labels), char(inferredState), numel(labels));
    end

    % ---- Build per-run trial-level structs (only 1 state per file)
    stRaw = inferredState;
    fld   = matlab.lang.makeValidName(stRaw);  % state field

    % === Broadband (trials x chan x bands), linear ===
    % original/mixed
    frB_orig = struct();
    frB_orig.label     = freq_orig.label;
    frB_orig.dimord    = 'rpt_chan_freq';
    frB_orig.freq      = 1:nBands;    % band index; human names in .bandnames
    frB_orig.bandnames = bandNames;
    frB_orig.powspctrm = absPow_orig;
    frB_orig.meta.subject  = rawSub;
    frB_orig.meta.task     = rawTask;
    frB_orig.meta.run      = rawRun;
    frB_orig.meta.filename = fname;

    % fractal
    frB_frac = frB_orig;
    frB_frac.powspctrm = absPow_frac;

    % oscillatory
    frB_osc = frB_orig;
    frB_osc.powspctrm = absPow_osc;

    % === Full-spectrum (trials x chan x freqs), linear ===
    frF_orig = struct();
    frF_orig.label     = freq_orig.label;
    frF_orig.dimord    = 'rpt_chan_freq';
    frF_orig.freq      = freqsFull;
    frF_orig.powspctrm = powFull_orig;
    frF_orig.meta.subject  = rawSub;
    frF_orig.meta.task     = rawTask;
    frF_orig.meta.run      = rawRun;
    frF_orig.meta.filename = fname;

    frF_frac           = frF_orig;
    frF_frac.powspctrm = powFull_frac;

    frF_osc            = frF_orig;
    frF_osc.powspctrm  = powFull_osc;

    % ---- Populate result struct
    R.hasData  = true;
    R.subjFld  = subjFld;
    R.stateFld = fld;
    R.frB_orig = frB_orig;
    R.frB_frac = frB_frac;
    R.frB_osc  = frB_osc;
    R.frF_orig = frF_orig;
    R.frF_frac = frF_frac;
    R.frF_osc  = frF_osc;
end

%% ----------------- helpers -----------------
function [d_band, d_full, d_band_trials, d_full_trials] = aggregate_irasa_component(runsMapBand, runsMapFull, epsP, VERBOSE, composites)
% Build subject×state averages + composites for one IRASA component.

  d_band        = struct();
  d_full        = struct();
  d_band_trials = struct();
  d_full_trials = struct();

  subs = fieldnames(runsMapBand);
  for i = 1:numel(subs)
    sfld  = subs{i};
    conds = fieldnames(runsMapBand.(sfld));
    for j = 1:numel(conds)
      cfld = conds{j};

      % ---- Broadband: manual append of trials across runs ----
      runsB = runsMapBand.(sfld).(cfld);
      [freqB_trials, runCountsB, runIDsB] = append_trials_freq(runsB, true);   % true => hasBandnames
      nTrialsB = size(freqB_trials.powspctrm,1);
      if VERBOSE
        fprintf('[BAND ] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
            sfld, pretty_state(cfld), char(strjoin(runIDsB,',')), mat2str(runCountsB), nTrialsB);
      end

      % mean of per-epoch dB  + linear mean
      avgB_db  = squeeze(mean(10*log10(freqB_trials.powspctrm + epsP), 1));  % chan x bands
      avgB_lin = squeeze(mean(freqB_trials.powspctrm, 1));                   % chan x bands (linear)

      freqB_avg                       = freqB_trials;
      freqB_avg.powspctrm             = reshape(avgB_db,  1, size(avgB_db,1),  size(avgB_db,2));
      freqB_avg.powspctrm_linear_mean = reshape(avgB_lin, 1, size(avgB_lin,1), size(avgB_lin,2));
      freqB_avg.dimord                = 'rpt_chan_freq';
      freqB_avg.units                 = 'dB';
      freqB_avg.avg_method            = 'mean_of_log10';
      freqB_avg.linear_units          = 'linear_power';
      freqB_avg.linear_avg_method     = 'mean_linear';
      freqB_avg.ntrials               = nTrialsB;

      d_band_trials.(sfld).(cfld) = freqB_trials;
      d_band.(sfld).(cfld)        = freqB_avg;

      % ---- Full-spectrum: manual append of trials across runs ----
      runsF = runsMapFull.(sfld).(cfld);
      [freqF_trials, runCountsF, runIDsF] = append_trials_freq(runsF, false);
      nTrialsF = size(freqF_trials.powspctrm,1);
      if VERBOSE
        fprintf('[FULL ] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
            sfld, pretty_state(cfld), char(strjoin(runIDsF,',')), mat2str(runCountsF), nTrialsF);
      end

      avgF_db  = squeeze(mean(10*log10(freqF_trials.powspctrm + epsP), 1));  % chan x freqs
      avgF_lin = squeeze(mean(freqF_trials.powspctrm, 1));                   % chan x freqs (linear)

      freqF_avg                       = freqF_trials;
      freqF_avg.powspctrm             = reshape(avgF_db,  1, size(avgF_db,1),  size(avgF_db,2));
      freqF_avg.powspctrm_linear_mean = reshape(avgF_lin, 1, size(avgF_lin,1), size(avgF_lin,2));
      freqF_avg.dimord                = 'rpt_chan_freq';
      freqF_avg.units                 = 'dB';
      freqF_avg.avg_method            = 'mean_of_log10';
      freqF_avg.linear_units          = 'linear_power';
      freqF_avg.linear_avg_method     = 'mean_linear';
      freqF_avg.ntrials               = nTrialsF;

      d_full_trials.(sfld).(cfld) = freqF_trials;
      d_full.(sfld).(cfld)        = freqF_avg;
    end
  end

  % ---- Composite jhana categories (trial-weighted via pooling) ----
  subs = fieldnames(d_band_trials);
  for i = 1:numel(subs)
    sfld = subs{i};
    for c = 1:size(composites,1)
      compName = composites{c,1};
      list     = composites{c,2};

      % Broadband composite
      poolB = {};
      for k = 1:numel(list)
        fld = matlab.lang.makeValidName(list{k});
        if isfield(d_band_trials.(sfld), fld), poolB{end+1} = d_band_trials.(sfld).(fld); end
      end
      if ~isempty(poolB)
        freqB_trials = append_trials_freq(poolB, true);   % pooled epochs
        nTrialsB     = size(freqB_trials.powspctrm,1);
        avgB_db      = squeeze(mean(10*log10(freqB_trials.powspctrm + epsP), 1));
        avgB_lin     = squeeze(mean(freqB_trials.powspctrm, 1));

        freqB_avg                       = freqB_trials;
        freqB_avg.powspctrm             = reshape(avgB_db,  1, size(avgB_db,1),  size(avgB_db,2));
        freqB_avg.powspctrm_linear_mean = reshape(avgB_lin, 1, size(avgB_lin,1), size(avgB_lin,2));
        freqB_avg.dimord                = 'rpt_chan_freq';
        freqB_avg.units                 = 'dB';
        freqB_avg.avg_method            = 'mean_of_log10';
        freqB_avg.linear_units          = 'linear_power';
        freqB_avg.linear_avg_method     = 'mean_linear';
        freqB_avg.ntrials               = nTrialsB;

        d_band.(sfld).(compName) = freqB_avg;
      end

      % Full-spectrum composite
      poolF = {};
      for k = 1:numel(list)
        fld = matlab.lang.makeValidName(list{k});
        if isfield(d_full_trials.(sfld), fld), poolF{end+1} = d_full_trials.(sfld).(fld); end
      end
      if ~isempty(poolF)
        freqF_trials = append_trials_freq(poolF, false);
        nTrialsF     = size(freqF_trials.powspctrm,1);
        avgF_db      = squeeze(mean(10*log10(freqF_trials.powspctrm + epsP), 1));
        avgF_lin     = squeeze(mean(freqF_trials.powspctrm, 1));

        freqF_avg                       = freqF_trials;
        freqF_avg.powspctrm             = reshape(avgF_db,  1, size(avgF_db,1),  size(avgF_db,2));
        freqF_avg.powspctrm_linear_mean = reshape(avgF_lin, 1, size(avgF_lin,1), size(avgF_lin,2));
        freqF_avg.dimord                = 'rpt_chan_freq';
        freqF_avg.units                 = 'dB';
        freqF_avg.avg_method            = 'mean_of_log10';
        freqF_avg.linear_units          = 'linear_power';
        freqF_avg.linear_avg_method     = 'mean_linear';
        freqF_avg.ntrials               = nTrialsF;

        d_full.(sfld).(compName) = freqF_avg;
      end
    end
  end
end

function [d_band, d_full] = add_ctrl_combined(d_band, d_full)
% Equal-weight average of counting + memory -> ctrl_combined

  subs = fieldnames(d_band);
  for i = 1:numel(subs)
    sfld = subs{i};

    if isfield(d_band.(sfld),'counting') && isfield(d_band.(sfld),'memory')
      A = d_band.(sfld).counting;  B = d_band.(sfld).memory;
      dB  = 0.5 * (squeeze(A.powspctrm) + squeeze(B.powspctrm));
      lin = 0.5 * (squeeze(A.powspctrm_linear_mean) + squeeze(B.powspctrm_linear_mean));
      outB                       = A;
      outB.powspctrm             = reshape(dB,  1, size(dB,1),  size(dB,2));
      outB.powspctrm_linear_mean = reshape(lin, 1, size(lin,1), size(lin,2));
      outB.avg_method            = 'mean_of_log10 (equal-weight across states)';
      outB.linear_avg_method     = 'mean_linear (equal-weight)';
      outB.ntrials               = A.ntrials + B.ntrials;
      outB.components            = {'counting','memory'};
      d_band.(sfld).ctrl_combined = outB;
    end

    if isfield(d_full.(sfld),'counting') && isfield(d_full.(sfld),'memory')
      A = d_full.(sfld).counting;  B = d_full.(sfld).memory;
      dB  = 0.5 * (squeeze(A.powspctrm) + squeeze(B.powspctrm));
      lin = 0.5 * (squeeze(A.powspctrm_linear_mean) + squeeze(B.powspctrm_linear_mean));
      outF                       = A;
      outF.powspctrm             = reshape(dB,  1, size(dB,1),  size(dB,2));
      outF.powspctrm_linear_mean = reshape(lin, 1, size(lin,1), size(lin,2));
      outF.avg_method            = 'mean_of_log10 (equal-weight across states)';
      outF.linear_avg_method     = 'mean_linear (equal-weight)';
      outF.ntrials               = A.ntrials + B.ntrials;
      outF.components            = {'counting','memory'};
      d_full.(sfld).ctrl_combined = outF;
    end
  end
end

function out = canonical_state(raw)
% Return the canonical state name as a string scalar.
    s = lower(strtrim(string(raw)));
    s = regexprep(s, '[-\s]+', '_');
    if s == "memory",   out = "memory";   return; end
    if s == "counting", out = "counting"; return; end
    tok = regexp(s, '^(?:jhana_)?([1-8])(st|nd|rd|th)?(?:_jhana)?$', 'tokens', 'once');
    if ~isempty(tok)
        n = str2double(tok{1}); ord = ["st","nd","rd","th","th","th","th","th"];
        out = string(sprintf('%d%s_jhana', n, ord(n))); return
    end
    tokw = regexp(s, '^(?:jhana_)?(first|second|third|fourth|fifth|sixth|seventh|eighth)(?:_jhana)?$', 'tokens', 'once');
    if ~isempty(tokw)
        names = ["first","second","third","fourth","fifth","sixth","seventh","eighth"];
        idx = find(strcmp(tokw{1}, names));
        ord = ["st","nd","rd","th","th","th","th","th"];
        out = string(sprintf('%d%s_jhana', idx, ord(idx))); return
    end
    out = "";
end

function out = task_to_state(taskStr)
% Map filename task (e.g., 'jhana5','memory','counting') -> canonical state string.
    t = lower(strtrim(taskStr));
    % jhana1..jhana8
    m = regexp(t,'^jhana([1-8])$','tokens','once');
    if ~isempty(m)
        out = canonical_state(m{1});   % canonical_state handles numeric -> "1st_jhana", etc.
        return
    end
    % skip plain 'jhana' (no number)
    if strcmp(t,'jhana')
        out = ""; return
    end
    % memory, counting
    out = canonical_state(t);
end

function lab = pretty_state(fld)
    base = regexprep(fld, '^x', '');
    lab  = strrep(base, '_', ' ');
end

function [out, runCounts, runIDs] = append_trials_freq(runList, hasBandnames)
% Manually concatenate trials across runs along rpt dim.
% - Validates identical channel labels and frequency axes.
% - Returns pooled trial-level struct (linear power).
% - Also returns run-wise trial counts and run IDs (strings).

    if nargin < 2, hasBandnames = false; end
    if isstruct(runList), runList = {runList}; end
    assert(~isempty(runList), 'append_trials_freq: empty run list');

    % validate compatibility and collect meta
    base    = runList{1};
    nRuns   = numel(runList);
    runCounts = zeros(1,nRuns);
    runIDs    = strings(1,nRuns);
    for k = 1:nRuns
        r = runList{k};
        assert(isequal(r.label, base.label), 'Channel labels differ across runs');
        assert(isequal(r.freq,  base.freq),  'Frequency axis differs across runs');
        runCounts(k) = size(r.powspctrm,1);
        if isfield(r,'meta') && isfield(r.meta,'run')
            runIDs(k) = string(r.meta.run);
        else
            runIDs(k) = string(k);
        end
    end

    % gather powspctrm into a cell, then cat
    P = cellfun(@(r) r.powspctrm, runList, 'UniformOutput', false);
    out             = base;                             % copy fields/axes from first run
    out.powspctrm   = cat(1, P{:});                     % rpt concatenation
    if isfield(base,'dimord'), out.dimord = base.dimord; else, out.dimord = 'rpt_chan_freq'; end
    if hasBandnames && isfield(base,'bandnames')
        out.bandnames = base.bandnames;
    end

    % optional: which run each trial came from
    try
        trlIdx = arrayfun(@(k) k*ones(runCounts(k),1), 1:nRuns, 'UniformOutput', false);
        out.trialinfo = vertcat(trlIdx{:});
    catch
        % non-fatal
    end

    % sanity check: pooled trials equal sum of per-run trials
    assert(size(out.powspctrm,1) == sum(runCounts), ...
        'Pooled ntrials mismatch: got %d, expected %d', size(out.powspctrm,1), sum(runCounts));
end

function qcTable = qc_ntrials_table(d_band, cutoff)
% Build a long table: subject | state | ntrials | flag_below_cutoff
    if nargin < 2 || isempty(cutoff), cutoff = 25; end

    % Preferred state order & a map from sanitized -> human label
    rawStates = {'counting','memory','1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
                 '5th_jhana','6th_jhana','7th_jhana','8th_jhana', ...
                 'form_jhana','formless_jhana','all_jhana','ctrl_combined'};
    sanStates = cellfun(@matlab.lang.makeValidName, rawStates, 'UniformOutput', false);
    S2H = containers.Map(sanStates, rawStates);  % e.g., 'x1st_jhana' -> '1st_jhana'

    subs = fieldnames(d_band);

    subj  = strings(0,1);
    state = strings(0,1);
    ntr   = zeros(0,1);

    for i = 1:numel(subs)
        sfld = subs{i};
        statesFlds = fieldnames(d_band.(sfld));
        for j = 1:numel(statesFlds)
            stFld = statesFlds{j};
            entry = d_band.(sfld).(stFld);

            % only averaged entries have ntrials
            if ~isfield(entry, 'ntrials'), continue; end

            % map struct field -> human-readable label
            if isKey(S2H, stFld)
                stHuman = string(S2H(stFld));
            else
                stHuman = string(pretty_state(stFld));   % fallback
            end

            subj(end+1,1)  = string(sfld);
            state(end+1,1) = stHuman;
            ntr(end+1,1)   = double(entry.ntrials);
        end
    end

    qcTable = table(subj, state, ntr, 'VariableNames', {'subject','state','ntrials'});

    % Sort: subject then preferred state order
    preferred = {'counting','memory','1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
                 '5th_jhana','6th_jhana','7th_jhana','8th_jhana', ...
                 'form_jhana','formless_jhana','all_jhana','ctrl_combined'};
    [~, ordState] = ismember(qcTable.state, preferred);
    ordState(ordState==0) = numel(preferred) + 1; % unknowns last
    [~, idx] = sortrows([double(categorical(qcTable.subject)), ordState]);
    qcTable   = qcTable(idx, :);

    % Flag anything below cutoff
    qcTable.flag_below_cutoff = qcTable.ntrials < cutoff;
end
