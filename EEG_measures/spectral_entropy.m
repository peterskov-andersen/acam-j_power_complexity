%==========================================================================
%
% Spectral Entropy per epoch×channel using MATLAB spectralEntropy
%
% Intended data:
%   - epoched EEG .set files
%   - 2.048 s epochs (2048 samples)
%   - Fs = 1000 Hz
%   - 1 Hz high-pass (cleaning), 35 Hz low-pass, surface Laplacian
%
% Feature (feat = 1):
% 1) SpecEn_1to35Hz : spectral entropy of the full epoch, computed with
%                     MATLAB spectralEntropy() on the power spectrum
%
% Outputs (FieldTrip-like):
% d_specen_trials : trial-level (rpt × chan × feat=1)
% d_specen        : averaged   (1 × chan × feat=1) + std + ntrials
%
% Notes:
% - Call used (MATLAB R2023b):
%       spectralEntropy(x, Fs, 'Window', fullEpochWindow, ...
%                           'OverlapLength', 0, ...
%                           'FFTLength', N, ...
%                           'SpectrumType', 'power', ...
%                           'Range', [1 35])
%
% - Because Window = full epoch and OverlapLength = 0, each epoch produces
%   exactly one spectral entropy value per channel.
%
% Dependencies: EEGLAB, FieldTrip-lite (optional), MATLAB spectralEntropy
%==========================================================================
clc; clear; close all;

%% ----- Paths -----
% EEGLAB
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab

% FieldTrip-lite (optional)
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0/plugins/Fieldtrip-lite250523');
if exist('ft_defaults','file'), ft_defaults; end

% Check spectralEntropy exists
assert(exist('spectralEntropy','file')==2, ...
    'spectralEntropy() not found on path/license.');

%% 1) User parameters
dataDir    = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_2s/35hz_laplacian';
saveSpecEn = 'specen_singlechan_broadband_spectralentropy.mat';

qcCutoff   = 15;     % flag subject×state cells below this epoch count (QC table only)
USE_PARFOR = true;
VERBOSE    = true;

% Spectral entropy settings
Fs           = 1000;
epochSec     = 2.048;
expectedN    = Fs * epochSec;      % 2048 samples
rangeHz      = [1 35];
spectrumType = 'power';            % power spectrum

% Selection controls
SELECT_BY_INDEX = true;
INDEX_RANGE = [1 1000];   % inclusive indices into sorted file list
resume = true;

% Checkpoint settings
checkpointDir = fullfile(dataDir, 'specen_checkpoints');
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end

% Composite categories (trial-weighted via pooling)
composites = {
  'form_jhana',     {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana'};
  'formless_jhana', {'5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
  'all_jhana',      {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
                     '5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
  'ctrl_combined',  {'counting','memory'};
};

%% 2) Prepare per-run map (subject -> state -> {runs}) and file selection
runsMapSE = struct();
fileList = dir(fullfile(dataDir,'*.set'));
if isempty(fileList)
  error('No .set files found under: %s', dataDir);
end

% Stable sort by filename so indices are reproducible
[~,ord] = sort({fileList.name});
fileList = fileList(ord);

if SELECT_BY_INDEX
  i1 = max(1, INDEX_RANGE(1));
  i2 = min(numel(fileList), INDEX_RANGE(2));
  selIdx = i1:i2;
else
  selIdx = 1:numel(fileList);
end

fprintf('Processing files [%d..%d] of %d (sorted by name)\n', ...
    selIdx(1), selIdx(end), numel(fileList));

%% 3) Loop files -> compute SpecEn metrics per epoch×channel
for iF = selIdx
  fname = fileList(iF).name;
  fpath = fileList(iF).folder;

  % ---- 3.1 Parse subject, task and run from the filename
  tok = regexp(fname, ...
    'sub-(\d+)_ses-[^_]+_task-([^_]+)_run-(\d+)_eeg.*\.set$', ...
    'tokens','once');

  if isempty(tok)
    warning('Skipping %s (name parse failed)', fname);
    continue;
  end

  rawSub  = ['sub-' tok{1}];
  rawTask = tok{2};
  rawRun  = tok{3};

  subjFld = matlab.lang.makeValidName(rawSub);

  % ---- 3.2 Infer state from task; skip plain "jhana"
  inferredState = task_to_state(rawTask);
  if inferredState == ""
    if VERBOSE
      fprintf('[SKIP] %s | task=%s (no usable state label)\n', rawSub, rawTask);
    end
    continue;
  end

  % ---- 3.2.1 Checkpoint path for this run (subject×state×run)
  ckName = sprintf('%s__%s__run-%s.mat', rawSub, char(inferredState), rawRun);
  ckPath = fullfile(checkpointDir, ckName);

  % Resume
  if resume && exist(ckPath,'file')
    if VERBOSE, fprintf('[RESUME] Using checkpoint: %s\n', ckName); end
    S = load(ckPath,'fr');
    fr = S.fr;

    if ~isfield(runsMapSE, subjFld), runsMapSE.(subjFld) = struct(); end
    fld = matlab.lang.makeValidName(inferredState);
    if ~isfield(runsMapSE.(subjFld), fld), runsMapSE.(subjFld).(fld) = {fr};
    else, runsMapSE.(subjFld).(fld){end+1} = fr; end
    continue;
  end

  % ---- 3.3 Load EEG & convert to FieldTrip struct
  EEG = pop_loadset('filename', fname, 'filepath', fpath);
  ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % epoched -> trial{}
  nE = numel(ft_data.trial);
  nC = numel(ft_data.label);

  if VERBOSE
    fprintf('[SE ] %s | %s run-%s : epochs=%d, chans=%d -> state=%s\n', ...
      rawSub, rawTask, rawRun, nE, nC, char(inferredState));
  end

  FEAT = 1;
  A = nan(nE, nC, FEAT, 'double'); % rpt × chan × feat

  usePar = USE_PARFOR && ~isempty(ver('parallel')) && nE >= 4;
  if usePar
    pool = gcp('nocreate');
    if isempty(pool)
      try
        parpool;
      catch
        usePar = false;
      end
    end
  end

  if usePar
    parfor e = 1:nE
      X = double(ft_data.trial{e}); % chan × time
      row = nan(nC, FEAT);

      for c = 1:nC
        ts = X(c,:);
        if any(~isfinite(ts)), continue; end
        row(c,1) = compute_specen_feature_spectralentropy( ...
            ts, Fs, rangeHz, spectrumType, expectedN);
      end

      A(e,:,:) = row;
    end
  else
    for e = 1:nE
      if mod(e,25)==0 && VERBOSE, fprintf(' epoch %d/%d\n', e, nE); end
      X = double(ft_data.trial{e}); % chan × time

      for c = 1:nC
        ts = X(c,:);
        if any(~isfinite(ts)), continue; end
        A(e,c,1) = compute_specen_feature_spectralentropy( ...
            ts, Fs, rangeHz, spectrumType, expectedN);
      end
    end
  end

  % ---- 3.4 Build FR struct for this run
  labels = ft_data.label(:);

  if ~isfield(runsMapSE, subjFld), runsMapSE.(subjFld) = struct(); end
  fld = matlab.lang.makeValidName(inferredState);

  fr = struct();
  fr.label = labels;
  fr.dimord = 'rpt_chan_feat';
  fr.freq = FEAT; % number of features (stored in the freq field)
  fr.featname = 'specen_singlechan';
  fr.featlabels = {'SpecEn_1to35Hz'};
  fr.units = 'unitless (spectral entropy)';
  fr.powspctrm = A;

  fr.settings = struct();
  fr.settings.method        = 'MATLAB spectralEntropy';
  fr.settings.Fs            = Fs;
  fr.settings.epochSec      = epochSec;
  fr.settings.expectedN     = expectedN;
  fr.settings.rangeHz       = rangeHz;
  fr.settings.spectrumType  = spectrumType;
  fr.settings.window        = 'hamming(full epoch)';
  fr.settings.overlapLength = 0;
  fr.settings.fftLength     = expectedN;

  fr.meta.subject  = rawSub;
  fr.meta.task     = rawTask;
  fr.meta.run      = rawRun;
  fr.meta.filename = fname;

  % ---- 3.5 Save checkpoint
  try
    save(ckPath,'fr','-v7.3');
    if VERBOSE, fprintf('[SAVE] Checkpointed: %s\n', ckName); end
  catch ME
    warning('Failed to save checkpoint %s: %s', ckName, ME.message);
  end

  % ---- 3.6 Stash in runsMap for pooling
  if ~isfield(runsMapSE.(subjFld), fld), runsMapSE.(subjFld).(fld) = {fr};
  else, runsMapSE.(subjFld).(fld){end+1} = fr; end
end

%% 3.9 Merge ALL checkpoints into runsMapSE (cumulative final)
have = containers.Map('KeyType','char','ValueType','logical');
subs0 = fieldnames(runsMapSE);
for ii = 1:numel(subs0)
  sfld = subs0{ii};
  states = fieldnames(runsMapSE.(sfld));
  for jj = 1:numel(states)
    st = states{jj};
    runs = runsMapSE.(sfld).(st);
    for kk = 1:numel(runs)
      fr0 = runs{kk};
      key = sprintf('%s|%s|%s', fr0.meta.subject, task_to_state(fr0.meta.task), fr0.meta.run);
      have(key) = true;
    end
  end
end

files = dir(fullfile(checkpointDir,'*.mat'));
for k = 1:numel(files)
  S = load(fullfile(files(k).folder, files(k).name), 'fr');
  fr = S.fr;
  st = task_to_state(fr.meta.task);
  key = sprintf('%s|%s|%s', fr.meta.subject, st, fr.meta.run);
  if isKey(have,key), continue; end

  sfld = matlab.lang.makeValidName(fr.meta.subject);
  fld  = matlab.lang.makeValidName(st);

  if ~isfield(runsMapSE, sfld), runsMapSE.(sfld) = struct(); end
  if ~isfield(runsMapSE.(sfld), fld), runsMapSE.(sfld).(fld) = {fr};
  else, runsMapSE.(sfld).(fld){end+1} = fr; end

  have(key) = true;
end

%% 4) Append runs and average per subject×state
d_specen_trials = struct();
d_specen        = struct();

subs = fieldnames(runsMapSE);
for i = 1:numel(subs)
  sfld = subs{i};
  d_specen_trials.(sfld) = struct();
  d_specen.(sfld)        = struct();

  conds = fieldnames(runsMapSE.(sfld));
  for j = 1:numel(conds)
    cfld = conds{j};
    runs = runsMapSE.(sfld).(cfld);

    [se_trials, runCounts, runIDs] = append_trials_features(runs);
    nTrials = size(se_trials.powspctrm,1);

    if VERBOSE
      fprintf('[POOL] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
        sfld, pretty_state(cfld), char(strjoin(runIDs,',')), mat2str(runCounts), nTrials);
    end

    meanFeat = squeeze(mean(se_trials.powspctrm, 1, 'omitnan')); % chan×feat
    stdFeat  = squeeze(std( se_trials.powspctrm, 0, 1, 'omitnan')); % chan×feat

    se_avg = se_trials;
    se_avg.powspctrm = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
    se_avg.powspctrm_std = reshape(stdFeat, [1, size(stdFeat,1), size(stdFeat,2)]);
    se_avg.avg_method = 'mean';
    se_avg.ntrials = nTrials;

    d_specen_trials.(sfld).(cfld) = se_trials;
    d_specen.(sfld).(cfld)        = se_avg;
  end
end

%% 5) Composite categories (trial-weighted via pooling)
subs = fieldnames(d_specen_trials);
for i = 1:numel(subs)
  sfld = subs{i};

  for c = 1:size(composites,1)
    compName = composites{c,1};
    list     = composites{c,2};

    poolList = {};
    for k = 1:numel(list)
      fld = matlab.lang.makeValidName(list{k});
      if isfield(d_specen_trials.(sfld), fld)
        poolList{end+1} = d_specen_trials.(sfld).(fld); %#ok<AGROW>
      end
    end

    if ~isempty(poolList)
      se_trials = append_trials_features(poolList);
      nTrials   = size(se_trials.powspctrm,1);

      meanFeat = squeeze(mean(se_trials.powspctrm, 1, 'omitnan'));
      stdFeat  = squeeze(std( se_trials.powspctrm, 0, 1, 'omitnan'));

      se_avg = se_trials;
      se_avg.powspctrm = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
      se_avg.powspctrm_std = reshape(stdFeat, [1, size(stdFeat,1), size(stdFeat,2)]);
      se_avg.avg_method = 'mean';
      se_avg.ntrials = nTrials;

      d_specen_trials.(sfld).(compName) = se_trials;
      d_specen.(sfld).(compName)        = se_avg;
    end
  end

  % ctrl_combined: equal-weight mean of the counting and memory averages
  % (replaces the pooled-epoch version computed above)
  if isfield(d_specen.(sfld),'counting') && isfield(d_specen.(sfld),'memory')
    A = d_specen.(sfld).counting;
    B = d_specen.(sfld).memory;

    mu = 0.5*(squeeze(A.powspctrm) + squeeze(B.powspctrm)); % chan×feat
    sd = sqrt(0.5*(squeeze(A.powspctrm_std).^2 + squeeze(B.powspctrm_std).^2));

    out = A;
    out.powspctrm = reshape(mu, [1, size(mu,1), size(mu,2)]);
    out.powspctrm_std = reshape(sd, [1, size(sd,1), size(sd,2)]);
    out.avg_method = 'mean (equal-weight across states)';
    out.ntrials = (A.ntrials + B.ntrials);
    out.components = {'counting','memory'};

    d_specen.(sfld).ctrl_combined = out;
  end
end

%% 6) Save outputs
save(fullfile(dataDir, saveSpecEn), 'd_specen', 'd_specen_trials', '-v7.3');
fprintf('Saved spectral entropy data to %s\n', fullfile(dataDir, saveSpecEn));

%% 7) QC: subject×state epoch counts (from d_specen)
qcTable = qc_ntrials_table(d_specen, qcCutoff);
disp(qcTable);

if any(qcTable.flag_below_cutoff)
  fprintf('\nStages below %d epochs:\n', qcCutoff);
  disp(qcTable(qcTable.flag_below_cutoff,:));
else
  fprintf('\nNo stages below %d epochs.\n', qcCutoff);
end

writetable(qcTable, fullfile(dataDir, sprintf('ntrials_qc_cutoff%d_specen.csv', qcCutoff)));

%% ============================= Helpers =============================

function val = compute_specen_feature_spectralentropy(ts, Fs, rangeHz, spectrumType, expectedN)
% Returns scalar spectral entropy for one time series.
%
% Uses full epoch as one analysis window, so output is one value per epoch.
% Call used (MATLAB R2023b):
%   spectralEntropy(ts, Fs, 'Window', hamming(N), 'OverlapLength', 0,
%                   'FFTLength', N, 'SpectrumType', 'power', 'Range', [1 35])

  ts = double(ts(:)); % column vector

  if any(~isfinite(ts))
    val = NaN;
    return
  end

  N = numel(ts);
  if N < 16
    val = NaN;
    return
  end

  if nargin >= 5 && ~isempty(expectedN) && N ~= expectedN
    % Unexpected epoch length: warn and use the actual length.
    warning('Epoch length is %d samples, expected %d. Using actual length.', N, expectedN);
  end

  win = hamming(N);

  % Restrict the spectrum to rangeHz
  try
    val = spectralEntropy(ts, Fs, ...
        'Window',        win, ...
        'OverlapLength', 0, ...
        'FFTLength',     N, ...
        'SpectrumType',  spectrumType, ...
        'Range',         rangeHz);
  catch
    % Fallback: if Range is not accepted in the local release, compute on
    % the already-filtered signal without explicit range restriction.
    val = spectralEntropy(ts, Fs, ...
        'Window',        win, ...
        'OverlapLength', 0, ...
        'FFTLength',     N, ...
        'SpectrumType',  spectrumType);
  end

  val = double(val);

  % With a full-epoch window the result is already a scalar.
  if ~isscalar(val)
    val = mean(val(:), 'omitnan');
  end
end

function out = canonical_state(raw)
  s = lower(strtrim(string(raw)));
  s = regexprep(s, '[-\s]+', '_');

  if s == "memory",   out = "memory";   return; end
  if s == "counting", out = "counting"; return; end

  tok = regexp(s, '^(?:jhana_)?([1-8])(st|nd|rd|th)?(?:_jhana)?$', 'tokens','once');
  if ~isempty(tok)
    n = str2double(tok{1});
    ord = ["st","nd","rd","th","th","th","th","th"];
    out = string(sprintf('%d%s_jhana', n, ord(n)));
    return
  end

  tokw = regexp(s, '^(?:jhana_)?(first|second|third|fourth|fifth|sixth|seventh|eighth)(?:_jhana)?$', 'tokens','once');
  if ~isempty(tokw)
    names = ["first","second","third","fourth","fifth","sixth","seventh","eighth"];
    idx = find(strcmp(tokw{1}, names));
    ord = ["st","nd","rd","th","th","th","th","th"];
    out = string(sprintf('%d%s_jhana', idx, ord(idx)));
    return
  end

  out = "";
end

function out = task_to_state(taskStr)
  t = lower(strtrim(taskStr));
  m = regexp(t,'^jhana([1-8])$','tokens','once');
  if ~isempty(m)
    out = canonical_state(m{1});
    return
  end
  if strcmp(t,'jhana')
    out = "";
    return
  end
  out = canonical_state(t);
end

function lab = pretty_state(fld)
  base = regexprep(fld, '^x', '');
  lab = strrep(base, '_', ' ');
end

function [out, runCounts, runIDs] = append_trials_features(runList)
% Concatenate trials across runs for multi-feature data (feat dim >=1).
  if isstruct(runList), runList = {runList}; end
  assert(~isempty(runList), 'append_trials_features: empty run list');

  base = runList{1};
  nRuns = numel(runList);

  runCounts = zeros(1,nRuns);
  runIDs    = strings(1,nRuns);

  for k = 1:nRuns
    r = runList{k};
    assert(isequal(r.label, base.label), 'Channel labels differ across runs');
    assert(isequal(r.featname, base.featname), 'Feature name mismatch across runs');

    runCounts(k) = size(r.powspctrm,1);
    if isfield(r,'meta') && isfield(r.meta,'run')
      runIDs(k) = string(r.meta.run);
    else
      runIDs(k) = string(k);
    end
  end

  P = cellfun(@(r) r.powspctrm, runList, 'UniformOutput', false);
  out = base;
  out.powspctrm = cat(1, P{:});
  out.dimord = 'rpt_chan_feat';

  assert(size(out.powspctrm,1) == sum(runCounts), ...
    'Pooled ntrials mismatch: got %d, expected %d', size(out.powspctrm,1), sum(runCounts));
end

function qcTable = qc_ntrials_table(d_struct, cutoff)
% Build long table: subject | state | ntrials | flag_below_cutoff
  if nargin < 2 || isempty(cutoff), cutoff = 25; end

  preferred = {'counting','memory','1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
               '5th_jhana','6th_jhana','7th_jhana','8th_jhana', ...
               'form_jhana','formless_jhana','all_jhana','ctrl_combined'};

  subs = fieldnames(d_struct);
  subj = strings(0,1);
  state = strings(0,1);
  ntr = zeros(0,1);

  for i = 1:numel(subs)
    sfld = subs{i};
    statesFlds = fieldnames(d_struct.(sfld));
    for j = 1:numel(statesFlds)
      stFld = statesFlds{j};
      entry = d_struct.(sfld).(stFld);
      if ~isfield(entry,'ntrials'), continue; end
      subj(end+1,1)  = string(sfld);
      state(end+1,1) = string(strrep(stFld,'_',' '));
      ntr(end+1,1)   = double(entry.ntrials);
    end
  end

  qcTable = table(subj, state, ntr, 'VariableNames', {'subject','state','ntrials'});

  [~, ordState] = ismember(qcTable.state, preferred);
  ordState(ordState==0) = numel(preferred) + 1;

  [~, idx] = sortrows([double(categorical(qcTable.subject)), ordState]);
  qcTable = qcTable(idx,:);

  qcTable.flag_below_cutoff = qcTable.ntrials < cutoff;
end
