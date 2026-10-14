%==========================================================================
%
% Higuchi Fractal Dimension (HFD) per epoch×channel, single-feature HFD
% using FileExchange package "Higuchi and Katz fractal dimension measures"
%
% Feature (feat = 1):
% 1) HFD_kmax200 : Higuchi_FD(ts, Kmax=200)
%
% Outputs (FieldTrip-like):
% d_hfd_trials : trial-level (rpt × chan × feat=1)
% d_hfd        : averaged   (1 × chan × feat=1) + std + ntrials
%
% Dependencies:
%   - EEGLAB
%   - FieldTrip-lite (optional)
%   - Higuchi_FD.m on MATLAB path (FileExchange 50290)
%==========================================================================

clc; clear; close all;

%% ----- Paths -----
% === Higuchi/Katz FD package (File Exchange 50290) ===
% Point this to the folder where Higuchi_FD.m lives
higuchiPkgDir = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/HiguchiKatzFD_50290';
addpath(genpath(higuchiPkgDir));
assert(exist('Higuchi_FD','file')==2, 'Higuchi_FD() not found. Check higuchiPkgDir.');

% EEGLAB
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab

% FieldTrip-lite (optional)
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0/plugins/Fieldtrip-lite250523');
if exist('ft_defaults','file'), ft_defaults; end

%% 1) User parameters
dataDir  = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_2s/35hz_laplacian';
saveHFD  = 'hfd_singlechan_kmax200_higuchi.mat';

qcCutoff   = 15;     % flag subject×state cells below this epoch count (QC table only)
USE_PARFOR = true;
VERBOSE    = true;

% HFD settings
Kmax = 200;
featLabel = sprintf('HFD_kmax%d', Kmax);

% Selection controls
SELECT_BY_INDEX = true;
INDEX_RANGE = [1 1000];   % inclusive indices into sorted file list
resume = true;

% Checkpoint settings
checkpointDir = fullfile(dataDir, 'hfd_checkpoints');
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
runsMapHFD = struct();
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

fprintf('Processing files [%d..%d] of %d (sorted by name)\n', selIdx(1), selIdx(end), numel(fileList));

%% 3) Loop files -> compute HFD metric per epoch×channel
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

    if ~isfield(runsMapHFD, subjFld), runsMapHFD.(subjFld) = struct(); end
    fld = matlab.lang.makeValidName(inferredState);
    if ~isfield(runsMapHFD.(subjFld), fld), runsMapHFD.(subjFld).(fld) = {fr};
    else, runsMapHFD.(subjFld).(fld){end+1} = fr; end
    continue;
  end

  % ---- 3.3 Load EEG & convert to FieldTrip struct
  EEG = pop_loadset('filename', fname, 'filepath', fpath);

  if isfield(EEG,'srate') && ~isempty(EEG.srate) && EEG.srate ~= 1000
    warning('EEG.srate = %.2f (not 1000) for %s', EEG.srate, fname);
  end

  ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % epoched -> trial{}
  nE = numel(ft_data.trial);
  nC = numel(ft_data.label);

  if VERBOSE
    fprintf('[HFD] %s | %s run-%s : epochs=%d, chans=%d -> state=%s | Kmax=%d\n', ...
      rawSub, rawTask, rawRun, nE, nC, char(inferredState), Kmax);
  end

  FEAT = 1;
  A = nan(nE, nC, FEAT, 'double'); % rpt × chan × feat

  usePar = USE_PARFOR && ~isempty(ver('parallel')) && nE >= 4;
  if usePar
    pool = gcp('nocreate');
    if isempty(pool)
      try parpool; catch, usePar=false; end
    end
  end

  if usePar
    parfor e = 1:nE
      X = double(ft_data.trial{e}); % chan × time
      row = nan(nC, 1);
      for c = 1:nC
        ts = X(c,:);
        if any(~isfinite(ts)), continue; end
        row(c) = compute_hfd_single(ts, Kmax);
      end
      A(e,:,1) = row;
    end
  else
    for e = 1:nE
      if mod(e,25)==0 && VERBOSE, fprintf(' epoch %d/%d\n', e, nE); end
      X = double(ft_data.trial{e}); % chan × time
      for c = 1:nC
        ts = X(c,:);
        if any(~isfinite(ts)), continue; end
        A(e,c,1) = compute_hfd_single(ts, Kmax);
      end
    end
  end

  % ---- 3.4 Build FR struct for this run
  labels = ft_data.label(:);

  if ~isfield(runsMapHFD, subjFld), runsMapHFD.(subjFld) = struct(); end
  fld = matlab.lang.makeValidName(inferredState);

  fr = struct();
  fr.label = labels;
  fr.dimord = 'rpt_chan_feat';
  fr.freq = FEAT; % number of features (stored in the freq field)
  fr.featname = 'hfd_singlechan';
  fr.featlabels = {featLabel};
  fr.units = 'unitless (fractal dimension)';
  fr.powspctrm = A;
  fr.meta.subject  = rawSub;
  fr.meta.task     = rawTask;
  fr.meta.run      = rawRun;
  fr.meta.filename = fname;
  fr.meta.kmax     = Kmax;

  % ---- 3.5 Save checkpoint
  try
    save(ckPath,'fr','-v7.3');
    if VERBOSE, fprintf('[SAVE] Checkpointed: %s\n', ckName); end
  catch ME
    warning('Failed to save checkpoint %s: %s', ckName, ME.message);
  end

  % ---- 3.6 Stash in runsMap for pooling
  if ~isfield(runsMapHFD.(subjFld), fld), runsMapHFD.(subjFld).(fld) = {fr};
  else, runsMapHFD.(subjFld).(fld){end+1} = fr; end
end

%% 3.9 Merge ALL checkpoints into runsMapHFD (cumulative final)
have = containers.Map('KeyType','char','ValueType','logical');
subs0 = fieldnames(runsMapHFD);
for ii = 1:numel(subs0)
  sfld = subs0{ii};
  states = fieldnames(runsMapHFD.(sfld));
  for jj = 1:numel(states)
    st = states{jj};
    runs = runsMapHFD.(sfld).(st);
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

  if ~isfield(runsMapHFD, sfld), runsMapHFD.(sfld) = struct(); end
  if ~isfield(runsMapHFD.(sfld), fld), runsMapHFD.(sfld).(fld) = {fr};
  else, runsMapHFD.(sfld).(fld){end+1} = fr; end

  have(key) = true;
end

%% 4) Append runs and average per subject×state
d_hfd_trials = struct();
d_hfd        = struct();

subs = fieldnames(runsMapHFD);
for i = 1:numel(subs)
  sfld = subs{i};
  d_hfd_trials.(sfld) = struct();
  d_hfd.(sfld)        = struct();

  conds = fieldnames(runsMapHFD.(sfld));
  for j = 1:numel(conds)
    cfld = conds{j};
    runs = runsMapHFD.(sfld).(cfld);

    [hfd_trials, runCounts, runIDs] = append_trials_features(runs);
    nTrials = size(hfd_trials.powspctrm,1);

    if VERBOSE
      fprintf('[POOL] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
        sfld, pretty_state(cfld), char(strjoin(runIDs,',')), mat2str(runCounts), nTrials);
    end

    meanFeat = squeeze(mean(hfd_trials.powspctrm, 1, 'omitnan')); % chan×feat(=1)
    stdFeat  = squeeze(std( hfd_trials.powspctrm, 0, 1, 'omitnan')); % chan×feat

    hfd_avg = hfd_trials;
    hfd_avg.powspctrm = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
    hfd_avg.powspctrm_std = reshape(stdFeat, [1, size(stdFeat,1), size(stdFeat,2)]);
    hfd_avg.avg_method = 'mean';
    hfd_avg.ntrials = nTrials;

    d_hfd_trials.(sfld).(cfld) = hfd_trials;
    d_hfd.(sfld).(cfld)        = hfd_avg;
  end
end

%% 5) Composite categories (trial-weighted via pooling)
subs = fieldnames(d_hfd_trials);
for i = 1:numel(subs)
  sfld = subs{i};

  for c = 1:size(composites,1)
    compName = composites{c,1};
    list     = composites{c,2};

    poolList = {};
    for k = 1:numel(list)
      fld = matlab.lang.makeValidName(list{k});
      if isfield(d_hfd_trials.(sfld), fld)
        poolList{end+1} = d_hfd_trials.(sfld).(fld); %#ok<AGROW>
      end
    end

    if ~isempty(poolList)
      hfd_trials = append_trials_features(poolList);
      nTrials   = size(hfd_trials.powspctrm,1);

      meanFeat = squeeze(mean(hfd_trials.powspctrm, 1, 'omitnan'));
      stdFeat  = squeeze(std( hfd_trials.powspctrm, 0, 1, 'omitnan'));

      hfd_avg = hfd_trials;
      hfd_avg.powspctrm = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
      hfd_avg.powspctrm_std = reshape(stdFeat, [1, size(stdFeat,1), size(stdFeat,2)]);
      hfd_avg.avg_method = 'mean';
      hfd_avg.ntrials = nTrials;

      d_hfd_trials.(sfld).(compName) = hfd_trials;
      d_hfd.(sfld).(compName)        = hfd_avg;
    end
  end

  % ctrl_combined: equal-weight mean of the counting and memory averages
  % (replaces the pooled-epoch version computed above)
  if isfield(d_hfd.(sfld),'counting') && isfield(d_hfd.(sfld),'memory')
    A = d_hfd.(sfld).counting;
    B = d_hfd.(sfld).memory;

    mu = 0.5*(squeeze(A.powspctrm) + squeeze(B.powspctrm)); % chan×feat
    sd = sqrt(0.5*(squeeze(A.powspctrm_std).^2 + squeeze(B.powspctrm_std).^2));

    out = A;
    out.powspctrm = reshape(mu, [1, size(mu,1), size(mu,2)]);
    out.powspctrm_std = reshape(sd, [1, size(sd,1), size(sd,2)]);
    out.avg_method = 'mean (equal-weight across states)';
    out.ntrials = (A.ntrials + B.ntrials);
    out.components = {'counting','memory'};

    d_hfd.(sfld).ctrl_combined = out;
  end
end

%% 6) Save outputs
save(fullfile(dataDir, saveHFD), 'd_hfd', 'd_hfd_trials', '-v7.3');
fprintf('Saved HFD data to %s\n', fullfile(dataDir, saveHFD));

%% 7) QC: subject×state epoch counts (from d_hfd)
qcTable = qc_ntrials_table(d_hfd, qcCutoff);
disp(qcTable);

if any(qcTable.flag_below_cutoff)
  fprintf('\nStages below %d epochs:\n', qcCutoff);
  disp(qcTable(qcTable.flag_below_cutoff,:));
else
  fprintf('\nNo stages below %d epochs.\n', qcCutoff);
end

writetable(qcTable, fullfile(dataDir, sprintf('ntrials_qc_cutoff%d_hfd.csv', qcCutoff)));

%% ============================= Helpers =============================

function val = compute_hfd_single(ts, Kmax)
% Compute Higuchi FD for a single channel time series.
  ts = double(ts(:))'; % row
  if any(~isfinite(ts))
    val = NaN;
    return;
  end
  try
    val = Higuchi_FD(ts, Kmax);
  catch
    % Some implementations prefer column input; try fallback
    val = Higuchi_FD(ts(:), Kmax);
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
    assert(isequal(r.featlabels, base.featlabels), 'Feature labels mismatch across runs');

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
