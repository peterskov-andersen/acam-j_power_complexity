%==========================================================================
%
% Permutation Entropy (EntropyHub) per epoch×channel, two-feature broadband PE
%
% Features (feat = 2):
% 1) PE_tau4_stride2 : tau=4 (samples @1000Hz), stride=2  -> implemented via
%                      subsampling by 2 and using tau'=2 in the subsampled series
% 2) PE_tau8_stride4 : tau=8 (samples @1000Hz), stride=4  -> subsampling by 4 and tau'=2
%
% Outputs (FieldTrip-like):
% d_pe_trials : trial-level (rpt × chan × feat=2)
% d_pe        : averaged   (1 × chan × feat=2) + std + ntrials
%
% Dependencies: EEGLAB, FieldTrip-lite (optional), EntropyHub (PermEn)
%==========================================================================
clc; clear; close all;

%% ----- Paths -----
% EntropyHub folder (contains PermEn.m)
entropyHubDir = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/EntropyHub_v2.0.0';

% EEGLAB
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab

% FieldTrip-lite (optional)
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0/plugins/Fieldtrip-lite250523');
if exist('ft_defaults','file'), ft_defaults; end

% Add EntropyHub to path
addpath(genpath(entropyHubDir));
assert(exist('PermEn','file')==2, 'EntropyHub PermEn() not found. Check entropyHubDir.');

%% 1) User parameters
dataDir  = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_2s/35hz_laplacian';
savePE   = 'pe_singlechan_m3_tau4s2_tau8s4_entropyhub.mat';

qcCutoff   = 15;     % flag subject×state cells below this epoch count (QC table only)
USE_PARFOR = true;
VERBOSE    = true;

% PE settings
m = 3;   % embedding dimension (n)

% Two PE settings; tau and stride in samples at 1000 Hz
cfgPE(1).tau    = 4;  cfgPE(1).stride = 2; cfgPE(1).name = 'PE_tau4_stride2';
cfgPE(2).tau    = 8;  cfgPE(2).stride = 4; cfgPE(2).name = 'PE_tau8_stride4';

% Selection controls
SELECT_BY_INDEX = true;
INDEX_RANGE = [1 1000];   % inclusive indices into sorted file list
resume = true;

% Checkpoint settings
checkpointDir = fullfile(dataDir, 'pe_checkpoints');
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
runsMapPE = struct();
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

%% 3) Loop files -> compute PE metrics per epoch×channel
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

    if ~isfield(runsMapPE, subjFld), runsMapPE.(subjFld) = struct(); end
    fld = matlab.lang.makeValidName(inferredState);
    if ~isfield(runsMapPE.(subjFld), fld), runsMapPE.(subjFld).(fld) = {fr};
    else, runsMapPE.(subjFld).(fld){end+1} = fr; end
    continue;
  end

  % ---- 3.3 Load EEG & convert to FieldTrip struct
  EEG = pop_loadset('filename', fname, 'filepath', fpath);
  ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % epoched -> trial{}
  nE = numel(ft_data.trial);
  nC = numel(ft_data.label);

  if VERBOSE
    fprintf('[PE ] %s | %s run-%s : epochs=%d, chans=%d -> state=%s\n', ...
      rawSub, rawTask, rawRun, nE, nC, char(inferredState));
  end

  FEAT = numel(cfgPE);
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
      row = nan(nC, FEAT);
      for c = 1:nC
        ts = X(c,:);
        if any(~isfinite(ts)), continue; end
        row(c,:) = compute_pe_features_entropyhub(ts, m, cfgPE);
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
        A(e,c,:) = compute_pe_features_entropyhub(ts, m, cfgPE);
      end
    end
  end

  % ---- 3.4 Build FR struct for this run
  labels = ft_data.label(:);

  if ~isfield(runsMapPE, subjFld), runsMapPE.(subjFld) = struct(); end
  fld = matlab.lang.makeValidName(inferredState);

  fr = struct();
  fr.label = labels;
  fr.dimord = 'rpt_chan_feat';
  fr.freq = FEAT; % number of features (stored in the freq field)
  fr.featname = 'pe_singlechan';
  fr.featlabels = {cfgPE.name};
  fr.units = 'unitless (normalized PE)';
  fr.powspctrm = A;
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
  if ~isfield(runsMapPE.(subjFld), fld), runsMapPE.(subjFld).(fld) = {fr};
  else, runsMapPE.(subjFld).(fld){end+1} = fr; end
end

%% 3.9 Merge ALL checkpoints into runsMapPE (cumulative final)
have = containers.Map('KeyType','char','ValueType','logical');
subs0 = fieldnames(runsMapPE);
for ii = 1:numel(subs0)
  sfld = subs0{ii};
  states = fieldnames(runsMapPE.(sfld));
  for jj = 1:numel(states)
    st = states{jj};
    runs = runsMapPE.(sfld).(st);
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

  if ~isfield(runsMapPE, sfld), runsMapPE.(sfld) = struct(); end
  if ~isfield(runsMapPE.(sfld), fld), runsMapPE.(sfld).(fld) = {fr};
  else, runsMapPE.(sfld).(fld){end+1} = fr; end

  have(key) = true;
end

%% 4) Append runs and average per subject×state
d_pe_trials = struct();
d_pe        = struct();

subs = fieldnames(runsMapPE);
for i = 1:numel(subs)
  sfld = subs{i};
  d_pe_trials.(sfld) = struct();
  d_pe.(sfld)        = struct();

  conds = fieldnames(runsMapPE.(sfld));
  for j = 1:numel(conds)
    cfld = conds{j};
    runs = runsMapPE.(sfld).(cfld);

    [pe_trials, runCounts, runIDs] = append_trials_features(runs);
    nTrials = size(pe_trials.powspctrm,1);

    if VERBOSE
      fprintf('[POOL] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
        sfld, pretty_state(cfld), char(strjoin(runIDs,',')), mat2str(runCounts), nTrials);
    end

    meanFeat = squeeze(mean(pe_trials.powspctrm, 1, 'omitnan')); % chan×feat
    stdFeat  = squeeze(std( pe_trials.powspctrm, 0, 1, 'omitnan')); % chan×feat

    pe_avg = pe_trials;
    pe_avg.powspctrm = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
    pe_avg.powspctrm_std = reshape(stdFeat, [1, size(stdFeat,1), size(stdFeat,2)]);
    pe_avg.avg_method = 'mean';
    pe_avg.ntrials = nTrials;

    d_pe_trials.(sfld).(cfld) = pe_trials;
    d_pe.(sfld).(cfld)        = pe_avg;
  end
end

%% 5) Composite categories (trial-weighted via pooling)
subs = fieldnames(d_pe_trials);
for i = 1:numel(subs)
  sfld = subs{i};

  for c = 1:size(composites,1)
    compName = composites{c,1};
    list     = composites{c,2};

    poolList = {};
    for k = 1:numel(list)
      fld = matlab.lang.makeValidName(list{k});
      if isfield(d_pe_trials.(sfld), fld)
        poolList{end+1} = d_pe_trials.(sfld).(fld); %#ok<AGROW>
      end
    end

    if ~isempty(poolList)
      pe_trials = append_trials_features(poolList);
      nTrials   = size(pe_trials.powspctrm,1);

      meanFeat = squeeze(mean(pe_trials.powspctrm, 1, 'omitnan'));
      stdFeat  = squeeze(std( pe_trials.powspctrm, 0, 1, 'omitnan'));

      pe_avg = pe_trials;
      pe_avg.powspctrm = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
      pe_avg.powspctrm_std = reshape(stdFeat, [1, size(stdFeat,1), size(stdFeat,2)]);
      pe_avg.avg_method = 'mean';
      pe_avg.ntrials = nTrials;

      d_pe_trials.(sfld).(compName) = pe_trials;
      d_pe.(sfld).(compName)        = pe_avg;
    end
  end

  % ctrl_combined: equal-weight mean of the counting and memory averages
  % (replaces the pooled-epoch version computed above)
  if isfield(d_pe.(sfld),'counting') && isfield(d_pe.(sfld),'memory')
    A = d_pe.(sfld).counting;
    B = d_pe.(sfld).memory;

    mu = 0.5*(squeeze(A.powspctrm) + squeeze(B.powspctrm)); % chan×feat
    sd = sqrt(0.5*(squeeze(A.powspctrm_std).^2 + squeeze(B.powspctrm_std).^2));

    out = A;
    out.powspctrm = reshape(mu, [1, size(mu,1), size(mu,2)]);
    out.powspctrm_std = reshape(sd, [1, size(sd,1), size(sd,2)]);
    out.avg_method = 'mean (equal-weight across states)';
    out.ntrials = (A.ntrials + B.ntrials);
    out.components = {'counting','memory'};

    d_pe.(sfld).ctrl_combined = out;
  end
end

%% 6) Save outputs
save(fullfile(dataDir, savePE), 'd_pe', 'd_pe_trials', '-v7.3');
fprintf('Saved PE data to %s\n', fullfile(dataDir, savePE));

%% 7) QC: subject×state epoch counts (from d_pe)
qcTable = qc_ntrials_table(d_pe, qcCutoff);
disp(qcTable);

if any(qcTable.flag_below_cutoff)
  fprintf('\nStages below %d epochs:\n', qcCutoff);
  disp(qcTable(qcTable.flag_below_cutoff,:));
else
  fprintf('\nNo stages below %d epochs.\n', qcCutoff);
end

writetable(qcTable, fullfile(dataDir, sprintf('ntrials_qc_cutoff%d_pe.csv', qcCutoff)));

%% ============================= Helpers =============================

function feats = compute_pe_features_entropyhub(ts, m, cfgPE)
% Returns 1×FEAT vector of permutation entropy values normalized to [0,1]
% by log(m!) using EntropyHub PermEn(...,'Norm',true).

  FEAT = numel(cfgPE);
  feats = nan(1, FEAT);

  ts = double(ts(:))'; % row

  for i = 1:FEAT
    tau    = cfgPE(i).tau;
    stride = cfgPE(i).stride;

    if mod(tau, stride) ~= 0
      error('tau (%d) must be divisible by stride (%d).', tau, stride);
    end

    ts_s  = ts(1:stride:end);
    tau_s = tau / stride;

    try
      [~, Pnorm, ~] = PermEn(ts_s, 'm', m, 'tau', tau_s, 'Norm', true, 'Logx', 2);
      feats(i) = Pnorm(end);
    catch
      % alternative call syntax (other EntropyHub versions)
      [~, Pnorm, ~] = PermEn(ts_s, m, 'tau', tau_s, 'Norm', true, 'Logx', 2);
      feats(i) = Pnorm(end);
    end

    % sanity check
    if feats(i) > 1.0001 || feats(i) < -0.0001
      error('Normalized PE out of range: %.6f (check PermEn version/call).', feats(i));
    end
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
