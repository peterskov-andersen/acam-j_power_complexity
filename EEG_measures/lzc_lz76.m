%==========================================================================
% Lempel-Ziv complexity (LZ76) per epoch and channel (2.048 s epochs).
%
% Features (feat = 6):
%   1) C0          = raw LZ76 count on median-thresholded binary signal
%   2) c0          = C0 * log2(T) / T       (Mediano entropy-rate estimate)
%   3) Cshuf_mean  = mean LZ76 over Nshuf time-shuffled versions
%                    of the same binary sequence
%   4) LZsN        = C0 / Cshuf_mean        (time-shuffle-normalized ratio)
%   5) Cphase_mean = mean LZ76 over Nphase FT phase-randomized
%                    surrogates of the continuous signal
%   6) LZphaseN    = C0 / Cphase_mean       (phase-surrogate-normalized ratio)
%
% Outputs (FieldTrip-like):
%   d_lz_trials : trial-level (rpt × chan × feat=6)
%   d_lz        : averaged   (1   × chan × feat=6)
%
% featname stored as: 'lz_singlechan'
%
% Dependencies: EEGLAB, FieldTrip (ft_defaults), EntRate (LZ76.c/mex)
% Mediano P, Rosas F, Luppi A. EntRate: entropy rate estimators
%==========================================================================

clc; clear; close all;

% ----- Paths -----
% EntRate path
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/EntRate-master');

% EEGLAB:
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab

% Fieldtrip:
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0/plugins/Fieldtrip-lite250523');

%% 0) Initialize FieldTrip
if exist('ft_defaults','file'), ft_defaults; end

%% 1) User parameters
dataDir     = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_2s/35hz_laplacian';
saveLZ      = 'lz76_singlechan2norm.mat';   % main output (d_lz, d_lz_trials)
qcCutoff    = 15;                           % flag subject×state cells below this epoch count (QC table only)
USE_PARFOR  = true;                        % parallelize across epochs if available
VERBOSE     = true;

% LZ settings
Nshuf       = 50;        % # time-shuffled surrogates for Cshuf_mean
Nphase      = 100;       % # phase-randomized surrogates for Cphase_mean
rngSeed     = 1;         % seeds the main MATLAB session only, not parfor workers
if ~isempty(rngSeed), rng(rngSeed); end

% T in c0 = C0*log2(T)/T:
%   []        -> epoch length in samples (used here)
%   scalar    -> fixed T for all epochs
T_override  = [];        % e.g. [] or 2048 or 5000

% ----- Selection controls -----
SELECT_BY_INDEX = true;              % set true to use numeric indices
INDEX_RANGE     = [1 1000];          % inclusive 1-based indices into sorted file list
resume          = true;              % if true, skip files that already have checkpoints

% ----- Checkpoint settings -----
checkpointDir = fullfile(dataDir, 'lz_checkpoints');
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end

% Composite categories (trial-weighted via pooling)
composites = {
  'form_jhana',     {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana'}; ...
  'formless_jhana', {'5th_jhana','6th_jhana','7th_jhana','8th_jhana'}; ...
  'all_jhana',      {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
                     '5th_jhana','6th_jhana','7th_jhana','8th_jhana'}; ...
  'ctrl_combined',  {'counting','memory'}
};

%% 2) Prepare per-run map (subject -> state -> {runs}) and file selection
runsMapLZ = struct();
fileList  = dir(fullfile(dataDir,'*.set'));
if isempty(fileList)
    error('No .set files found under: %s', dataDir);
end

% Stable sort by filename so indices are reproducible
[~,ord] = sort({fileList.name});
fileList = fileList(ord);

% Compute the indices to run
if SELECT_BY_INDEX
    i1 = max(1, INDEX_RANGE(1));
    i2 = min(numel(fileList), INDEX_RANGE(2));
    selIdx = i1:i2;
else
    selIdx = 1:numel(fileList);
end
fprintf('Processing files [%d..%d] of %d (sorted by name)\n', selIdx(1), selIdx(end), numel(fileList));

%% 3) Loop files -> compute LZ metrics per epoch×channel
useFixedT = ~isempty(T_override);   % for parfor compatibility

for iF = selIdx
    fname = fileList(iF).name;
    fpath = fileList(iF).folder;

    % ---- 3.1 Parse subject, task and run from the filename
    tok = regexp(fname, ...
        'sub-(\d+)_ses-[^_]+_task-([^_]+)_run-(\d+)_eeg.*\.set$', 'tokens','once');
    if isempty(tok)
        warning('Skipping %s (name parse failed)', fname);
        continue;
    end
    rawSub  = ['sub-' tok{1}];
    rawTask = tok{2};
    rawRun  = tok{3};
    subjFld = matlab.lang.makeValidName(rawSub);

    % ---- 3.2 Infer state from task; skip plain "jhana"
    inferredState = task_to_state(rawTask);  % "" means skip
    if inferredState == ""
        if VERBOSE
            fprintf('[SKIP] %s | task=%s (no usable state label)\n', rawSub, rawTask);
        end
        continue;
    end

    % ---- 3.2.1 Checkpoint path for this run (subject×state×run)
    ckName = sprintf('%s__%s__run-%s.mat', rawSub, char(inferredState), rawRun);
    ckPath = fullfile(checkpointDir, ckName);

    % If resuming, skip compute and load the checkpoint directly
    if resume && exist(ckPath, 'file')
        if VERBOSE, fprintf('[RESUME] Using checkpoint: %s\n', ckName); end
        S = load(ckPath, 'fr');
        fr = S.fr;
        if ~isfield(runsMapLZ,subjFld), runsMapLZ.(subjFld) = struct(); end
        fld = matlab.lang.makeValidName(inferredState);
        if ~isfield(runsMapLZ.(subjFld),fld), runsMapLZ.(subjFld).(fld) = {fr};
        else,  runsMapLZ.(subjFld).(fld){end+1} = fr; end
        continue;
    end

    % ---- 3.3 Load EEG & convert to FieldTrip struct
    EEG     = pop_loadset('filename', fname, 'filepath', fpath);
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');  % epoched -> trial{}

    nE = numel(ft_data.trial);
    nC = numel(ft_data.label);
    if VERBOSE
        fprintf('[LZ  ] %s | %s run-%s : epochs=%d, chans=%d -> state=%s\n', ...
            rawSub, rawTask, rawRun, nE, nC, char(inferredState));
    end

    % ---- 3.4 LZ per epoch×channel (features: C0, c0, Cshuf_mean, LZsN, Cphase_mean, LZphaseN)
    FEAT = 6;
    A = nan(nE, nC, FEAT, 'double');  % rpt × chan × feat

    usePar = USE_PARFOR && ~isempty(ver('parallel')) && nE >= 4;
    if usePar
        pool = gcp('nocreate'); if isempty(pool), try parpool; catch, usePar=false; end, end
    end

    if usePar
        parfor e = 1:nE
            X = double(ft_data.trial{e}); % chan × time
            row = nan(nC, FEAT);
            for c = 1:nC
                ts = X(c,:);
                if any(isfinite(ts))
                    % ----- binarize around median -----
                    thr = median(ts, 'omitnan');
                    v   = ts > thr;

                    % ----- raw LZ76 -----
                    C0 = LZ76(v);

                    % ----- Mediano-normalized entropy rate c0 -----
                    if useFixedT
                        T  = T_override;
                    else
                        T  = numel(v);
                    end
                    c0 = C0 * log2(T) / T;

                    % ----- time-shuffled surrogate mean -----
                    if Nshuf > 0
                        Cshuf_mean = meanLZ76_over_shuffles(v, Nshuf);
                    else
                        Cshuf_mean = NaN;
                    end

                    % time-shuffle-normalized ratio
                    if ~isfinite(Cshuf_mean) || Cshuf_mean <= eps
                        LZsN = NaN;
                    else
                        LZsN = C0 / Cshuf_mean;
                    end

                    % ----- phase-randomized surrogate mean (FT) -----
                    if Nphase > 0
                        Cphase_mean = meanLZ76_over_phase_surrogates(ts, Nphase);
                    else
                        Cphase_mean = NaN;
                    end

                    % phase-surrogate-normalized ratio
                    if ~isfinite(Cphase_mean) || Cphase_mean <= eps
                        LZphaseN = NaN;
                    else
                        LZphaseN = C0 / Cphase_mean;
                    end

                    row(c,1) = C0;
                    row(c,2) = c0;
                    row(c,3) = Cshuf_mean;
                    row(c,4) = LZsN;
                    row(c,5) = Cphase_mean;
                    row(c,6) = LZphaseN;
                end
            end
            A(e,:,:) = row;
        end
    else
        for e = 1:nE
            if mod(e,25)==0 && VERBOSE, fprintf('  epoch %d/%d\n', e, nE); end
            X = double(ft_data.trial{e}); % chan × time
            for c = 1:nC
                ts = X(c,:);
                if any(isfinite(ts))
                    thr = median(ts, 'omitnan');
                    v   = ts > thr;

                    % raw LZ76
                    C0 = LZ76(v);

                    % Mediano-normalized entropy rate c0
                    if useFixedT
                        T  = T_override;
                    else
                        T  = numel(v);
                    end
                    c0 = C0 * log2(T) / T;

                    % time-shuffled surrogate mean
                    if Nshuf > 0
                        Cshuf_mean = meanLZ76_over_shuffles(v, Nshuf);
                    else
                        Cshuf_mean = NaN;
                    end

                    if ~isfinite(Cshuf_mean) || Cshuf_mean <= eps
                        LZsN = NaN;
                    else
                        LZsN = C0 / Cshuf_mean;
                    end

                    % phase-randomized surrogate mean (FT)
                    if Nphase > 0
                        Cphase_mean = meanLZ76_over_phase_surrogates(ts, Nphase);
                    else
                        Cphase_mean = NaN;
                    end

                    if ~isfinite(Cphase_mean) || Cphase_mean <= eps
                        LZphaseN = NaN;
                    else
                        LZphaseN = C0 / Cphase_mean;
                    end

                    A(e,c,1) = C0;
                    A(e,c,2) = c0;
                    A(e,c,3) = Cshuf_mean;
                    A(e,c,4) = LZsN;
                    A(e,c,5) = Cphase_mean;
                    A(e,c,6) = LZphaseN;
                end
            end
        end
    end

    % ---- 3.5 Channel labels
    labels = ft_data.label(:);

    % ---- 3.6 Build FR struct for this run
    if ~isfield(runsMapLZ,subjFld), runsMapLZ.(subjFld) = struct(); end
    stRaw = inferredState;                         % single state per file
    fld   = matlab.lang.makeValidName(stRaw);

    fr = struct();
    fr.label        = labels;
    fr.dimord       = 'rpt_chan_feat';
    fr.freq         = FEAT;                        % number of features (stored in the freq field)
    fr.featname     = 'lz_singlechan';
    fr.featlabels   = {'C0','c0','Cshuf_mean','LZsN','Cphase_mean','LZphaseN'};
    fr.units        = 'unitless';                  % C0 and surrogate means are counts; c0 is bits/sample; LZsN, LZphaseN are ratios
    fr.powspctrm    = A;                           % rpt×chan×feat
    fr.meta.subject = rawSub;
    fr.meta.task    = rawTask;
    fr.meta.run     = rawRun;
    fr.meta.filename= fname;

    % ---- 3.7 Save checkpoint for this run immediately
    try
        save(ckPath, 'fr', '-v7.3');
        if VERBOSE, fprintf('[SAVE] Checkpointed: %s\n', ckName); end
    catch ME
        warning('Failed to save checkpoint %s: %s', ckName, ME.message);
    end

    % ---- 3.8 Stash in runsMap for pooling
    if ~isfield(runsMapLZ.(subjFld),fld), runsMapLZ.(subjFld).(fld) = {fr};
    else,  runsMapLZ.(subjFld).(fld){end+1} = fr; end
end

%% 3.9 Merge ALL checkpoints into runsMapLZ (cumulative final)
% Build a set of keys already in runsMapLZ to avoid duplicates
have = containers.Map('KeyType','char','ValueType','logical');
subs0 = fieldnames(runsMapLZ);
for ii = 1:numel(subs0)
    sfld = subs0{ii};
    states = fieldnames(runsMapLZ.(sfld));
    for jj = 1:numel(states)
        st = states{jj};
        runs = runsMapLZ.(sfld).(st);
        for kk = 1:numel(runs)
            fr0 = runs{kk};
            key = sprintf('%s|%s|%s', fr0.meta.subject, task_to_state(fr0.meta.task), fr0.meta.run);
            have(key) = true;
        end
    end
end

% Load every checkpoint and add missing ones
files = dir(fullfile(checkpointDir, '*.mat'));
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'fr');
    fr = S.fr;
    st = task_to_state(fr.meta.task);
    key = sprintf('%s|%s|%s', fr.meta.subject, st, fr.meta.run);
    if isKey(have, key), continue; end
    sfld = matlab.lang.makeValidName(fr.meta.subject);
    fld  = matlab.lang.makeValidName(st);
    if ~isfield(runsMapLZ, sfld), runsMapLZ.(sfld) = struct(); end
    if ~isfield(runsMapLZ.(sfld), fld), runsMapLZ.(sfld).(fld) = {fr};
    else, runsMapLZ.(sfld).(fld){end+1} = fr; end
    have(key) = true;
end

%% 4) Append runs and average per subject×state
d_lz_trials = struct();   % pooled trial-level (rpt×chan×feat)
d_lz        = struct();   % averaged (1×chan×feat), with std + ntrials

subs = fieldnames(runsMapLZ);
for i = 1:numel(subs)
    sfld  = subs{i};
    d_lz_trials.(sfld) = struct();
    d_lz.(sfld)        = struct();

    conds = fieldnames(runsMapLZ.(sfld));
    for j = 1:numel(conds)
        cfld = conds{j};                  % sanitized state field
        runs = runsMapLZ.(sfld).(cfld);

        [lz_trials, runCounts, runIDs] = append_trials_features(runs);
        nTrials = size(lz_trials.powspctrm,1);

        if VERBOSE
            fprintf('[POOL] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
                sfld, pretty_state(cfld), char(strjoin(runIDs,',')), mat2str(runCounts), nTrials);
        end

        % Mean & SD across epochs (preserve features dimension)
        meanFeat = squeeze(mean(lz_trials.powspctrm, 1, 'omitnan'));   % chan×feat
        stdFeat  = squeeze(std( lz_trials.powspctrm, 0, 1, 'omitnan'));% chan×feat

        lz_avg                    = lz_trials;     % copy axes/meta
        lz_avg.powspctrm          = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]); % 1×chan×feat
        lz_avg.powspctrm_std      = reshape(stdFeat,  [1, size(stdFeat,1),  size(stdFeat,2)]);
        lz_avg.avg_method         = 'mean';
        lz_avg.ntrials            = nTrials;

        d_lz_trials.(sfld).(cfld) = lz_trials;
        d_lz.(sfld).(cfld)        = lz_avg;
    end
end

%% 5) Composite jhana categories (trial-weighted via pooling)
subs = fieldnames(d_lz_trials);
for i = 1:numel(subs)
    sfld = subs{i};

    for c = 1:size(composites,1)
        compName = composites{c,1};
        list     = composites{c,2};

        % Pool trials across listed states if present
        poolList = {};
        for k = 1:numel(list)
            fld = matlab.lang.makeValidName(list{k});
            if isfield(d_lz_trials.(sfld), fld)
                poolList{end+1} = d_lz_trials.(sfld).(fld); %#ok<AGROW>
            end
        end

        if ~isempty(poolList)
            lz_trials = append_trials_features(poolList); % pooled trial-level
            nTrials   = size(lz_trials.powspctrm,1);

            meanFeat = squeeze(mean(lz_trials.powspctrm, 1, 'omitnan'));
            stdFeat  = squeeze(std( lz_trials.powspctrm, 0, 1, 'omitnan'));

            lz_avg                    = lz_trials;
            lz_avg.powspctrm          = reshape(meanFeat, [1, size(meanFeat,1), size(meanFeat,2)]);
            lz_avg.powspctrm_std      = reshape(stdFeat,  [1, size(stdFeat,1), size(stdFeat,2)]);
            lz_avg.avg_method         = 'mean';
            lz_avg.ntrials            = nTrials;

            d_lz_trials.(sfld).(compName) = lz_trials;
            d_lz.(sfld).(compName)        = lz_avg;
        end
    end

    % ctrl_combined: equal-weight mean of the counting and memory averages
    % (replaces the pooled-epoch version computed above)
    if isfield(d_lz.(sfld),'counting') && isfield(d_lz.(sfld),'memory')
        A = d_lz.(sfld).counting;    % 1×chan×feat
        B = d_lz.(sfld).memory;      % 1×chan×feat

        mu = 0.5*(squeeze(A.powspctrm) + squeeze(B.powspctrm));             % chan×feat
        sd = sqrt(0.5*(squeeze(A.powspctrm_std).^2 + squeeze(B.powspctrm_std).^2));

        out               = A;   % copy axes/meta
        out.powspctrm     = reshape(mu, [1, size(mu,1), size(mu,2)]);
        out.powspctrm_std = reshape(sd, [1, size(sd,1), size(sd,2)]);
        out.avg_method    = 'mean (equal-weight across states)';
        out.ntrials       = (A.ntrials + B.ntrials);  % total epochs (counting + memory)
        out.components    = {'counting','memory'};
        d_lz.(sfld).ctrl_combined = out;
    end
end

%% 6) Save outputs
save(fullfile(dataDir, saveLZ), 'd_lz', 'd_lz_trials', '-v7.3');
fprintf('Saved LZ data to %s\n', fullfile(dataDir,saveLZ));

%% 7) QC: subject×state epoch counts (from d_lz)
qcTable = qc_ntrials_table(d_lz, qcCutoff);
disp(qcTable)
if any(qcTable.flag_below_cutoff)
    fprintf('\nStages below %d epochs:\n', qcCutoff);
    disp(qcTable(qcTable.flag_below_cutoff,:));
else
    fprintf('\nNo stages below %d epochs.\n', qcCutoff);
end
writetable(qcTable, fullfile(dataDir, sprintf('ntrials_qc_cutoff%d.csv', qcCutoff)));

% ----------------- helpers -----------------
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
% Map filename task (e.g., 'jhana5','memory','counting') -> canonical state.
    t = lower(strtrim(taskStr));
    % Match jhana1..jhana8
    m = regexp(t,'^jhana([1-8])$','tokens','once');
    if ~isempty(m)
        out = canonical_state(m{1});   % numeric -> "1st_jhana", etc.
        return
    end
    % Skip plain 'jhana'
    if strcmp(t,'jhana')
        out = ""; return
    end
    % memory/counting or already-canonical text
    out = canonical_state(t);
end

function lab = pretty_state(fld)
    base = regexprep(fld, '^x', '');
    lab  = strrep(base, '_', ' ');
end

function [out, runCounts, runIDs] = append_trials_features(runList)
% Concatenate trials across runs for multi-feature data (feat dim >=1).
    if nargin < 1, error('append_trials_features: missing runList'); end
    if isstruct(runList), runList = {runList}; end
    assert(~isempty(runList), 'append_trials_features: empty run list');

    base      = runList{1};
    nRuns     = numel(runList);
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
    out            = base;
    out.powspctrm  = cat(1, P{:});                 % rpt concatenation
    out.dimord     = 'rpt_chan_feat';
    try
        trlIdx = arrayfun(@(k) k*ones(runCounts(k),1), 1:nRuns, 'UniformOutput', false);
        out.trialinfo = vertcat(trlIdx{:});
    catch
    end
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

    subj  = strings(0,1);
    state = strings(0,1);
    ntr   = zeros(0,1);

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

    % Sort: subject then preferred state order
    [~, ordState] = ismember(qcTable.state, preferred);
    ordState(ordState==0) = numel(preferred) + 1;
    [~, idx] = sortrows([double(categorical(qcTable.subject)), ordState]);
    qcTable   = qcTable(idx, :);

    % Flag anything below cutoff
    qcTable.flag_below_cutoff = qcTable.ntrials < cutoff;
end

% ----------------- LZ utilities -----------------
function mu = meanLZ76_over_shuffles(v, N)
% v: logical row vector; N: number of shuffles
% Returns the mean LZ76 (EntRate) over N independent random permutations of v.
    L = numel(v);
    vals = zeros(1,N);
    for n = 1:N
        idx = randperm(L);
        vals(n) = LZ76(v(idx));
    end
    mu = mean(vals);
end

function mu = meanLZ76_over_phase_surrogates(ts, N)
% ts: 1 × L continuous time series (double)
% N : number of FT phase-randomized surrogates
% Returns mean LZ76 over N surrogates:
%   - preserve power spectrum (FT surrogate)
%   - randomize phases
%   - binarize (median threshold)
%   - compute LZ76 on binary surrogate

    L = numel(ts);
    vals = zeros(1,N);

    % FFT of the original signal
    X     = fft(ts);
    amp   = abs(X);
    phase = angle(X);

    % Indices for positive frequencies (excluding DC and Nyquist)
    if rem(L,2) == 0
        % even length: 1=DC, L/2+1=Nyquist
        posIdx = 2:(L/2);
        nyqIdx = L/2 + 1;
    else
        % odd length: 1=DC, no exact Nyquist
        posIdx = 2:((L+1)/2);
        nyqIdx = [];
    end

    for n = 1:N
        newPhase = phase;

        % Random phases for positive frequencies
        randPh = 2*pi*rand(size(posIdx));
        newPhase(posIdx) = randPh;

        % Keep Nyquist purely real if present
        if ~isempty(nyqIdx)
            newPhase(nyqIdx) = 0;
        end

        % Enforce Hermitian symmetry for negative frequencies
        for k = posIdx
            k_conj = L + 2 - k;      % index of conjugate partner
            newPhase(k_conj) = -newPhase(k);
        end

        % Reconstruct complex spectrum with preserved amplitudes
        Xs = amp .* exp(1i * newPhase);

        % Back to time domain
        xs = real(ifft(Xs));

        % Binarize surrogate around its median
        thr_s = median(xs, 'omitnan');
        v_s   = xs > thr_s;

        vals(n) = LZ76(v_s);
    end

    mu = mean(vals);
end
