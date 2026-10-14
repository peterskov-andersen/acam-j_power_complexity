% DFA exponent per epoch and channel (10 s epochs).
% The state is taken from the task label in the filename; files labelled
% task-jhana with no number are skipped.
%
% Outputs:
%   d_dfa_trials (trial-level, rpt×chan×1)
%   d_dfa        (averaged 1×chan×1)
%
% Dependencies: EEGLAB, FieldTrip (ft_defaults), dfa_custom_box_size.m (on path)

clc; clear; close all;

%% 0) Initialize FieldTrip
if exist('ft_defaults','file'), ft_defaults; end

%% 1) User parameters
dataDir     = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoch_10s/Epoched/10s_35hz_laplacian';
saveDFA     = 'dfaData.mat';       % main output (d_dfa, d_dfa_trials)
qcCutoff    = 6;                   % flag subject×state cells below this epoch count (QC table only)
USE_PARFOR  = true;                % parallelize across epochs if available
VERBOSE     = true;

% DFA function (DFluC with fixed box sizes, see dfa_custom_box_size.m)
DFA_FUNC    = @dfa_custom_box_size;

% Composite categories (trial-weighted via pooling)
composites = {
  'form_jhana',     {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana'};
  'formless_jhana', {'5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
  'all_jhana',      {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
                     '5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
  'ctrl_combined',  {'counting','memory'}
};

%% 2) Prepare per-run map (subject -> state -> {runs})
runsMapDFA = struct();
fileList   = dir(fullfile(dataDir,'*.set'));
if isempty(fileList)
    error('No .set files found under: %s', dataDir);
end

%% 3) Loop files -> compute DFA per epoch×channel, stash by subject×state
for iF = 1:numel(fileList)
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

    % ---- 3.3 Load EEG & convert to FieldTrip struct
    EEG     = pop_loadset('filename', fname, 'filepath', fpath);
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');  % epoched -> trial{}

    nE = numel(ft_data.trial);
    nC = numel(ft_data.label);

    if VERBOSE
        fprintf('[DFA ] %s | %s run-%s : epochs=%d, chans=%d -> state=%s\n', ...
            rawSub, rawTask, rawRun, nE, nC, char(inferredState));
    end

    % ---- 3.4 DFA per epoch×channel (alpha exponent)
    A = nan(nE, nC, 'double');  % epochs × channels

    usePar = USE_PARFOR && ~isempty(ver('parallel')) && nE >= 4;
    if usePar
        pool = gcp('nocreate'); if isempty(pool), try parpool; catch, usePar=false; end, end
    end

    if usePar
        parfor e = 1:nE
            X = double(ft_data.trial{e}); % chan × time
            row = nan(1,nC);
            for c = 1:nC
                ts = X(c,:);
                if any(isfinite(ts))
                    ts = ts - mean(ts,'omitnan');   % center
                    y  = cumsum(ts);                % integrate (DFA profile)
                    try
                        row(c) = DFA_FUNC(y);
                    catch
                        row(c) = nan;
                    end
                end
            end
            A(e,:) = row;
        end
    else
        for e = 1:nE
            if mod(e,25)==0 && VERBOSE, fprintf('  epoch %d/%d\n', e, nE); end
            X = double(ft_data.trial{e}); % chan × time
            for c = 1:nC
                ts = X(c,:);
                if any(isfinite(ts))
                    ts = ts - mean(ts,'omitnan');
                    y  = cumsum(ts);
                    try
                        A(e,c) = DFA_FUNC(y);
                    catch
                        A(e,c) = nan;
                    end
                end
            end
        end
    end

    % ---- 3.5 State label (one state per file)
    labels = repmat(inferredState, nE, 1); %#ok<NASGU>

    % ---- 3.6 Stash by subject×state
    if ~isfield(runsMapDFA,subjFld), runsMapDFA.(subjFld) = struct(); end
    stRaw = inferredState;                         % single state per file
    fld   = matlab.lang.makeValidName(stRaw);

    fr = struct();
    fr.label      = ft_data.label(:);
    fr.dimord     = 'rpt_chan_freq';
    fr.freq       = 1;                      % single "feature" slot
    fr.featname   = 'dfa_alpha';
    fr.units      = 'alpha';
    fr.powspctrm  = reshape(A, [nE, nC, 1]);      % rpt×chan×1
    fr.meta.subject  = rawSub;
    fr.meta.task     = rawTask;
    fr.meta.run      = rawRun;
    fr.meta.filename = fname;

    if ~isfield(runsMapDFA.(subjFld),fld), runsMapDFA.(subjFld).(fld) = {fr};
    else,  runsMapDFA.(subjFld).(fld){end+1} = fr; end
end

%% 4) Append runs and average per subject×state
d_dfa_trials = struct();   % pooled trial-level (rpt×chan×1)
d_dfa        = struct();   % averaged (1×chan×1), with std + ntrials

subs = fieldnames(runsMapDFA);
for i = 1:numel(subs)
    sfld  = subs{i};
    d_dfa_trials.(sfld) = struct();
    d_dfa.(sfld)        = struct();

    conds = fieldnames(runsMapDFA.(sfld));
    for j = 1:numel(conds)
        cfld = conds{j};                  % sanitized state field
        runs = runsMapDFA.(sfld).(cfld);

        [dfa_trials, runCounts, runIDs] = append_trials_scalar(runs);
        nTrials = size(dfa_trials.powspctrm,1);

        if VERBOSE
            fprintf('[POOL] %s | %s : runs=%s, per-run trials=%s -> pooled ntrials=%d\n', ...
                sfld, pretty_state(cfld), char(strjoin(runIDs,',')), mat2str(runCounts), nTrials);
        end

        % Mean & SD across epochs (alpha)
        meanAlpha = squeeze(mean(dfa_trials.powspctrm, 1, 'omitnan'));   % chan×1
        stdAlpha  = squeeze(std( dfa_trials.powspctrm, 0, 1, 'omitnan'));% chan×1

        dfa_avg                    = dfa_trials;
        dfa_avg.powspctrm          = reshape(meanAlpha, [1, numel(meanAlpha), 1]); % 1×chan×1
        dfa_avg.powspctrm_std      = reshape(stdAlpha,  [1, numel(stdAlpha),  1]);
        dfa_avg.avg_method         = 'mean';
        dfa_avg.ntrials            = nTrials;

        d_dfa_trials.(sfld).(cfld) = dfa_trials;
        d_dfa.(sfld).(cfld)        = dfa_avg;
    end
end

%% 5) Composite jhana categories (trial-weighted via pooling)
subs = fieldnames(d_dfa_trials);
for i = 1:numel(subs)
    sfld = subs{i};

    for c = 1:size(composites,1)
        compName = composites{c,1};
        list     = composites{c,2};

        % Pool trials across listed states if present
        poolList = {};
        for k = 1:numel(list)
            fld = matlab.lang.makeValidName(list{k});
            if isfield(d_dfa_trials.(sfld), fld)
                poolList{end+1} = d_dfa_trials.(sfld).(fld); %#ok<AGROW>
            end
        end

        if ~isempty(poolList)
            dfa_trials = append_trials_scalar(poolList); % pooled trial-level
            nTrials    = size(dfa_trials.powspctrm,1);

            meanAlpha = squeeze(mean(dfa_trials.powspctrm, 1, 'omitnan'));
            stdAlpha  = squeeze(std( dfa_trials.powspctrm, 0, 1, 'omitnan'));

            dfa_avg                    = dfa_trials;
            dfa_avg.powspctrm          = reshape(meanAlpha, [1, numel(meanAlpha), 1]);
            dfa_avg.powspctrm_std      = reshape(stdAlpha, [1, numel(stdAlpha),  1]);
            dfa_avg.avg_method         = 'mean';
            dfa_avg.ntrials            = nTrials;

            d_dfa_trials.(sfld).(compName) = dfa_trials;
            d_dfa.(sfld).(compName)        = dfa_avg;
        end
    end

    % ctrl_combined: equal-weight mean of the counting and memory averages
    % (replaces the pooled-epoch version computed above)
    if isfield(d_dfa.(sfld),'counting') && isfield(d_dfa.(sfld),'memory')
        A = d_dfa.(sfld).counting;
        B = d_dfa.(sfld).memory;

        mu = 0.5*(squeeze(A.powspctrm) + squeeze(B.powspctrm));
        sd = sqrt(0.5*(squeeze(A.powspctrm_std).^2 + squeeze(B.powspctrm_std).^2)); % approximate pooled SD

        out                 = A;   % copy axes/meta
        out.powspctrm       = reshape(mu, [1, numel(mu), 1]);
        out.powspctrm_std   = reshape(sd, [1, numel(sd),  1]);
        out.avg_method      = 'mean (equal-weight across states)';
        out.ntrials         = (A.ntrials + B.ntrials);  % total epochs (counting + memory)
        out.components      = {'counting','memory'};
        d_dfa.(sfld).ctrl_combined = out;
    end
end

%% 6) Save outputs
save(fullfile(dataDir, saveDFA), 'd_dfa', 'd_dfa_trials', '-v7.3');
fprintf('Saved DFA data to %s\n', fullfile(dataDir,saveDFA));

%% 7) QC: subject×state epoch counts (from d_dfa)
qcTable = qc_ntrials_table(d_dfa, qcCutoff);
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

function [out, runCounts, runIDs] = append_trials_scalar(runList)
% Concatenate trials across runs for single-feature data (DFA alpha).
    if nargin < 1, error('append_trials_scalar: missing runList'); end
    if isstruct(runList), runList = {runList}; end
    assert(~isempty(runList), 'append_trials_scalar: empty run list');

    base      = runList{1};
    nRuns     = numel(runList);
    runCounts = zeros(1,nRuns);
    runIDs    = strings(1,nRuns);

    for k = 1:nRuns
        r = runList{k};
        assert(isequal(r.label, base.label), 'Channel labels differ across runs');
        assert(isfield(r,'freq') && isequal(r.freq, 1), 'Expected r.freq == 1');
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
    out.dimord     = 'rpt_chan_freq';
    out.featname   = 'dfa_alpha';
    out.units      = 'alpha';
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
