%==========================================================================
% CBPT_example.m
%
% Example script: spectral entropy. The same settings were used for all
% complexity and power measures; only the input file and feature differ.
%
% Cluster-based permutation tests (CBPT):
%   ctrl_combined  vs  1st–8th jhāna
%
% Using Spectral Entropy outputs from:
%   spectral_entropy.m
%
% Expected input:
%   d_specen.(subj).(state).powspctrm  ->  1 x chan x feat(=1)
%   d_specen.(subj).(state).featlabels ->  {'SpecEn_1to35Hz'}
%
% Correction across the 8 contrasts:
%   - raw_p  : min cluster p per jhāna (uncorrected)
%   - holm_p : Holm-Bonferroni across 8 jhānas
%   - bh_q   : BH-FDR across 8 jhānas
%
% Sign convention:
%   If jhanaMinusCtrl = true, A=jhāna and B=ctrl => t > 0 means jhāna > ctrl.
%
% Dependencies:
%   - FieldTrip (ft_defaults, ft_freqstatistics, ft_prepare_neighbours, ft_prepare_layout)
%   - EEGLAB (for reading chanlocs from a .set file)
%==========================================================================

clear; clc;

% ---- FieldTrip ----
if exist('ft_defaults','file'), ft_defaults; end

% ---- EEGLAB (adjust path if needed) ----
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab;

%% ---------- 0) Paths and settings ----------
dataDir = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_2s/35hz_laplacian';

% Output of spectral_entropy.m
matPath = fullfile(dataDir, 'specen_singlechan_broadband_spectralentropy.mat');

% Conditions
ctrlName   = 'ctrl_combined';
jh_names_h = {'1st_jhana','2nd_jhana','3rd_jhana','4th_jhana', ...
              '5th_jhana','6th_jhana','7th_jhana','8th_jhana'};
jh_names_f = cellfun(@matlab.lang.makeValidName, jh_names_h, 'UniformOutput', false);

% Feature(s) to use
% Leave empty {} to use all available features in the file
requested_feats = {'SpecEn_1to35Hz'};

% CBPT parameters
nPerm         = 'all';    % 'all' or numeric
clusterAlpha  = 0.05;     % cfg.clusteralpha
alphaCluster  = 0.025;    % cfg.alpha
minNbChan     = 2;        % cfg.minnbchan

% Sign convention: if true, A=jhāna and B=ctrl => t>0 means jhāna>ctrl
jhanaMinusCtrl = true;

% Parallelization
USE_PARFOR = true;
if USE_PARFOR && isempty(ver('parallel'))
    warning('Parallel Computing Toolbox not found; running serial.');
    USE_PARFOR = false;
end
if USE_PARFOR
    p = gcp('nocreate');
    if isempty(p)
        try
            parpool;
        catch
            warning('Could not start parpool; running serial.');
            USE_PARFOR = false;
        end
    end
end

% Output
outDir = fullfile(dataDir, 'cbpt_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
tag     = datestr(now,'yyyymmdd_HHMM');
safe    = @(s) matlab.lang.makeValidName(s);
outBase = sprintf('cbpt_SpecEn_%s_vs_all_jhanas_%s', safe(ctrlName), tag);
outMat  = fullfile(outDir, [outBase '.mat']);

%% ---------- 1) Load d_specen; get subjects & features ----------
S        = load(matPath, 'd_specen');
d_specen = S.d_specen;
subsAll  = fieldnames(d_specen);

assert(~isempty(subsAll), 'No subjects found in d_specen.');

% Find a reference subject with ctrl and any jhāna
refSub = "";
for i = 1:numel(subsAll)
    sf = subsAll{i};
    if isfield(d_specen.(sf), ctrlName)
        hasAnyJh = any(cellfun(@(f) isfield(d_specen.(sf), f), jh_names_f));
        if hasAnyJh
            refSub = string(sf);
            break;
        end
    end
end
assert(refSub~="", 'No subject found with %s and any jhāna.', ctrlName);

ref = d_specen.(refSub).(ctrlName);
labels = ref.label(:);
nChan  = numel(labels);

% Feature labels
assert(isfield(ref, 'featlabels') && ~isempty(ref.featlabels), ...
    'featlabels missing in d_specen reference entry.');
featlabels_all = ref.featlabels(:);

% Select requested features
if isempty(requested_feats)
    featIdxUse = 1:numel(featlabels_all);
else
    featIdxUse = zeros(1, numel(requested_feats));
    for r = 1:numel(requested_feats)
        idx = find(strcmp(featlabels_all, requested_feats{r}), 1, 'first');
        assert(~isempty(idx), ...
            'Requested feature "%s" not found in featlabels: %s', ...
            requested_feats{r}, strjoin(featlabels_all.', ', '));
        featIdxUse(r) = idx;
    end
end

featlabels = featlabels_all(featIdxUse);
nFeat      = numel(featlabels);
Kc         = numel(jh_names_h);

fprintf('Reference subject: %s | nChan=%d\n', refSub, nChan);
fprintf('Available featlabels: %s\n', strjoin(featlabels_all.', ', '));
fprintf('Using SpecEn features: %s\n', strjoin(featlabels.', ', '));

%% ---------- 2) Build neighbours & layout from first .set ----------
setList = dir(fullfile(dataDir, '*.set'));
assert(~isempty(setList), 'No .set files found in %s', dataDir);

EEG = pop_loadset('filename', setList(1).name, 'filepath', setList(1).folder);
assert(isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs), 'chanlocs missing in the .set');

[elecAll, haveXYZ] = eeglab_chanlocs_to_ft_elec(EEG.chanlocs);
assert(haveXYZ, 'Some channels lack X/Y/Z; fix chanlocs and re-run.');

labels_cell = cellstr(labels);
[idxOK, idxLabInAll] = ismember(labels_cell, elecAll.label);
assert(all(idxOK), 'Electrode coordinates missing for: %s', strjoin(labels_cell(~idxOK), ', '));
ord = idxLabInAll;

elec = [];
elec.label   = elecAll.label(ord);
elec.elecpos = elecAll.elecpos(ord,:);
elec.chanpos = elecAll.chanpos(ord,:);
elec.unit    = elecAll.unit;

cfgN = [];
cfgN.method = 'triangulation';
cfgN.elec   = elec;
neighbours  = ft_prepare_neighbours(cfgN);

cfgL = [];
cfgL.elec   = elec;
layout      = ft_prepare_layout(cfgL);

deg = cellfun(@numel, {neighbours.neighblabel});
fprintf('Neighbours (triangulation): mean degree=%.1f (min=%d, max=%d) over %d channels\n', ...
        mean(deg), min(deg), max(deg), numel(deg));

%% ---------- 3) Run the 8 contrasts per feature ----------
features = repmat(struct('featlabel',[], 'contrast_names',[], 'stats',[], ...
                         'raw_p',[], 'holm_p',[], 'bh_q',[], ...
                         'n_subj_used',[], 'subjects_used',[]), 1, nFeat);

for fIdx = 1:nFeat
    featName   = featlabels{fIdx};
    featIdxPow = featIdxUse(fIdx);

    fprintf('\n===== FEATURE: %s =====\n', featName);

    stats_cell  = cell(1,Kc);
    raw_p       = nan(Kc,1);
    n_subj_used = nan(Kc,1);
    subs_lists  = cell(1,Kc);

    if USE_PARFOR
        parfor k = 1:Kc
            [stat_k, pmin_k, nSubj_k, subs_k] = run_one_contrast_specen( ...
                k, d_specen, subsAll, ctrlName, jh_names_h, jh_names_f, ...
                featName, featIdxPow, nChan, elec, neighbours, layout, ...
                jhanaMinusCtrl, clusterAlpha, alphaCluster, nPerm, minNbChan);

            stats_cell{k}  = stat_k;
            raw_p(k)       = pmin_k;
            n_subj_used(k) = nSubj_k;
            subs_lists{k}  = subs_k;
        end
    else
        for k = 1:Kc
            [stat_k, pmin_k, nSubj_k, subs_k] = run_one_contrast_specen( ...
                k, d_specen, subsAll, ctrlName, jh_names_h, jh_names_f, ...
                featName, featIdxPow, nChan, elec, neighbours, layout, ...
                jhanaMinusCtrl, clusterAlpha, alphaCluster, nPerm, minNbChan);

            stats_cell{k}  = stat_k;
            raw_p(k)       = pmin_k;
            n_subj_used(k) = nSubj_k;
            subs_lists{k}  = subs_k;
        end
    end

    % Holm & BH within this feature (1×8 family)
    [p_holm, ~] = holm_stepdown(raw_p);
    q_bh        = bh_fdr(raw_p);

    features(fIdx).featlabel      = featName;
    features(fIdx).contrast_names = jh_names_h;
    features(fIdx).stats          = stats_cell;
    features(fIdx).raw_p          = raw_p;
    features(fIdx).holm_p         = p_holm;
    features(fIdx).bh_q           = q_bh;
    features(fIdx).n_subj_used    = n_subj_used;
    features(fIdx).subjects_used  = subs_lists;
end

%% ---------- 4) Summary table ----------
summary = table();
for fIdx = 1:nFeat
    for j = 1:Kc
        st = features(fIdx).stats{j};
        if isempty(st)
            ppos = NaN; pneg = NaN; pmin = NaN;
        else
            ppos = NaN; pneg = NaN;
            if isfield(st,'posclusters') && ~isempty(st.posclusters)
                ppos = min([st.posclusters.prob]);
            end
            if isfield(st,'negclusters') && ~isempty(st.negclusters)
                pneg = min([st.negclusters.prob]);
            end
            pmin = min([ppos pneg], [], 'omitnan');
        end

        summary = [summary; table( ...
            string(featlabels{fIdx}), string(jh_names_h{j}), ppos, pneg, pmin, ...
            'VariableNames', {'feature','jhana','min_pos_p','min_neg_p','min_two_sided_p'})]; %#ok<AGROW>
    end
end

%% ---------- 5) Subject bookkeeping ----------
allUsed = strings(0,1);
for fIdx = 1:numel(features)
    L = features(fIdx).subjects_used;
    for k = 1:numel(L)
        if ~isempty(L{k})
            allUsed = [allUsed; L{k}(:)]; %#ok<AGROW>
        end
    end
end
subjects_unique = unique(allUsed, 'stable');

%% ---------- 6) Save ----------
results = struct();
results.features    = features;
results.labels      = elec.label;
results.featlabels  = featlabels;
results.neighbours  = neighbours;
results.layout      = layout;

results.subjects_unique      = subjects_unique;
results.subjects_per_feature = arrayfun(@(f) {f.subjects_used}, features);
results.summary              = summary;

results.params = struct( ...
    'family',            'ctrl_vs_each_jhana_per_feature_specen', ...
    'ctrl',              ctrlName, ...
    'jhanas',            {jh_names_h}, ...
    'features',          {featlabels}, ...
    'nPerm',             nPerm, ...
    'clusterAlpha',      clusterAlpha, ...
    'alphaCluster',      alphaCluster, ...
    'minNbChan',         minNbChan, ...
    'neighbour_method',  'triangulation', ...
    'jhanaMinusCtrl',    jhanaMinusCtrl, ...
    'use_parfor',        USE_PARFOR, ...
    'matPath',           matPath);

save(outMat, '-struct', 'results', '-v7.3');
fprintf('\nSaved SpecEn CBPT results to: %s\n', outMat);

%% ================== Helpers ==================

function [stat, pmin, nSubj, subs_here] = run_one_contrast_specen( ...
    k, d_specen, subsAll, ctrlName, jh_names_h, jh_names_f, ...
    featName, featIdxPow, nChan, elec, neighbours, layout, ...
    jhanaMinusCtrl, clusterAlpha, alphaCluster, nPerm, minNbChan)

    jh_h = jh_names_h{k};
    jh_f = jh_names_f{k};

    % Subjects that have BOTH ctrl and this jhāna
    subs_here = strings(0,1);
    for i = 1:numel(subsAll)
        sf = subsAll{i};
        if isfield(d_specen.(sf), ctrlName) && isfield(d_specen.(sf), jh_f)
            subs_here(end+1,1) = string(sf);
        end
    end
    nSubj = numel(subs_here);

    if nSubj < 4
        fprintf('  %s vs %s: skipped (only %d subjects)\n', ctrlName, jh_h, nSubj);
        stat = [];
        pmin = NaN;
        return;
    end

    % Build subj × chan matrices for this feature
    A = nan(nSubj, nChan);  % jhāna (if jhanaMinusCtrl)
    B = nan(nSubj, nChan);  % ctrl

    for i = 1:nSubj
        sf = subs_here(i);

        Pctrl = squeeze(d_specen.(sf).(ctrlName).powspctrm); % chan x feat
        Pjh   = squeeze(d_specen.(sf).(jh_f).powspctrm);     % chan x feat

        if isvector(Pctrl), Pctrl = Pctrl(:); end
        if isvector(Pjh),   Pjh   = Pjh(:);   end

        if jhanaMinusCtrl
            A(i,:) = double(Pjh(:,featIdxPow)).';
            B(i,:) = double(Pctrl(:,featIdxPow)).';
        else
            A(i,:) = double(Pctrl(:,featIdxPow)).';
            B(i,:) = double(Pjh(:,featIdxPow)).';
        end
    end

    % Wrap as FieldTrip-like "freq" data with singleton freq dim
    freqA = [];
    freqA.label     = elec.label;
    freqA.dimord    = 'subj_chan_freq';
    freqA.freq      = 1;
    freqA.powspctrm = reshape(double(A), [nSubj, nChan, 1]);
    freqA.units     = 'unitless (spectral entropy)';
    freqA.featlabel = featName;

    freqB = freqA;
    freqB.powspctrm = reshape(double(B), [nSubj, nChan, 1]);

    % Paired design
    design = zeros(2, 2*nSubj);
    design(1,:) = [1:nSubj 1:nSubj];
    design(2,:) = [ones(1,nSubj) 2*ones(1,nSubj)];

    % CBPT config
    cfg = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.parameter        = 'powspctrm';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = clusterAlpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = alphaCluster;
    cfg.numrandomization = nPerm;
    cfg.uvar             = 1;
    cfg.ivar             = 2;
    cfg.design           = design;
    cfg.channel          = 'all';
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = minNbChan;

    if jhanaMinusCtrl
        fprintf('  CBPT: %s minus %s — %s (n=%d)\n', jh_h, ctrlName, featName, nSubj);
    else
        fprintf('  CBPT: %s minus %s — %s (n=%d)\n', ctrlName, jh_h, featName, nSubj);
    end

    stat = ft_freqstatistics(cfg, freqA, freqB);   % sign = A − B
    stat.featlabel = featName;
    stat.layout    = layout;

    % One p per contrast: min cluster p across pos/neg
    pmin = 1;
    if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
        pmin = min(pmin, min([stat.posclusters.prob]));
    end
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
        pmin = min(pmin, min([stat.negclusters.prob]));
    end
end

function [elec, haveXYZ] = eeglab_chanlocs_to_ft_elec(chanlocs)
    labels = {chanlocs.labels}';
    hasX   = arrayfun(@(c) isfield(c,'X') && ~isempty(c.X), chanlocs);
    hasY   = arrayfun(@(c) isfield(c,'Y') && ~isempty(c.Y), chanlocs);
    hasZ   = arrayfun(@(c) isfield(c,'Z') && ~isempty(c.Z), chanlocs);
    keep   = hasX & hasY & hasZ;

    haveXYZ = all(keep);
    labels  = labels(keep);
    pos     = [[chanlocs(keep).X]' [chanlocs(keep).Y]' [chanlocs(keep).Z]'];

    elec = [];
    elec.label   = labels;
    elec.elecpos = double(pos);
    elec.chanpos = double(pos);
    elec.unit    = 'mm';
end

function [p_adj, order] = holm_stepdown(p)
    p = p(:);
    p(isnan(p)) = 1;
    [p_sorted, order] = sort(p, 'ascend');
    m = numel(p_sorted);
    holm_vals = (m - (1:m)' + 1) .* p_sorted;
    p_adj_sorted = cummax_local(holm_vals);
    p_adj_sorted(p_adj_sorted>1) = 1;
    p_adj = zeros(size(p));
    p_adj(order) = p_adj_sorted;
end

function q = bh_fdr(p)
    p = p(:);
    p(isnan(p)) = 1;
    [p_sorted, order] = sort(p, 'ascend');
    m = numel(p_sorted);
    bh_vals = (m ./ (1:m)') .* p_sorted;
    q_sorted = cummin_local(bh_vals(end:-1:1));
    q_sorted = q_sorted(end:-1:1);
    q_sorted(q_sorted>1) = 1;
    q = zeros(size(p));
    q(order) = q_sorted;
end

function y = cummax_local(x)
    y = x;
    for i = 2:numel(x)
        y(i) = max(y(i-1), x(i));
    end
end

function y = cummin_local(x)
    y = x;
    for i = 2:numel(x)
        y(i) = min(y(i-1), x(i));
    end
end