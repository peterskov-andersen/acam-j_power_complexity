% =========================================================================
% Compute channel-wise aperiodic exponent / intercept from saved IRASA
% fractal spectra.
%
% INPUT:
%   powerData_full_spectrum_irasa.mat
%       -> uses d_full_frac
%
% OUTPUT:
%   aperiodic_params_irasa.mat
%
% SAVED STRUCTURE:
%   d_aperiodic.(subject).(state)
%       .label
%       .dimord      = 'rpt_chan_freq'
%       .freq        = 1:nFeat
%       .featlabels  = {'beta_3_35','offset_3_35','beta_3_13','beta_13_30'}
%       .powspctrm   = 1 x chan x feat
%       .ntrials
%
% NOTES:
% - Fits are done on the IRASA FRACTAL spectrum only
% - Fits use log10(power) ~ log10(freq)
% - beta = -slope
% - offset/intercept is the regression intercept in log10(power) units
%   for the 3-35 Hz fit
% - Uses powspctrm_linear_mean (the script stops if it is missing)
% =========================================================================

clear; clc;

%% ---------------- USER SETTINGS ----------------
dataDir = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_5s/80Hz_laplacian';
inFile  = fullfile(dataDir, 'powerData_full_spectrum_irasa.mat');
outFile = fullfile(dataDir, 'aperiodic_params_irasa.mat');

epsP    = 1e-20;   % guard against log10(0)
VERBOSE = true;

% Features to save in powspctrm
featlabels = { ...
    'beta_3_35', ...
    'offset_3_35', ...
    'beta_3_13', ...
    'beta_13_30' ...
    };

%% ---------------- LOAD ----------------
S = load(inFile, 'd_full_frac');

if ~isfield(S, 'd_full_frac')
    error('Could not find d_full_frac in %s', inFile);
end

d_full_frac = S.d_full_frac;
subs = fieldnames(d_full_frac);

if isempty(subs)
    error('d_full_frac is empty.');
end

%% ---------------- MAIN LOOP ----------------
d_aperiodic = struct();

for iS = 1:numel(subs)
    subjFld = subs{iS};
    states  = fieldnames(d_full_frac.(subjFld));

    for iSt = 1:numel(states)
        stateFld = states{iSt};
        A = d_full_frac.(subjFld).(stateFld);

        % ---- Get linear mean spectrum: chan x freq
        if isfield(A, 'powspctrm_linear_mean')
            X = squeeze(A.powspctrm_linear_mean);   % expected chan x freq
        else
            error('Field powspctrm_linear_mean missing for %s | %s', subjFld, stateFld);
        end

        freq = double(A.freq(:))';

        % enforce chan x freq
        if size(X,2) ~= numel(freq)
            X = X';
        end

        if size(X,2) ~= numel(freq)
            error('Could not align channels x frequencies for %s | %s', subjFld, stateFld);
        end

        nChan = size(X,1);
        nFeat = numel(featlabels);

        % Feature matrix: 1 x chan x feat (same layout as the other measures)
        P = nan(1, nChan, nFeat);

        % Raw fit details per channel
        rawSlope_3_35   = nan(nChan,1);
        rawOffset_3_35  = nan(nChan,1);
        rawR2_3_35      = nan(nChan,1);
        rawSlope_3_13   = nan(nChan,1);
        rawR2_3_13      = nan(nChan,1);
        rawSlope_13_30  = nan(nChan,1);
        rawR2_13_30     = nan(nChan,1);

        for ch = 1:nChan
            yLin = double(X(ch,:));

            fit_3_35  = fit_loglog(freq, yLin, [3 35], epsP);
            fit_3_13  = fit_loglog(freq, yLin, [3 13], epsP);
            fit_13_30 = fit_loglog(freq, yLin, [13 30], epsP);

            % Save requested features in powspctrm
            P(1,ch,1) = fit_3_35.beta;       % beta_3_35
            P(1,ch,2) = fit_3_35.intercept;  % offset_3_35 (log10 power intercept)
            P(1,ch,3) = fit_3_13.beta;       % beta_3_13
            P(1,ch,4) = fit_13_30.beta;      % beta_13_30

            % Save extra raw fit details
            rawSlope_3_35(ch)  = fit_3_35.slope;
            rawOffset_3_35(ch) = fit_3_35.intercept;
            rawR2_3_35(ch)     = fit_3_35.r2;

            rawSlope_3_13(ch)  = fit_3_13.slope;
            rawR2_3_13(ch)     = fit_3_13.r2;

            rawSlope_13_30(ch) = fit_13_30.slope;
            rawR2_13_30(ch)    = fit_13_30.r2;
        end

        % ---- Build output struct (same layout as the other measures)
        out = struct();
        out.label      = A.label;
        out.dimord     = 'rpt_chan_freq';
        out.freq       = 1:nFeat;              % feature index
        out.featlabels = featlabels;
        out.powspctrm  = P;
        out.ntrials    = get_ntrials_safe(A);

        % Fit metadata
        out.fit_ranges_hz = { ...
            [3 35], ...
            [3 35], ...
            [3 13], ...
            [13 30] ...
            };
        out.feature_notes = { ...
            'positive spectral exponent beta = -slope from 3-35 Hz', ...
            'intercept from log10(power) ~ log10(freq) fit over 3-35 Hz', ...
            'positive spectral exponent beta = -slope from 3-13 Hz', ...
            'positive spectral exponent beta = -slope from 13-30 Hz' ...
            };

        % Raw slopes, offsets and R^2 per channel
        out.slope_3_35_raw    = rawSlope_3_35;
        out.offset_3_35_raw   = rawOffset_3_35;
        out.r2_3_35           = rawR2_3_35;

        out.slope_3_13_raw    = rawSlope_3_13;
        out.r2_3_13           = rawR2_3_13;

        out.slope_13_30_raw   = rawSlope_13_30;
        out.r2_13_30          = rawR2_13_30;

        d_aperiodic.(subjFld).(stateFld) = out;

        if VERBOSE
            fprintf('[APERIODIC] %s | %s | ntrials=%g | channels=%d\n', ...
                subjFld, stateFld, out.ntrials, nChan);
        end
    end
end

%% ---------------- SAVE ----------------
save(outFile, 'd_aperiodic', '-v7.3');
fprintf('\nSaved aperiodic parameters to:\n%s\n', outFile);

%% ================= HELPER FUNCTIONS =================

function fitOut = fit_loglog(freq, powLin, fitRange, epsP)
% Fit log10(power) ~ log10(freq)
% Returns:
%   slope      raw regression slope (usually negative)
%   beta       positive exponent = -slope
%   intercept  regression intercept in log10(power) units
%   r2         coefficient of determination

    idx = freq >= fitRange(1) & freq <= fitRange(2);
    idx = idx & isfinite(freq) & isfinite(powLin) & powLin > 0;

    f = double(freq(idx));
    p = double(powLin(idx));

    if numel(f) < 3
        fitOut.slope     = NaN;
        fitOut.beta      = NaN;
        fitOut.intercept = NaN;
        fitOut.r2        = NaN;
        return;
    end

    x = log10(f(:));
    y = log10(p(:) + epsP);

    coeff = polyfit(x, y, 1);
    yhat  = polyval(coeff, x);

    ss_res = sum((y - yhat).^2);
    ss_tot = sum((y - mean(y)).^2);

    if ss_tot == 0
        r2 = NaN;
    else
        r2 = 1 - (ss_res / ss_tot);
    end

    fitOut.slope     = coeff(1);
    fitOut.beta      = -coeff(1);
    fitOut.intercept = coeff(2);
    fitOut.r2        = r2;
end

function n = get_ntrials_safe(A)
    if isfield(A, 'ntrials') && ~isempty(A.ntrials)
        n = double(A.ntrials);
    else
        n = NaN;
    end
end