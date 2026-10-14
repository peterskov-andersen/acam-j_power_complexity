%----------------------------------------------------------------------
% laplacian_35hz.m
% 35 Hz low-pass filter, then spline surface Laplacian, for all .set files
% in one folder. Settings shown are for the 2.048 s epochs. The script was
% run once per epoch length; the README lists the low-pass used for each.
%----------------------------------------------------------------------
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
% 0) Start EEGLAB (loads plugins, including currentdensity)
eeglab; close;  % start EEGLAB, then close the GUI

% 1) Input folder (epoched .set files)
inputDir  = '/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data/Epoched_2s';

% 2) Output folder
outputDir = fullfile(inputDir, '35hz_laplacian');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 3) List .set files
fileList = dir(fullfile(inputDir, '*.set'));
if isempty(fileList)
    error('No .set files found in %s', inputDir);
end

% 4) Process each file
for iFile = 1:numel(fileList)
    fn = fileList(iFile).name;
    fprintf('Processing [%d/%d]: %s\n', iFile, numel(fileList), fn);

    % 4a) Load
    EEG = pop_loadset('filename', fn, 'filepath', inputDir);

    % 4b) Low-pass @35 Hz
    EEG = pop_eegfiltnew(EEG, [], 35);   % 35 Hz low-pass

    % 4c) Spline surface Laplacian (lambda = 1e-5)
    EEG = pop_currentdensity(EEG, 'method', 'spline', 'lambda', 1e-05);

    % 4d) Rename and save
    [~, nm] = fileparts(fn);
    EEG.setname = [nm '35hz_lp'];
    outName = [EEG.setname '.set'];
    pop_saveset(EEG, 'filename', outName, 'filepath', outputDir);

    clear EEG;
end

disp('All files filtered at 35 Hz and Laplacian applied.');
