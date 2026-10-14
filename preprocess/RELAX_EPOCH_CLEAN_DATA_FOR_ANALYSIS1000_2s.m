%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2022) Neil Bailey

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see https://www.gnu.org/licenses/.

%% RELAX_EPOCH_CLEAN_DATA_FOR_ANALYSIS:

clear all; close all; clc;
%% DEPENDENCIES (toolboxes you need to install, and cite if you use this script):

% MATLAB signal processing toolbox (from MATLAB website)

% EEGLAB:
% https://sccn.ucsd.edu/eeglab/index.php
% Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of neuroscience methods, 134(1), 9-21.
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0');
eeglab('nogui');

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0/plugins/Fieldtrip-lite250523');

% Specify  RELAX folder location (this toolbox):
addpath('/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/software/eeglab2025.1.0/plugins/RELAX2.0.1');

% Specify rejection parameters:
RELAX_epoching_cfg.InterpolateRejectedChannels='yes';
RELAX_epoching_cfg.SingleChannelImprobableDataThreshold=5; %MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
RELAX_epoching_cfg.AllChannelImprobableDataThreshold=3; %SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
RELAX_epoching_cfg.SingleChannelKurtosisThreshold=5; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis 
RELAX_epoching_cfg.AllChannelKurtosisThreshold=3; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis 
RELAX_epoching_cfg.reject_amp=120; % Absolute voltage amplitude threshold - if an epoch has voltages that deviate from 0 by more than this value the epoch is rejected

% The following muscle affected epoch rejection may 
% be useful if you want to ensure no epochs included in
% the analysis are still affected by muscle activity:
RELAX_epoching_cfg.MuscleSlopeThreshold=-0.31; % log frequency/ log power slope threshold for muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or more stringent = -0.72, more negative thresholds remove more muscle. 
RELAX_epoching_cfg.MaxProportionOfMuscleEpochsToClean=0.50; % Maximum proportion of muscle epochs to remove (if more epochs than this are affected, only the worst effected up to this proportion will be removed)
RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity='no'; % 'yes' or 'no'

% set to reject epochs if HEOG movements larger than a specified threshold were detected using the step function described by Luck (2014):
% Luck, S. J. (2014). An introduction to the event-related potential technique. MIT press.
% (this setting is only helpful if participants were required to focus on a specific point, and to not move their eyes during the task)
RELAX_epoching_cfg.Reject_HEOG_movements=0; % set to 1 to reject epochs based on horizontal eye movements during the epoch
RELAX_epoching_cfg.HEOG_electrode_labels={}; % set the labels used for the HEOG electrodes
RELAX_epoching_cfg.HEOG_rejection_period=[-200 300]; % in ms, the period within which to search for HEOG movements and reject epochs if HEOG movements are found, relative to stimuli onset
RELAX_epoching_cfg.HEOG_rejection_threshold=32; % HEOG rejection threshold (16 microvolts = 1 degree of eye movement)

RELAX_epoching_cfg.Reject_VEOG_movements=0; % reject based on probable blinks?
RELAX_epoching_cfg.VEOG_electrode_labels={}; % set the labels used for the VEOG electrodes
RELAX_epoching_cfg.VEOG_rejection_period=[-200 200]; % in ms, the period within which to search for VEOG movements and reject epochs if HEOG movements are found, relative to stimuli onset
RELAX_epoching_cfg.VEOG_rejection_threshold=60; % VEOG rejection threshold
RELAX_epoching_cfg.compute_VEOG_difference=0; % compute difference between a lower and upper VEOG electrode? (if set to 1, make sure both are listed in the 1st and 2nd position for VEOG_electrode_labels)

RELAX_epoching_cfg.DataType='Resting'; % 'Task' for cognitive tasks, 'Resting' for data without stimuli presented
RELAX_epoching_cfg.restingdatatriggerinterval=2.048; % (seconds) sets how often an 'X' trigger is inserted into resting data for epoching. 3.5 with [-2.5 2.5] PeriodToEpoch provides 5s epochs with 1.5s overlaps
RELAX_epoching_cfg.TriggersToEpoch={ }; % triggers to use for epoching
RELAX_epoching_cfg.PeriodToEpoch=[0 2.048]; % [start end] surrounding trigger
RELAX_epoching_cfg.RemoveOtherTriggers='no'; % 'yes'|'no' - removes all other triggers except the trigger used for timelocking the epoch

% Baseline correct the data?
RELAX_epoching_cfg.BL_correction_method='none'; % Set type of baseline correction to perform ('regression' = recommended, 'subtraction' = traditional but not recommended)
RELAX_epoching_cfg.BLperiod=[-200 0]; % Set baseline period for baseline correction of data to reduce potential influence of drift
RELAX_epoching_cfg.NumberOfFactors=2;
RELAX_epoching_cfg.BL_correction_Factor_1_Level_1={ }; % triggers to include in factor 1's level 1 for regression BL correction (all other triggers are included in factor 1's level 2)
RELAX_epoching_cfg.BL_correction_Factor_2_Level_1={ }; % triggers to include in factor 2's level 1 for regression BL correction (all other triggers are included in factor 2's level 2)

% Specify whether all data is in a single folder or data are in BIDS format
% (each EEG file within its own separate folder):
RELAX_epoching_cfg.all_data_in_1_folder_or_BIDS_format='folder'; % set to: 'BIDS' or 'folder'

% Specify the 'to be epoched' file locations: 
% If including all data in a single folder, specify the folder containing
% the cleaned data, otherwise if using BIDS format (one sub-folder per
% participant), then specify the folder containing all the sub-folders.
RELAX_epoching_cfg.CleanedPath=['/data/mgh-meditation/adept_group_study/data_processed/eeg_relax/all_raw/all_set/RELAXProcessed_final/Cleaned_Data'];

% Load pre-processing statistics file for these participants if it already
% exists:
RELAX_epoching_cfg.OutputPath=[RELAX_epoching_cfg.CleanedPath filesep 'Epoched'];
mkdir(RELAX_epoching_cfg.OutputPath);

%% List all files in directory
cd(RELAX_epoching_cfg.CleanedPath);
if strcmpi(RELAX_epoching_cfg.all_data_in_1_folder_or_BIDS_format,'folder')
    RELAX_epoching_cfg.dirList=dir('*.set');
    RELAX_epoching_cfg.files={RELAX_epoching_cfg.dirList.name};
elseif strcmpi(RELAX_epoching_cfg.all_data_in_1_folder_or_BIDS_format,'BIDS')
    RELAX_epoching_cfg.dirList=dir([RELAX_epoching_cfg.CleanedPath '\**\*_RELAX.set']); 
    % If you have cleaned your data with RELAX, then the suffix "_RELAX" will ensure only the clean data is found here.
    % If you used a different method to clean your data, include the suffix that indicates the cleaned data 
    % (otherwise the script will epoch the raw data as well as the epoched data)
    already_processed=[];
    for f=1:size(RELAX_epoching_cfg.dirList,1)
        if contains(RELAX_epoching_cfg.dirList(f).name,'Epoched')
            already_processed(f)=1;
        end
    end
    RELAX_epoching_cfg.dirList(already_processed==1,:)=[]; % remove filenames from this list if they have already been epoched
    RELAX_epoching_cfg.files={RELAX_epoching_cfg.dirList.name};
    RELAX_epoching_cfg.folders={RELAX_epoching_cfg.dirList.folder};
end
if isempty(RELAX_epoching_cfg.files)
    disp('No files found..')
end

%% RUN SCRIPT BELOW:
RELAX_epoching_cfg.FilesToProcess=1:numel(RELAX_epoching_cfg.files); % Set which files to process

[OutlierParticipantsToManuallyCheck,EpochRejections,RELAX_epoching_cfg] = RELAX_epoch_the_clean_data_Wrapper(RELAX_epoching_cfg);
    
        