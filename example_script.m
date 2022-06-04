%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

% clear all
% clc
%% Preprocessing

% Parameters initialization
init_path = pwd;
choice_struct.path_data = [init_path, '\example_signal\AbleC EO 01.EDF'];
choice_struct.path_eeglab = [init_path, '\EEGlab\eeglab2021.0'];
choice_struct.begin_cropp = 385; choice_struct.end_cropp = 385;
choice_struct.line_noise = "Spectrum Interpolation"; choice_struct.line_noise_f = 60;
choice_struct.path_fieldtrip = [init_path, '\fieldtrip']; choice_struct.high_pass = 1;
choice_struct.low_pass = 100; choice_struct.rereference = "Median";
choice_struct.clean_data = "y"; choice_struct.options = load('bad_chans_options');
choice_struct.clean_artifacts = "y"; choice_struct.wind = 4; choice_struct.overlap = 10;
choice_struct.Kth = 1.15; choice_struct.ArtefThreshold = 6;
choice_struct.resample = "n"; choice_struct.drop_channs = "y";
choice_struct.channels2drop = ["AA","LABEL"];

return_struct = preprocessing(choice_struct);

%% Swarm Decomposition

% Parameters initialization
param_struct.rhythm_names = ["delta","theta","alpha","lbeta","hbeta","lgamma","hgamma"];
param_struct.rhythm_lims.delta = [1 4]; param_struct.rhythm_lims.theta = [4 8];
param_struct.rhythm_lims.alpha = [8 13]; param_struct.rhythm_lims.lbeta = [13 18];
param_struct.rhythm_lims.hbeta = [18 25]; param_struct.rhythm_lims.lgamma = [25 57];
param_struct.rhythm_lims.hgamma = [63 100];
param_struct.SwD_parameters.min_peak = 0.02; param_struct.SwD_parameters.component_std = 0.01;
param_struct.SwD_parameters.welch_wind = 0.5; param_struct.SwD_parameters.to_fullfillSpectrum = 0;
param_struct.SwD_parameters.clustering_factor = 0.1;

SwDResults = SwD_preprocessed_EEGs(return_struct,param_struct);

%% Connectivity Analsysis & Feature Extraction
conn_struct.feature_name = "coher"; conn_struct.chans = "all";
conn_struct.reg = "intra-regional";
feature = ConnectivityAnalysis(SwDResults, conn_struct);




