%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function [data_wICA,icaEEG2]= wICA(data,Fs,chan_names,flag,Kthr,ArtefThreshold)
% In this function the wICA algorithm is implemented to remove artifacts
% from EEG signals, as proposed in N.P. Castellanos and V.A. Makarov (2006), 
% "Recovering EEG brain signals: Artifact suppression with wavelet enhanced
% independent component analysis", J. Neurosci. Methods, 158, 300-312.
% The original script and some examples can be found http://www.mat.ucm.es/~vmakarov/downloads.php

%% Inputs: 
% data            -double matrix. NxM matrix, N the number of electrodes
%                  and M the data length
% Fs              -double. The sampling rate in Hz
% chan_names      -string array. Contains the names of all electrodes.
%                  Useful only if flag ~= 0 (for plots)
% flag            -double. If flag == 0 no plots, if flag == 1 only the
%                  cleaned EEGs will be plotted, if flag == 2 only the
%                  cleaned ICs will be plotted and if flag == 3 everything
%                  will be plotted.
% Kthr            -Tolerance for cleaning artifacts. try 1, 1.15, 1.2 etc
% ArtefThreshold  -Threshold for detection of ICs with artefacts
%                  Set lower values if you manually select ICs with 
%                  artifacts. Otherwise increase.
%% Outputs:
% data_wICA       -double matrix. The wICA cleaned EEGs
% icaEEG2         -double matrix. The wavelet filters ICs

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------
    data = double(data);
    
    %% Remove mean values from channels 
    data  = detrend(data,'constant');
    %% Plot raw data
    if flag==3
        figure('color','w');
        PlotEEG(data, Fs, chan_names, 200, 'Raw data (19 channels 10-20 system)');
    end
    %% Find independent components
    [icaEEG, A, W] = fastica(data,'stabilization','off','verbose','off'); 
    
    %Plot independent components
    if flag==3
        figure('color','w');
        PlotEEG(icaEEG, Fs, [], [], 'Independent Components');
        xlabel('Time (s)')
    end
    %% wICA Artifact Rejection
    
    nICs = 1:size(icaEEG,1); % Components to be processed, e.g. [1, 4:7]
    verbose = 'on';                         

    icaEEG2 = RemoveStrongArtifacts(icaEEG, nICs, Kthr, ArtefThreshold, Fs, verbose); 
    if (flag==3 || flag == 2)
        figure('color','w');
        PlotEEG(icaEEG2, Fs, [], [], 'wavelet filtered ICs');
    end
    data_wICA = A*icaEEG2;
    if (flag == 3 || flag == 1)
        figure('color','w');
        PlotEEG(data_wICA, Fs, chan_names, 100, 'wICA cleanned EEG');
    end
end