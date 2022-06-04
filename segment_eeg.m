%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function [segmented_EEG] = segment_eeg(swarm_struct,paramStruct)
%%Segments preprocessed EEG according to channel and window criteria
%
%% Inputs: 
%        swarm_struct   - struct of structs. Contains swarm info: 
%                         new_Fs: sampling rate, scalar
%
%                         info: activity_info, struct
%
%                         chan_names: channel names, string
%
%                         rhythms: struct of structs, MUST BE N x 1 struct,
%                         with C fields. (N: number of windows, C:
%                         channels)
%
%                         activity_info - struct with fields the names of each rhythm. To
%                                         be used for window selection
%
%                         win_length - Duration (in sec) of EEG epochs
%
%         paramStruct   - struct. Contains parameter info:
%                         wanted_dur    - Duration (in min) after EEG segmenting
%                         rhythm_names  - Names of rhythms
%                         channel_act_thresh - Percentage (0-1) of the specified new time 
%                                               that must be active in a channel 
%                         regional_activity_thresh - Thresholds for regional activity
%                                   One threshold for each group of electrodes: regions = ["F","P","C","O","T"]
%                                   Example for 10-20 system: [3,2,2,1,2] (Out of 7,3,3,2,4 electrodes respectively
%
%                         window_activity_thresh - Threshold for discarding a window
%                                                   Must be equal to or lower than 0.25
%
%% Outputs: 
%         segmented_EEG - struct of structs, with fields:
%
%                rhythms      -  struct that contains rhythms from the remaining
%                          windows after segmenting
%
%                info         -  struct that contains for each rhythm, two
%                          tables, exists and how_many_times. These tables contain binary
%                          and numeric info about existence of rhythm in each window
%
%                window_regions - Brain regions active in each window
%
%                chan_names   -  Same as input
%
%                Fs           -  Same as input
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------
%% Extract swarm info
    Fs = swarm_struct.new_Fs;
    chan_names = swarm_struct.chan_names;
    activity_info = swarm_struct.info;
    EEG_rhythms = swarm_struct.rhythms;
    win_length = swarm_struct.wind_dur;
    name = swarm_struct.name;
    
%% Extract parameter info
    wanted_dur = paramStruct.wanted_length;
    rhythm_names = paramStruct.wanted_rhythm_names;
    channel_act_thresh = paramStruct.channel_act_thresh;
    regional_activity_thresh = paramStruct.regional_activity_thresh;
    window_activity_thresh = paramStruct.window_activity_thresh;
    
%%
    wanted_dur = wanted_dur*60; %In seconds
    keep_windows = wanted_dur/win_length; %Number of windows after segmenting
    num_ry = length(rhythm_names);
    num_ch = length(chan_names);
    num_win = height(activity_info.delta.exists);
    
    window_flag = strings(num_ry,num_win);
    window_points = NaN(num_ry,num_win);
    window_sum_activity = NaN(1,num_win);
    
%% First discard channels based on channel activity criterion
%% then segment EEG based on Window Criterion

%% For each rhythm
    for k = 1:num_ry
        ry_info = activity_info.(rhythm_names(k));% Extract activity info
        exists = ry_info.exists;% Keep exists table
%% Channel criterion, for each channel
        for m = 1:num_ch
             [decision(m),~,~] = channel_criterion(exists,...
                 chan_names(m),wanted_dur,win_length,channel_act_thresh);
             if strcmp(decision(m),"discard")
                %Discard rhythm from all windows of channel m
                for n = 1:num_win
                     EEG_rhythms(n).(chan_names(m)).(rhythm_names(k)) = [];
                end
             end            
        end
%% Window Criterion, for each window  
% Assign points to segment according to total score of a window
        points = [0,1,2,3,4]; % [discard, local ,medium, high, very high]
        n_points = [points; 1.5*points; points; 1.2*points; 1.2*points; 1.2*points; points];
        for n = 1:num_win
            window = exists(n,:); % (1 x num_ch) table
            [window_flag(k,n),window_regions{k,n}] = window_criterion(window,chan_names,...
                regional_activity_thresh,window_activity_thresh);
            if strcmp(window_flag(k,n),"discard")
                window_points(k,n) = n_points(k,1);
            elseif strcmp(window_flag(k,n),"local")
                window_points(k,n) = n_points(k,2);
            elseif strcmp(window_flag(k,n),"medium")
                window_points(k,n) = n_points(k,3);
            elseif strcmp(window_flag(k,n),"high")
                window_points(k,n) = n_points(k,4);
            elseif strcmp(window_flag(k,n),"vhigh")
                window_points(k,n) = n_points(k,5);
            end
            
        end
    end
    
    window_sum_activity = sum(window_points,1,'omitnan'); %Sum of each column, a score for each window is calculated
 
    %Windows are sorted based on activity
    [sorted_activity,indexes] = sort(window_sum_activity,2,'descend');
    
    best_indexes = indexes(1:keep_windows);
    
    %Transform rhythms struct to cell, to segment EEG
    fnames = fieldnames(EEG_rhythms);
    size_x = length(fnames);
    size_y = length(EEG_rhythms);
    best_EEG_rhythms = struct2cell(EEG_rhythms);    
    best_EEG_rhythms = reshape(best_EEG_rhythms,[size_x, size_y]);
    best_EEG_rhythms = transpose(best_EEG_rhythms);
    best_EEG_rhythms = best_EEG_rhythms(best_indexes,:);
    
    best_EEG_rhythms = cell2struct(best_EEG_rhythms,fnames,2);%Keep best windows
    window_regions = transpose(window_regions);
    best_window_regions = window_regions(best_indexes,:); %Keep regions of best windows
    
    for k = 1:num_ry
        ry_info = activity_info.(rhythm_names(k));% Extract activity info
        exists = ry_info.exists;% Keep exists table
        how_many_times = ry_info.how_many_times;
        best_exists = how_many_times(best_indexes,:); %Keep info for best windows
        best_how_many_times = exists(best_indexes,:); %Keep info for best windows
        best_activity_info.(rhythm_names(k)).exists = best_exists; 
        best_activity_info.(rhythm_names(k)).how_many_times = best_how_many_times;
    end

%% Outputs
    % Reconstruct input struct
    segmented_EEG.rhythms = best_EEG_rhythms;
    segmented_EEG.window_regions = best_window_regions;
    segmented_EEG.info = best_activity_info;
    segmented_EEG.new_Fs = Fs;
    segmented_EEG.chan_names = chan_names;
    segmented_EEG.name = name;
    segmented_EEG.wind_dur = win_length;%in sec
    segmented_EEG.new_length = wanted_dur/60; %in min
    %Add a new field for time indexes of windows
    segmented_EEG.time_indexes = best_indexes;

end

