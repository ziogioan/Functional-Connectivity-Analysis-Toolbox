%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function [results] = SwD_brain_activity_criteria(SwDdata,paramStruct)

% Description: This function segments a SwD-processed EEG and keeps only its best parts. 
% The part of the EEG that is kept are the windows that demonstrate the 
% higher activity according to some brain activity criteria. The criteria are mainly based on the
% quantity (not the quality) of brain activity in each channel, window and major brain
% regions. Windows with "high quantitative activity" are considered more
% important and thus more likely to be kept in the shorter version of the
% input EEG. An option to segment the EEG after this process into k equal
% parts is given.

%% Inputs:
%          SwDdata - struct array or path of a .mat file or path of the folder where the 
%                            SwD-processed datasets (.mat files) have been saved (passed on from SwD_preprocessed_EEGs.m).
%
%               Must contain fields:
%
%                       rhythms - 1xN struct array with K fields, where N the number of windows and K the number of scalp electrodes. 
%                       
%                       chan_names - 1xK string array. Names of scalp channels
%               
%                       new_Fs - int scalar. New sampling rate after possible changes within the SwD function
%
%                       info - struct array. One field for each rhythm calculated, as specified by 'which_rhythms' variable.
%                              Each field is a struct, with fields:
%                                   'how_many_times' -> NxK table, how many times this rhythm existed
%                                                       in each electrode and in each time window
%                                   'exists' -> NxK table, whether a rhythm appeared or not
%                                                in a certain electrode and certain time window
%                       wind_dur - scalar. Window duration of input data in seconds
%
%                       new_length - scalar. Duration of EEG input in minutes
%
%                       num_wind - int scalar. Number of windows in EEG signal
%
%                       name - char array or string. Name of EEG file
%                          (e.g. patient name)
%
%          paramStruct   - struct array. Contains necessary parameters
%                               for Swarm Decomposition function
%
%               Must contain fields:
%
%                         wanted_length - desired length of EEG after segmenting
%
%                         wanted_rhythm_names - string array, names of
%                                   rhythms to be included in the segmenting process criteria
%
%                         saveRes - boolean. Save results in designated path
%
%                         savePath - char array or string scalar. Path(folder) to save results
%
%                         channel_act_thresh - Percentage (0-1) of the specified new time 
%                                               that must be active in a channel 
%
%                         regional_activity_thresh - Thresholds for regional activity
%                                   One threshold for each group of electrodes: regions = ["F","P","C","O","T"]
%                                   Example for 10-20 system: [3,2,2,1,2] (Out of 7,3,3,2,4 electrodes respectively)
%                         window_activity_thresh - Threshold for discarding a window
%                                                   Must be equal to or
%                                                   lower than 0.25 (0-1)
%
%                         divide_in_k - int scalar, number of parts to divide the EEG after application of all criteria
%% Outputs: 
%          results - cell array with k cells, one for each part.
%                    Each cell is a SwD struct with identical format as the
%                    input struct SwDdata.
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------
%% Determine input type

if isfile(SwDdata)
    warning('You have provided the path to a single .mat SwD-processed EEG file');
elseif isfolder(SwDdata)
    warning('You have provided the path to a folder of .mat SwD-processed EEG files');
elseif isstruct(SwDdata)
    warning('You have provided an SwD-processed EEG data struct');
end

%% Check universal input parameters
if ~isfield(paramStruct,'wanted_rhythm_names') || (isfield(paramStruct,'wanted_rhythm_names') && isempty(paramStruct.wanted_rhythm_names))
    paramStruct.wanted_rhythm_names = ["delta","theta","alpha","lbeta","hbeta","lgamma","hgamma"];
end
if ~isfield(paramStruct,'channel_act_thresh') || (isfield(paramStruct,'channel_act_thresh') && isempty(paramStruct.channel_act_thresh))
    paramStruct.channel_act_thresh = 0.25; %At least 1/4 of the specified new time must be active in a channel  
end
if ~isfield(paramStruct,'regional_activity_thresh') || (isfield(paramStruct,'regional_activity_thresh') && isempty(paramStruct.regional_activity_thresh))
    paramStruct.regional_activity_thresh = [3,2,2,1,2]; %Thresholds for regional activity
        %One threshold for each group of electrodes: regions = ["F","P","C","O","T"]
end
if ~isfield(paramStruct,'window_activity_thresh') || (isfield(paramStruct,'window_activity_thresh') && isempty(paramStruct.window_activity_thresh))
    paramStruct.window_activity_thresh = 0; %Threshold for discarding a window. Must be equal to or lower than 0.25
end
if ~isfield(paramStruct,'divide_in_k') || (isfield(paramStruct,'divide_in_k') && isempty(paramStruct.divide_in_k))
    paramStruct.divide_in_k = 1; %default value for dividing EEG into parts
end
if ~isfield(paramStruct,'saveRes') || (isfield(paramStruct,'saveRes') && isempty(paramStruct.saveRes))
        paramStruct.saveRes = 0; %default value is no save results
end

%% Case single file or data
if isfile(SwDdata) || isstruct(SwDdata)
    if isfile(SwDdata)
        pathDataDir = dir(SwDdata);
        name = pathDataDir.name;
        SwDstruct = load([pathDataDir.folder, '/',name]);
    elseif isstruct(SwDdata)
        SwDstruct = SwDdata;       
    end
    %% Check if SwD struct contains all necessary fields
    if ~isfield(SwDstruct,'rhythms')
        error('SwD struct must have a rhythms field')
    end
    if isfield(SwDstruct,'name')
        error('SwD struct must have a name field')
    end
    if ~isfield(SwDstruct,'new_Fs')
        error('SwD struct must have a new_Fs field')
    end
    if ~isfield(SwDstruct,'info')
        error('SwD struct must have an info field')
    end
    if ~isfield(SwDstruct,'wind_dur')
        error('SwD struct must have a wind_dur field')
    end
    if ~isfield(SwDstruct,'new_length')
        error('SwD struct must have a new_length field')
    end
    if ~isfield(SwDstruct,'chan_names')
        error('SwD struct must have a chan_names field')
    end
    
    %% Check case input parameters
    if ~isfield(paramStruct,'wanted_length') || (isfield(paramStruct,'wanted_length') && isempty(paramStruct.wanted_length))
        paramStruct.wanted_length = SwDstruct.new_length/4; %default value for new EEG length
    end
    
    [segmented_eeg] = segment_eeg(SwDstruct,paramStruct.wanted_length,paramStruct.wanted_rhythm_names);
        
    %% Divide the remaining EEG in k parts
    [results] = kparts_rhythms(segmented_eeg,paramStruct.divide,paramStruct.wanted_rhythm_names);
    
    %% Return results OR/AND Save k parts of segmented EEG separately
       
    if paramStruct.saveRes
        for i = 1:length(results)
            segmented_eeg = results{i};
            try
                save([paramStruct.savePath, name,'_part',int2str(i), '_rhythms_',int2str(new_length),...
                    'min','.mat'], 'segmented_eeg')
            catch
                warning('Results were not saved since you havent provided a valid folder path')
            end
        end
    end
    %Function returns
    
%% Case folder with multiple files
elseif isfolder(SwDdata)
    pathDataDir = dir(SwDdata);
    %counter = 1; 
    for j = 1:length(pathDataDir) - 2   
        name = pathDataDir(j+2).name;
        SwDstruct = load([SwDdata, '/',name]);
        %% Check if SwD struct contains all necessary fields
        if ~isfield(SwDstruct,'rhythms')
            error('SwD struct must have a rhythms field')
        end
        if isfield(SwDstruct,'name')
            error('SwD struct must have a name field')
        end
        if ~isfield(SwDstruct,'new_Fs')
            error('SwD struct must have a new_Fs field')
        end
        if ~isfield(SwDstruct,'info')
            error('SwD struct must have an info field')
        end
        if ~isfield(SwDstruct,'wind_dur')
            error('SwD struct must have a wind_dur field')
        end
        if ~isfield(SwDstruct,'new_length')
            error('SwD struct must have a new_length field')
        end
        if ~isfield(SwDstruct,'chan_names')
            error('SwD struct must have a chan_names field')
        end
        %% Check case input parameters
        if ~isfield(paramStruct,'wanted_length') || (isfield(paramStruct,'wanted_length') && isempty(paramStruct.wanted_length))
            paramStruct.wanted_length = SwDstruct.new_length/4; %default value for new EEG length
        end

        [segmented_eeg] = segment_eeg(SwDstruct,paramStruct.wanted_length,paramStruct.wanted_rhythm_names);

        %% Divide the remaining EEG in k parts
        [results] = kparts_rhythms(segmented_eeg,paramStruct.divide,paramStruct.wanted_rhythm_names);

        %% Return results OR/AND Save k parts of segmented EEG separately

        if paramStruct.saveRes
            for i = 1:length(results)
                segmented_eeg = results{i};
                try
                    save([paramStruct.savePath, name,'_part',int2str(i), '_rhythms_',int2str(new_length),...
                        'min','.mat'], 'segmented_eeg')
                catch
                    warning('Results were not saved since you havent provided a valid folder path')
                end
            end
        end
        clearvars SwDstruct segmented_eeg results 
%         disp(counter)
%         counter = counter + 1;
    end

end