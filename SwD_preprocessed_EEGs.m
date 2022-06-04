%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function SwDResults = SwD_preprocessed_EEGs(EEGdata,paramStruct)

% Description: This function decomposes the EEG into its oscillatory modes
% as per the Swarm Decomposition algorithm. The SwD process takes place for
% each EEG channel and for each EEG epoch/window if available.

%% Inputs:  
%          EEGdata       - struct array or path of a .mat file or path of the folder where the 
%                            preprocessed datasets (.mat files) have been saved (passed on from preprocessing).
%
%               Must contain fields:
%
%                        EEG_data - cell array, each cell is a window and contains all scalp electrodes.  
%
%                        Fs - int scalar. Sampling rate of input data
%
%                        wind_dur - scalar. Window duration of input data in seconds
%
%                        new_length - scalar. Duration of EEG input in minutes
%
%                        num_wind - int scalar. Number of windows in EEG signal
%
%                        name - char array or string. Name of EEG file
%                          (e.g. patient name)
%
%          paramStruct   - struct array. Contains necessary parameters
%                               for Swarm Decomposition function
%
%               Must contain fields:
%
%                          which_rhythms - string array. Names of rhythms
%                                          to be saved from the Swarm Decomposition procedure
%                                          Format of names must be as follows:
%                                          rhythm_names = ["delta","theta","alpha","lbeta","hbeta","lgamma","hgamma"];
%
%                          rhythm_lims - struct array. Fields of struct are
%                                        the names of the rhythms to save in the SwD
%                                        procedure. Format must be identical to 'which_rhythms'
%
%                          saveRes - boolean. Save results in designated path
%
%                          savePath - char array or string scalar. Path(folder) to save results
%
%                          pathSwDFiles  - char array. Path to location of SwD function
%
%% Outputs:
%          SwDResults - struct array. If a path is provided (savePath
%                       variable) the results will be also saved in specified path
%
%              Should contain fields:
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
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------
%% Determine input type

if ischar(EEGdata)
    if isfile(EEGdata)
        warning('You have provided the path to a single .mat EEG file');
    elseif isfolder(EEGdata)
        warning('You have provided the path to a folder of .mat EEG files');
    end
end
if isstruct(EEGdata)
    warning('You have provided an EEG data struct');
end

%% Check input parameters and apply default values where needed
if ~isfield(paramStruct,'rhythm_names') || (isfield(paramStruct,'rhythm_names') && isempty(paramStruct.rhythm_names))
        paramStruct.rhythm_names = ["delta","theta","alpha","lbeta","hbeta","lgamma","hgamma"];
end
rhythm_names = paramStruct.rhythm_names;
if ~isfield(paramStruct,'rhythm_lims') || (isfield(paramStruct,'rhythm_lims') && isempty(paramStruct.rhythm_lims))
        paramStruct.rhythm_lims.delta = [1 4]; paramStruct.rhythm_lims.theta = [4 8];
        paramStruct.rhythm_lims.alpha = [8 13]; paramStruct.rhythm_lims.lbeta = [13 18];
        paramStruct.rhythm_lims.hbeta = [18 25]; paramStruct.rhythm_lims.lgamma = [25 57];
        paramStruct.rhythm_lims.hgamma = [63 100]; 
end
if ~isfield(paramStruct,'which_rhythms') || (isfield(paramStruct,'which_rhythms') && isempty(paramStruct.which_rhythms))
        paramStruct.which_rhythms = ["delta","theta","alpha","lbeta","hbeta","lgamma","hgamma"]; %default value is all rhythms
end
if ~isfield(paramStruct,'saveRes') || (isfield(paramStruct,'saveRes') && isempty(paramStruct.saveRes))
        paramStruct.saveRes = 0; %default value is no save results
end

%% Case single file or data
if ischar(EEGdata) || isstruct(EEGdata)
    if ischar(EEGdata) && contains(EEGdata,'.mat')
        pathDataDir = dir(EEGdata);
        name = pathDataDir.name;
        EEGstruct = load([pathDataDir.folder, '/',name]);
        data = EEGstruct.EEG_data;    chan_names = EEGstruct.chan_names;
        chan_names_all = EEGstruct.chan_names_all;
    elseif isstruct(EEGdata)
        EEGstruct = EEGdata;
        data = EEGstruct.EEG_data;    chan_names = EEGstruct.chan_names;
        chan_names_all = EEGstruct.chan_names_all;
        if isfield(EEGstruct,'EEG_data')
            data = EEGstruct.EEG_data;
        else
            error('EEG struct must have an EEG_data field')
        end
        if isfield(EEGstruct,'name')
            name = EEGstruct.name;
        else
            error('EEG struct must have a name field')
        end
        if isfield(EEGstruct,'chan_names')
            chan_names = EEGstruct.chan_names;
        else
            error('EEG struct must have a chan_names field')
        end
        if isfield(EEGstruct,'chan_names_all')
            chan_names_all = EEGstruct.chan_names_all;
        else
            error('EEG struct must have a chan_names_all field')
        end
    end
    
    wind_dur = EEGstruct.wind_dur; new_length = EEGstruct.new_length;
    num_wind = EEGstruct.num_wind;
    %Fill possible spaces in chan_names
    for i = 1:length(chan_names)
        idx = strfind(chan_names(i),' ');
        if ~isempty(idx)
            temp = convertStringsToChars(chan_names(i));
            temp(idx) = '_';
            chan_names(i) = convertCharsToStrings(temp);
        end
        idx = [];
    end
    
    %% Rhythm Extraction
    SwD_par = paramStruct.SwD_parameters;
    counter = 1;
    for i = 1:length(data)           
%         try
            
            [rhythms(counter), new_Fs] = rhythm_extraction(table2array(data{i}), EEGstruct.Fs, chan_names, paramStruct.rhythm_lims, paramStruct.which_rhythms,SwD_par);
%         catch
%             cd(paramStruct.pathSwDFiles);
%             [rhythms(counter), new_Fs] = rhythm_extraction(table2array(data{i}), EEGstruct.Fs, chan_names, paramStruct.rhythm_lims, paramStruct.which_rhythms,SwD_par);
%         end
        counter = counter + 1;
        fprintf("\n")
        disp(i)
        fprintf("\n")
    end

    % Counter for how many times a rhythm exists.
    rhythm_exists = zeros(counter-1,length(chan_names),length(rhythm_names));
    how_many_times = zeros(counter-1,length(chan_names),length(rhythm_names));
    for i = 1:(counter-1) %For each window
        for m = 1:length(chan_names) % For each channel
            for k = 1:length(rhythm_names) % For each rhythm
                if ~isempty(rhythms(i).(chan_names(m)).(rhythm_names(k)))
                    temp = rhythms(i).(chan_names(m)).(rhythm_names(k));
                    rhythm_exists(i,m,k) = 1;
                    how_many_times(i,m,k) = length(temp(1,:));
                end
            end
        end 
    end

    for k = 1:length(rhythm_names)
        how_many_times_tbl = array2table(how_many_times(:,:,k),'VariableNames',cellstr(chan_names));
        rhythm_exists_tbl = array2table(rhythm_exists(:,:,k),'VariableNames',cellstr(chan_names));
        info.(rhythm_names(k)).how_many_times = how_many_times_tbl;
        info.(rhythm_names(k)).exists = rhythm_exists_tbl;
    end
    
    %% Return results
    SwDResults.rhythms = rhythms; SwDResults.chan_names = chan_names;
    SwDResults.new_Fs = new_Fs; SwDResults.info = info;
    SwDResults.wind_dur = EEGstruct.wind_dur; SwDResults.new_length = EEGstruct.new_length;
    SwDResults.name = name; SwDResults.num_wind = EEGstruct.num_wind;
    SwDResults.chan_names_all = chan_names_all;
    if isfield(paramStruct,'saveRes')
        if paramStruct.saveRes
            try
                save([paramStruct.savePath, name, '_rhythms.mat'],...
                'rhythms', 'chan_names', 'chan_names_all', 'new_Fs', 'info','wind_dur','new_length','num_wind','name')
            catch
                warning('SwD results were not saved since you havent provided a valid folder path')
            end
        end
    end
    clearvars rhythms info chan_names new_Fs

%% Case folder with multiple files 
elseif isfolder(EEGdata)
    pathDataDir = dir(EEGdata);
    
    for j = 1:length(pathDataDir) - 2
        name = pathDataDir(j+2).name;
        EEGstruct = load([EEGdata, '/',name]);

        %% Fetch EEG parameters from structs
        data = EEGstruct.EEG_data; chan_names = EEGstruct.chan_names;
        wind_dur = EEGstruct.wind_dur; new_length = EEGstruct.new_length;
        num_wind = EEGstruct.num_wind; chan_names_all = EEGstruct.chan_names_all;
        %Fill possible spaces in chan_names
        for i = 1:length(chan_names)
            idx = strfind(chan_names(i),' ');
            if ~isempty(idx)
                temp = convertStringsToChars(chan_names(i));
                temp(idx) = '_';
                chan_names(i) = convertCharsToStrings(temp);
            end
            idx = [];
        end
%         Fs = paramStruct.Fs; %Resampling was made, 256
%         wind_dur = paramStruct.wind_dur;%4 sec
%         new_length = paramStruct.new_length;%2 min
        %% Rhythm Extraction
        counter = 1;
        for i = 1:length(data)  %Segmenting 30*wind_dur:wind_dur:30*wind_dur% 
%             try
                [rhythms(counter), new_Fs] = rhythm_extraction(table2array(data{i}) , EEGstruct.Fs, chan_names, paramStruct.rhythm_lims, paramStruct.which_rhythms);
%             catch
%                 cd(paramStruct.pathSwDFiles);
%                 [rhythms(counter), new_Fs] = rhythm_extraction(table2array(data{i}) , EEGStruct.Fs, chan_names, paramStruct.rhythm_lims, paramStruct.which_rhythms);
%             end
            counter = counter + 1;
            fprintf("\n")
            disp(i)
            fprintf("\n")
        end

        % Counter for how many times a rhythm exists.
        rhythm_exists = zeros(counter-1,length(chan_names),length(rhythm_names));
        how_many_times = zeros(counter-1,length(chan_names),length(rhythm_names));
        for i = 1:(counter-1) %For each window
            for m = 1:length(chan_names) % For each channel
                for k = 1:length(rhythm_names) % For each rhythm
                    if ~isempty(rhythms(i).(chan_names(m)).(rhythm_names(k)))
                        temp = rhythms(i).(chan_names(m)).(rhythm_names(k));
                        rhythm_exists(i,m,k) = 1;
                        how_many_times(i,m,k) = length(temp(1,:));
                    end
                end
            end 
        end

        for k = 1:length(rhythm_names)
            how_many_times_tbl = array2table(how_many_times(:,:,k),'VariableNames',cellstr(chan_names));
            rhythm_exists_tbl = array2table(rhythm_exists(:,:,k),'VariableNames',cellstr(chan_names));
            info.(rhythm_names(k)).how_many_times = how_many_times_tbl;
            info.(rhythm_names(k)).exists = rhythm_exists_tbl;
        end
    %% Return results
        SwDResults.rhythms = rhythms; SwDResults.chan_names = chan_names;
        SwDResults.new_Fs = new_Fs; SwDResults.info = info;
        SwDResults.wind_dur = EEGstruct.wind_dur; SwDResults.new_length = EEGstruct.new_length;
        SwDResults.name = name; SwDResults.num_wind = EEGstruct.num_wind;
        SwDResults.chan_names_all = chan_names_all;
        if isfield(paramStruct,'saveRes')
            if paramStruct.saveRes
                try
                    save([paramStruct.savePath, name, '_rhythms.mat'],...
                    'rhythms', 'chan_names', 'chan_names_all', 'new_Fs', 'info','wind_dur','new_length','num_wind','name')
                catch
                    warning('SwD results were not saved since you havent provided a valid folder path')
                end
            end
        end 
    end
end

end

