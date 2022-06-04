%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function return_struct = preprocessing(choice_struct)

% In this function the preprocessing of the EEG takes place.
% The data must be in EDF format.
% 
%% Inputs:
%              choice_struct   -Struct that contains the following
%                                information
% Obligatory:  path_data       -char. It should be the path of the folder where the 
%                               datasets are saved.
%              path_eeglab     -char. Path of eeglab.
% 
% Optional:    path_save       -char. The path where the data will be saved
%              begin_cropp     -double. length in seconds that will be
%                               cropped from the start [default = []]
%              end_cropp       -double. length in seconds that will be
%                               cropped from the end [default = []]
%              line_noise      -string. Line noise removal algorithm.
%                               Spectrum Interpolation or Clean Line or
%                               None [Default = None]
%              line_noise_f    -If line_noise ~= "None"
%                               the frequency of the noise should be
%                               specified in Hz.
%              path_fielftrip  -char. Path of fieldtrip (obligatory, if spectrum
%                               interpolation is used)
%              high_pass       -double. Frequency for the high pass filter
%                               in Hz [default = 1]
%              low_pass        -double. Frequency for the low pass filter
%                               in Hz [default = 100]
%              rereference     -string. Median or Average or None 
%                               [default =  Median]
%              clean_data      -string. y (for yes) or n (for no). [default
%                               = "n"]
%              options         -struct. if clean_data == "y" a struct with
%                               the options for the algorithm should be
%                               included. For more info chech the eeglab
%                               documentation and the bad_chans_options.mat
%                               file.
%              clean_artifacts -string. y (for yes) or n (for no).[default
%                               = "y"]
%              wind            -double. If clean_artifacts == "y" the
%                               length of the window, in which the algorithm
%                               will be implemented, should be included in
%                               seconds [default = 4]
%              overlap         -double. If clean_artifacts == "y" the
%                               overlap between the windows should be
%                               specified in percentage [default = 10]
%              Kth             -double. Tolerance for cleaning artifacts, 
%                               [default = 1.15]
%              ArtefThreshold  -double. Threshold for detection of ICs with artefacts
%                               [default = 6]
%              resample        -string. y (for yes) or n (for no). [default
%                               = "n"
%              Fs              -double. If resample == "y" the new sampling
%                               rate should be included in Hz
%              drop_channs    -string. y (for yes) or n (for no). [default
%                               = "n"
%              channels2drop  -string array. If dropp_channs == "y", an array
%                               with the labels of the channels that will
%                               be dropped should be included e.g
%                               ["Fp1","Fp2"]
%% Outputs:     
%              return_struct  -Struct with the necessary outputs (only in
%                              case of full path and not folder)
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas  
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

%% Fetch input parameters - Default values
path_data = choice_struct.path_data;
path_eeglab = choice_struct.path_eeglab;
if isfield(choice_struct,'path_save')
    path_save = choice_struct.path_save;
end
if ~isfield(choice_struct,'line_noise')
    choice_struct.line_noise = "None";
end
if choice_struct.line_noise ~= "None" && ~isfield(choice_struct,'line_noise_f')
    error('The frequency of the line noise should be given')
end
if choice_struct.line_noise == "Spectrum interpolation" && ~isfield(choice_struct,'path_fieldtrip')
    error('The path of fieldtrip is necessary to implement Spectrum Interpolation')
end
if ~isfield(choice_struct,'high_pass')
    choice_struct.high_pass = 1;
end
if ~isfield(choice_struct,'low_pass')
    choice_struct.low_pass = 100;
end
if ~isfield(choice_struct,'rereference')
    choice_struct.rereference = "Median";
end
if ~isfield(choice_struct,'clean_data')
    choice_struct.clean_data = "n";
end
if choice_struct.clean_data == "y" && ~isfield(choice_struct,'options')
    error('A struct with the options of the clean data algorithm should be included')
end
if ~isfield(choice_struct,'clean_artifacts')
    choice_struct.clean_artifacts = "y";
    choice_struct.wind = 4;
    choice_struct.overlap = 10;
    choice_struct.Kth = 1.15;
    choice_struct.ArtefThreshold = 6;
end

if ~isfield(choice_struct,'resample')
    choice_struct.resample = "n";
end

if ~isfield(choice_struct,'drop_channs')
    choice_struct.drop_channs = "n";
end

%% Find paths to necessary folders
path_clean_rawdata = [path_eeglab,'/plugins/clean_rawdata'];
if isfield(choice_struct, 'path_fieldtrip')
    path_preproc = [choice_struct.path_fieldtrip,'/preproc'];
end

init_path = pwd;
cd(init_path);
cd(path_eeglab);
eeglab;

%% Multiple EEG - serial preprocessing case
if ischar(path_data) && ~contains(path_data,'.EDF')
    dir_data = dir(path_data);
    for j = 1:length(dir_data) - 2
        %% Load EEG
        name = dir_data(j+2).name;
        EEG = pop_biosig([path_data,'/',name]);
        
        %% Line noise removal
        if choice_struct.line_noise == "Spectrum Interpolation"
            filepath = path_data;
            cd(path_data);
            filename = [filepath,'/',name(1:end-4) ,'.set'];
            if ~isfile([name(1:end-4) ,'.set'])
                pop_saveset(EEG, 'filepath',filepath, 'filename' ,[name(1:end-4),'.set']);
            end
            cd(choice_struct.path_fieldtrip);
            hdr = ft_read_header(filename);
            data = ft_read_data(filename, 'header', hdr);
            events = ft_read_event(filename, 'header', hdr);
            cd(path_preproc);
            [filt] = ft_preproc_dftfilter(data, EEG.srate, choice_struct.line_noise_f, 'dftreplace', 'neighbour');
            EEG = pop_fileio(hdr, filt, events);
        elseif choice_struct.line_noise == "Clean Line"
            EEG = pop_cleanline(EEG,choice_struct.line_noise_f);
        end
        %% Lowpass & Highpass filters
        EEG = pop_eegfiltnew(EEG, 'locutoff', choice_struct.high_pass, 'hicutoff', [], 'minphase', 1);
        EEG = pop_eegfiltnew(EEG, 'locutoff', [], 'hicutoff', choice_struct.low_pass, 'minphase', 1);
        
        %% Crop data
        if isfield(choice_struct,'begin_cropp')
            EEG.data(:,1:choice_struct.begin_cropp*EEG.srate) = [];
        end
        if isfield(choice_struct,'end_cropp')
            EEG.data(:,end-choice_struct.end_cropp*EEG.srate:end) = [];
        end
        EEG.pnts = length(EEG.data(1,:));
        %% Resample data
        if choice_struct.resample ~= "n"
            EEG = pop_resample(EEG,choice_struct.Fs);
        end
        %% Re-reference data to new reference
        if choice_struct.rereference == "Median"
            med = median(EEG.data);
            data_z = EEG.data - med;
            EEG.data = data_z;
            EEG.ref = 'Median';
        elseif choice_struct.rereference == "Average"
            EEG.data = reref(EEG.data, []);
        end
        
        %% Save channel names of all wanted electrodes 1
        chan_names_all = strings(length(EEG.data(:,1)),1);
        for k = 1:length(EEG.data(:,1))
            chan_names_all(k) = EEG.chanlocs(k).labels; % Get channel names
        end
        
        %% Apply artifact removal criteria (EEGLAB function)
        if choice_struct.clean_data == "y"
            cd(path_clean_rawdata);
            [EEG,~,~,~] = clean_artifacts(EEG,choice_struct.options);
        end
        
        %% Drop unwanted channels (specified by input parameter)
        temp = struct2table(EEG.chanlocs);
        temp = temp(1:end,1);
        if choice_struct.drop_channs == "y"
            drop = [];
            for d = 1:length(choice_struct.channels2drop)
                for t = 1:height(temp)
                    if contains(temp{t,1}{1,1},choice_struct.channels2drop(d))
                        drop = [drop,t];
                    end
                end
            end
            EEG.data(drop,:) = [];
            temp = struct2table(EEG.chanlocs);
            temp(drop,:) = [];
            EEG.chanlocs = table2struct(temp);
        end
        
        %% Save channel names of all wanted electrodes 2
        chan_names_all = strings(length(EEG.data(:,1)),1);
        for k = 1:length(EEG.data(:,1))
            for n = 1: length(dropped)
                if convertCharsToStrings(EEG.chanlocs(k).labels)~= convertCharsToStrings(dropped{n})
                    chan_names_all(k) = EEG.chanlocs(k).labels; % Get channel names
                end
            end
        end
        %Fill possible spaces in chan_names_all
        for i = 1:length(chan_names_all)
            idx = strfind(chan_names_all(i),' ');
            if ~isempty(idx)
                temp = convertStringsToChars(chan_names_all(i));
                temp(idx) = '_';
                chan_names_all(i) = convertCharsToStrings(temp);
            end
            idx = [];
        end
        
        %% Save channel names that survived
        chan_names = strings(length(EEG.data(:,1)),1);
        for k = 1:length(EEG.data(:,1))
            chan_names(k) = EEG.chanlocs(k).labels; % Get channel names
            chan_names_cell{k} = EEG.chanlocs(k).labels; % Get channel names
        end
        %% Artifact removal with wavelet-Independent Component Analysis (wICA)
        if choice_struct.clean_artifacts == "y"
            wind_len = EEG.srate*choice_struct.wind;
            step = round(wind_len*(1-choice_struct.overlap/100));
            counter = 1;
            for i = 1:step:length(EEG.data(1,:)) - wind_len 
                EEG_win = EEG.data(:,i:i+wind_len-1);
                cd(init_path);
                eeg_wICA = wICA(EEG_win, EEG.srate, chan_names, 0,choice_struct.Kth,choice_struct.ArtefThreshold); % wICA in 10s not overlapping window
                EEG_wICA = EEG;
                EEG_wICA.data = eeg_wICA;

                EEG_wICA.pnts = length(EEG_wICA.data(1,:));
                EEG_wICA.nbchan = length(EEG_wICA.data(:,2));

                EEG_wICA.data = (EEG_wICA.data - mean(EEG_wICA.data,2))./std(EEG_wICA.data,0,2);
                EEG_data{counter} = array2table(EEG_wICA.data);
                EEG_data{counter}.Properties.RowNames = chan_names_cell;
                counter = counter + 1;
            end
            %% Save results and parameters
            Fs = EEG.srate;
            num_wind = counter - 1;
            wind_dur = wind_len/EEG.srate;
            new_length = (counter-1)*choice_struct.wind/60;
            name = name(1:end-4);
            save([path_save,'/' name, '_preproc.mat'], 'EEG_data',...
                'chan_names','chan_names_all','wind_dur','new_length','name','num_wind','Fs')
            clearvars EEG_data EEG chan_names chan_names_cell EEG_wICA
            return_struct = NaN;
        
        else
            %% No wICA case
            EEG_data{1} = EEG.data;
            %% Save results and parameters
            Fs = EEG.srate;
            num_wind = 1;
            wind_dur = length(EEG.data(1,:))/EEG.srate;
            new_length = wind_dur/60;
            name = name(1:end-4);
            save([path_save,'/' name, '_preproc.mat'], 'EEG_data',...
                'chan_names','chan_names_all','wind_dur','new_length','name','num_wind','Fs')
            clearvars EEG_data EEG chan_names chan_names_cell EEG_wICA
            return_struct = NaN;
        end      
    end
%% Single EEG case
elseif ischar(path_data) && contains(path_data,'.EDF')
    %% Load EEG
    name = dir(path_data).name;
    EEG = pop_biosig(path_data);

    %% Line noise removal
    if choice_struct.line_noise == "Spectrum Interpolation"
        pathSep = strfind(path_data,filesep);
        filepath = path_data(1:pathSep(end)-1);
        cd(filepath);
        filename = [filepath,'/',name(1:end-4) ,'.set'];
        if ~isfile([name(1:end-4) ,'.set'])
            pop_saveset(EEG, 'filepath',filepath, 'filename' ,[name(1:end-4),'.set']);
        end
        cd(choice_struct.path_fieldtrip);
        hdr = ft_read_header(filename);
        data = ft_read_data(filename, 'header', hdr);
        events = ft_read_event(filename, 'header', hdr);
        cd(path_preproc);
        [filt] = ft_preproc_dftfilter(data, EEG.srate, choice_struct.line_noise_f, 'dftreplace', 'neighbour');
        EEG = pop_fileio(hdr, filt, events);
    elseif choice_struct.line_noise == "Clean Line"
        EEG = pop_cleanline(EEG,choice_struct.line_noise_f);
    end
    %% Lowpass & Highpass filters
    EEG = pop_eegfiltnew(EEG, 'locutoff', choice_struct.high_pass, 'hicutoff', [], 'minphase', 1);
    EEG = pop_eegfiltnew(EEG, 'locutoff', [], 'hicutoff', choice_struct.low_pass, 'minphase', 1);
    %% Crop data
    if isfield(choice_struct,'begin_cropp')
        EEG.data(:,1:choice_struct.begin_cropp*EEG.srate) = [];
    end
    if isfield(choice_struct,'end_cropp')
        EEG.data(:,end-choice_struct.end_cropp*EEG.srate:end) = [];
    end
    EEG.pnts = length(EEG.data(1,:));
    %% Resample data
    if choice_struct.resample ~= "n"
        EEG = pop_resample(EEG,choice_struct.Fs);
    end
    %% Re-reference data to new reference
    if choice_struct.rereference == "Median"
        med = median(EEG.data);
        data_z = EEG.data - med;
        EEG.data = data_z;
        EEG.ref = 'Median';
    elseif choice_struct.rereference == "Average"
        EEG.data = reref(EEG.data, []);
    end
    
    %% Save channel names of all wanted electrodes 1
    chan_names_all = strings(length(EEG.data(:,1)),1);
    for k = 1:length(EEG.data(:,1))
        chan_names_all(k) = EEG.chanlocs(k).labels; % Get channel names
    end
    
    %% Apply artifact removal criteria (EEGLAB function)
    if choice_struct.clean_data == "y"
        cd(path_clean_rawdata);
        [EEG,~,~,~] = clean_artifacts(EEG,choice_struct.options);
    end
   
    %% Drop unwanted channels (specified by input parameter)
    temp = struct2table(EEG.chanlocs);
    temp = temp(1:end,1);
    if choice_struct.drop_channs == "y"
        drop = [];
        for d = 1:length(choice_struct.channels2drop)
            for t = 1:height(temp)
                if contains(temp{t,1}{1,1},choice_struct.channels2drop(d))
                    drop = [drop,t];
                end
            end
        end
        EEG.data(drop,:) = [];
        temp = struct2table(EEG.chanlocs);
        dropped = temp{drop,1};
        temp(drop,:) = [];
        EEG.chanlocs = table2struct(temp);
    end
    
    %% Save channel names of all wanted electrodes 2
    chan_names_all = strings(length(EEG.data(:,1)),1);
    for k = 1:length(EEG.data(:,1))
        for n = 1: length(dropped)
            if convertCharsToStrings(EEG.chanlocs(k).labels)~= convertCharsToStrings(dropped{n})
                chan_names_all(k) = EEG.chanlocs(k).labels; % Get channel names
            end
        end
    end
    
    %Fill possible spaces in chan_names_all
    for i = 1:length(chan_names_all)
        idx = strfind(chan_names_all(i),' ');
        if ~isempty(idx)
            temp = convertStringsToChars(chan_names_all(i));
            temp(idx) = '_';
            chan_names_all(i) = convertCharsToStrings(temp);
        end
        idx = [];
    end
    
    %% Save channel names that survived
    chan_names = strings(length(EEG.data(:,1)),1);
    for k = 1:length(EEG.data(:,1))
        chan_names(k) = EEG.chanlocs(k).labels; % Get channel names
        chan_names_cell{k} = EEG.chanlocs(k).labels; % Get channel names
    end
    %% Artifact removal with wavelet-Independent Component Analysis (wICA)
    if choice_struct.clean_artifacts == "y"
        wind_len = EEG.srate*choice_struct.wind;
        step = round(wind_len*(1-choice_struct.overlap/100));
        counter = 1;
        for i = 1:step:length(EEG.data(1,:)) - wind_len 
            EEG_win = EEG.data(:,i:i+wind_len-1);
            cd(init_path);
            eeg_wICA = wICA(EEG_win, EEG.srate, chan_names, 0,choice_struct.Kth,choice_struct.ArtefThreshold); % wICA in 10s not overlapping window
            EEG_wICA = EEG;
            EEG_wICA.data = eeg_wICA;

            EEG_wICA.pnts = length(EEG_wICA.data(1,:));
            EEG_wICA.nbchan = length(EEG_wICA.data(:,2));

            EEG_wICA.data = (EEG_wICA.data - mean(EEG_wICA.data,2))./std(EEG_wICA.data,0,2);
            EEG_data{counter} = array2table(EEG_wICA.data);
            EEG_data{counter}.Properties.RowNames = chan_names_cell;
            counter = counter + 1;
        end
        %% Save results and parameters
        Fs = EEG.srate;
        num_wind = counter - 1;
        wind_dur = wind_len/EEG.srate;
        new_length = (counter-1)*choice_struct.wind/60;
        name = name(1:end-4);
        if isfield(choice_struct,'path_save')
            save([path_save,'/' name, '_preproc.mat'], 'EEG_data',...
                'chan_names','chan_names_all','wind_dur','new_length','name','num_wind','Fs')
        end
        return_struct.EEG_data = EEG_data; return_struct.chan_names = chan_names;
        return_struct.wind_dur = wind_dur; return_struct.name = name;
        return_struct.new_length = new_length; return_struct.num_wind = num_wind;
        return_struct.Fs = Fs; return_struct.chan_names_all = chan_names_all;
    else
        %% No wICA case
        EEG_data{1} = EEG.data;
        %% Save results and parameters
        Fs = EEG.srate;
        num_wind = 1;
        wind_dur = length(EEG.data(1,:))/EEG.srate;
        new_length = wind_dur/60;
        name = name(1:end-4);
        if isfield(choice_struct,'path_save')
            save([path_save,'/' name, '_preproc.mat'], 'EEG_data',...
                'chan_names','chan_names_all','wind_dur','new_length','name','num_wind','Fs')
        end
        return_struct.EEG_data = EEG_data; return_struct.chan_names = chan_names;
        return_struct.wind_dur = wind_dur; return_struct.name = name;
        return_struct.new_length = new_length; return_struct.num_wind = num_wind;
        return_struct.Fs = Fs; return_struct.chan_names_all = chan_names_all;
    end
end
        
end