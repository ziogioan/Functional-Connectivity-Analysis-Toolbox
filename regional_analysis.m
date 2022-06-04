%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE

function [regional_feat] = regional_analysis(func_name,param_struct)

% Calculates regional connectivity measures among several region choices: a)
% inter-regional coherence, b) intra-regional coherence, c)left-right
% hemisphere coherence, d) anti-symmetric coherence (only anti-lateral
% coherence)
% 
% Inputs: func_name     - name of the wanted feature function
%                       (available atm are:"coher_swarm_calc",
%                       "plv_swarm_calc","bspec_swarm_calc","cross_bspec_calc","pac_swarm_calc"
%         param_struct  - struct that contains all necessary input arguments.
%               These can be:
%                       regional  - string specifying regional analysis type:
%                           "inter-regional","intra-regional","left-right"
%                       sumup     - binary flag, to sum or not to sum electrodes in inter and left-right
%                       rhythms_struct   - struct of eeg data
%                       chans     - string specifying desired channels:
%                                        can be "all" or specific channels
%                       Fs        - sampling rate
%                       bands     - vector containing start and end of 
%                                           non-conventional bands
%                       new_length - duration of eeg in min
%                       wind_dur  - duration of windows in seconds
%                       /////////////////////////////////////////
%                       For surrogate calculation, depending on feature
%                       measured: 
%                                 1. Coherence
%                                  window - window for mscohere calculation
%                                  overlap - overlap for mscohere
%                                  method - surrogate data threshold calculation method,
%                                           available choices are:
%                                           "Koopmans","permutation","constant"
%                                  alpha  - 1 - confidence for "Koopmans" method
%                                 2. Phase Locking Value
%                                  method - surrogate threshold calculation
%                                       method, choose from: "block-resampling"
%                                       and "permutation"
%                                  alpha  - 1 - confidence 
%                                  N      - number of surrogate samples
%                                  blocksize - argument for
%                                  "block-resampling method"
%                                 3. Bispectrum & Bicoherence
%
%                                 4. Anti-Symmetric Cross Bispectrum OR
%                                 Cross Bispectrum and Cross Bicoherence
%
%                                 5. Phase-Amplitude Coupling

% Outputs: regional_feat          - structure containing calculated feature
%                                   for all regional combinations for each window 
%
%% Written by: Ioannis Ziogas && Charalampos Lamprou, August 2021

%% Initializations
chan_names = param_struct.chan_names;
chans = param_struct.chans;
if chans == "all"
    chans = param_struct.chan_names_all;
%     chans = ["Fp1", "Fp2", "F7", "F3","Fz", "F4", "F8", "T3","C3","Cz","C4","T4","T5","P3","Pz","P4","T6","O1","O2"];
end
wind_dur = param_struct.wind_dur; Fs = param_struct.Fs; num_wind = param_struct.num_wind;
wind_length = wind_dur*Fs; regional = param_struct.regional; bands_length = param_struct.bands_length;

%% Load signal 
data = param_struct.rhythms;
tempfnames = fieldnames(data(1));
rhythm_names = fieldnames(data(1).(tempfnames{1})); param_struct.rhythm_names = rhythm_names;

switch regional
    case "inter-regional"
        reg = ["F","P","C","T","O"];
        for j = 1:num_wind   
            for r = 1:length(reg)     
                count = 1; % counter for the electrode pairs (1:171)
                for k = 1:length(chans)
                    if contains(chans(k),reg(r))
                        flag_k = elec_exists(chan_names, chans(k)); %1 if electrode k exists
                        if flag_k % if electrode k exists in region r   
                            sigrhy = struct2cell(data(j).(chans(k)));
                            if ~isempty(sigrhy)
                                reg_table.(chans(k)) = sigrhy;  
                                reg_chan_names(count) = chans(k);    
                                count = count + 1;                       
                            end
                        end
                    end
                end
                if exist('reg_table','var')
                    regional_eeg(j).(reg(r)).data = reg_table;
                    regional_eeg(j).(reg(r)).name_info = reg_chan_names;
                    clearvars reg_table reg_chan_names
                else
                    temp{1,1} = NaN(wind_length,1);
                    regional_eeg(j).(reg(r)).data.NoChannels = temp;
                end
            end
            counter = 1;
            for g = 1:r
                for m = counter+1:r
                    region_a = regional_eeg(j).(reg(g));
                    region_b = regional_eeg(j).(reg(m));
                    param_struct.a = region_a.data;
                    param_struct.b = region_b.data;
                    [conn_feat_mean,new_names] = feval(func_name,param_struct);              
                    if ~isempty(new_names)
                        names = new_names;
                    end
                    if contains(func_name,'pac') || contains(func_name,'cross_bspec')
                        comb_name = [convertStringsToChars(reg(g)),'_',...
                            convertStringsToChars(reg(m))];
                        regional_feat.(comb_name)(j,:) = conn_feat_mean(1:end/2);
                        comb_name = [convertStringsToChars(reg(m)),'_',...
                            convertStringsToChars(reg(g))];
                        regional_feat.(comb_name)(j,:) = conn_feat_mean(end/2+1:end);
                    else                     
                        comb_name = [convertStringsToChars(reg(g)),'_',...
                            convertStringsToChars(reg(m))];
                        regional_feat.(comb_name)(j,:) = conn_feat_mean;  
                    end
                end
                counter = counter + 1; 
                if counter == r
                    counter = 1000; 
                end
            end
        end
    case "intra-regional"
        reg = ["F","P","C","T","O"];
        for j = 1:num_wind   
            for r = 1:length(reg)     
                count = 1; % counter for the electrode pairs (1:171)
                for k = 1:length(chans)
                    if contains(chans(k),reg(r))
                        flag_k = elec_exists(chan_names, chans(k)); %1 if electrode k exists
                        if flag_k % if a pair of electrodes exists calculate mean & std of coherence  
                            sigrhy = struct2cell(data(j).(chans(k)));
                            if ~isempty(sigrhy)
                                reg_table.(chans(k)) = sigrhy;  
                                reg_chan_names(count) = chans(k);    
                                count = count + 1;                       
                            end
                        end
                    end
                end
                if exist('reg_table','var')
                    param_struct.a = reg_table;
                    param_struct.b = reg_table;
                    [conn_feat_mean,new_names] = feval(func_name,param_struct);
                    if ~isempty(new_names)
                        names = new_names;
                    end
                    regional_feat.(reg(r))(j,:) = conn_feat_mean;   
                    clearvars reg_table reg_chan_names
                else 
                    % This if statement is executed only for the intra
                    % regional case: Because here region electrodes here are
                    % gathered and instantly passed to swarm_calc functions,
                    % a "NoChannels" electrode is not necessary here. In
                    % other regional cases, the regional electrodes are
                    % first gathered into respective subregions, and then a
                    % new for loop is executed, so "NoChannels" electrodes 
                    % are necessary 
                    if contains(func_name,"bspec") || contains(func_name,"pac")
                        
                        regional_feat.(reg(r))(j,:) = NaN(1,param_struct.feat_length*bands_length);
                    else
                        regional_feat.(reg(r))(j,:) = NaN(1,bands_length);
                    end
                end
            end 
        end 
    case "left-right"
        reg = ["L1","L2","L3","R1","R2","R3"];
        reg1 = ["Fp1F7F3","T3C3","T5P3O1","Fp2F8F4","T4C4","T6P4O2"];
        for j = 1:num_wind
            for r = 1:length(reg)     
                count = 1; % counter for the electrode pairs (1:171)
                for k = 1:length(chans)
                    if contains(reg1(r),chans(k))
                        flag_k = elec_exists(chan_names, chans(k)); %1 if electrode k exists
                        if flag_k % if a pair of electrodes exists calculate mean & std of coherence  
                            sigrhy = struct2cell(data(j).(chans(k)));
                            if ~isempty(sigrhy)
                                reg_table.(chans(k)) = sigrhy;  
                                reg_chan_names(count) = chans(k);    
                                count = count + 1;                       
                            end
                        end
                    end
                end
                if exist('reg_table','var')
                    regional_eeg(j).(reg(r)).data = reg_table;
                    regional_eeg(j).(reg(r)).name_info = reg_chan_names;
                    clearvars reg_table reg_chan_names
                else
                    temp{1,1} = NaN(wind_length,1);
                    regional_eeg(j).(reg(r)).data.NoChannels = temp;
                end
            end

            for g = 1:length(fieldnames(regional_eeg))/2
                for m = length(fieldnames(regional_eeg))/2+1:length(fieldnames(regional_eeg))
                    region_a = regional_eeg(j).(reg(g));
                    region_b = regional_eeg(j).(reg(m));
                    param_struct.a = region_a.data;
                    param_struct.b = region_b.data;
                    [conn_feat_mean,new_names] = feval(func_name,param_struct);               
                    if ~isempty(new_names)
                        names = new_names;
                    end
                    if contains(func_name,'pac') || contains(func_name,'cross_bspec')
                        %PAC names are handled differently from other
                        %features, because e.g. the F - P regional pair may
                        %have different PAC value from the P - F pair 
                        comb_name = [convertStringsToChars(reg(g)),'_',...
                            convertStringsToChars(reg(m))];
                        regional_feat.(comb_name)(j,:) = conn_feat_mean(1:end/2);
                        comb_name = [convertStringsToChars(reg(m)),'_',...
                            convertStringsToChars(reg(g))];
                        regional_feat.(comb_name)(j,:) = conn_feat_mean(end/2+1:end);
                    else                     
                        comb_name = [convertStringsToChars(reg(g)),'_',...
                            convertStringsToChars(reg(m))];
                        regional_feat.(comb_name)(j,:) = conn_feat_mean;  
                    end     
                end
            end
        end
    case "anti-symmetric"
        reg = ["UL","CL","DL","OL","IL","TL","UR","CR","DR","OR","IR","TR"];
        reg1 = ["Fp1F7F3T3C3","T5P3F7T3C3F3","T5P3O1T3C3","F7T3T5","F3C3P3","Fp1F7F3T3C3T5P3O1"...
            ,"Fp2F8F4T4C4","T6P4F8T4C4F4","T6P4O2T4C4","F8T4T6","F4C4P4","Fp2F8F4T4C4T6P4O2"];
%         reg = ["UL","DL","UR","DR"];
%         reg1 = ["Fp1F7F3T3C3","T5P3O1T3C3","Fp2F8F4T4C4","T6P4O2T4C4"];
        for j = 1:num_wind
            for r = 1:length(reg)     
                count = 1; % counter for the electrode pairs (1:171)
                for k = 1:length(chans)
                    if contains(reg1(r),chans(k))
                        flag_k = elec_exists(chan_names, chans(k)); %1 if electrode k exists
                        if flag_k % if a pair of electrodes exists calculate mean & std of coherence  
                            sigrhy = struct2cell(data(j).(chans(k)));
                            if ~isempty(sigrhy)
                                reg_table.(chans(k)) = sigrhy;  
                                reg_chan_names(count) = chans(k);    
                                count = count + 1;                       
                            end
                        end
                    end
                end
                if exist('reg_table','var')
                    regional_eeg(j).(reg(r)).data = reg_table;
                    regional_eeg(j).(reg(r)).name_info = reg_chan_names;
                    clearvars reg_table reg_chan_names
                else
                    temp{1,1} = NaN(wind_length,1);
                    regional_eeg(j).(reg(r)).data.NoChannels = temp;
                end
            end

            for g = 1:length(fieldnames(regional_eeg))/2
                region_a = regional_eeg(j).(reg(g));
                region_b = regional_eeg(j).(reg(g+length(fieldnames(regional_eeg))/2));
                param_struct.a = region_a.data;
                param_struct.b = region_a.data;
                [conn_feat_mean_a,~] = feval(func_name,param_struct); 
                param_struct.a = region_b.data;
                param_struct.b = region_b.data;
                [conn_feat_mean_b,new_names] = feval(func_name,param_struct); 
                if ~isempty(new_names)
                        names = new_names;
                end
                comb_name_minus = [convertStringsToChars(reg(g)),'_',...
                    convertStringsToChars(reg(g+length(fieldnames(regional_eeg))/2)),'_minus'];
%                 comb_name_div = [convertStringsToChars(reg(g)),'_',...
%                     convertStringsToChars(reg(g+length(fieldnames(regional_eeg))/2)),'_div'];
%                 comb_name = [convertStringsToChars(reg(g)),'-',...
%                     convertStringsToChars(reg(g+length(fieldnames(regional_eeg))/2))];
%                 regional_feat.(comb_name)(j,:) = ...
%                     abs((conn_feat_mean_a - conn_feat_mean_b));
                regional_feat.(comb_name_minus)(j,:) = conn_feat_mean_a - conn_feat_mean_b;
                %regional_feat.(comb_name_div)(j,:) = conn_feat_mean_a./conn_feat_mean_b;
            end
        end
end
regional_feat = structfun(@mean,regional_feat,'UniformOutput',false);
regional_feat.names = names;
end