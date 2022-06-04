%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function [window_flag,window_regions] = window_criterion(window,chan_names,regional_activity_thresh,window_activity_thresh)
%Window Criterion for evaluating windows of a rhythm based on activity

%% Inputs: 
%          window - 1 x C binary table, C is number of channels, 
%                 Is 0 if current rhythm of channel N of current window has activity 
%          chan_names - string array, names of channels being examined
%          regional_activity_thresh - 1x5 vector, contains activity thresholds of regions
%                   For a region to be considered as active, a rule of thumb is considered: 
%                   half the number of electrodes in a region should be active. e.g. for
%                   10-20 system with 19 electrodes, Frontal region has 7 electrodes.
%                   A good choice would be 3/7.
%
%          window_activity_thresh - scalar, 0-1, Threshold for discarding a window
%                                   Must be equal to or lower than 0.25
%
%% Outputs: 
%          window_flag - Flag that characterizes window activity
%          window_regions - string array, names of regions being active in
%                           current window
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------
    regions = ["F","P","C","O","T"];
    regional_elecs_ideal = [7,3,3,2,4];
    regional_activity = NaN(1,length(regions));
    regional_index = NaN(1,length(regions));
    for r = 1:length(regions)
        region_chans = [];
        for m = 1:length(chan_names)
            if contains(chan_names(m),regions(r))
               region_chans = [region_chans chan_names(m)]; 
            end
        end
        activity_subtable = window(1,region_chans);
        activity_subtable = table2array(activity_subtable);
        %regional_elecs = length(activity_subtable);
        regional_activity(r) = sum(activity_subtable)/regional_elecs_ideal(r);
        
        regional_activity_thresh(r) = regional_activity_thresh(r)/regional_elecs_ideal(r); 
        
        if regional_activity(r) < regional_activity_thresh(r)
            regional_index(r) = 0;
        else
            regional_index(r) = 1;
        end
    end

    if istable(window)
        window = table2array(window);
    end

    num_ch = length(window);
    window_activity_per = sum(window)/num_ch;
    if window_activity_per < window_activity_thresh && sum(regional_index) < 1
        window_flag = "discard"; %Discard
        window_regions = "None";
    elseif window_activity_per < window_activity_thresh && sum(regional_index) >= 1
        window_flag = "local"; %Keep, low / local activity
        window_regions = regions(find(regional_index == 1));
    elseif window_activity_per < 2*window_activity_thresh && ...
            window_activity_per >= window_activity_thresh
        window_flag = "medium"; %Keep, medium activity
        window_regions = regions(find(regional_index == 1));
    elseif window_activity_per < 3*window_activity_thresh && ...
            window_activity_per >= 2*window_activity_thresh 
        window_flag = "high"; %Keep, high activity
        window_regions = regions(find(regional_index == 1));
    elseif window_activity_per >= 3*window_activity_thresh 
        window_flag = "vhigh"; %Keep, very high activity
        window_regions = regions(find(regional_index == 1));
    end
    
end

