%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function [decision,channel_activity_per,activity_thresh] = channel_criterion(exists,chan_name,wanted_dur,win_length,activity_thresh)
%Channel Criterion for evaluating channels of a rhythm based on activity
%% Inputs: 
%         exists - N x 1 binary table, N is number of windows, 
%                 Is 0 if current rhythm of window N of current channel has activity 
%         chan_name - string scalar, name of channel being examined
%         wanted_dur - integer scalar, in seconds, desired time length of EEG after segmenting
%         win_length - integer scalar, in seconds, EEG epoch length
%
%% Outputs: 
%         decision - string scalar, discard or keep according to criterion
%         channel_activity_per - float scalar, percentage (0-1) of active windows across whole channel
%         activity_thresh - float scalar, percentage, specified as a value between 0-1 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

num_win = height(exists);

channel_activity = exists(:,chan_name);
channel_activity = table2array(channel_activity);
channel_activity_per = sum(channel_activity)/num_win;

%If wanted duration of signal is e.g. 2 min, then threshold is set to half,
%in this case, 1 min. So if percentage of channel activity was equal to
%desired new duration e.g. for 4s windows, 120/4 = 30 windows -> 30/153 =
%0.196. This percentage corresponds to 2 minutes of EEG, so the threshold
%is set to half of this percentage to match 1 minute of EEG. 
activity_thresh = ((wanted_dur/win_length)/num_win)*activity_thresh;

if channel_activity_per < activity_thresh
    decision = "discard";
else
    decision = "keep";
end

end

