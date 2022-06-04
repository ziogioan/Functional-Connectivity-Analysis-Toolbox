%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function [tf_PAC_all,tf_PAC] = tf_PAC_surr(param_struct)
% This function calculates the tf_PAC comodulogram, given two signals, based
% on the work of Tamanna T. K. Munia & SelinAviyente (https://doi.org/10.1038/s41598-019-48870-2).
% A big part of the script was taken by the original scripts provided by
% the authors.
% Three different PAC measures are available, Mean Vector Length (MVL),
% Modulation Index (MI) and General Linear Model (GLM). Four different
% surrogate methods are available. 
%
%% Inputs: 
% param_struct   -struct. Contains the necessary parameters
%
%% Outputs: 
% tf_PAC_all     -double matrix. Comodulogram matrix
% tf_PAC         -double. Max value of comodulogram matrix

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


%% Fetch necessary parameters from param_struct

[tf_PAC_all,~] = tfPAC_surrogate(param_struct);
tf_PAC = abs(max(max(tf_PAC_all))); % Computed tf-MVL value
% Surr_tf = abs(mean2(tfsurr_PAC)); 
% [high_in, low_in] = find((tf_PAC_all==tf_PAC));
% 
% high_pacf = highfreq(high_in); % Detected amplitude providing Frequency
% low_pacf = lowfreq(low_in); % Detected phase providing Frequency
% pacfreq = [low_pacf, high_pacf];
% plot_comodulogram(tf_MVL_all,high,low) %plot comodulogram
 
 function [tf_canolty,PAC_surr] = tfPAC_surrogate(param_struct)
% This function computes the phase amplitude coupling using TF-MVL method.
% Choice over three different surrogate data methods is included.
% NoSurrogate, Constant Thresholding, Block-Swapping and Block-Resampling
% methods are available.
% Two signals enter the pipeline through 'param_struct'. A Min-Max Scaling
% option is given, to minimize the effect of amplitude on the resulting PAC
% value. The Rid-Rihadzek distribution is used to extract Amplitude and Phase signals.
% Surrogate threshold is computed and then real PAC value is calculated
% through calc_pac.m function.
%
% Input:   x            : input signal 
%          high_freq    : Amplitude Frequency range 
%          low_freq     : Phase Frequency range 
%          Fs           : Sampling frequency
% Output:  tf_canolty   : Computed PAC using TF-MVL method
%          PAC_surr     : mean surrogate values using one of various
%                         surrogate methods
%

% Written by: Tamanna T. K. Munia, January 2019

% These scripts have been optimised for the Windows operating systm  
% MATLAB version used 2018a.

%% Modified by: Charalampos Lamprou && Ioannis Ziogas, December 2021

low = param_struct.low;
high = param_struct.high;
step_low = param_struct.step_low;
step_high = param_struct.step_high;


%% Compute TF_MVL 
highfreq = high(1):step_high:high(2); %Divide high frequency vector into small intervals
amp_length = length(highfreq);
lowfreq = low(1):step_low:low(2); %Divide low frequency vector into small intervals
phase_length = length(lowfreq);
tf_canolty = zeros(amp_length-1,phase_length-1);
tfsurr_PAC = zeros(amp_length-1,phase_length-1);
%   thsur = zeros(amp_length,phase_length);


%% Amplitude and Phase calculation
x_low = param_struct.x_low;
x_high = param_struct.x_high;
if param_struct.normalise
    %Min-Max Scaling for Amplitude
    x_high = (x_high - min(x_high))./ (max(x_high) - min(x_high));
end
Fs = param_struct.Fs;
method = param_struct.surr_method;
[tfr] = rid_rihaczek4(x_low+x_high,Fs);
W = tfr;
W2 = W(2:end,:);
%Iterate through frequency blocks in the comodulogram plane
for i = 1:phase_length
  for j = 1:amp_length
        l_freq = lowfreq(i);
        h_freq = highfreq(j);
        Amp = abs(mean(W2(h_freq:h_freq,:),1));
        tfd_low = mean(W2(l_freq:l_freq,:),1);
        angle_low = angle(tfd_low);
        Phase = (angle_low);
        param_struct.Phase = Phase;
        param_struct.Amp = Amp;
        tf_canolty_r = calc_pac(param_struct); %A single PAC value is returned here
        if method == "block-resampling"
            [PAC_surr] = PAC_surrogate_block_resampling(param_struct);
        elseif method == "permutation"
            [PAC_surr] = PAC_surrogate_block_swapping(param_struct);
        elseif method == "noSurrogate"
            [PAC_surr] = calc_pac(param_struct);
        elseif method == "Constant"
            PAC_surr = param_struct.per_thresh;
        end

        tf_canolty_r(tf_canolty_r<PAC_surr) = 0; 
        tf_canolty(j,i) = tf_canolty_r;        
  end
end
end
function [PAC_surr] = PAC_surrogate_block_resampling(param_struct)

% This function generates the surrogate data from given phase and amplitude using block method.
% Input:   Amp: detected amplitude envelope
%          Phase: detected phase envelope
% Output:  MVL_surr: mean surrogate values using block resampling approach

%% ***Modified to contain all inputs in 'param_struct'

% Written by: Tamanna T. K. Munia, January 2019

% These scripts have been optimised for the Windows operating systm  
% MATLAB version used 2018a.

%% Surrogate data calculation
 Amp = param_struct.Amp;
 Phase = param_struct.Phase;
 block_size = param_struct.block_size;
 N = param_struct.N;
 L_a = length(Amp) - mod(length(Amp),block_size); % If Amp cannot be perfectly divided into N blocks, make it  
 Amp_block = reshape(Amp(1:L_a), block_size, []); % divide the amp into blocks
 L_p = length(Phase) - mod(length(Phase),block_size); % If Phase cannot be perfectly divided into N blocks, make it
 ph_block = reshape(Phase(1:L_p), block_size, []); % divide the phase into blocks

 for surr = 1:N  
     random_Phblock = randperm(size(ph_block,1),1); %randomly select phase block number
     ph_surr = ph_block(random_Phblock,:); %extract phase signal for random block
     random_Ampblock = randperm(size(Amp_block,1),1); %randomly select Amp block number
     Amp_ran = Amp_block(random_Ampblock,:); %extract Amp signal for random block
     Amp_sur = Amp_ran(randperm(length(Amp_ran))); % Shuffle the random Amp block only
     param_struct.Amp = Amp_sur;
     param_struct.Phase = ph_surr;
     [M_sur] = calc_pac(param_struct); %Calc PAC for Phase block and shuffled Amp block
     MVL_surra(surr) = M_sur;
 end

  PAC_surr = mean(MVL_surra); %Surrogate Threshold value is the Mean Value of N surrogate PAC values
 
end

function [PAC_surr] = PAC_surrogate_block_swapping(param_struct)
%% ***Modified to contain all inputs in 'param_struct'

Amp = param_struct.Amp;
Phase = param_struct.Phase;
numsurrogate = param_struct.N;
srate = param_struct.Fs;
phase = Phase;
amplitude = Amp;
numpoints = length(Amp);   %% number of sample points in raw signal 
if numpoints ~= length(phase)
    error('Phase signal and Amplitude signal must have the same length')
end
minskip = srate;   %% time lag must be at least this big 
maxskip = numpoints-srate; %% time lag must be smaller than this 
skip = ceil(numpoints.*rand(numsurrogate*2,1)); % Generate N random breakpoints 
skip(find(skip>maxskip)) = []; % Keep only breakpoints that fall into the minskip-maxskip interval
skip(find(skip<minskip)) = []; 
mn = length(skip);
skip = skip(1:mn,1); 
surrogate_m = zeros(mn,1);
for s = 1:mn
    %Reverse amplitude vector, around breakpoint s - Phase remains intact
    %Calculate PAC for this pair
    surrogate_amplitude = [amplitude(skip(s):end) amplitude(1:skip(s)-1)];
    param_struct.Amp = surrogate_amplitude;
    param_struct.Phase = phase;
    surrogate_m(s) = calc_pac(param_struct);        
end 
PAC_surr = mean(surrogate_m); %Surrogate Threshold is the mean of all surrogate PAC values
  
end
end

