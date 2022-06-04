%% Thesis - Ioannis Ziogas 9132 & Charalampos Lamprou 9114- AUTh ECE

function [hil_PAC_all,tf_PAC] = hil_PAC_surr(param_struct)
% This function employs existing TF-MVL method for PAC calculation 
% plus different surrogate methods to determine computed PAC statistical significance.
% A normalisation option in the form of Min-Max Scaling is given. 
% Three different PAC measures are available: Mean Vector Length (MVL), 
% Modulation Index (MI) and General Linear Model (GLM).

%% Written by: Ioannis Ziogas && Charalampos Lamprou, December 2021

[hil_PAC_all,hil_surr_PAC] = hilPAC_surrogate(param_struct);

tf_PAC = abs(max(max(hil_PAC_all))); % Computed tf-MVL value

 
 function [PAC_true,PAC_surr] = hilPAC_surrogate(param_struct)
% This function computes the phase amplitude coupling.
% Choice over three different PAC measures and three different surrogate
% data methods is included. MVL, MI and GLM measures are available. 
% NoSurrogate, Constant Thresholding, Block-Swapping and Block-Resampling
% methods are available.
% Two signals enter the pipeline through 'param_struct'. A Min-Max Scaling
% option is given, to minimize the effect of amplitude on the resulting PAC
% value. The envelope of the Amplitude - High Frequency signal is
% extracted. To ensure that the envelope extracted is close to ideal, the
% frequency corresponding to the maximum spectral concentration is found
% through Welch's method, and envelope is extracted according to this
% frequency. Phase is extracted through Hilbert Transform of the Phase -
% Low Frequency signal. Various surrogate methods are provided. PAC
% surrogate threshold is computed and then real PAC value is calculated
% through calc_pac.m function.
% Input:   param_struct - contains:
%                                  x            : input signal 
%                                  high_freq    : Amplitude Frequency range 
%                                  low_freq     : Phase Frequency range 
%                                  Fs           : Sampling frequency
%
% Output:  tf_canolty   : Computed PAC using TF-MVL method
%          PAC_surr     : mean surrogate values using block-swapping method


%% Written by: Ioannis Ziogas && Charalampos Lamprou, December 2021

%% Amplitude and Phase calculation
x_low = param_struct.x_low;
x_high = param_struct.x_high;
if param_struct.normalise
    %Min-Max Scaling for Amplitude
    x_high = (x_high - min(x_high))./ (max(x_high) - min(x_high));
    
end
Fs = param_struct.Fs;
method = param_struct.st_sig_method;
[p,f] = pwelch(x_low,length(x_low)/4,50,0:0.5:15,Fs);
[~,m] = max(p);
peaks = f(m);

% if round(peaks) == 0
%     Amp = envelope(x_high,1,'peak'); %Calculate for 1 Hz
% else
%     Amp = envelope(x_high,round(peaks),'peak');
% end

Amp = envelope(x_high);

Amp = ensurePositiveAmp(Amp);

Phase = angle(hilbert(x_low));

% figure
% subplot(3,1,1);
% plot(x_low,'LineWidth',1.2)
% xlabel('????????')
% ylabel('???????????????? ??????')
% title('x_f_p(k): ???? "?????" ??????? ?????????? ????? f_p')
% subplot(3,1,3);
% plot(x_high); hold on; plot(Amp,'LineWidth',1.2)
% xlabel('????????')
% ylabel('???????????????? ??????')
% title('x_f_a(k) ??? A_f_a(k): ???? "???????" ?????? ?????????? ????? f_a ??? ??????? ???????')
% % legend(['???? "???????"'],['??????? ???????'])
% subplot(3,1,2);
% plot(Phase,'LineWidth',1.2)
% ylabel('????? (???????)')
% xlabel('????????')
% title('?_f_p(k): ???? ??????? ?????????? ?????')

param_struct.Phase = Phase;
param_struct.Amp = Amp;
PAC_true = calc_pac(param_struct);
if method == "block-resampling"
    [PAC_surr] = PAC_surrogate_block_resampling(param_struct);
elseif method == "block-swapping"
    [PAC_surr] = PAC_surrogate_block_swapping(param_struct);
elseif method == "noSurrogate"
    [PAC_surr] = calc_pac(param_struct);
elseif method == "Constant"
    PAC_surr = param_struct.per_thresh;
end

PAC_true(PAC_true<PAC_surr) = 0; 

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
numpoints = length(Amp); %% number of sample points in raw signal 
if numpoints ~= length(phase)
    error('Phase signal and Amplitude signal must have the same length')
end
minskip = srate; %% time lag must be at least this big 
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
    surrogate_amplitude = [amplitude(skip(s):end)' amplitude(1:skip(s)-1)'];
    param_struct.Amp = surrogate_amplitude;
    param_struct.Phase = phase;
    surrogate_m(s) = calc_pac(param_struct);        
end 
PAC_surr = mean(surrogate_m); %Surrogate Threshold is the mean of all surrogate PAC values
  
end

%--------------------------------------------------------------------------

function Amp = ensurePositiveAmp(Amp)
    %In case that envelope contains negative values:
    negEnv = [1];
    counter = 2;
    while ~isempty(negEnv) 
        negEnv = find(Amp < 0); %Find negative values
        serial = diff(negEnv)';
        try
            cross = [1 find(serial > 1) length(serial)];
        catch
            cross = [1 [find(serial > 1)]' length(serial)];
        end
        maxlen = 0;
        for cr = 1:length(cross)-1
            sublen = length(serial(cross(cr):cross(cr+1)-1));
            if sublen > maxlen
                maxlen = sublen;
            end
        end
    %If they are more than 10 in a row, find another frequency peak to extract envelope
        if maxlen > 10  
            [~,I] = sort(p,'descend');
            if counter <= 3 %Search until 3rd higher peak - First has been already examined above             
                I = I(counter);
                counter = counter + 1;
                pk = round(f(I));
            else
            %If other peaks also fail, then go again to higher peak and search
            %around the peak +-2 frequency bins
                range = [2,-2,4,-4];
                len = 1001;
                for r = 1:length(range)
                    if I(1) + range(r) < length(f) && I(1) + range(r) > 0
                        pk = round(f(round(I(1) + range(r))));
                        Amp = envelope(x_high,pk,'peak');
                        negEnv2 = find(Amp < 0);
                        if length(negEnv2) < len
                            pos = range(r);
                            len = length(negEnv2);
                            negEnv = negEnv2;
                        end
                    end
                end
                pk = round(f(round(I(1) + pos)));
            end
            Amp = envelope(x_high,pk,'peak');
        else
            %If only a few values are under 0, just set those to 0
            Amp(negEnv) = 0;
            negEnv = find(Amp < 0);
        end
        if ~isempty(negEnv) 
            disp(-1)
        end
    end
end
end

