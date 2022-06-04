%Thesis Project: Classification and Characterization of the Effect of Migraine 
%through Functional Connectivity Characteristics: Application to EEG 
%Recordings from a Multimorbid Clinical Sample

function [rhythms,new_Fs] = rhythm_extraction(data,Fs,chan_names, rhythm_lims, which_rhythms, SwD_par)
%RHYTHM_EXTRACTION - Takes as input data in matrix format(channels X points)
%   then applies SwD. Then the residual is discarded from each channel, and 
%   and components are stored in struct according to which channel they
%   belong to. For each component the PSD is computed and normalised, and 
%   the power percentage of each rhythm is calculated. Then each component
%   is classified to a band and finally a struct with the rhythms of every 
%   channel is returned.
%   
%%   Inputs:  
%            data - MxN matrix (M channels, N data length)
%            Fs   - sampling rate
%            chan_names - list of strings with channel labels
%            rhythm_lims - struct array, each field contains frequency domain limits for
%                          each rhythm provided
%            which_rhythms - string array, contains names of rhythms to be
%                            returned
%
%%   Outputs: 
%            rhythms - struct of structs
%                       Each struct contains the channels rhythms
%            new_Fs  - new sampling rate
%
%   Examples of rhythm frequency limits:
%     delta = [1 4]; theta = [4 8]; alpha = [8 13]; lbeta = [13 18];
%     hbeta = [18 25]; lgamma = [25 57]; hgamma = [63 100]; 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

%% Check rhythms to save and their provided frequency limits
    finddelta = sum(contains(which_rhythms,"delta"));findtheta = sum(contains(which_rhythms,"theta"));
    findalpha = sum(contains(which_rhythms,"alpha"));findlbeta = sum(contains(which_rhythms,"lbeta"));
    findhbeta = sum(contains(which_rhythms,"hbeta"));findlgamma = sum(contains(which_rhythms,"lgamma"));
    findhgamma = sum(contains(which_rhythms,"hgamma"));
    if finddelta
        deltaLims = rhythm_lims.delta; 
    end
    if findtheta
        thetaLims = rhythm_lims.theta; 
    end
    if findalpha 
        alphaLims = rhythm_lims.alpha;
    end
    if findlbeta
        lbetaLims = rhythm_lims.lbeta; 
    end
    if findhbeta
        hbetaLims = rhythm_lims.hbeta; 
    end
    if findlgamma 
        lgammaLims = rhythm_lims.lgamma;
    end
    if findhgamma
        hgammaLims = rhythm_lims.hgamma;
    end

%% Calculate Swarm Decomposition components
%   [SwD_all,chanlocs,chan_ind,new_Fs]= our_SwD(data,Fs,chan_names,0,1); %4th arg is flag for plots, 5th is flag for filt
    [SwD_all,chanlocs,chan_ind,new_Fs] = prepare_for_SwD(data,SwD_par,chan_names,Fs,0);
    
    ind = cell(length(chan_names),1);
    for i = 1:length(chan_names)
        ind{i} = find(chan_ind==i);
    end
%% Store components in their electrodes & Discard the residual from each electrode
    elec_comp = struct();
    for i = 1:length(chan_names)
        elec_comp.(chan_names(i)) = SwD_all(ind{i},:);
    end
  
    for i = 1:length(chan_names)
        rhythms.(chan_names(i)) = struct('delta',[],'theta',[],'alpha',[],'lbeta',[],'hbeta',[],'lgamma',[],'hgamma',[]);
    end
%% Calculate FFT power spectrum and find power percentage contained in each frequency rhythm band
    L = length(SwD_all(1,:));
    nfft = 2^nextpow2(L);
    ff = new_Fs/2*linspace(0,1,nfft/2);
    step = ff(2)-ff(1); %step = Fs/nfft;
    %For each channel
    for i = 1:length(chan_names) 
        x = elec_comp.(chan_names(i));
        X = (abs(fft(x',nfft))/L).^2;
        X = X(1:end/2,:);
%       [vals,peaks] = max(X);
%       peaks = peaks*step;
%       avg = mean(X);
        X_norm = X./sum(X);
        
        x_delta = []; 
        if finddelta
            delta_per = sum(X_norm(floor(deltaLims(1)/step):floor(deltaLims(2)/step),:));
        end
        x_theta = []; 
        if findtheta
            theta_per = sum(X_norm(ceil(thetaLims(1)/step):floor(thetaLims(2)/step),:));
        end
        x_alpha = []; 
        if findalpha
            alpha_per = sum(X_norm(ceil(alphaLims(1)/step):floor(alphaLims(2)/step),:));
        end
        x_lo_beta = []; 
        if findlbeta
            lo_beta_per = sum(X_norm(ceil(lbetaLims(1)/step):floor(lbetaLims(2)/step),:));
        end
        x_hi_beta = []; 
        if findhbeta
            hi_beta_per = sum(X_norm(ceil(hbetaLims(1)/step):floor(hbetaLims(2)/step),:));
        end
        x_lo_gamma = []; 
        if findlgamma
            lo_gamma_per = sum(X_norm(ceil(lgammaLims(1)/step):floor(lgammaLims(2)/step),:));
        end
        x_hi_gamma = []; 
        if findhgamma
            hi_gamma_per = sum(X_norm(ceil(hgammaLims(1)/step):floor(hgammaLims(2)/step),:));
        end
        %For each component, sort components to a rhythm according to their power distribution
        for j = 1:length(X(1,:))
            if delta_per(j) > 0.4 && delta_per(j) > theta_per(j) && finddelta
                    x_delta(:,j) = x(j,:);
            elseif theta_per(j) > 0.4 && theta_per(j) > delta_per(j) && theta_per(j) > alpha_per(j) && findtheta
                    x_theta(:,j) = x(j,:);
            elseif alpha_per(j) > 0.4 && alpha_per(j) > theta_per(j) && alpha_per(j) > lo_beta_per(j) && findalpha
                    x_alpha(:,j) = x(j,:);
            elseif lo_beta_per(j) > 0.4 && lo_beta_per(j) > hi_beta_per(j) && lo_beta_per(j) > alpha_per(j) && findlbeta
                    x_lo_beta(:,j) = x(j,:);
            elseif hi_beta_per(j) > 0.4 && hi_beta_per(j) > lo_beta_per(j) && hi_beta_per(j) > lo_gamma_per(j) && findhbeta
                    x_hi_beta(:,j) = x(j,:);
            elseif lo_gamma_per(j) > 0.3 && lo_gamma_per(j) > hi_beta_per(j) && findlgamma 
                    x_lo_gamma(:,j) = x(j,:);
            elseif hi_gamma_per(j) > 0.3 && findhgamma
                    x_hi_gamma(:,j) = x(j,:);        
            else
                    warning("Rhythm was found, but it doesn't belong in any band")
                    fprintf("delta: %f \n, theta: %f \n, alpha: %f \n, low beta: %f \n",...
                        delta_per(j),theta_per(j),alpha_per(j),lo_beta_per(j))
                    fprintf("high beta: %f \n, low gamma: %f \n,high gamma: %f \n",hi_beta_per(j),...
                        lo_gamma_per(j),hi_gamma_per(j))
            end      
        end
 
        if ~isempty(x_delta)
            x_delta(:,all(x_delta == 0)) = [];
        end
        if ~isempty(x_theta)
            x_theta(:,all(x_theta == 0)) = [];
        end
        if ~isempty(x_alpha)
            x_alpha(:,all(x_alpha == 0)) = [];
        end
        if ~isempty(x_lo_beta)
            x_lo_beta(:,all(x_lo_beta == 0)) = [];
        end
        if ~isempty(x_hi_beta)
            x_hi_beta(:,all(x_hi_beta == 0)) = [];
        end
        if ~isempty(x_lo_gamma)
            x_lo_gamma(:,all(x_lo_gamma == 0)) = [];
        end
        if ~isempty(x_hi_gamma)
            x_hi_gamma(:,all(x_hi_gamma == 0)) = [];
        end
        rhythms.(chan_names(i)).delta = x_delta;
        rhythms.(chan_names(i)).theta = x_theta;
        rhythms.(chan_names(i)).alpha = x_alpha;
        rhythms.(chan_names(i)).lbeta = x_lo_beta;
        rhythms.(chan_names(i)).hbeta = x_hi_beta;
        rhythms.(chan_names(i)).lgamma = x_lo_gamma;
        rhythms.(chan_names(i)).hgamma = x_hi_gamma;

    end
end

