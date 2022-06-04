%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE

function [varargout] = calc_swarm_feature(func_name,varargin)
%%Calls function of a swarm feature and returns flattened outputs
% Inputs: func_name - string, name of function
%         varargin  - arguments required by the function called
%
% Commonly used arguments: segmented_eeg - struct of segmented eeg data
%                          chans         - array of desired channels or 'all'
%                          Fs
%                          bandwidth - BW of the new, non - conventional bands
%                          overlap - overlap of the new, non - conventional bands
%                          N - number of surrogates for coherence and plv calc
%                          thres - threshold of minimum psd percentage to determine if an
%                                   SwD belongs to a band.
%                          new_length - new eeg time length after segmenting ( in minutes )
%                          wind_dur - duration of windows
%
% Outputs: varargout - Cell array, contains outputs of called function 
    
    if length(varargin) == 1
        varargin = varargin{1,1};
    end
    %Find how many output arguments func has
    nout = nargout(func_name);
    %Create cell to store outputs of function
    results = cell(1,nout);
    %Call function and get outputs
    [results{:}] = feval(func_name,varargin{:});


% mean_mean_feature is a B x P matrix, where B is number of bands,and P is 
% number of electrode pairs or number of electrodes. To store it, it needs 
% to be flattened to a vector 1 x B*P, where 1st element is the 1st 
% electrode pair of the 1st band,2nd element is 2nd electrode pair of the 
% 1st band, etc...
 
%% Only for running swd_features
%     if strcmp(func_name,"power_spectrum_features") == 0
%         mean_mean_feature = results{1};
%         std_mean_feature = results{2};
%         [size_x,size_y] = size(mean_mean_feature);
% 
%     %Flattening
%         flat_mmf = reshape(mean_mean_feature',[1,size_x*size_y]);
%         flat_smf = reshape(std_mean_feature',[1,size_x*size_y]);
%         results{1} = flat_mmf;
%         results{2} = flat_smf;
%     end
    
    varargout = results;
end