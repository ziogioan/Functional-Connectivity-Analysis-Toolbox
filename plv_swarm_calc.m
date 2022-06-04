function [mean_PLV, names] = plv_swarm_calc(param_struct)

% Calculates the phase locking value between SwDs of different electrodes
% when both of them belong to the same rhythm.
% Inputs: elec_case    - string specifying method used: can be "single" for
%                        electrode to electrode coherence, or any of the regional options
%                        specified in regional_analysis.m for region to region coherence
%         x            - data of the one electrode (struct)
%         y            - data of the other electrode (struct)
%         surr_struct  - struct containing  parameters for
%                        surrogate method used for determining significance threshold in
%                        coherence calculation: method and method parameters
%         For more details, see regional_analysis.m

% Outputs: mean_plv - mean of plv at each band and pair of electrodes
%                           across windows
%          std_plv  - std of plv at each band and pair of electrodes
%                           across windows

plvThresh = param_struct.plvThresh;

%------------------- Check if input is in correct form --------------------
x = param_struct.a;
y = param_struct.b;

if ~isstruct(x) || ~isstruct(y)
    error("rhythms should be in struct form")
end
%x-region: struct, each field is an electrode
electrodesx = fieldnames(x);
electrodesy = fieldnames(y);

elec_case = param_struct.regional;
if ~strcmp(electrodesx{1,1}, "NoChannels") && ~strcmp(electrodesy{1,1}, "NoChannels")
    for ex = 1:length(electrodesx)
        xx = x.(electrodesx{ex}); %xx-electrode: struct, each field is a rhythm
        xx = cell2struct(xx, param_struct.rhythm_names);
        rhythms = param_struct.conv_bands;
        for ey = 1:length(electrodesy)
            if elec_case == "intra-regional" || elec_case == "anti-symmetric" || elec_case == "single"
                ifcond = (ey > ex);
            elseif elec_case == "inter-regional" || elec_case == "left-right"
                ifcond = 1;
            end
                    
            if ifcond               
                yy = y.(electrodesy{ey});
                yy = cell2struct(yy,param_struct.rhythm_names);
                for r = 1:length(rhythms)
                    for l = 1:length(rhythms)
                        xxx = xx.(rhythms{r});%xxx-rhythm: matrix,each col is an oscillation
                        yyy = yy.(rhythms{l});
                        if ~isempty(xxx) && ~isempty(yyy)
                            lenx = size(xxx,2);
                            leny = size(yyy,2);
                            for lx = 1:lenx
                                for ly = 1:leny
                                    if ~exist('mPLV','var')
                                        [mPLV.(join(["mPLV_",rhythms{r},'_',rhythms{l}],'')),~] ...    
                                            = PLV_calc(xxx(:,lx),yyy(:,ly),param_struct); % calculate PLV                                 
                                    elseif exist('mPLV','var') && isfield(mPLV,join(["mPLV_",rhythms{r},'_',rhythms{l}],''))
                                        [mPLV.(join(["mPLV_",rhythms{r},'_',rhythms{l}],''))(end+1), ~] ...
                                            = PLV_calc(xxx(:,lx),yyy(:,ly),param_struct); % calculate PLV
                                    elseif exist('mPLV','var') && ~isfield(mPLV,join(["mPLV_",rhythms{r},'_',rhythms{l}],''))...
                                            && isfield(mPLV,join(["mPLV_",rhythms{l},'_',rhythms{r}],''))
                                        [mPLV.(join(["mPLV_",rhythms{l},'_',rhythms{r}],''))(end+1),~] ...    
                                            = PLV_calc(xxx(:,lx),yyy(:,ly),param_struct); % calculate PLV
                                    else
                                        [mPLV.(join(["mPLV_",rhythms{r},'_',rhythms{l}],'')), ~] ...    
                                            = PLV_calc(xxx(:,lx),yyy(:,ly),param_struct); % calculate PLV
                                    end  
                                end
                            end
                        else
                            if ~exist('mPLV','var')
                                mPLV.(join(["mPLV_",rhythms{r},'_',rhythms{l}],'')) = NaN;                              
                            elseif exist('mPLV','var') && isfield(mPLV,join(["mPLV_",rhythms{r},'_',rhythms{l}],''))
                                mPLV.(join(["mPLV_",rhythms{r},'_',rhythms{l}],''))(end+1) = NaN;
                            elseif exist('mPLV','var') && ~isfield(mPLV,join(["mPLV_",rhythms{r},'_',rhythms{l}],''))...
                                        && isfield(mPLV,join(["mPLV_",rhythms{l},'_',rhythms{r}],''))
                                    mPLV.(join(["mPLV_",rhythms{l},'_',rhythms{r}],''))(end+1) = NaN;
                            else
                                mPLV.(join(["mPLV_",rhythms{r},'_',rhythms{l}],'')) = NaN;
                            end  
                        end
                    end
                end
            end         
        end
    end
end
mean_PLV = NaN(1,param_struct.bands_length);
if exist('mPLV','var')
    names = fieldnames(mPLV);
    for i = 1:length(names)
        name = names{i};
        if all(isnan(mPLV.(name)))
            mPLV.(name) = 0;
        else
            %oldPLV(i) = mean(mPLV.(name));
            inds = find(mPLV.(name) > plvThresh);               
            plvBoost = length(inds) / length(mPLV.(name))*mean(mPLV.(name)(inds));             
            if isempty(inds)
                mPLV.(name) = 0;
            else
                mPLV.(name) = mean(mPLV.(name)(inds)) + plvBoost;
    
            end
        end
        mean_PLV(i) = mPLV.(name);
    end
elseif ~exist('mPLV','var') 
    mean_PLV = NaN(1,param_struct.bands_length);
    names = [];
end
end