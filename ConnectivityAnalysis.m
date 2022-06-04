function feature = ConnectivityAnalysis(data, param_struct)

% In this function Functional Connectivity featres are calculated
% 
% Inputs:      data            -struct or char.
%              param_struct    -struct containing the following info
% Info:
%              feature_name    -string. "coher" or "plv" or "cross_bspec"
%                               or "pac"
%              chans           -string or string array. "all" or array of
%                               strings e.g ["Fp1","Fp2",...]
%              name            -char. For saving the results 
%              reg             -string. Regional analysis type. "intra-regional"
%                               or "inter-regional" or "left-right" or 
%                               "anti-symmetric"
%%                         ***If feature_name == "coher"***
%              window          -double. Window for the calculation of
%                               coherence via the welch method
%              overlap         -double. Overlap for the calculation of
%                               coherence via the welch method
%              st_sig_method   -string. method for assessing the
%                               statistical significance. "Koopmans" or 
%                               "permutation" or "constant"
%              N               -double. If st_sig_method == "permutation" N
%                               is the number of permutations
%              alpha           -double. If st_sig_method == "permutation"
%                               or "Koopmans", alpha is the level of
%                               statistical significance
%              thres           -double. If st_sig_method == "constant",
%                               thres is the threshold of statistical
%                               significance
%              band_method     -string. "conventional" or
%                               "non-conventional"
%%                          If band_method == "non-conventional"
%              bandwidth       -double array. Bandwidth for segmenting the
%                               coherence function
%              untouched       -double array. limits where the different bandwidths
%                               will be applied
%              overlap         -double. Overlap between bands. If 
%                               st_sig_method == "Koopmans" overlap must be
%                               0
%              lower_limit     -double. lower frequency of interset
%              upper_limit     -double. higher frequency of interset
%%                          If band_method == "conventional"
%              bands           -double array. The limits of the desired
%                               bands. e.g [1,4,4,10,...]
%              conv_bands      -string array. The names of the bands in
%                               "bands", e.g ["delta","theta",...]
%%                         ***If feature_name == "plv"***
%              conv_bands      -string array. The names of the bands of
%                               interest as given in SwD_preprocessed_EEGs.m
%              st_sig_method   -string. method for assessing the
%                               statistical significance. "block-resampling"
%                               or "permutation" or "noSurrogate"
%              N               -double. If st_sig_method == "permutation" or
%                               "block-resampling" N is the number of permutations
%              alpha           -double. If st_sig_method == "permutation"
%                               or "block-resampling", alpha is the level of
%                               statistical significance
%              blocksize       -double. if st_sig_method == "block-resampling"
%                               this parameter decides the size of the
%                               blocks
%              plvThresh       -double. Threshold for plv statistical
%                               significance algorithm. For more info check
%                               plv_swarm_calc
%%                         ***If feature_name == "cross_bspec"***
%              nfft            -double. Nfft for the calculation of cross
%                               Bispectrum. 
%              cBS_wind        -double. Window for the calculation of
%                               cross Bispectrum
%              cBS_overlap     -double. Overlap for the calculation of
%                               cross Bispectrum
%              st_sig_method   -string. method for assessing the
%                               statistical significance. "normalization"
%                               or "constant" or "nothing"
%              per_thres       -double. If st_sig_method == "constant",
%                               thres is the threshold of statistical
%                               significance as a percentage of the max
%                               value of the Bispectrum
%              features        -string array. Features to be calculated in 
%                               each bispectral block ["TotalBIC","RelPow",
%                               "BSEnt","BICEnt","BSnPks","BICnPks"]
%              band_method     -string. "conventional" or
%                               "non-conventional"
%%                          If band_method == "non-conventional"
%              bandwidth       -double array. Bandwidth for segmenting the
%                               Bispectrum
%              untouched       -double array. limits where the different bandwidths
%                               will be applied
%              overlap         -double. Overlap between bands
%              lower_limit     -double. lower frequency of interset
%              upper_limit     -double. higher frequency of interset
%%                          If band_method == "conventional"
%              bands           -double array. The limits of the desired
%                               bands. e.g [1,4,4,10,...]
%              conv_bands      -string array. The names of the bands in
%                               "bands", e.g ["delta","theta",...]
%%                         ***If feature_name == "pac"***
%              measure         -string. Algorithm for pac calculation "MI"
%                               or "GLM" or "MVL" 
%              nBins           -double. If measure == "MI". For more info
%                               check the modulationIndex.m function.
%              flag_AIC        -double. If measure == "GLM". For more info
%                               check the GLM_CFC.m function.
%              normalise       -double. 1 or 0. Min-Max scaling for the
%                               envelope (amplitude) of the high frequency
%                               signal.
%              st_sig_method   -string. "block-resampling" or "block-swapping"
%                               or "noSurrogate"
%              N               -double. If st_sig_method == "block-resampling" or
%                               "block-swapping" N is the number of permutations
%              alpha           -double. If st_sig_method == "permutation"
%                               or "block-resampling", alpha is the level of
%                               statistical significance
%              block_size      -double. if st_sig_method == "block-resampling"
%                               or "block-swapping" this parameter defines
%                               the size of the blocks
%              phAmpExtr       -string. This parameter defines the method
%                               that will be used for the phase and the 
%                               amplitude extraction. "tf" or "hilbert"
%%                          If band_method == "tf"
%              step_low        -double. step in Hz for calculating pac. For
%                               more info check tf_PAC_surr.m
%              step_high       -double. step in Hz for calculating pac. For
%                               more info check tf_PAC_surr.m
%              band_method     -string. "conventional" or
%                               "non-conventional"
%%                          If band_method == "conventional"
%              bandsLow        -double array. Limits of the desired bands
%                               e.g [1,5,5,11] for the low frequencies. For the
%                               bandwidth of each band the following must
%                               be true: bandwidth mod step_low = 0
%              bandsHigh       -double array. Limits of the desired bands
%                               e.g [25,55,65,100] for the high frequencies. For
%                               the bandwidth of each band the following must
%                               be true: bandwidth mod step_high = 0
%              conv_bandsLow   -string array. The names of the given low bands
%                               e.g ["delta","theta"]
%              conv_bandsHigh  -string array. The names of the given high bands
%                               e.g ["lo-gamma","hi-gamma"]
%              untouched       -double array. limits where the different bandwidths
%                               will be applied
%              overlap         -double. Overlap between bands
%              lower_limit     -double. lower frequency of interest
%              upper_limit     -double. higher frequency of interest
%%                          If band_method == "non-conventional"
%              bandwidthLow    -double array. The bandwidths of the desired
%                               low frequency bands. e.g [4,6]
%              bandwidthHigh   -double array. The bandwidths of the desired
%                               high frequency bands. e.g [10,15]
%              untouchedLow    -double array. It's length should be 
%                               length(bandwidthLow)-1. e.g [10]. So up to
%                               10Hz the bandwidth will be 4 and for
%                               frequencies higher than 10 the bandwidth
%                               will be 6
%              untouchedHigh   -double array. It's length should be 
%                               length(bandwidthHigh)-1. e.g [40]. So up to
%                               40Hz the bandwidth will be 10 and for
%                               frequencies higher than 40 the bandwidth
%                               will be 15
%              overlapLow      -double. The overlap of the low frequency
%                               bands
%              overlapHigh     -double. The overlap of the high frequency
%                               bands
%              lower_limitLow  -double. The lower limit of the low frequency
%                               bands in Hz
%              upper_limitLow  -double. The upper limit of the low frequency
%                               bands in Hz
%              lower_limitHigh -double. The lower limit of the high frequency
%                               bands in Hz
%              upper_limitHigh -double. The upper limit of the high frequency
%                               bands in Hz
%              features        -string array. The features that will be calculated
%                               in each frequency block defined by the created bands
%                               ["totalPac","Ent1",","RelPow"]
%%                          If band_method == "hilbert"
%              pacThresh       -double. Threshold for pac statistical
%                               significance algorithm. For more info check
%                               pac_swarm_calc
%              low             -cell of chars. The bands that will be
%                               considered as low for the phase extraction
%                               e.g {'delta','theta'}
%              high            -cell of chars. The bands that will be
%                               considered as high for the envelope extraction
%                               e.g {'hi_beta','hi_gamma'}


if ischar(data)
    if isfile(data)
        warning('You have provided the path to a single .mat file');
    elseif isfolder(data)
        warning('You have provided the path to a folder of .mat files');
    end
end
if isstruct(data)
    warning('You have provided a data struct');
end

if ~isfield(param_struct,'feature_name')
    error('Feature name is necessary')
end
if param_struct.feature_name ~= "coher" && param_struct.feature_name ~= "plv"...
    && param_struct.feature_name ~= "cross_bspec" && param_struct.feature_name ~= "pac"
    error('Enter a valid feature name')
end

if ~isfield(param_struct,'channs')
    param_struct.channs = "all";
end

if ~isfield(param_struct,'name')
    warning('Results will not be saved without the name parameter')
end

if ~isfield(param_struct,'reg')
    error('reg is necessary')
end
if param_struct.reg ~= "intra-regional" && param_struct.feature_name ~= "inter-regional"...
    && param_struct.feature_name ~= "left-right" && param_struct.feature_name ~= "anti-symmetric"
    error('Enter a valid reg name')
end

warning('If a parameter misses it will be automatically be set to default\nand it might not be compatible. So check that you entered all the nessecary parameters')

if param_struct.feature_name == "coher"
    str = 'COH';
    if ~isfield(param_struct,'band_method')
        param_struct.band_method = "non-conventional";
    end
    if ~isfield(param_struct,'bandwidth')
        param_struct.bandwidth = 3;
    end
    if ~isfield(param_struct,'untouched')
        param_struct.untouched = [];
    end
    if ~isfield(param_struct,'overlap')
        param_struct.overlap = 0;
    end
    if ~isfield(param_struct,'lower_limit')
        param_struct.lower_limit = 1;
    end
    if ~isfield(param_struct,'upper_limit')
        param_struct.upper_limit = 50;
    end
    if ~isfield(param_struct,'window')
        param_struct.window = 250;
    end
    if ~isfield(param_struct,'overlap')
        param_struct.overlap = 0;
    end
    if ~isfield(param_struct,'st_sig_method')
        param_struct.st_sig_method = "Koopmans";
    end
    if ~isfield(param_struct,'alpha')
        param_struct.alpha = 0.05;
    end
    
    if param_struct.band_method == "conventional" && (~isfield(param_struct,'bands') || ~isfield(param_struct,'conv_bands'))
        error('bands and conv_bands are necessary when conventional bands are used\n\n')
    end
    if param_struct.st_sig_method == "Koopmans" && param_struct.overlap ~= 0
        error('Overlap should be 0 when the Koopamns method is used')
    end
elseif param_struct.feature_name == "plv"
    str = 'PLV';
    if ~isfield(param_struct,'conv_bands')
        param_struct.conv_bands = ["delta","theta","alpha","lo_beta","hi_beta","lo_gamma","hi_gamma"];
    end
    if ~isfield(param_struct,'st_sig_method')
        param_struct.st_sig_method = "noSurrogate";
    end
    if ~isfield(param_struct,'plvThresh')
        param_struct.plvThresh = 0.35;
    end


elseif param_struct.feature_name == "cross_bspec"
    str = 'crossBspec';
    if ~isfield(param_struct,'band_method')
        param_struct.band_method = "non-conventional";
    end
    if ~isfield(param_struct,'bandwidth')
        param_struct.bandwidth = 6;
    end
    if ~isfield(param_struct,'untouched')
        param_struct.untouched = [];
    end
    if ~isfield(param_struct,'overlap')
        param_struct.overlap = 25;
    end
    if ~isfield(param_struct,'lower_limit')
        param_struct.lower_limit = 1;
    end
    if ~isfield(param_struct,'upper_limit')
        param_struct.upper_limit = 29;
    end
    if ~isfield(param_struct,'nfft')
        param_struct.nfft = 256;
    end
    if ~isfield(param_struct,'cBS_wind')
        param_struct.cBS_wind = 5;
    end
    if ~isfield(param_struct,'cBS_overlap')
        param_struct.cBS_overlap = 50;
    end
    if ~isfield(param_struct,'st_sig_method')
        param_struct.st_sig_method = "normalization";
    end
    if ~isfield(param_struct,'features')
        param_struct.features = ["TotalBIC","RelPow","BSEnt","BICEnt","BSnPks","BICnPks"];
    end
    
    if param_struct.band_method == "conventional" && (~isfield(param_struct,'bands') || ~isfield(param_struct,'conv_bands'))
        error('bands and conv_bands are necessary when conventional bands are used\n\n')
    end
elseif param_struct.feature_name == "pac"
    str = 'PAC';
    if ~isfield(param_struct,'measure')
        param_struct.measure = "MI";
    end
    if ~isfield(param_struct,'nBins')
        param_struct.nBins = 18;
    end
    if ~isfield(param_struct,'normalise')
        param_struct.normalise = 0;
    end
    if ~isfield(param_struct,'st_sig_method')
        param_struct.st_sig_method = "noSurrogate";
    end
    if ~isfield(param_struct,'phAmpExtr')
        param_struct.phAmpExtr = "hilbert";
    end
    if ~isfield(param_struct,'pacThresh')
        param_struct.pacThresh = 0;
    end
    if ~isfield(param_struct,'low')
        param_struct.low = {'delta','theta','alpha'};
    end
    if ~isfield(param_struct,'high')
        param_struct.high = {'hi_beta','lo_gamma','hi_gamma'};
    end
end

if param_struct.reg == "intra-regional"
    regions = ["F","P","C","T","O"];
    combinations = length(regions);
    param_struct.regional = param_struct.reg;
    str2 = 'intra';
elseif param_struct.reg == "inter-regional"
    regions = ["F","P","C","T","O"];
    combinations = sum(1:length(regions)-1);
    param_struct.regional = param_struct.reg;
    str2 = 'inter';
elseif param_struct.reg == "left-right"
    regions = ["L1","L2","L3","R1","R2","R3"];
    combinations = (length(regions)/2)^2;
    param_struct.regional = param_struct.reg;
    str2 = 'left_right';
elseif param_struct.reg == "anti-symmetric"
    regions = ["L1","L2","L3","R1","R2","R3"];
    combinations = (length(regions)/2)*2;
    param_struct.regional = param_struct.reg;
    str2 = 'anti-symmetric';
end

param_struct.regional = param_struct.reg;
feature_name = join([param_struct.feature_name,"_swarm_calc"],'');
if ~contains(feature_name,'pac')
    if contains(feature_name,'coh') && param_struct.band_method == "non-conventional"
        bands = create_bands(param_struct.lower_limit,param_struct.upper_limit,...
            param_struct.bandwidth,param_struct.overlap,param_struct.untouched);
        bands_length = length(bands)/2;
        param_struct.feat_length = 1;
        param_struct.bands_length = bands_length;
        param_struct.bands = bands;
    elseif contains(feature_name,'coh') && param_struct.band_method == "conventional"
        bands = param_struct.bands;
        bands_length = length(bands)/2;
        param_struct.bands_length = bands_length;
        param_struct.bands = bands;
    elseif contains(feature_name,'plv')
        bands_length = sum(1:length(param_struct.conv_bands));
        param_struct.bands_length = bands_length;
        param_struct.feat_length = 1;
    elseif contains(feature_name,'cross_bspec') && param_struct.band_method == "non-conventional"
        bands = create_bands(param_struct.lower_limit,param_struct.upper_limit,...
            param_struct.bandwidth,param_struct.overlap,param_struct.untouched);
        bands_length = (length(bands)/2)^2; %Scan the whole bi-frequency plane
        param_struct.feat_length = length(param_struct.features);
        param_struct.bands_length = bands_length;
        param_struct.bands = bands;
    elseif contains(feature_name,'cross_bspec') && param_struct.band_method == "conventional"
        bands = param_struct.bands;
        bands_length = length(bands)/2;
        param_struct.bands_length = bands_length;
        param_struct.bands = bands;
    end
else
    if param_struct.phAmpExtr == "tf"
        warning('bandwidth mod step must be 0')
        bandsLow = create_bands(lower_limitLow,upper_limitLow,bandwidthLow...
            ,overlapLow,untouchedLow);
        bandsHigh = create_bands(lower_limitHigh,upper_limitHigh,bandwidthHigh...
            ,overlapHigh,untouchedHigh);
        bands_length = (length(bandsLow)/2)*(length(bandsHigh)/2);
        param_struct.feat_length = length(param_struct.features);
        param_struct.bandsLow = bandsLow;
        param_struct.bandsHigh = bandsHigh;
    else
        bands_length = length(param_struct.low)*length(param_struct.high);
        param_struct.feat_length = 1; %No features are extracted for Swarm strategy, only PAC value
    end
    param_struct.bands_length = bands_length;
end


if ischar(data) && ~contains(data,'.mat')
    dir_data = dir(data);
    if contains(feature_name,'bspec') || contains(feature_name,'pac')
        feature = NaN((length(dir_rhythm_mat_files_cond)-2)/2,...
            combinations*bands_length*param_struct.feat_length);
    else
        feature = NaN((length(dir_rhythm_mat_files_cond)-2)/2,combinations*bands_length);
    end
    count = 1;
    for i = 1:length(dir_data) - 2
        name = dir_rhythm_mat_files_cond(i+2).name;
        load(fullfile(path_rhythm_mat_files_cond,name));
        param_struct.new_length = new_length; param_struct.wind_dur = wind_dur;
        param_struct.num_wind = num_wind; param_struct.Fs = new_Fs;
        param_struct.chan_names = chan_names; param_struct.rhythms = rhythms;
        param_struct.chan_names_all = chan_names_all;
        [~,feat] = feature_extraction_pipeline("regional", feature_name,param_struct);
        [feature(count,:),names] = unstruct(feat,bands_length,param_struct.feat_length);
        count = count + 1;
    end
    
elseif ischar(data) && contains(data,'.mat')
    load(data);
    param_struct.new_length = new_length; param_struct.wind_dur = wind_dur;
    param_struct.num_wind = num_wind; param_struct.Fs = new_Fs;
    param_struct.chan_names = chan_names; param_struct.rhythms = rhythms;
    param_struct.chan_names_all = chan_names_all;
    [~,feat] = feature_extraction_pipeline("regional", feature_name,param_struct);
    [feature,names] = unstruct(feat,bands_length,param_struct.feat_length);
elseif isstruct(data)
    param_struct.new_length = data.new_length; param_struct.wind_dur = data.wind_dur;
    param_struct.num_wind = data.num_wind; param_struct.Fs = data.new_Fs;
    param_struct.chan_names = data.chan_names; param_struct.rhythms = data.rhythms;
    param_struct.chan_names_all = data.chan_names_all;
    [~,feat] = feature_extraction_pipeline("regional", feature_name,param_struct);
    [feature,names] = unstruct(feat,bands_length,param_struct.feat_length);
elseif ~ischar(data) && ~isstruct(data)
    error('data should be char or struct\n\n')
end
names = convertStringsToChars(names);
feature = array2table(feature,'VariableNames',names);
if isfield(param_struct,'name')
    writetable(feature, [str, '_', convertStringsToChars(name), '_',str2,'.csv'])
end

end