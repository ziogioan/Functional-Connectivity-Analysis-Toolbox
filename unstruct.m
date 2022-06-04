%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample
function [feature,out_names] = unstruct(feature_struct,bands_length,feat_length)
% This function is used to turn the struct returned from
% regional_analysis.m to a double array.
% 
%% Inputs:   
% feature_struct    -struct. Contains the results of the connectivity
%                    analysis. It is returned by regional_analysis.m
% bands_length      -double. The number of bands that used for feature
%                    extraction
% feat_length       -double. The number of features that were extracted
%% Outputs:
% feature:          -double array. Contains the values of the
%                    feature_struct
% out_names:        -string array. Contains the name of each feature value

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

feat_names = feature_struct.names;
feature_struct = rmfield(feature_struct, 'names');
names = fieldnames(feature_struct);
feature = NaN(1,bands_length*length(names)*feat_length);
out_names = strings(1,bands_length*length(names)*feat_length);
count = 1;
for i = 1:bands_length*feat_length
    for j = 1:length(names)
        name1 = names{j};
        name2 = feat_names{i};
        out_names(count) = convertCharsToStrings([name2,'_',name1]);
        feature(count) = feature_struct.(name1)(i);
        count = count + 1;
    end
end

end