%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function x_new = unwrap_SwDs(x)
% This function is used to unwrap the cell of the SwDs. 
%% Inputs: 
% x       -An electrode in form of cell vector, containing SwD oscillations, grouped as
%          rhythms. Each cell of the cell vector can contain one or more oscillations
%          This function separates oscillations that belong to the same
%          rhythm.
% 
%% OUTPUTS: 
% x_new   -An electrode, containing SwD oscillations, in form of a
%          cell vector. Each cell contains one oscillation.

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

[l,h] = size(x);
if l < h
    x = transpose(x);
end
    
for i = 1:length(x)
    [~,c] = size(x{i});
    col = 1:c;
    if c > 1
        for j = 2:c
            x{end+1} = x{i}(:,col(j));
        end
        x{i} = x{i}(:,col(1));
    end
end

x_new = x;
end