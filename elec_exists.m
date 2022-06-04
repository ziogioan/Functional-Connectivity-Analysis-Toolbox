function [flag,pos] = elec_exists(elec_mat, elec)

flag = 0;
counter = 1;
pos = NaN;
for i = 1:length(elec_mat)
    if elec_mat{i} == elec
        flag = 1;
        pos = counter;
    end
    counter = counter + 1;
end