function bands = create_bands(lower_limit,upper_limit,bandwidth,overlap,untouched)
% Function to create the desired non conventional bands
% with specified bandwidth and overlap.

if length(bandwidth) == 1
    if length(untouched) >= 1
        error('utouched should be empty')
    end
    bands = create_bands_pre(lower_limit,upper_limit,bandwidth,overlap);
else
    if length(bandwidth) - 1 ~= length(untouched)
        error('The length of the untouched must be equal to length of bandwidth - 1\n')
    end
      
    untouched(end + 1) = upper_limit;
    bands = create_bands_pre(lower_limit,untouched(1),bandwidth(1),overlap);
    for i = 2:length(bandwidth)
        temp = bands(end) - bandwidth(i)*overlap/100;
        bands = [bands create_bands_pre(temp,untouched(i),bandwidth(i),overlap)];
    end

end

if bands(end) > upper_limit && (bands(end-2) + round(0.5*bandwidth) < upper_limit)
    bands(end) = upper_limit;
elseif bands(end) == upper_limit
    
else
    bands = bands(1:end-2);
end
end