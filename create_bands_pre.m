function bands = create_bands_pre(lower_limit,upper_limit,bandwidth,overlap)
% Function to create the desired non conventional bands
% with specified bandwidth and overlap.

bands(1) = lower_limit;
bands(2) = lower_limit + bandwidth;

counter = 3;
while bands(end) < upper_limit 
    bands(counter) = bands(end) - bandwidth*(overlap/100);
    counter = counter + 1;
    bands(counter) = bands(end) + bandwidth;
    counter = counter + 1;
end

end