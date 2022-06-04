function coher = COH_calc(x,y,Fs,surr_struct)
%Calculate coherence with Welch method and find significance threshold
%using a surrogate data method. Available methods: a) random permutation,
% b) Koopmans method, c) constant threshold

surr_method = surr_struct.st_sig_method;
f = surr_struct.f;
window = surr_struct.window;
overlap = surr_struct.overlap;
coher = mscohere(x,y,window,overlap,f,Fs);


if surr_method == "permutation"
    N = surr_struct.N;
    sur_coher = NaN(length(f),N);
    alpha = surr_struct.alpha;
    if alpha > 1
        error("Alpha must be between 0-1")
    end
    if N < 1
        error("Number of surrogates(N) must be positive integer")
    end
    for i = 1:N
        ind = randperm(length(x));
        z = y(ind);
        sur_coher(:,i) = mscohere(x,z,window,overlap,f,Fs);
    end
    coher_thresh = prctile(sur_coher',(1-alpha)*100);

elseif surr_method == "Koopmans" && overlap == 0
    alpha = surr_struct.alpha;
    lenx = length(x);
    noseg = lenx/window;
    K = noseg;
    coher_thresh = 1 - alpha^(1/(K-1));
    
elseif surr_method == "constant" 
    coher_thresh = surr_struct.thres;
end

coher(find(coher < coher_thresh)) = 0;

end