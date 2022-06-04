function V = RelPow_calc(bspec_cut,bspec,symmetric)

bspec = abs(bspec);
if symmetric
    bspec = triu(bspec);
    bspec_cut = triu(abs(bspec_cut));
end

Vtotal = trapz(trapz(bspec.^2,2));
V = trapz(trapz(bspec_cut.^2,2));
V = V/Vtotal*100;

end
