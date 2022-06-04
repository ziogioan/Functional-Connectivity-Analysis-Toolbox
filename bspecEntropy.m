function Ent = bspecEntropy(bspec,upLim1,upLim2,lowLim1,lowLim2,q)

lowLim1 = ceil(lowLim1); lowLim2 = ceil(lowLim2); 
upLim1 = ceil(upLim1); upLim2 = ceil(upLim2);
%try
bspec = abs(bspec(lowLim1:upLim1,lowLim2:upLim2)).^q;
%catch
%    a = 1;
%end
p = bspec./sum(sum(bspec));
if isnan(p)
    p = zeros(size(p,1),size(p,2));   
end
p = p + 10e-10;
Ent = -sum(sum(p.*log(p)));

end