function totalBic = sBic_calc(bic,upLim1,upLim2,lowLim1,lowLim2)

lowLim1 = ceil(lowLim1); lowLim2 = ceil(lowLim2); 
upLim1 = ceil(upLim1); upLim2 = ceil(upLim2);
bic = bic(lowLim1:upLim1,lowLim2:upLim2);
totalBic = sum(sum(bic));

end