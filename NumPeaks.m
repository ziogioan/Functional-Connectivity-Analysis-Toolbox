function noPeaks = NumPeaks(Peaks,upLim1,upLim2,lowLim1,lowLim2)

noPeaks = 0;
for i = 1:length(Peaks(1,:))
    if Peaks(1,i) > lowLim1 && Peaks(1,i) < upLim1...
            && Peaks(2,i) > lowLim2 && Peaks(2,i) < upLim2
        noPeaks = noPeaks + 1;
    end
end

end