%Advanced DSP Project
%Effect of Deep Brain Stimulation on Parkinsonian Tremor

% Algorithm that finds Bispectral Area Peaks. The process is initialized
% by setting as starting point the leftmost element in the coord matrix.
% Then all neighboring elements of the startPoint are checked. If
% neighbours with magnitude (z coord) less than the startPoint's magnitude 
% are found, the algorithm deletes those elements. Else if one or more 
% neighbours with magnitude greater than the startPoint's are found, then
% the neighbour with the max magnitude survives and becomes the startPoint
% while the others are deleted. 
function Peaks = findPeaks(xcoord, ycoord, Bispec)

[r,c] = size(xcoord);
if r > c
    xcoord = xcoord';
    ycoord = ycoord';
end
    
%Leftmost point is startPoint 
startPoint = [xcoord(1); ycoord(1)];
coord = [xcoord; ycoord];
[~,n] = size(coord);
Peaks = zeros(2,n);
oldNei = [];
k = 1;
if n == 1
    Peaks = startPoint;
    n = 0;
end

while n>=1
    if n ~= 1
        Nei = zeros(2,8);%Possible neighbors are all neighboring boxes in a radius of 1 (8 neighbors)
        for i = 2:n %%%Locating neighbors based on distance
             posNei = coord(:,i);
            if sqrt((startPoint(1)-posNei(1))^2 + (startPoint(2)-posNei(2))^2)<=sqrt(2) &&...
                    sqrt((startPoint(1)-posNei(1))^2 + (startPoint(2)-posNei(2))^2)>0
                Nei(:,i-1) = posNei;  
            end
        end %%%
        if isempty(oldNei) == 0
            for i = 1:length(oldNei)
                posNei = oldNei(:,i);
                if sqrt((startPoint(1)-posNei(1))^2 + (startPoint(2)-posNei(2))^2)<=sqrt(2) &&...
                    sqrt((startPoint(1)-posNei(1))^2 + (startPoint(2)-posNei(2))^2)>0
                    Nei(:,i) = posNei;  
                end
            end
        end
        %%%Finding and deleting found neighbors from the coord matrix, to
        %%%avoid picking them again
        for i = 1:length(Nei(1,:))
                for j =1:n
                    if Nei(:,i) == coord(:,j)
                        coord(:,j) = [0,0];
                    end
                end    
        end
        coord(:,all(coord == 0 )) = [];
        Nei(:,all(Nei == 0 )) = []; %Because Nei is initialized with 0's
        Nei = unique(Nei','rows')';
        oldNei = [oldNei Nei];
        oldNei = unique(oldNei','rows')';
        %%%
        if isempty(Nei) == 0
            NeiBispec = zeros(1,length(Nei(1,:)));
            for i = 1:length(Nei(1,:))
                 NeiBispec(:,i) = Bispec(Nei(1,i),Nei(2,i)); %Finding z values in Bispec for each Nei 
            end
            %%%Finds max Nei and compares it to z value of startPoint
            %%%If maxNei>startPoint, maxNei is new startPoint and
            %%%startPoint becomes the 1st element
            if max(NeiBispec) > Bispec(startPoint(1),startPoint(2)) + 10^(-7)
                oldNei = [oldNei startPoint];
                [~,s] = max(NeiBispec);
                startPoint = Nei(:,s);
                coord(:,1) = startPoint;
            else
                % If maxNei < startPoint, then startPoint is a Peak,
                % is deleted from coord and the next element becomes startPoint
                Peaks(:,k) = startPoint;
                oldNei = [oldNei startPoint];
                if isempty(coord) == 0
                    coord(:,1) = [];
                end
                if (isempty(coord) == 0)
                    startPoint = coord(:,1);
                    k = k + 1; 
                end
            end
            %%%If Nei empty, startPoint is Peak by default
        elseif isempty(Nei) == 1 
            Peaks(:,k) = startPoint;
            k = k+1;
            coord(:,1) = [];
            startPoint = coord(:,1);
        end 
        n = length(coord(1,:)); 
    else %%%If coord has only 1 element, then it is Peak by default
        Nei = zeros(2,8);
        if isempty(oldNei) == 0
            for i = 1:length(oldNei)
                posNei = oldNei(:,i);
                if sqrt((startPoint(1)-posNei(1))^2 + (startPoint(2)-posNei(2))^2)<=sqrt(2) &&...
                    sqrt((startPoint(1)-posNei(1))^2 + (startPoint(2)-posNei(2))^2) > 0
                    Nei(:,i) = posNei;  
                end
            end
        end
        Nei(:,all(Nei == 0 )) = []; %Because Nei is initialized with 0's
        if isempty(Nei) == 0
            NeiBispec = zeros(1,length(Nei(1,:)));
            for i = 1:length(Nei(1,:))
                 NeiBispec(:,i) = Bispec(Nei(1,i),Nei(2,i)); %Finding z values in Bispec for each Nei 
            end
            %%%Finds max Nei and compares it to z value of startPoint
            %%%If maxNei>startPoint, maxNei is new startPoint and
            %%%startPoint becomes the 1st element
            if max(NeiBispec) <= Bispec(startPoint(1),startPoint(2))
                Peaks(:,k) = startPoint;
                k = k+1;
            end
        end
        n = -1;
    end
end
Peaks(:,all(Peaks == 0 )) = [];
Peaks = unique(Peaks','rows')';
