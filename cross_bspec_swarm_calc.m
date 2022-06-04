function [mean_feat,out_names] = cross_bspec_swarm_calc(param_struct)

% Calculates the coherence between SwDs of different electrodes when
% both of them belong to the same non - conventional band.
%
% Inputs: elec_case    - string specifying method used: can be "single" for
%                        electrode to electrode coherence, or any of the regional options
%                        specified in regional_analysis.m for region to region coherence
%         x            - data of the one electrode (struct)
%         y            - data of the other electrode (struct)
%         bands        - non-conventional bands of interest
%         Fs           - sampling rate
%         surr_struct  - struct containing various parameters 
%  For more details, see regional_analysis.m


% Outputs: mean_coher - mean of coherence at each band and pair of electrodes
%                       across windows
%          std_coher  - std of coherence at each band and pair of electrodes
%                       across windows

%Collect signals from input struct and sum all components of each electrode


a = param_struct.a;
b = param_struct.b;
a = gather_sigs(a);
b = gather_sigs(b);

%------------------- Check if input is in correct form --------------------
if isstruct(a)
    a = struct2cell(a);
    a = unwrap_SwDs(a); %if there are two SwDs in a band then separate them
elseif ismatrix(a)
    [l,h] = size(a);
    if l < h
        a = transpose(a);
    end
    for t = 1:length(a(1,:))
        new_a{t,1} = a(:,t);
    end
    a = new_a;
end
if isstruct(b)
    b = struct2cell(b);
    b = unwrap_SwDs(b); %if there are two SwDs in a band then separate them
elseif ismatrix(b)
    [l,h] = size(b);
    if l < h
        b = transpose(b);
    end
    for t= 1:length(b(1,:))
        new_b{t,1} = b(:,t);
    end
    b = new_b;
end

elec_type = param_struct.regional;
% sumup = param_struct.sumup;
% if sumup == 1 && (elec_type == "single" || elec_type == "intra-regional" || elec_type == "anti-symmetric")
%     error('Sum up can only be used for inter-regional & left-right')
% end
% if sumup == 1
%     a = sum(cell2mat(a'),2);
%     b = sum(cell2mat(b'),2);
%     a = mat2cell(a,length(a),1);
%     b = mat2cell(b,length(b),1);
% end

if (elec_type == "single" || elec_type == "intra-regional" || elec_type == "anti-symmetric")
    param_struct.x = a;
    param_struct.y = b;
    [mean_feat,out_names] = cBisCalc(param_struct);
elseif elec_type == "left-right" || elec_type == "inter-regional"
    param_struct.x = a;
    param_struct.y = b;
    [mean_feat1,out_names1] = cBisCalc(param_struct);

    param_struct.x = b;
    param_struct.y = a;
    [mean_feat2,out_names2] = cBisCalc(param_struct);
    if all(convertCharsToStrings(out_names1) == convertCharsToStrings(out_names2)) ~= 1
        error('Different number of features for a/b & b/a is not possible here')
    else
        out_names = out_names1;
    end

    mean_feat = [mean_feat1 mean_feat2];
end
    

function [mean_feat,out_names] = cBisCalc(param_struct)
    elec_case = param_struct.regional;
    features = param_struct.features;
    if ischar(features)
        features = convertCharsToStrings(features);
    end
    symmetric = 0;
    x = param_struct.x;
    y = param_struct.y;
    band_method = param_struct.band_method;
    if band_method == "conventional"
        conv_bands = param_struct.conv_bands;
    end
    bands = param_struct.bands;
    Fs = param_struct.Fs;
    nfft = param_struct.nfft;
    method = param_struct.st_sig_method;
    wind = param_struct.cBS_wind;
    nsamp = param_struct.nfft;
    overlap = param_struct.cBS_overlap;
    for i = 1:length(x)
        for j = 1:length(y)
            if elec_case == "single" || elec_case == "intra-regional" || elec_case == "anti-symmetric"
                ifcond = (j ~= i); ifcondAS = (j ~= i);
            elseif elec_case == "inter-regional"  || elec_case == "left-right"
                ifcond = 1; ifcondAS = 1;
            end
            
            %Anti-Symmetric Cross-Bispectrum
            % j > i is needed because Anti-Symmetric cBS is
            % Bi-Directional - Fp1-Fp2 pair is enough to describe
            % bi-directional interactions
            if ifcond && ~isempty(x{i}) && ~isempty(y{j}) && ~(all(isnan(x{i})) || all(isnan(y{j})))
                xin = (x{i} - mean(x{i})) / std(x{i}); %Normalize
                yin = (y{j} - mean(y{j})) / std(y{j});
                
                cBic = bicoherx(xin,yin,yin,nfft,wind,nsamp,overlap);

                if ifcondAS
                    [bspecxxy,nrecs] = bispecdx(xin,xin,yin,nfft,wind,nsamp,overlap,0);
                    bspecxyx = bispecdx(xin,yin,xin,nfft,wind,nsamp,overlap,0);
                
                    %A normalization procedure to remove artificial interactions
%                     if method == "Normalization"
%                         bspecxyx = sqrt(nrecs)*real(bspecxyx)/std(reshape(real(bspecxyx),1,[])) + sqrt(nrecs)*1i*imag(bspecxyx)/std(reshape(imag(bspecxyx),1,[]));
%                         bspecxxy = sqrt(nrecs)*real(bspecxxy)/std(reshape(real(bspecxxy),1,[])) + sqrt(nrecs)*1i*imag(bspecxxy)/std(reshape(imag(bspecxxy),1,[]));                        
%                     end
                    
%                     antiScBspec = abs(bspecxxy - bspecxyx);
                    antiScBspec = bspecxxy - bspecxyx;
                    if method == "Normalization"
                        antiScBspec = sqrt(nrecs).*real(antiScBspec)./std(reshape(real(antiScBspec),1,[])) + sqrt(nrecs)*1i.*imag(antiScBspec)./std(reshape(imag(antiScBspec),1,[]));
                    end
                    antiScBspec = abs(antiScBspec);
                    antiScBic = bicoherx_temp(xin,xin,yin,nfft,wind,nsamp,overlap);
                    
                    if method == "Constant"
                        per_thres = param_struct.per_thres;
                        antiScBspec(antiScBspec < per_thres*max(max(antiScBspec))) = 0;               
                    end 
                end
                
                step = Fs/nfft; %Frequency bin size
                start = floor((nfft/2)+1 + bands(1)/step);
                ende = floor(nfft/2+1+bands(end)/step);
                range = start:ende;
                %Crop the cBS and cBIC in the range of the specified bands
                if ifcondAS
                    antiScBspec = antiScBspec(range,range);
                    antiScBic = antiScBic(range,range);
                    %Find Bispectrum Peaks
                    [~,c] = contour(antiScBspec,4);
                    level_list = c.LevelList;
                    if level_list == 0
                        bspecNon0 = 0; %If Bspec is 0
                    else
                        bspecNon0 = 1;
                    end
                    if sum(contains(features,"BSnPks")) > 0 && bspecNon0
                        yellow = level_list(end);
                        [xy,yy] = find(antiScBspec>=yellow);
                        BSPeaks = findPeaks(xy, yy, antiScBspec);
                    end
                     %Find Bicoherence Peaks
                    [~,c] = contour(antiScBic,4);
                    level_list = c.LevelList;
                    if level_list == 0
                        bicNon0 = 0; %If Bspec is 0
                    else
                        bicNon0 = 1;
                    end
                    if sum(contains(features,"BICnPks")) > 0 && bicNon0
                        yellow = level_list(end);
                        [xy,yy] = find(antiScBic>=yellow);
                        BICPeaks = findPeaks(xy, yy, antiScBic);
                    end
                end
                %Find Bicoherence Peaks
%                 cBic = cBic(range,range);
%                 [~,c] = contour(cBic,4);
%                 level_list = c.LevelList;
%                 if level_list == 0
%                     bicNon0 = 0; %If Bspec is 0
%                 else
%                     bicNon0 = 1;
%                 end                
%                 if sum(contains(features,"BICnPks")) > 0 && bicNon0
%                     yellow = level_list(end);
%                     [xy,yy] = find(cBic>=yellow);
%                     BICPeaks = findPeaks(xy, yy, cBic);
%                 end
                for k = 1:2:length(bands) - 1
                    %Find frequencies that correspond to f2 = [bands(i),bands(i+1)]
                    ind1 = (k+1)/2;
                    if band_method == "conventional"
                        indn1 = conv_bands(ind1);
                    else
                        indn1 = string(ind1);
                    end
                    lowLim1 = bands(k)/step - bands(1)/step + 1;
                    upLim1 = bands(k+1)/step - bands(1)/step + 1;
                    lowLim1 = floor(lowLim1);
                    upLim1 = floor(upLim1);
                    if ifcondAS
                        bspecTemp1 = zeros(size(antiScBspec,1),size(antiScBspec,2));
                        bspecTemp1(lowLim1:upLim1,:) = antiScBspec(lowLim1:upLim1,:);
                        bicTemp1 = zeros(size(antiScBic,1),size(antiScBic,2));
                        bicTemp1(lowLim1:upLim1,:) = antiScBic(lowLim1:upLim1,:);
                    end
%                     bicTemp1 = zeros(size(cBic,1),size(cBic,2));
%                     bicTemp1(lowLim1:upLim1,:) = cBic(lowLim1:upLim1,:);
            %         figure('Visible','off')
            %         contour(waxis,waxis,abs(bspecTemp1),4)

                    for n = 1:2:length(bands) - 1
                        %Find frequencies that correspond to f1 = [bands(j),bands(j+1)]
                        ind2 = (n+1)/2;
                        if band_method == "conventional"
                            indn2 = conv_bands(ind2);
                        else
                            indn2 = string(ind2);
                        end
                        lowLim2 = bands(n)/step - bands(1)/step + 1;
                        upLim2 = bands(n+1)/step - bands(1)/step + 1;
                        lowLim2 = floor(lowLim2);
                        upLim2 = floor(upLim2);
                        if ifcondAS
                            bspecTemp2 = zeros(size(antiScBspec,1),size(antiScBspec,2));                        
                            %bspecTemp2 now is the region that corresponds to bandX -bandY
                            %coupling
                            bspecTemp2(:,lowLim2:upLim2) = bspecTemp1(:,lowLim2:upLim2);
                            bicTemp2 = zeros(size(antiScBic,1),size(antiScBic,2));
                            bicTemp2(:,lowLim2:upLim2) = bicTemp1(:,lowLim2:upLim2);
                        end
%                         bicTemp2 = zeros(size(cBic,1),size(cBic,2));
%                         bicTemp2(:,lowLim2:upLim2) = bicTemp1(:,lowLim2:upLim2);
        %                 figure('Visible','off')
        %                 contour(waxis,waxis,abs(bspecTemp2),4)

                        if sum(contains(features,"TotalBIC")) > 0 && ifcondAS
                            if ~exist('TotalBIC','var')
                                if bicNon0 == 1
                                    TotalBIC.(join(["TotalBIC_band",indn1,"_",indn2],'')) = sBic_calc(antiScBic,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    TotalBIC.(join(["TotalBIC_band",indn1,"_",indn2],'')) = 0;
                                end
                            elseif exist('TotalBIC','var') && isfield(TotalBIC,join(["TotalBIC_band",indn1,"_",indn2],''))
                                if bicNon0 == 1
                                    TotalBIC.(join(["TotalBIC_band",indn1,"_",indn2],''))(end+1) = sBic_calc(antiScBic,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    TotalBIC.(join(["TotalBIC_band",indn1,"_",indn2],''))(end+1) = 0;
                                end
                            else
                                if bicNon0 == 1
                                    TotalBIC.(join(["TotalBIC_band",indn1,"_",indn2],'')) = sBic_calc(antiScBic,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    TotalBIC.(join(["TotalBIC_band",indn1,"_",indn2],'')) = 0;
                                end
                            end
%                         else 
%                             warning("Bicoherence feature was not calculated. You should specify a feature name that contains Bic if you want to calculate it")
                        end

                        if sum(contains(features,"RelPow")) > 0 && ifcondAS
                            if ~exist('RelPow','var')
                                if bspecNon0 == 1
                                    RelPow.(join(["RelPow_band",indn1,"_",indn2],'')) = RelPow_calc(bspecTemp2,antiScBspec,symmetric);
                                else
                                    RelPow.(join(["RelPow_band",indn1,"_",indn2],'')) = 0;
                                end
                            elseif exist('RelPow','var') && isfield(RelPow,join(["RelPow_band",indn1,"_",indn2],''))
                                if bspecNon0 == 1
                                    RelPow.(join(["RelPow_band",indn1,"_",indn2],''))(end+1) = RelPow_calc(bspecTemp2,antiScBspec,symmetric);
                                else
                                    RelPow.(join(["RelPow_band",indn1,"_",indn2],''))(end+1) = 0;
                                end
                            else
                                if bspecNon0 == 1
                                    RelPow.(join(["RelPow_band",indn1,"_",indn2],'')) = RelPow_calc(bspecTemp2,antiScBspec,symmetric);
                                else
                                    RelPow.(join(["RelPow_band",indn1,"_",indn2],'')) = 0;
                                end
                            end
%                         else 
%                             warning("RelativePower feature was not calculated. You should specify a feature name that contains RelPow if you want to calculate it")
                        end

                        if sum(contains(features,"BSEnt")) > 0 && ifcondAS
                            if ~exist('BSEnt','var')
                                if bspecNon0 == 1
                                    BSEnt.(join(["BSEnt_band",indn1,"_",indn2],'')) = bspecEntropy(bspecTemp2,upLim1,upLim2,lowLim1,lowLim2,1);
                                else
                                    BSEnt.(join(["BSEnt_band",indn1,"_",indn2],'')) = 0;
                                end
                            elseif exist('BSEnt','var') && isfield(BSEnt,join(["BSEnt_band",indn1,"_",indn2],''))
                                if bspecNon0 == 1
                                    BSEnt.(join(["BSEnt_band",indn1,"_",indn2],''))(end+1) = bspecEntropy(bspecTemp2,upLim1,upLim2,lowLim1,lowLim2,1);
                                else
                                    BSEnt.(join(["BSEnt_band",indn1,"_",indn2],''))(end+1) = 0;
                                end
                            else
                                if bspecNon0 == 1
                                    BSEnt.(join(["BSEnt_band",indn1,"_",indn2],'')) = bspecEntropy(bspecTemp2,upLim1,upLim2,lowLim1,lowLim2,1);
                                else
                                    BSEnt.(join(["BSEnt_band",indn1,"_",indn2],'')) = 0;
                                end
                            end
%                         else 
%                             warning("BspecEntropy1 feature was not calculated. You should specify a feature name that contains Ent1 if you want to calculate it")
                        end

                        if sum(contains(features,"BICEnt")) > 0 && ifcondAS
                            if ~exist('BICEnt','var')
                                if bicNon0 == 1
                                    BICEnt.(join(["BICEnt_band",indn1,"_",indn2],'')) = bspecEntropy(bicTemp2,upLim1,upLim2,lowLim1,lowLim2,1);
                                else
                                    BICEnt.(join(["BICEnt_band",indn1,"_",indn2],'')) = 0;
                                end
                            elseif exist('BICEnt','var') && isfield(BICEnt,join(["BICEnt_band",indn1,"_",indn2],''))
                                if bicNon0 == 1
                                    BICEnt.(join(["BICEnt_band",indn1,"_",indn2],''))(end+1) = bspecEntropy(bicTemp2,upLim1,upLim2,lowLim1,lowLim2,1);
                                else
                                    BICEnt.(join(["BICEnt_band",indn1,"_",indn2],''))(end+1) = 0;
                                end
                            else
                                if bicNon0 == 1
                                    BICEnt.(join(["BICEnt_band",indn1,"_",indn2],'')) = bspecEntropy(bicTemp2,upLim1,upLim2,lowLim1,lowLim2,1);
                                else
                                    BICEnt.(join(["BICEnt_band",indn1,"_",indn2],'')) = 0;
                                end
                            end
%                         else 
%                             warning("BicEntropy feature was not calculated. You should specify a feature name that contains Ent2 if you want to calculate it")
                        end

                        if sum(contains(features,"BSnPks")) > 0 && ifcondAS
                            if ~exist('BSnPks','var')
                                if bspecNon0 == 1
                                    BSnPks.(join(["BSnPks_band",indn1,"_",indn2],'')) = NumPeaks(BSPeaks,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    BSnPks.(join(["BSnPks_band",indn1,"_",indn2],'')) = 0;
                                end
                            elseif exist('BSnPks','var') && isfield(BSnPks,join(["BSnPks_band",indn1,"_",indn2],''))
                                if bspecNon0 == 1
                                    BSnPks.(join(["BSnPks_band",indn1,"_",indn2],''))(end+1) = NumPeaks(BSPeaks,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    BSnPks.(join(["BSnPks_band",indn1,"_",indn2],''))(end+1) = 0;
                                end
                            else
                                if bspecNon0 == 1
                                    BSnPks.(join(["BSnPks_band",indn1,"_",indn2],'')) = NumPeaks(BSPeaks,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    BSnPks.(join(["BSnPks_band",indn1,"_",indn2],'')) = 0;
                                end
                            end
%                         else 
%                             warning("Bispectrum Number of peaks feature was not calculated. You should specify a feature name that contains Bic if you want to calculate it")
                        end

                        if sum(contains(features,"BICnPks")) > 0 && ifcondAS
                            if ~exist('BICnPks','var')
                                if bicNon0 == 1
                                    BICnPks.(join(["BICnPks_band",indn1,"_",indn2],'')) = NumPeaks(BICPeaks,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    BICnPks.(join(["BICnPks_band",indn1,"_",indn2],'')) = 0;
                                end
                            elseif exist('BICnPks','var') && isfield(BICnPks,join(["BICnPks_band",indn1,"_",indn2],''))
                                if bicNon0 == 1
                                    BICnPks.(join(["BICnPks_band",indn1,"_",indn2],''))(end+1) = NumPeaks(BICPeaks,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    BICnPks.(join(["BICnPks_band",indn1,"_",indn2],''))(end+1) = 0;
                                end
                            else
                                if bicNon0 == 1
                                    BICnPks.(join(["BICnPks_band",indn1,"_",indn2],'')) = NumPeaks(BICPeaks,upLim1,upLim2,lowLim1,lowLim2);
                                else
                                    BICnPks.(join(["BICnPks_band",indn1,"_",indn2],'')) = 0;
                                end
                            end
%                         else 
%                             warning("Bicoherence Number of peaks feature was not calculated. You should specify a feature name that contains Bic if you want to calculate it")
                        end
                    end
                end
            end
        end
    end

    out_names = [];
    if exist('TotalBIC','var')
        names = fieldnames(TotalBIC);
        out_names = [out_names;names];
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(TotalBIC.(name)))) >= 1
                error("TotalBIC contains NaN")
            else
                TotalBIC.(name) = mean(TotalBIC.(name));
            end
            mean_TotalBIC(i) = TotalBIC.(name);
        end
    elseif ~exist('mean_TotalBIC','var')
        mean_TotalBIC = [];
    end

    if exist('RelPow','var')
        names = fieldnames(RelPow);
        out_names = [out_names;names];
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(RelPow.(name)))) >= 1
                error("RelPow contains NaN")
            else
                RelPow.(name) = mean(RelPow.(name));
            end
            mean_RelPow(i) = RelPow.(name);
        end
    elseif ~exist('RelPow','var')
        mean_RelPow = [];
    end

    if exist('BSEnt','var')
        names = fieldnames(BSEnt);
        out_names = [out_names;names];
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(BSEnt.(name)))) >= 1
                error("BSEnt contains NaN")
            else
                BSEnt.(name) = mean(BSEnt.(name));
            end
            mean_BSEnt(i) = BSEnt.(name);
        end
    elseif ~exist('BSEnt','var')
        mean_BSEnt = [];
    end

    if exist('BICEnt','var')
        names = fieldnames(BICEnt);
        out_names = [out_names;names];  
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(BICEnt.(name)))) >= 1
                error("BICEnt contains NaN")
            else
                BICEnt.(name) = mean(BICEnt.(name));
            end
            mean_BICEnt(i) = BICEnt.(name);
        end
    elseif ~exist('BICEnt','var')
        mean_BICEnt = [];
    end

    if exist('BSnPks','var')
        names = fieldnames(BSnPks);
        out_names = [out_names;names];  
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(BSnPks.(name)))) >= 1
                error("BSnPks contains NaN")
            else
                BSnPks.(name) = mean(BSnPks.(name));
            end
            mean_BSnPks(i) = BSnPks.(name);
        end
    elseif ~exist('BSnPks','var')
        mean_BSnPks = [];
    end

    if exist('BICnPks','var')
        names = fieldnames(BICnPks);
        out_names = [out_names;names];  
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(BICnPks.(name)))) >= 1
                error("BICnPks contains NaN")
            else
                BICnPks.(name) = mean(BICnPks.(name));
            end
            mean_BICnPks(i) = BICnPks.(name);
        end
    elseif ~exist('BICnPks','var')
        mean_BICnPks = [];
    end

    if exist('mean_TotalBIC','var') || exist('mean_RelPow','var') || ...
            exist('mean_BSEnt','var') || exist('mean_BICEnt','var') || ...
            exist('mean_BSnPks','var') || exist('mean_BICnPks','var')
        mean_feat = [mean_TotalBIC, mean_RelPow, mean_BSEnt,...
                    mean_BICEnt,mean_BSnPks,mean_BICnPks];
    end
    if isempty(mean_feat)
        mean_feat = NaN(1,param_struct.bands_length*length(param_struct.features));
    end
end
    
end

function [reg_table] = gather_sigs(x)
    chans = fieldnames(x);
    count = 1;
    for t = 1:length(chans)    
        sig = unwrap_SwDs(x.(chans{t}));
        sig = cell2table(sig);
        idx = all(cellfun(@isempty,sig{:,:}),2);
        sig(idx,:)=[];
        sig = table2array(sig);
        sig = horzcat(sig{:});
        if ~isempty(sig)
            [lx,hx] = size(sig);
            if hx > lx
                sig = transpose(sig);
                warning("Rhythms should be in shape [samples,realizations]")
            end                              
            reg_table(:,count) = sum(sig,2);  
%             reg_chan_names(count) = chans{k};    
            count = count + 1;                       
        end
    end
end