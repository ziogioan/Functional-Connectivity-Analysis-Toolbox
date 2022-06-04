function EEG_reref = rereference(EEG, chancoords,type)
AA = EEG.data(20,:);
if length(EEG.data(:,1)) == 21
    Label = EEG.data(21,:);
end
if type == "REST"
    for nc = 1:19
        xyz_elec(nc,1) = chancoords(nc).X;
        xyz_elec(nc,2) = chancoords(nc).Y;
        xyz_elec(nc,3) = chancoords(nc).Z;
    end
    xyz_dipoles = load(['corti869-3000dipoles.dat']);
    % Calculate the dipole orientations.
    xyz_dipOri = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
    xyz_dipOri ( 2601: 3000, 1 ) = 0;
    xyz_dipOri ( 2601: 3000, 2 ) = 0;
    xyz_dipOri ( 2601: 3000, 3 ) = 1;
    % define headmodel
    headmodel        = [];
    headmodel.type   = 'concentricspheres';
    headmodel.o      = [ 0.0000 0.0000 0.0000 ];
    headmodel.r      = [ 0.8700,0.9200,1];
    headmodel.cond   = [ 1.0000,0.0125,1];
    headmodel.tissue = { 'brain' 'skull' 'scalp' };
    [G,~] = dong_calc_leadfield3(xyz_elec,xyz_dipoles,xyz_dipOri,headmodel);
    G = G';
    data_z = rest_refer(EEG.data(1:19,:),G);
    data_z = [data_z; AA];
    if length(EEG.data(:,1)) == 21
        [data_z] = [data_z; Label];
    end
    EEG_reref = EEG;
    EEG_reref.data = data_z;
    EEG_reref.ref = 'REST';
elseif type == "Median"
    med = median(EEG.data(1:19,:));
    data_z = EEG.data(1:19,:) - med;
    data_z = [data_z; AA];
    if length(EEG.data(:,1)) == 21
        [data_z] = [data_z; Label];
    end
    EEG_reref = EEG;
    EEG_reref.data = data_z;
    EEG_reref.ref = 'Median';
end        

end

