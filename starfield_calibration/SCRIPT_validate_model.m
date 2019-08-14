function err =SCRIPT_validate_model(set, extrinsic_calibration_mode, is_intrinsic_calibration, is_distortion_calibration)
% Returns error on validation set.

% Args:
% extrinsic_calibration_mode: yes, no, image. 
% is_intrinsic_calibration: true / false.
% is_distortion_calibration: true / false.
addpath(genpath('../libraries'));
image_size = [2048 2048];

% initialize spice
cspice_furnsh(set.spice);

%% read stars 
starSummary = readtable(set.inlierStarSummary);
[XX(:,1), XX(:,2), XX(:,3)] = ...
raDec2XYZ(deg2rad(starSummary.ra), deg2rad(starSummary.dec));    
xx(:,1) = starSummary.x;
xx(:,2) = starSummary.y; 
nb_points = height(starSummary);

%% read extirinsic
% git
% Q = [extrinsic.Q_1 extrinsic.Q_2 extrinsic.Q_3 extrinsic.Q_4];
% for npoint = 1:nb_points
%     idx = find(cassis_time2num(starSummary.time(npoint)) == cassis_time2num(extrinsic.time));
%     q(:,npoint) = quaternion(Q(idx,:));
%     R(:,:,npoint) = RotationMatrix(q(:,npoint));
% end
if strcmp(extrinsic_calibration_mode, 'image')
    % Use image base extrinsics.
    extrinsic = readtable(set.extrinsic_ba); 
    Q = [extrinsic.Q_1 extrinsic.Q_2 extrinsic.Q_3 extrinsic.Q_4];
    for npoint = 1:nb_points
        idx = find(cassis_time2num(starSummary.time(npoint)) == cassis_time2num(extrinsic.time));
        q(:,npoint) = quaternion(Q(idx,:));
        R{npoint} = RotationMatrix(q(:,npoint));
    end
else
    if  strcmp(extrinsic_calibration_mode, 'no')
        % Use improved corrected spice kernel.
        time_int = cassis_time2num( starSummary.time(1) );
        time_str= datestr(time_int,'HH:MM:SS.FFF yyyy dd mmm');
        time_ET = cspice_str2et(time_str); 
        R_sc_2_cru =  cspice_pxform('TGO_SPACECRAFT', 'TGO_CASSIS_CRU',  time_ET);
        R_tel_2_fsa = cspice_pxform('TGO_CASSIS_TEL', 'TGO_CASSIS_FSA',  time_ET);
    else
        % Use original spice kernels.
        extrinsic = readtable(set.extrinsic_final);
        R_tel_2_fsa = reshape(table2array(extrinsic(1,:)), 3, 3);
        R_sc_2_cru = reshape(table2array(extrinsic(2,:)), 3, 3);
    end
    for npoint = 1:nb_points
        time_int = cassis_time2num( starSummary.time(npoint) );
        time_str= datestr(time_int,'HH:MM:SS.FFF yyyy dd mmm');
        time_ET = cspice_str2et(time_str); %any valid time will do
        R_cru_2_tel = cspice_pxform('TGO_CASSIS_CRU', 'TGO_CASSIS_TEL',  time_ET);
        R_j2000_2_sc = cspice_pxform('J2000', 'TGO_SPACECRAFT',  time_ET);
        R{npoint}= R_tel_2_fsa*R_cru_2_tel*R_sc_2_cru*R_j2000_2_sc;
    end
end

%% read intrinsics
if  ~is_intrinsic_calibration
    intrinsic = readtable(set.intrinsic0);
else
    intrinsic = readtable(set.intrinsic_final);
end
x0 = intrinsic.x0;
y0 = intrinsic.y0;
pixSize = intrinsic.pixSize;
f = intrinsic.f;
K = f_x0_y0_2K(f, x0, y0, pixSize);

%% read lens distortion model
if ~is_distortion_calibration 
    lensCorrection = readtable(set.lensCorrection0);
else
    lensCorrection = readtable(set.lensCorrection_final);
end
A_corr = [lensCorrection.A_corr_1 lensCorrection.A_corr_2 lensCorrection.A_corr_3 lensCorrection.A_corr_4 lensCorrection.A_corr_5 lensCorrection.A_corr_6];
[xx_norm(:,1), xx_norm(:,2)] = cassis_detector2focalplane(xx(:,1), xx(:,2), image_size(2), image_size(1), pixSize*1000);
chi = lift2D_to_6D(xx_norm); 
xx_corr_norm = chi*A_corr';
xx_corr_norm(:,1) = xx_corr_norm(:,1)./xx_corr_norm(:,3);
xx_corr_norm(:,2) = xx_corr_norm(:,2)./xx_corr_norm(:,3);
xx_corr_norm = xx_corr_norm(:,[1 2]);
[xx_corr(:,1), xx_corr(:,2)] = cassis_focalplane2detector(xx_corr_norm(:,1), xx_corr_norm(:,2), image_size(2), image_size(1), pixSize*1000);

%% compute error
for npoint = 1:nb_points 
    res_err(npoint,:) = stars2image_error( XX(npoint,:), xx_corr(npoint,:),  R{npoint}, K);
end


err = mean(sqrt(sum(res_err.^2, 2)));
    
end



