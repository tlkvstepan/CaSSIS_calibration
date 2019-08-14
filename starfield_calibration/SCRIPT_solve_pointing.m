function SCRIPT_solve_pointing(set)
% Given rotation of FSA with respect to J2000 estimated from images we find
% two correction angles:
%
% (1) R_cru_2_cru_real 
% (2) R_fsa_2_fsa_real 
%
% Overall transformatin from J2000 to FSA looks like this:
%
% R_fsa_2_fsa_real*R_tel_2_fsa*R_cru_2_tel*R_cru_2_cru_real*R_sp_2_cru*R_j2000_2_sp

% clear all
% addpath(genpath('../libraries'));
% dataset_path =  '/HDD1/Data/CASSIS_STARFIELD_FIX';
% set = DATASET_starfields(dataset_path, 'combined');

cspice_furnsh(set.spice);

% read intrinsics
intrinsic = readtable(set.intrinsic_ba);
f  = intrinsic.f;
x0 = intrinsic.x0;
y0 = intrinsic.y0;
pixSize = intrinsic.pixSize;

% read real extrinsics
extrinsic_ba = readtable(set.extrinsic_ba);

for nexp = 1:height(extrinsic_ba)

    time_int = cassis_time2num( extrinsic_ba.time(nexp) );
    time_str= datestr(time_int,'HH:MM:SS.FFF yyyy dd mmm');
    time_ET = cspice_str2et(time_str); %any valid time will do
     
    % get R_j2000_2_fsa_real estimated from images
    q_j2000_2_fsa_real =quaternion( [extrinsic_ba.Q_1(nexp) extrinsic_ba.Q_2(nexp) extrinsic_ba.Q_3(nexp) extrinsic_ba.Q_4(nexp)] );
    R_j2000_2_fsa_real{nexp} = RotationMatrix( q_j2000_2_fsa_real );

    % get corresponding spice extrinsics:
    % R_j2000_2_cru = R_sc_2_cru * R_j2000_2_sc
    % R_cru_2_fsa 
    R_j2000_2_cru{nexp} = cspice_pxform('J2000', 'TGO_CASSIS_CRU',  time_ET);
    R_cru_2_fsa{nexp}     = cspice_pxform('TGO_CASSIS_CRU', 'TGO_CASSIS_FSA',  time_ET);
        
end

% get sequence id
seq_id = time2seqIndex( cassis_time2num(extrinsic_ba.time) );
unique_id = unique(seq_id);

% RANSAC
nransac = 200;
train_set_size = 3;
inlier_th = 0.26; % deg
best_inlier_prc =  -1;
for i = 1:nransac
    
    % select subset of sequences
    unique_id = unique_id( randperm(length(unique_id)) );
    train_ids = unique_id(1:train_set_size);
    
    % estimate parameters
    q_fsa_2_fsa_real = quaternion.rotationmatrix( eye(3,3) );
    q_cru_2_cru_real = quaternion.rotationmatrix( eye(3,3) );
    sol0 = [q_fsa_2_fsa_real.e(:); q_cru_2_cru_real.e(:)];
    fun = @(sol) clc_matrix_res(sol, seq_id, train_ids, R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa);
    options = optimoptions('lsqnonlin', 'Algorithm',  'levenberg-marquardt', 'Display', 'off');
    [sol, ~, res] = lsqnonlin(fun, sol0, [], [], options);
    R_fsa_2_fsa_real  = RotationMatrix( quaternion( sol(1:4) ));
    R_cru_2_cru_real  = RotationMatrix( quaternion( sol(5:end) ));
    
    % find outliers
     [inlier_mask, inlier_nb] = find_inliers_nb(sol,  R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa, inlier_th);
     inlier_prc = inlier_nb/length(seq_id)*100;
     if best_inlier_prc < inlier_prc
         best_inlier_prc = inlier_prc;
         best_inlier_mask = inlier_mask;
         best_R_fsa_2_fsa_real = R_fsa_2_fsa_real;
         best_R_cru_2_cru_real = R_cru_2_cru_real;
     end    
     fprintf('RANSAC iteraton %i, current inliers %2.2f , best inliers %2.2f \n', i, inlier_prc, best_inlier_prc);
         
end

% fit again using all inliers
fprintf('Perform final fit using all inliers \n');
train_ids = unique(seq_id(best_inlier_mask));
q_fsa_2_fsa_real = quaternion.rotationmatrix( best_R_fsa_2_fsa_real );
q_cru_2_cru_real = quaternion.rotationmatrix( best_R_cru_2_cru_real );
sol0 = [q_fsa_2_fsa_real.e(:); q_cru_2_cru_real.e(:)];
fun = @(sol) clc_matrix_res(sol, seq_id, train_ids, R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa);
options = optimoptions('lsqnonlin', 'Algorithm',  'levenberg-marquardt', 'Display', 'off');
[sol, ~, res] = lsqnonlin(fun, sol0, [], [], options);
R_fsa_2_fsa_real  = RotationMatrix( quaternion( sol(1:4) ));
R_cru_2_cru_real  = RotationMatrix( quaternion( sol(5:end) ));
    
% compute angular error on inliers before and after optimization
angle_err_after = clc_angular_res(sol, seq_id, train_ids, R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa);
q_fsa_2_fsa_real = quaternion.rotationmatrix( eye(3,3) );
q_cru_2_cru_real = quaternion.rotationmatrix( eye(3,3) );
sol = [q_fsa_2_fsa_real.e(:); q_cru_2_cru_real.e(:)];
angle_err_before = clc_angular_res(sol, seq_id, train_ids, R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa);
     
 
fprintf('Absolute angular error BEFORE OPTIMIZATION mean = %0.3f [deg] \n', mean(angle_err_before));
fprintf('Absolute angular error AFTER OPTIMIZATION mean = %0.3f [deg] \n', mean(angle_err_after));

% show angles before and after calibration
time_str_spice = datestr(time_int,'HH:MM:SS.FFF yyyy dd mmm');
time_ET = cspice_str2et(time_str_spice); %any valid time will do
R_tel_2_fsa_real= R_fsa_2_fsa_real*cspice_pxform('TGO_CASSIS_TEL', 'TGO_CASSIS_FSA',  time_ET);
axis3 =2; axis2 =1 ; axis1 =3; 
[angle3, angle2, angle1] = cspice_m2eul(R_tel_2_fsa_real', axis3, axis2, axis1);

fprintf('TEL2FSA axis= (%i, %i, %i), angles = (%f, %f, %f) \n', axis3,axis2,axis1, rad2deg(angle3), rad2deg(angle2), rad2deg(angle1));

R_sc_2_cru_real = R_cru_2_cru_real*cspice_pxform('TGO_SPACECRAFT', 'TGO_CASSIS_CRU',  time_ET);
axis3 =1; axis2 =2 ; axis1 =3; 
[angle3, angle2, angle1] = cspice_m2eul(R_sc_2_cru_real', axis3, axis2, axis1);

fprintf('SC2CRU axis= (%i, %i, %i), angles = (%f, %f, %f) \n', axis3,axis2,axis1, rad2deg(angle3), rad2deg(angle2), rad2deg(angle1));

% save calibaration results
extrinsic=table( [R_tel_2_fsa_real(:)'; R_sc_2_cru_real(:)']);
writetable(extrinsic, set.extrinsic); 

end


function [inlier_mask, inlier_nb] = find_inliers_nb(sol,  R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa, inlier_th)
  
    nb = length(R_j2000_2_fsa_real);
    inlier_mask = true(nb,1);
    angle_err = clc_angular_res(sol, ones(nb,1), ones(nb,1), R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa);
    inlier_mask = (angle_err < inlier_th);
    inlier_nb = sum(inlier_mask);

end

function err = clc_angular_res(sol, seq_id, train_ids, R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa)

      R_fsa_2_fsa_real  = RotationMatrix( quaternion( sol(1:4) ));
      R_cru_2_cru_real  = RotationMatrix( quaternion( sol(5:end) ));
      
      err = [];
      for n = 1:length(R_j2000_2_fsa_real)
            if any(train_ids == seq_id(n)) 
                R_j2000_2_fsa = R_fsa_2_fsa_real*R_cru_2_fsa{n}*R_cru_2_cru_real*R_j2000_2_cru{n}; 
                q_spice = quaternion.rotationmatrix( R_j2000_2_fsa  );
                q_real = quaternion.rotationmatrix( R_j2000_2_fsa_real{n} ); 
                q_err = ldivide(q_real, q_spice);
                [angle_error, ~] = AngleAxis(q_err); 
                angle_error = rad2deg(angle_error);   
                angle_error = min( min((angle_error), (360-angle_error)), (360+angle_error));
                err = [err; angle_error];
            end
        end
  

end

function err = clc_matrix_res(sol, seq_id, train_ids, R_j2000_2_fsa_real, R_j2000_2_cru, R_cru_2_fsa)
  
    R_fsa_2_fsa_real  = RotationMatrix( quaternion( sol(1:4) ));
    R_cru_2_cru_real  = RotationMatrix( quaternion( sol(5:end) ));
    
    err = [];
    for n = 1:length(R_j2000_2_fsa_real)
        if any(train_ids == seq_id(n)) 
            R_j2000_2_fsa = R_fsa_2_fsa_real*R_cru_2_fsa{n}*R_cru_2_cru_real*R_j2000_2_cru{n}; 
            tmp= (eye(3) - R_j2000_2_fsa_real{n}'*R_j2000_2_fsa);
            err = [err; tmp(:)];
        end
        
    end
        
end
