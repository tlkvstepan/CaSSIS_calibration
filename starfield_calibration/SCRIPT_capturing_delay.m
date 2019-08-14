% Given star matches, camera intrinsics and extrinsics 

% 1. estimates pointing error mean and standard deviation in angles 
% 2. estimates pointing error deviation inside a sequence.

clear all;


addpath(genpath('../libraries'));
dataset_path =  '/HDD1/Data/cassis_starfield';
set = DATASET_starfields(dataset_path, 'stellar_cal_combined');

%% read all

cspice_furnsh(set.spice);

% rotation commands
rotCommands = readtable(set.rotCommand);

% predicted rotation
extrinsic_spice = readtable(set.extrinsic0_spice);

% read intrinsics
intrinsic = readtable(set.intrinsic_ba);
f  = intrinsic.f;
x0 = intrinsic.x0;
y0 = intrinsic.y0;
pixSize = intrinsic.pixSize;

% read extrinsics
extrinsic_ba = readtable(set.extrinsic_ba);

%% clean up
% these are outliers
%xtrinsic_spice([ 73,    74,    75,    76,    80],:)=[] ;

%% analyse pointing error
% find angular pointing error
nb_exp = height(extrinsic_ba);
for nexp = 1:nb_exp

    time_str = extrinsic_ba.time(nexp);
    time{nexp} =  time_str{1};
    idx = cassis_time2num(extrinsic_spice.time) ==  cassis_time2num(time_str);
    
    Q_real = [extrinsic_ba.Q_1(nexp) extrinsic_ba.Q_2(nexp) extrinsic_ba.Q_3(nexp) extrinsic_ba.Q_4(nexp)];
    q_real = quaternion( Q_real );
    R_real = RotationMatrix(q_real);
        
    Q_com = [extrinsic_spice.Q_1(idx) extrinsic_spice.Q_2(idx) extrinsic_spice.Q_3(idx) extrinsic_spice.Q_4(idx)];
    q_com = quaternion( Q_com);
    R_com = RotationMatrix(q_com);
    
    
    % axis-angle representation of the error
    q_err = ldivide(q_real,q_com);
    [angle_diff(nexp), axis(:,nexp)] = AngleAxis(q_err); 
    angle_diff(nexp) = rad2deg(angle_diff(nexp));   
    angle_diff(nexp) = min( min((angle_diff(nexp)), (360-angle_diff(nexp))), (360+angle_diff(nexp)));
  
    
    if nexp == 1
     time_str_spice = datestr(cassis_time2num(time_str),0);
     time_ET = cspice_str2et(time_str_spice); % time in seconds from J2000
     R_fsa_2_cru = cspice_pxform('TGO_CASSIS_FSA', 'TGO_CASSIS_TEL', time_ET);
    end
    axis_rot(:,nexp) = R_fsa_2_cru*axis(:,nexp)

     
    % 
    %[state, ~]=cspice_spkezr ( 'TGO_SPACECRAFT', time_ET, 'J2000', 'NONE', 'EARTH');
    
    %$q_err = ldivide(q_com, q_real);
    %[angle_diff(nexp), axis] = AngleAxis(q_err); 
    
    % to rotation matrix
    %R = RotationMatrix(q_err);
    %R_from_j2000_2_cru*R
    
    % commanded angle
    idx = cassis_time2num(rotCommands.time) == cassis_time2num(time_str);
    angle(nexp) =  rotCommands.angle(idx);    
    
end

% take care of outliers
outlier_th = 2;
time = time(angle_diff < outlier_th);
angle = angle(angle_diff < outlier_th);
angle_diff = angle_diff(angle_diff < outlier_th);


pointing = table(time', angle', angle_diff');
writetable(pointing, 'pointing.csv'); 

for i = 1:length(time)
    fprintf('%i %s\n',i, time{i});
end


seqids= time2seqIndex( cassis_time2num(time) );
for seqid = unique(seqids)
      mask = seqids == seqid;
      if nnz(mask) >= 5
           seq_angle_sig(seqid)  = std( angle_diff(mask) )
      end
end

% robust angular mean and std
angle_mu = mean(angle_diff);
angle_sig = std(angle_diff);

% mean and std of pointing error 
fprintf('Absolute angular error mean = %0.3f [deg], and deviation = %0.3f [deg] \n', angle_mu, angle_sig)

% std of pointing error inside sequence
fprintf('Mean deviation of pointing error inside sequences = %0.5f [deg] \n', mean(seq_angle_sig))

% plot error for sequence
% that have more than 5 measurments
colors = lines(length(unique(seqids)))
indexes = 1:length(angle_diff);
plot(angle_diff);
hold on;
for seqid = unique(seqids)
      mask = seqids == seqid;
      if nnz(mask) >=5
           seq_angle_sig(seqid)  = robustcov( angle_diff(mask) )
           plot( indexes(mask), angle_diff(mask), 'o', 'MarkerFaceColor', colors(seqid,:));
           hold on;
      end
end
ylabel('Pointing error, [deg]');

