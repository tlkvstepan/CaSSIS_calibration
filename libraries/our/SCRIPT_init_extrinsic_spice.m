function SCRIPT_init_extrinsic_spice(set)

%%

%dataset_path = '/home/tulyakov/Desktop/espace-server';
%if ~exist('dataset_name','var')
%    dataset_name = 'pointing_cassis';
%end
addpath(genpath('../libraries'));
min_points_per_image = 10;

%%
clc
fprintf('Initializing rotation from SPICE kernels\n');

% read exposures summary
expSummary = readtable(set.exposuresSummary);
nb_exp = height(expSummary);

% read inlier stars
inlierStarSummary = readtable(set.inlierStarSummary);
nb_points = height(inlierStarSummary);

% initialize spice
cspice_furnsh(set.spice);

% get angles
f = figure;
nexp_valid = 1;
for nexp = 1:nb_exp

    fprintf('%s\n', expSummary.exp_time{nexp});
            
    % convert time 
    t_str_cassis = expSummary.exp_time{nexp};
    t = cassis_time2num(t_str_cassis);
    time_str_spice = datestr(t,'HH:MM:SS.FFF yyyy dd mmm');
    time_ET = cspice_str2et(time_str_spice); % time in seconds from J2000
     
    % extract rotational matrices from SPICE kernel
    R_fsa_2_j2000 = cspice_pxform('J2000', 'TGO_CASSIS_FSA', time_ET);
    q = quaternion.rotationmatrix( R_fsa_2_j2000 );
    Q_fsa_2_j2000(nexp_valid,:) = q.e;
    time{nexp_valid} = t_str_cassis;
    
    nexp_valid = nexp_valid + 1; 
    
end

time = time';
Q = Q_fsa_2_j2000;
extrinsic0 = table(Q, time);

% remove frames that do not contain enougth stars
unique_time = unique(cassis_time2num(extrinsic0.time)) ;
nb_time = length(unique_time);
for ntime = 1:nb_time
    mask_extrinsic = unique_time(ntime) == cassis_time2num(extrinsic0.time); 
    mask_stars = unique_time(ntime) == cassis_time2num(inlierStarSummary.time);
    %nnz(mask_stars)
    if( nnz(mask_stars) < min_points_per_image )
        inlierStarSummary= inlierStarSummary(~mask_stars,:);
        extrinsic0 = extrinsic0(~mask_extrinsic,:);
    end
end

writetable(extrinsic0, set.extrinsic0_spice); 

writetable(inlierStarSummary, set.inlierStarSummary);


