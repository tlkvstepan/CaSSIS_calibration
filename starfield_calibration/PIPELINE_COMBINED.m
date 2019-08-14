clear all;
diary('hist')

%% PARAMETERS ESTIMATION
addpath(genpath('../libraries'));
dataset_path = '/HDD1/Data/CASSIS_STARFIELD';
dataset = 'training';

set = DATASET_starfields(dataset_path, dataset);
SCRIPT_search_folder(set);
SCRIPT_collect_sequences(set);

prm.skip_first_exp = 0;                   
prm.adjust_subExp_on = false;       
prm.virtualImage_on = true;            
prm.min_exposure_filter  = 1;
prm.max_exposure_filter = inf;    

SCRIPT_save_rawExp(setdataset = 'commissioning_2'
set = DATASET_starfields(dataset_path, dataset);
, prm);    
SCRIPT_denoise(set);
SCRIPT_recognize(set);
SCRIPT_collect_star(set);
SCRIPT_filter_outliers(set);
SCRIPT_init_intrinsic(set);
SCRIPT_init_lensDistortion(set);
SCRIPT_init_extrinsic_spice(set);
SCRIPT_find_rotCommand_spice(set);
SCRIPT_init_extrinsic_local(set);

iter = 1;
 while SCRIPT_bundle_adjustment(set,iter) > 0
     iter = iter + 1;
     pause(0.1);
 end

SCRIPT_solve_distortion(set);
SCRIPT_solve_pointing(set)

copyfile(set.intrinsic_ba, set.intrinsic_final)
copyfile(set.extrinsic, set.extrinsic_final)
copyfile(set.lensDistortion, set.lensDistortion_final)
copyfile(set.lensCorrection, set.lensCorrection_final )

%% VALIDATION
dataset = 'commissioning_2'
set = DATASET_starfields(dataset_path, dataset);

SCRIPT_search_folder(set);
SCRIPT_collect_sequences(set);

prm.skip_first_exp = 0;                   % how many first exposures to skip
prm.adjust_subExp_on = false;       % adjust intensity of sub exposures
prm.virtualImage_on = true;            % save virtual image
prm.min_exposure_filter  = 1;
prm.max_exposure_filter = inf;    

SCRIPT_save_rawExp(set, prm);    

SCRIPT_denoise(set);

SCRIPT_recognize(set);

SCRIPT_collect_star(set);

SCRIPT_filter_outliers(set);

SCRIPT_init_intrinsic(set);
       
SCRIPT_init_lensDistortion(set);
 
SCRIPT_init_extrinsic_spice(set);
 
SCRIPT_find_rotCommand_spice(set);
 
SCRIPT_init_extrinsic_local(set);

iter = 1;
 while SCRIPT_bundle_adjustment(set,iter) > 0
     iter = iter + 1;
     pause(0.1);
 end
ste
