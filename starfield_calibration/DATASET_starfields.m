function set = DATASET_starfields(dataset_path, subset_name)

    set.name = subset_name;
        
    % Function returns all folder pertained to CaSSIS starfield dataset
	 warning ('off','all'); % to prevent warnings on folder creation
    
	if( strcmp('commissioning_2', subset_name) )
        
        subsetpath = [dataset_path '/cruise/160407_commissioning_2'];
        
    elseif( strcmp('pointing_cassis', subset_name) )
      
        subsetpath = [dataset_path '/cruise/160413_pointing_cassis'];
   
    elseif( strcmp('mcc_motor', subset_name) )
   
        subsetpath = [dataset_path '/cruise/160614_mcc_motor'];
   
    elseif( strcmp('stellar_cal_orbit10', subset_name) )
   
        subsetpath = [dataset_path '/aerobraking/161124_stellar_cal_orbit10'];
   
    elseif( strcmp('stellar_cal_orbit09', subset_name) )
   
        subsetpath = [dataset_path '/aerobraking/161120_stellar_cal_orbit09'];
   
    elseif( strcmp('training', subset_name) )
           
        subsetpath = [dataset_path '/tests/training'];
   
    else
        error('No set with such name');
    end
    
    % root
    set.root = subsetpath;
    
        
    % ----------- input -----------------
   
    % framelets and xmls
    set.level0 = [subsetpath '/level0b'];    
    
    % kernel
    set.spice = '/HDD1/Data/KERNELS/mk/em16_ops.tm';
        
    % ----------- output --------------
    
    mkdir(subsetpath,  'OUTPUT');
        
    % combined frames
    set.sequences = [subsetpath '/OUTPUT/sequences'];
    mkdir([subsetpath '/OUTPUT/'], 'sequences');
    
    set.raw_exposures = [subsetpath '/OUTPUT/raw_exposures'];
    mkdir([subsetpath '/OUTPUT/'], 'raw_exposures');
    
    % denoised exposures
    set.denoise_exposure = [subsetpath '/OUTPUT/denoise_exposures'];
    mkdir([subsetpath '/OUTPUT/'],  'denoise_exposures');
    
    % recognized stars
    set.recognize = [subsetpath '/OUTPUT/recognize']; 
    mkdir([subsetpath '/OUTPUT/'], 'recognize');
    
    % ----------- summaries -------------

    % factory parameters
    set.intrinsic0 = [subsetpath '/OUTPUT/intrinsic0.csv'];
    set.lensCorrection0 = [subsetpath '/OUTPUT/lensCorrection0.csv'];
    set.extrinsic0_spice = [subsetpath '/OUTPUT/extrinsic0_spice.csv'];
    
    % individual rotation angle tuning
    set.extrinsic0_local = [subsetpath '/OUTPUT/extrinsic0_local.csv'];
    set.local_vs_spice_angel_diff_IMG = [subsetpath '/OUTPUT/local_vs_spice_angle_diff.png'];
    set.local_vs_spice_proj_err_IMstellar_cal_orbit09G = [subsetpath '/OUTPUT/local_vs_spice_proj_err.png'];
    set.ange_vs_angel_diff_IMG = [subsetpath '/OUTPUT/ange_vs_angel_diff.png'];
        
    % BA parameters
    set.extrinsic_ba = [subsetpath '/OUTPUT/extrinsic_ba.csv'];
    set.intrinsic_ba = [subsetpath '/OUTPUT/intrinsic_ba.csv'];
    set.ba_residuals_IMG = [subsetpath '/OUTPUT/ba_residuals%i.png']; 
    
        
    % commanded rotation
    set.extrinsic = [subsetpath '/OUTPUT/extrinsic.csv'];
    set.rotCommand = [subsetpath '/OUTPUT/rotCommand.csv'];
        
    % lens distortion estimation
    set.lensDistortion = [subsetpath '/OUTPUT/lensDistortion.csv'];
    set.lensCorrection = [subsetpath '/OUTPUT/lensCorrection.csv'];
    set.lensCorrection_field_IMG = [subsetpath '/OUTPUT/lensCorrection_field.png'];
    set.lensCorrection_empirical_field_IMG = [subsetpath '/OUTPUT/lensCorrection_empirical_field.png'];
    set.lensCorrection_residual_field_IMG = [subsetpath '/OUTPUT/lensCorrection_residual_field.png'];
    set.lensCorrection_residuals_IMG = [subsetpath '/OUTPUT/lensCorrection_residuals.png'];
    
    
    % pointining tuining
    set.pointing_error_IMG = [subsetpath '/OUTPUT/pointing_error.png'];
    set.pointing_correction_IMG = [subsetpath '/OUTPUT/pointing_error.png'];
    
    % final parameters
    set.intrinsic_final = 'OUTPUT/intrinsic_final.csv';
    set.lensCorrection_final ='OUTPUT/lensCorrection_final.csv';
    set.lensDistortion_final ='OUTPUT/lensDistortion_final.csv';
    set.extrinsic_final = 'OUTPUT/extrinsic_final.csv';
    mkdir('OUTPUT');
    
    % folder content summary
    set.folderContent = [subsetpath '/OUTPUT/folderContent.csv']; 
    
    % sequences summary
    set.sequencesSummary = [subsetpath '/OUTPUT/sequencesSummary.csv']; 
        
    % exposure summary
    set.exposuresSummary = [subsetpath '/OUTPUT/exposuresSummary.csv']; 
    
    % matched stars summary
    set.allStarSummary = [subsetpath '/OUTPUT/allStarSummary.csv']; 
    
    % star filtering (using brightness, 
    set.inlierStarSummary = [subsetpath '/OUTPUT/inlierStarSummary.csv']; 
    set.inlierStarSummary_IMG = [subsetpath '/OUTPUT/allStarSummary.png']; 
    
    
     warning ('on','all');
end
