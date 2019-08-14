% Script recogizes stars stars

function SCRIPT_recognize(set)

%%

addpath(genpath('../libraries'));


%%
clc
fprintf('Detecting and recognising stars\n');

% read exposures summary
expSummary = readtable(set.exposuresSummary);
nb_exp = height(expSummary);

%-------------------------------------------------------------------------
for nexp = 1:nb_exp
    
    fprintf('%s...\n', expSummary.fname_exp{nexp});
        
    % copy image to temp folder
    delete('work/1st/*'); 
    delete('work/2nd/*'); 
    fname = [set.denoise_exposure '/' expSummary.fname_exp{nexp}];
    copyfile(fname, 'work/1st/tmp.tif');
    copyfile(fname, 'work/2nd/tmp.tif');
    I = imread(fname);

    imwrite(I, 'work/1st/tmp.tif');
    
    % detect stars
    sys_command = ['solve-field --parity pos   --obj 100 --scale-low 1  --pixel-error 15 --scale-high 1.5 --no-plots --cpulimit  300 --downsample 2  work/1st/tmp.tif'];
    system(sys_command);
    
    % rerun to find more matches
    if exist('work/1st/tmp.wcs','file' )== 2
        copyfile('work/1st/tmp.wcs', 'work/2nd/1st_attemt.wcs');
        sys_command = ['solve-field --parity pos    --no-plots --pixel-error 15  --verify work/2nd/1st_attemt.wcs work/2nd/tmp.tif'];
        system(sys_command);
    end
       
    % parse output
    match_list = [];
    x = [];
    y = []; % since we flipped image
    ra = [];
    dec = [];
    flux=[];
    
    if exist('work/2nd/tmp.corr', 'file') == 2
    
        info = fitsinfo('work/2nd/tmp.corr');
        tableData = fitsread('work/2nd/tmp.corr','binarytable');
        flux = tableData{12};
        x = tableData{1};
        y = 2048 - tableData{2}; % since we flipped image
        ra = tableData{7};
        dec = tableData{8};
        
        shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor', 'white')
        I = step(shapeInserter, cat(3,I,I,I), int32([x(:) 2048-y(:) 5*ones(numel(x),1)])); 
        
        %valid_exp_list = [valid_exp_list; nexp];
    
        imwrite(I,[set.recognize '/' expSummary.fname_exp{nexp}])
        fprintf(' %i stars matched \n', size(match_list,1));
        pause(0.1)
        
        % save stars
        fname = [set.recognize '/'  expSummary.fname_exp{nexp} '.csv'];
        starSummary = table(x, y, ra, dec, flux);
        writetable(starSummary, fname); 
        
    end
end


