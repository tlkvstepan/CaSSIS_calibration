function SCRIPT_solve_distortion(set)
% Given star matches, camera intrinsics and extrinsics SCRIPT estimates rational
% distortion and undistortion models.

%dataset_path = '/home/tulyakov/Desktop/espace-server';
%dataset_name = 'pointing_cassis';
addpath(genpath('../libraries'));
image_size = [2048 2048];

%%
clc;
fprintf('Estimating rational lens distortion parameters\n');

% read folders structure
%set = DATASET_starfields(dataset_path, dataset_name);

% read stars 
starSummary = readtable(set.inlierStarSummary);
nb_points = height(starSummary);

% read lens distortions
lensCorrection0 = readtable(set.lensCorrection0);
A0_corr = [lensCorrection0.A_corr_1, lensCorrection0.A_corr_2, lensCorrection0.A_corr_3, lensCorrection0.A_corr_4, lensCorrection0.A_corr_5, lensCorrection0.A_corr_6];

% read intrinsics
intrinsic = readtable(set.intrinsic_ba);
x0 = intrinsic.x0;
y0 = intrinsic.y0;
pixSize = intrinsic.pixSize;
f  = intrinsic.f;
K = f_x0_y0_2K(f, x0, y0, pixSize);

% read extirinsic
extrinsic = readtable(set.extrinsic_ba);
nb_exp = height(extrinsic);

%% Apply estimated camera model

% prepare points
[XX(:,1), XX(:,2), XX(:,3)] = ...
raDec2XYZ(deg2rad(starSummary.ra), deg2rad(starSummary.dec));    
xx(:,1) = starSummary.x;
xx(:,2) = starSummary.y; 
nb_points = height(starSummary);

% prepare rotation matrices
Qvec = [extrinsic.Q_1 extrinsic.Q_2 extrinsic.Q_3 extrinsic.Q_4];
Q = quaternion(Qvec);
for n = 1:nb_exp
    R(:,:,n) = RotationMatrix(Q(n));
end

% prepare point-rotation matrix indexing for speed
for n = 1:height(starSummary)
    rotIdx(n) = find(cassis_time2num(starSummary.time(n)) == cassis_time2num(extrinsic.time));
end

% map stars to sensor plane
for n = 1:nb_points
    [res(n,:), xx_pred(n,:)] = stars2image_error( XX(n,:), xx(n,:), R(:,:,rotIdx(n)), K);
end

% show residual error vectorfield
[x_query, y_query] = meshgrid(0:50:2048, 0:50:2048);
dx_r2i = xx_pred(:,1) - xx(:,1);
dy_r2i = xx_pred(:,2) - xx(:,2);
dx_r2i_grid = griddata(xx(:,1), xx(:,2), dx_r2i, x_query(:), y_query(:));
dy_r2i_grid = griddata(xx(:,1), xx(:,2), dy_r2i, x_query(:), y_query(:));
dx_r2i_grid(dx_r2i_grid == NaN) = 0;
dy_r2i_grid(dy_r2i_grid == NaN) = 0;
dx_r2i_grid = medfilt2(reshape(dx_r2i_grid, size(x_query)) , [7 7]);
dy_r2i_grid = medfilt2(reshape(dy_r2i_grid, size(x_query)), [7 7]);
f=visualize_vector_field(x_query, y_query, reshape(x_query(:)+dx_r2i_grid(:), size(x_query)), reshape(y_query(:)+dy_r2i_grid(:), size(x_query)));
hgexport(f, set.lensCorrection_empirical_field_IMG,  ...
     hgexport('factorystyle'), 'Format', 'png'); 


% convert detector coordinates to focal plane coordinates (see ISIS)
[xx_pred_fp(:,1), xx_pred_fp(:,2)]= cassis_detector2focalplane( xx_pred(:,1), xx_pred(:,2), image_size(2), image_size(1), pixSize*1000); % we want focal plane coord in mm  
[xx_fp(:,1), xx_fp(:,2)] = cassis_detector2focalplane( xx(:,1), xx(:,2), image_size(2), image_size(1), pixSize*1000);

% compute lifted coordinates
chi_fp = lift2D_to_6D(xx_fp); 

% solve for correction model ( xx -> xx_pred ) 
fun = @(param) reshape(rational_model_image_side_error(param/param(end), xx_pred_fp, chi_fp), [], 1)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
sol0 = [A0_corr(1,:)'; A0_corr(2,:)'; A0_corr(3,:)'];
[sol,~,res] = lsqnonlin(fun, sol0, [], [], options);
res0 = reshape(rational_model_image_side_error(sol0/sol0(end), xx_pred_fp, chi_fp), [], 1);
A_corr = [sol(1:6)'; sol(7:12)'; sol(13:end)'];   
A_corr = A_corr/A_corr(end);

err = sqrt(sum(reshape(res/pixSize/1000, nb_points, 2).^2,2));
avgErr0 = mean(sqrt(sum(reshape(res0/pixSize/1000, nb_points, 2).^2,2)));
avgErr  = mean(err);

fprintf('Average error without correction model %d [pix] \n', avgErr0);
fprintf('Average error with correction model %d [pix] \n', avgErr);

% save correcton matrix
lensCorrection = table(A_corr);
writetable(lensCorrection, set.lensCorrection);

% show training residuals
f1 = figure;%figure('units','normalized','outerposition',[0 0 1 1]);;
C = err;
R = 100*err / 10;
scatter(xx(:,1), xx(:,2), R, C, 'filled'); hold on
caxis([0 10])
colorbar;
axis([0 2048 0 2048]);
grid on;
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'top'
colorbar;
hgexport(f1, set.lensCorrection_residuals_IMG,  ...
     hgexport('factorystyle'), 'Format', 'png'); 
 
xlabel('sample, [pix]')
ylabel('line, [pix]')

% show residual vector field
xx_corr_fp = undistort_rational_function(sol, xx_fp);
 [xx_corr(:,1), xx_corr(:,2)] = cassis_focalplane2detector(xx_corr_fp(:,1), xx_corr_fp(:,2),  image_size(2), image_size(1), pixSize*1000);
 [x_query, y_query] = meshgrid(0:10:2048, 0:10:2048);
dx_r2i = xx_pred(:,1) - xx_corr(:,1)
dy_r2i = xx_pred(:,2) - xx_corr(:,2)
dx_r2i_grid = griddata(xx_corr(:,1), xx_corr(:,2), dx_r2i, x_query(:), y_query(:));
dy_r2i_grid = griddata(xx_corr(:,1), xx_corr(:,2), dy_r2i, x_query(:), y_query(:));
dx_r2i_grid(dx_r2i_grid == NaN) = 0;
dy_r2i_grid(dy_r2i_grid == NaN) = 0;
dx_r2i_grid = medfilt2(reshape(dx_r2i_grid, size(x_query)) , [7 7]);
dy_r2i_grid = medfilt2(reshape(dy_r2i_grid, size(x_query)), [7 7]);
f=visualize_vector_field(x_query, y_query, reshape(x_query(:)+dx_r2i_grid(:), size(x_query)), reshape(y_query(:)+dy_r2i_grid(:), size(x_query)),[0.1, 0.2, 0.3, 0.4, 0.5]);
hgexport(f, set.lensCorrection_residual_field_IMG,  ...
     hgexport('factorystyle'), 'Format', 'png'); 

 % show distortion field
[x, y, i, j] = simulate_distortion_field(@undistort_rational_function, sol, image_size, pixSize*1000);
f = visualize_vector_field(i, j, x, y);
axis([0 2048 0 2048]);
hgexport(f, set.lensCorrection_field_IMG,  ...
     hgexport('factorystyle'), 'Format', 'png'); 


xlabel('sample, [pix]')
ylabel('line, [pix]')
 
% estimate distortion model
[A_dist, maxErr] = inverse_rational_model(A_corr, image_size, pixSize*1000)
fprintf('Maximum error of the distortion model with the respect to correction model %d [pix] \n', maxErr);
lensDistortion = table(A_dist);
writetable(lensDistortion, set.lensDistortion);


% 
% 
% figure;
% [x, y, i, j] = simulate_distortion_field(@undistort_rational_function, [invA(1,:)'; invA(2,:)'; invA(3,:)'], image_size);
% visualize_vector_field(i, j, x, y);
% axis([0 2048 0 2048]);

end
