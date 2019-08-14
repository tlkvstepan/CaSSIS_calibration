%function compare_distortion(trg_model_fn, src_model_fn)
src_model_fn = '/HDD1/Data/CASSIS_STARFIELD/tests/combined/OUTPUT/lensCorrection.csv';
trg_model_fn ='/HDD1/Data/CASSIS_STARFIELD_FIX/tests/combined/OUTPUT/lensCorrection.csv';

trg_model = readtable(trg_model_fn);
src_model = readtable(src_model_fn);

trg_A = [trg_model.A_corr_1  trg_model.A_corr_2 trg_model.A_corr_3 trg_model.A_corr_4 trg_model.A_corr_5 trg_model.A_corr_6];
trg_vec = reshape(trg_A',[],1);
src_A = [src_model.A_corr_1  src_model.A_corr_2 src_model.A_corr_3 src_model.A_corr_4 src_model.A_corr_5 src_model.A_corr_6];
src_vec = reshape(src_A',[],1);

pixSize = 1e-05;
image_size = [2048, 2048];
[trg_x, trg_y, i, j] = simulate_distortion_field(@undistort_rational_function, trg_vec, image_size, pixSize*1000);
[src_x, src_y, ~,~] = simulate_distortion_field(@undistort_rational_function, src_vec, image_size, pixSize*1000);

dx = trg_x - src_x; 
dy = trg_y - src_y;

f=visualize_vector_field(i, j, i+dx, j+dy, [0.5, 1, 1.5, 2]);
%hgexport(f, set.lensCorrection_residual_field_IMG,  ...
 %    hgexport('factorystyle'), 'Format', 'png'); 
axis([0 2048 0 2048]);
grid on;
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'top'




