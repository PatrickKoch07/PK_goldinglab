<<<<<<< HEAD
%%
ip = exp_20211029_InitializeExp();
gated_folder = [ip.exp.path '\spots_quantify\Run-Cy5-atrous100-thresh0-13-Dec-2021-gated'];
% gated_folder = [ip.exp.path 'spots_quantify\Test-TXRED-50atroustotalnorm50power19-Oct-2021'];

%%
%(ip, n_frame, n_channel, gated_folder)
reconstructimage_v7_w_channel(ip, 1, 2, gated_folder)

%%
reconstructimage_v7_w_phase(ip, 1, 2, gated_folder)

%%
%Variables: ip, frame#, channel directory, channel#, LUTS (can set to [] if not sure)
recon_ty_v2(ip, 1, gated_folder, 2, [])

%%
%Variables: ip, frame#, channel directory, channel#, focus channel#, LUTS (can set to [] if not sure)
recon_bestcv(ip, 1, gated_folder, 2, 1, [230 340]);

%%
recon_bestcv_max3proj(ip, 1, gated_folder, 2, 1, [230 340]);

%%
recon_bestcv_max3proj_norm(ip, 33, gated_folder, 2, 1, []);
=======
%%
ip = exp_20211029_InitializeExp();
gated_folder = [ip.exp.path '\spots_quantify\Run-Cy5-atrous100-thresh0-13-Dec-2021-gated'];
% gated_folder = [ip.exp.path 'spots_quantify\Test-TXRED-50atroustotalnorm50power19-Oct-2021'];

%%
%(ip, n_frame, n_channel, gated_folder)
reconstructimage_v7_w_channel(ip, 1, 2, gated_folder)

%%
reconstructimage_v7_w_phase(ip, 1, 2, gated_folder)

%%
%Variables: ip, frame#, channel directory, channel#, LUTS (can set to [] if not sure)
recon_ty_v2(ip, 1, gated_folder, 2, [])

%%
%Variables: ip, frame#, channel directory, channel#, focus channel#, LUTS (can set to [] if not sure)
recon_bestcv(ip, 1, gated_folder, 2, 1, [230 340]);

%%
recon_bestcv_max3proj(ip, 1, gated_folder, 2, 1, [230 340]);

%%
recon_bestcv_max3proj_norm(ip, 33, gated_folder, 2, 1, []);
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
