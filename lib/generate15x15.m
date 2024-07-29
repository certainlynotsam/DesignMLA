% Generate 15x15 grid for Powerphotonics
% Boya

lambda = 605e-9;
k=2*pi/lambda; % wavenumber
n_imm=1.406; % refractive index of immersion medium
n_s=1.33; % refractive index of sample (water)
NA=1.3; % objective numerical aperture
f_obj=3e-3; % focal length of objective lens
D_cam=10e-6; % pixel size of grid

D=2*NA*f_obj; % diameter of objective lens

%% 1. 3 Hex, 1.99 pitch, no tilt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_tube=200e-3; % focal length of tube lens
f_fl=125e-3; % focal length of fourier lens
f_mla=80e-3; % focal length of MLA
D_mla=1.99683e-3; % diameter of u-lens, between flat edges

M_relay=(f_fl/f_tube); % magnification of the relay system
size_bfp=D*M_relay*(n_s/n_imm); % diameter of relayed BFP
N_pixels_obj=round(size_bfp/D_cam); % number of pixels in objective plane (one dimension)
rho_max=size_bfp/2;
N_mla=size_bfp/D_mla;
N_pixels_mla=round(D_mla/D_cam); % number of pixels per microlens 
[x_mla,y_mla]=get_MLAcentres(N_mla,N_pixels_mla); % xy-coordinate of u-lens
curvature =36.662; % unit in mm
thickness1 = get_physical_mask(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);
thickness1 = thickness1*1e6;
figure(1); imagesc(thickness1); title('1. 3 Hex, 1.99 pitch, no tilt')

%% 7. 7 Hex, 1.53 pitch, no tilt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_tube=200e-3; % focal length of tube lens
f_fl=125e-3; % focal length of fourier lens
f_mla=61.63e-3; % focal length of MLA
D_mla=1.53e-3; % diameter of u-lens, between flat edges

M_relay=(f_fl/f_tube); % magnification of the relay system
size_bfp=D*M_relay*(n_s/n_imm); % diameter of relayed BFP
N_pixels_obj=round(size_bfp/D_cam); % number of pixels in objective plane (one dimension)
rho_max=size_bfp/2;
N_mla=size_bfp/D_mla;
N_pixels_mla=round(D_mla/D_cam); % number of pixels per microlens 
[x_mla,y_mla]=get_MLAcentres(N_mla,N_pixels_mla); % xy-coordinate of u-lens
curvature = 28.219;
thickness7 = get_physical_mask_7hex(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);
thickness7 = thickness7*1e6;
figure(7); imagesc(thickness7); title('7. 7 Hex, 1.53 pitch, no tilt')



%% Full 15x15 layout
full_map = NaN(1500,1500);
spacing1 = floor((1500-461*3)/4);
full_map((1:461)+spacing1,(1:461)+spacing1) = thickness1;
full_map((1:461)+spacing1,(1:461)+spacing1*3+461*2) = thickness7;

figure(10)
imagesc(full_map)
title('Full 15x15 map')

%% Generating LightForge format
size_map = size(full_map,1);
spacing = 10; 
coordinate_map = [];

XYZ_grid = zeros(size_map+1,size_map+1);
XYZ_grid(1,:) = 0:10:(size_map*10);
XYZ_grid(:,1) = 0:10:(size_map*10);
XYZ_grid(2:(size_map+1),2:(size_map+1)) = full_map;
XYZ_grid = round(XYZ_grid*1000)/1000;
writematrix(XYZ_grid,'spiramp2_XYZ_grid.csv')