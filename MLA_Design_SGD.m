%% MLA_Design_SGD.m
% by Sam Daly (sgd46@cam.ac.uk), University of Cambridge
% based on the work of Boya Zhang and Olivia Dovernor

% This code can be used to fabricate four MLA patterns compatible with
% LightForge by PowerPhotonics
% (https://www.powerphotonic.com/products/lightforge/). View the readme for
% futher information and instructions for use. 

clear all; close all; clc; addpath lib\
resultdir = fullfile(pwd, 'results', [datestr(now, 'yyyymmdd_HHMM') '_MLA_Design']);

%% Preliminary inputs and calculation of physical parameters
% Optical parameters of the microscope system that are independent of the
% light field emission path.

% input
lambda = 640e-9; % wavelength (m)
NA = 1.49; % objective numerical aperture
obj_mag = 100; % magnification of objective 
f_tube = 200e-3; % focal length of tube lens (m)
D_cam = 1e-6; % xy grid, will be downsampled from 1x1 to 10x10 later... why is this called D_cam? this is essentially *pixel size*... 
n_imm = 1.518; % refractive index of immersion medium
n_s = 1.33; % refractive index of sample

% calculations
k = 2 * pi / lambda; % calculate wavenumber 
f_obj = f_tube/obj_mag; % focal length of objective lens 
D = 2 * NA * f_obj; % back focal plane diameter (m)

%% Design 1: 3 Hex, 2.0 mm pitch, no tilt (what is tilt?)

% MLA inputs
f_fl = 175e-3; % focal length of fourier lens (m)
f_mla = 100e-3; % focal length of MLA (m)
D_mla = 2.0e-3; % pitch, aka in-diameter, of u-lens (m)
n_lens = 1.453; % refractive index of lens material

Mag = obj_mag * ( f_mla / f_fl ); % total magnification
M_relay = f_fl / f_tube; % magnification of the relay system
size_bfp = D * M_relay * ( n_s / n_imm ); % diameter of conj. BFP (m)

N_pixels_obj = round( size_bfp / D_cam ); % number of pixels along diameter of conj. BFP (see grid calculation, D_cam) 
rho_max = size_bfp/2; % maximum value for microlens coords
N_mla = size_bfp / D_mla; % number of lenses along the conj. BFP diameter
N_pixels_mla = round( D_mla / D_cam ); % number of 'pixels' per microlens 
[x_mla, y_mla] = get_MLAcentres( N_mla, N_pixels_mla ); % get xy-coordinates of MLA

curvature = calc_radius(f_mla, n_lens, 'plano'); % options: 'plano' - same units a f_mla (planoconvex only currently)

thickness1 = get_physical_mask(D_mla, D_cam, N_pixels_obj, N_pixels_mla, x_mla, y_mla, curvature); % compute MLA sagitta and pattern
thickness1 = thickness1*1e6; % convert from m into um

figure(1); imagesc(thickness1); title(['3 Hex, ' num2str(D_mla*10^3, '%.2f') ' mm pitch, ' num2str(f_mla*10^3) ' mm focal length']);
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');

T1Correction = abs(min(thickness1(:)));
T1Max = abs(max(thickness1(:)));

%% Design 2: 3 Hex, 1.9 pitch (slight adaptation)

% MLA inputs
f_fl = 175e-3; % focal length of fourier lens (m)
f_mla = 100e-3; % focal length of MLA (m)rho
D_mla = 1.9e-3; % pitch, aka in-diameter, of u-lens (m)
n_lens = 1.453; % refractive index of lens material

Mag = obj_mag * ( f_mla / f_fl ); % total magnification
M_relay = f_fl / f_tube; % magnification of the relay system
size_bfp = D * M_relay * ( n_s / n_imm ); % diameter of conj. BFP (m)

N_pixels_obj = round( size_bfp / D_cam ); % number of pixels along diameter of conj. BFP (see grid calculation, D_cam) 
rho_max = size_bfp/2; % maximum value for microlens coords
N_mla = size_bfp / D_mla; % number of lenses along the conj. BFP diameter
N_pixels_mla = round( D_mla / D_cam ); % number of 'pixels' per microlens 
[x_mla, y_mla] = get_MLAcentres( N_mla, N_pixels_mla ); % get xy-coordinates of MLA

curvature = calc_radius(f_mla, n_lens, 'plano'); % options: 'plano' - same units a f_mla (planoconvex only currently)

thickness2 = get_physical_mask(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature); % compute MLA sagitta and pattern (3 hex)
thickness2 = thickness2*1e6; % convert from m into um

figure(2); imagesc(thickness2); title(['3 Hex, ' num2str(D_mla*10^3, '%.2f') ' mm pitch, ' num2str(f_mla*10^3) ' mm focal length']); axis equal;
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; xlabel('x (μm)'); ylabel('y (μm)');

T3Correction = abs(min(thickness2(:)));
T3Max = abs(max(thickness2(:)));


%% Design 3: 7 Hex, 1.5 pitch, no tilt

% MLA inputs
f_fl = 175e-3; % focal length of fourier lens (m)
f_mla = 100e-3; % focal length of MLA (m)
D_mla = 1.5e-3; % pitch, aka in-diameter, of u-lens (m)
n_lens = 1.453; % refractive index of lens material

Mag = obj_mag * ( f_mla / f_fl ); % total magnification
M_relay = f_fl / f_tube; % magnification of the relay system
size_bfp = D * M_relay * ( n_s / n_imm ); % diameter of conj. BFP (m)

N_pixels_obj = round( size_bfp / D_cam ); % number of pixels along diameter of conj. BFP (see grid calculation, D_cam) 
rho_max = size_bfp/2; % maximum value for microlens coords
N_mla = size_bfp / D_mla; % number of lenses along the conj. BFP diameter
N_pixels_mla = round( D_mla / D_cam ); % number of 'pixels' per microlens 
[x_mla, y_mla] = get_MLAcentres( N_mla, N_pixels_mla ); % get xy-coordinates of MLA

curvature = calc_radius(f_mla, n_lens, 'plano'); % options: 'plano' - same units a f_mla (planoconvex only currently)

thickness3 = get_physical_mask_7hex(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature); % compute MLA sagitta and pattern (7 hex)
thickness3 = thickness3 * 1e6;% convert from m into um

figure(3); imagesc(thickness3); title(['7 Hex, ' num2str(D_mla*10^3, '%.2f') ' mm pitch, ' num2str(f_mla*10^3) ' mm focal length'])
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');

T7Correction = abs(min(thickness3(:)));
T7Max = abs(max(thickness3(:)));

%% Design 4: 7 Hex, 1.4 pitch, slight adaptation

% MLA inputs
f_fl = 175e-3; % focal length of fourier lens (m)
f_mla = 100e-3; % focal length of MLA (m)
D_mla = 1.4e-3; % pitch, aka in-diameter, of u-lens (m)
n_lens = 1.453; % refractive index of lens material

Mag = obj_mag * ( f_mla / f_fl ); % total magnification
M_relay = f_fl / f_tube; % magnification of the relay system
size_bfp = D * M_relay * ( n_s / n_imm ); % diameter of conj. BFP (m)

N_pixels_obj = round( size_bfp / D_cam ); % number of pixels along diameter of conj. BFP (see grid calculation, D_cam) 
rho_max = size_bfp/2; % maximum value for microlens coords
N_mla = size_bfp / D_mla; % number of lenses along the conj. BFP diameter
N_pixels_mla = round( D_mla / D_cam ); % number of 'pixels' per microlens 
[x_mla, y_mla] = get_MLAcentres( N_mla, N_pixels_mla ); % get xy-coordinates of MLA

curvature = calc_radius(f_mla, n_lens, 'plano'); % options: 'plano' - same units a f_mla (planoconvex only currently)

thickness4 = get_physical_mask_7hex(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);  % compute MLA sagitta and pattern (7 hex)
thickness4 = thickness4 * 1e6; % convert from m into um

figure(4); imagesc(thickness4); title(['7 Hex, ' num2str(D_mla*10^3, '%.2f') ' mm pitch, ' num2str(f_mla*10^3) ' mm focal length'])
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');

T4Correction = abs(min(thickness4(:)));
T4Max = abs(max(thickness4(:)));


%% Downsample
% Downsample the xy-grid from 1 um to 10 um to PowerPhotonics specifications 

design1_10um = thickness1(6:10:end, 6:10:end); 
figure(5); imagesc(design1_10um); title('Design 1 (downsampled grid 10 μm)')
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');
DownsampledMin1 = abs(min(design1_10um(:))); 
DownsampledMax1 = abs(max(design1_10um(:))); 

design2_10um = thickness2(6:10:end, 6:10:end); 
figure(6); imagesc(design2_10um); title('Design 2 (downsampled grid 10 μm)')
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');
DownsampledMin3 = abs(min(design2_10um(:)));
DownsampledMax3 = abs(max(design2_10um(:))); 

design3_10um = thickness3(6:10:end, 6:10:end); 
figure(7); imagesc(design3_10um); title('Design 3 (downsampled grid 10 μm)')
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');
DownsampledMin7 = abs(min(design3_10um(:))); 
DownsampledMax7 = abs(max(design3_10um(:))); 

design4_10um = thickness4(6:10:end, 6:10:end); 
figure(8); imagesc(design4_10um); title('Design 4 (downsampled grid 10 μm)')
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal; xlabel('x (μm)'); ylabel('y (μm)');
DownsampledMin4 = abs(min(design4_10um(:))); 
DownsampledMax4 = abs(max(design4_10um(:)));

%% Zemax modelling 
% Note: needs to be in mm not um

mkdir(resultdir)

ExtractedLens1 = design1_10um(:); % units in um
ExtractedLens1 = ExtractedLens1/1000; % units in mm
MLA1_Zemax_mm = reshape(ExtractedLens1',457^2,1); 
writematrix(MLA1_Zemax_mm, [resultdir '\MLA1_Zemax_mm']);
max_zemax1 = abs(max(ExtractedLens1(:)));
min_zemax1 = abs(min(ExtractedLens1(:)));

ExtractedLens2 = design2_10um(:); % units in um
ExtractedLens2 = ExtractedLens2/1000; % units in mm
MLA2_Zemax_mm = reshape(ExtractedLens2',457^2,1); 
writematrix(MLA2_Zemax_mm, [resultdir '\MLA2_Zemax_mm']);
max_zemax2 = abs(max(ExtractedLens2(:)));
min_zemax2 = abs(min(ExtractedLens2(:)));

ExtractedLens3 = design3_10um(:); % units in um
ExtractedLens3 = ExtractedLens3/1000; % units in mm
MLA3_Zemax_mm = reshape(ExtractedLens3',457^2,1); 
writematrix(MLA3_Zemax_mm, [resultdir '\MLA3_Zemax_mm']);
max_zemax3 = abs(max(ExtractedLens3(:)));
min_zemax3 = abs(min(ExtractedLens3(:)));

ExtractedLens4 = design4_10um(:); % units in um
ExtractedLens4 = ExtractedLens4/1000; % units in mm
MLA4_Zemax_mm = reshape(ExtractedLens4',457^2,1); 
writematrix(MLA4_Zemax_mm, [resultdir '\MLA4_Zemax_mm']);
max_zemax4 = abs(max(ExtractedLens4(:)));
min_zemax4 = abs(min(ExtractedLens4(:)));

%% Full 15x15 layout

full_map = NaN(1500,1500);
spacing1 = floor((1500-457*2)/4);
full_map((1:457)+(spacing1),(1:457)+(spacing1)) = design1_10um;
full_map((1:457)+(spacing1),(1:457)+(spacing1*3)+457) = design3_10um;
full_map((1:457)+(spacing1*3)+457,(1:457)+(spacing1)) = design4_10um;
full_map((1:457)+(spacing1*3)+457,(1:457)+(spacing1*3)+457) = design2_10um;

figure(9)
imagesc(full_map)
title('Full 15 \times 15 map')
cbar = colorbar; cbar.Label.String = 'Z height (μm)'; axis equal;
xlabel('x (μm)'); ylabel('y (μm)');

%% LightForge format 
% Produces the .csv file for Light Forge 
% (x y z in um, first line counting in 10s)

size_map = size(full_map,1);
spacing = 10; 
coordinate_map = [];

XYZ_grid = zeros(size_map+1,size_map+1); % add zeros to first row and first column 
XYZ_grid(1,:) = 0:10:(size_map*10); % adds array of 1500 in x 
XYZ_grid(:,1) = 0:10:(size_map*10); % and y to XYZ_grid
XYZ_grid(2:(size_map+1),2:(size_map+1)) = full_map; % data populated from the second row and second column of XYZ_grid
XYZ_grid = round(XYZ_grid*1000)/1000;% rounds z values to 0.001 um (z = 1nm resolution) 
DesignMLA = XYZ_grid;

figure(10)
visible_data = DesignMLA(2:end, 2:end);
imagesc(visible_data)
title('Powerphotonics - MLA Design')
cbar = colorbar; cbar.Label.String = 'Z height (um)'; 
caxis([min((visible_data(:))), max((visible_data(:)))]);

writematrix(DesignMLA, [resultdir '\MLADesign.csv'])
writematrix(DesignMLA, [resultdir '\MLADesign.dat'])
copyfile("MLA_Design_SGD.m", [resultdir '\script.txt']);
