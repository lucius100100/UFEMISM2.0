clc;
clear all;
close all;

%define filename
filename = "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\results_linear_time_no_ghf_higher_resolution_enable_jourdain_10000_test_implicit\main_output_ANT_00001.nc";

%read mesh from file
mesh = read_mesh_from_file(filename);

%read time (time)
time = ncread(filename, 'time');
ti = length(time);

%plot mesh
%plot_mesh(mesh);

%plot ocean temperature difference
depth_level = 1;
T_ocean_t1 = ncread(filename, 'T_ocean', [1, depth_level, 1], [Inf, 1, 1]);
T_ocean_t2 = ncread(filename, 'T_ocean', [1, depth_level, ti], [Inf, 1, 1]);
T_diff = T_ocean_t2 - T_ocean_t1;
H = plot_mesh_data(mesh, T_diff);

%colours for plot and colorbar
set(H.Ax, 'CLim');
colormap(H.Ax, 'jet'); 
colorbar(H.Ax, 'location', 'eastoutside');