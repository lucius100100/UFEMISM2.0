clc
clear all
close all

filename = "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\results_anomaly_field\main_output_ANT_00001.nc";

mesh = read_mesh_from_file(filename);
%plot_mesh(mesh);

time = ncread(filename, 'time');
ti = length(time);

%Plot ice thickness difference (LGM-PI, positive = ice growth)
Hi = ncread(filename, 'Hi');
%plot_mesh_data(mesh, Hi(:, 2)-Hi(:, 1));

%Plot ocean temperature difference
%T_ocean = ncread(filename, 'T_ocean', [1,1,ti],[Inf,Inf,1]);
%T_ocean = ncread(filename, 'T_ocean');
%T_ocean1 = ncread(filename, 'T_ocean', [1,1,1],[Inf,Inf,1]);
%ocean_T = T_ocean - T_ocean1;
%plot_mesh_data(mesh, ocean_T(:,1));
depth_level = 1;
T_ocean_t1 = ncread(filename, 'T_ocean', [1, depth_level, 1], [Inf, 1, 1]);
T_ocean_t2 = ncread(filename, 'T_ocean', [1, depth_level, ti], [Inf, 1, 1]);
T_diff = T_ocean_t2 - T_ocean_t1;
%plot_mesh_data(mesh, T_diff);
plot_mesh_data(mesh, T_ocean_t2);