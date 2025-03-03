clc;
clear all;
close all;

%---load in and prepare file---
filename = 'C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Results_realistic_climate_RACMO_matrix_ocean\main_output_ANT_00001.nc';

%mesh
mesh = read_mesh_from_file(filename);

plot_mesh(mesh);