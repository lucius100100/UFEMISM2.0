function plot_mesh_test

    clc
    clear all
    close all
    
    filename = "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\results_ant_template_test_1\main_output_ANT_00001.nc";
    
    mesh = read_mesh_from_file(filename);
    % plot_mesh(mesh);
    
    Hs = ncread(filename, 'T_ocean');
    
    plot_mesh_data(mesh, Hs(:,1));
    plot_mesh_data(mesh, Hs(:,11));
    
    % diff = Hs(:,11) - Hs(:,1);
    
    % plot_mesh_data(mesh, diff);
    
    % bed_roughness = ncread(filename, 'bed_roughness');
    % plot_mesh_data(mesh, bed_roughness(:,11));

end