clc
clear all
close all

%%

foldername = '../../results_DIVA';

mesh = read_mesh_from_file( [foldername '/main_output_ANT_00001.nc']);

fieldnames = {
  'u_3D_b'
  'v_3D_b'
  };

for fi = 1: length( fieldnames)
  fieldname = fieldnames{ fi};
  filename = [foldername '/' fieldname '.nc'];
  ice.(fieldname) = ncread( filename, fieldname);
%   plot_mesh_data( mesh, mean( ice.(fieldname),2));
%   title( strrep( fieldname,'_','\_'));
end

plot_mesh_data( mesh, sqrt( ice.u_3D_b( :,1).^2 + ice.v_3D_b( :,1).^2))
set( gca,'colorscale','log','clim',[0.05,5000]);

%%

foldername = '../../results_BPA';

mesh = read_mesh_from_file( [foldername '/main_output_ANT_00001.nc']);

fieldnames = {
  'u_bk'
  'v_bk'
  };

for fi = 1: length( fieldnames)
  fieldname = fieldnames{ fi};
  filename = [foldername '/' fieldname '.nc'];
  ice.(fieldname) = ncread( filename, fieldname);
%   plot_mesh_data( mesh, mean( ice.(fieldname),2));
%   title( strrep( fieldname,'_','\_'));
end

plot_mesh_data( mesh, sqrt( ice.u_bk( :,1).^2 + ice.v_bk( :,1).^2))
set( gca,'colorscale','log','clim',[0.05,5000]);