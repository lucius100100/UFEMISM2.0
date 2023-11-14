clc
clear all
close all

filename = '../../results_ant_spinup_part1/main_output_ANT_00001.nc';
% filename = '../../results_ant_spinup_part2/main_output_ANT_00001.nc';
% filename = '../../results_ant_spinup_part3/main_output_ANT_00001.nc';
% filename = '../../results_ant_spinup_part4/main_output_ANT_00001.nc';
% filename = '../../results_ant_spinup_part5/main_output_ANT_00001.nc';
% filename = '../../results_ant_spinup_part6/main_output_ANT_00001.nc';
% filename = '../../results_ant_spinup_part6_finalise/main_output_ANT_00001.nc';
% filename = '../../results_ant_retreat_DIVA/main_output_ANT_00001.nc';

mesh = read_mesh_from_file( filename);
time = ncread( filename,'time');

%%
close all

wa = 300;
ha = 300;
margins_hor = [25,40,40,25];
margins_ver = [40,100];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: size( H.Ax,1)
  for j = 1: size( H.Ax,2)
    set( H.Ax{ i,j},'xtick',[],'ytick',[],'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
    H.Patch{ i,j} = patch('parent',H.Ax{ i,j},'vertices',mesh.V,'faces',mesh.Tri,'facecolor','none',...
      'facevertexcdata',zeros( mesh.nV,1)+0.1,'edgecolor','interp');
    
    % Color bars
    pos = get( H.Ax{ i,j},'position');
    H.Cbar{i,j} = colorbar( H.Ax{ i,j},'location','southoutside');
    set( H.Ax{ i,j},'position',pos);
  end
end

% u
colormap( H.Ax{ 1,1}, turbo( 256))
set( H.Ax{ 1,1},'colorscale','log','clim',[0.02,5000]);
ylabel( H.Cbar{ 1,1},'Velocity [m yr^{-1}]')
set( H.Cbar{ 1,1},'ticks',[0.1,1,10,100,1000])

% dHi
colormap( H.Ax{ 1,2}, bluewhiteredmap( 33))
set( H.Ax{ 1,2},'clim',[-200,200]);
ylabel( H.Cbar{ 1,2},'Thickness error [m]')

% pore water fraction
colormap( H.Ax{ 1,3}, turbo( 256))
set( H.Ax{ 1,3},'clim',[0,1]);
ylabel( H.Cbar{ 1,3},'Pore water fraction [0-1]')

Hi0  = ncread( filename,'Hi',[1,1],[Inf,1]);

for ti = 1: length( time)
  
  title( H.Ax{ 1,2},['Time: ' num2str( time( ti)) ' yr'])
  
  Hi  = ncread( filename,'Hi'                 ,[1,ti],[Inf,1]);
  u   = ncread( filename,'uabs_surf'          ,[1,ti],[Inf,1]);
  pwf = ncread( filename,'pore_water_fraction',[1,ti],[Inf,1]);
  dHi = Hi - Hi0;
  
  u(   Hi<0.1) = NaN;
  dHi( Hi<0.1) = NaN;
  pwf( Hi<0.1) = NaN;
  
  set( H.Patch{ 1,1},'facevertexcdata',u)
  set( H.Patch{ 1,2},'facevertexcdata',dHi)
  set( H.Patch{ 1,3},'facevertexcdata',pwf)
  
  drawnow('update')
%   pause(0.1)
%   pause
  
end