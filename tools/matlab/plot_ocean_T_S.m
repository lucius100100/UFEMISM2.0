clc;
clear all;
close all;

%---load in and prepare file---
filename = 'E:\Master\Guided_research\Results_low_resolution\Results_realistic_climate_RACMO_matrix_ocean\main_output_ANT_00001.nc';

%mesh
mesh = read_mesh_from_file(filename);

%time
time = ncread(filename,'time');
ti   = length(time);

%---plot ocean T---
depth_level = 1;
T_ocean_t1 = ncread(filename, 'T_ocean', [1, depth_level, 1],  [Inf, 1, 1]);
T_ocean_t2 = ncread(filename, 'T_ocean', [1, depth_level, ti], [Inf, 1, 1]);
T_diff     = T_ocean_t2 - T_ocean_t1;

%plot T_diff
H = plot_mesh_data(mesh, T_ocean_t2);
set(H.Ax, 'CLim');
colormap(H.Ax, 'jet');
title(H.Ax, ['Ocean temperature difference']);
%title(H.Ax, ['Ocean temperature difference (depth=', num2str(depth_level), ')']);
ylabel(H.Cbar, 'Temperature (^\circ C)');

%position colorbar labels
H.Cbar.Label.Units = 'normalized';
posLabel = H.Cbar.Label.Position;
posLabel(1) = posLabel(1) - 1;
H.Cbar.Label.Position = posLabel;

%position colorbar
pos = get(H.Cbar, 'Position');
pos(1) = pos(1) - 0.005;
set(H.Cbar, 'Position', pos);

H.Cbar.FontSize = 12;
H.Cbar.Label.FontSize = 12;

%---overlay basins and shelves---
ax = H.Ax;
hold(ax, 'on');

%polar grid lines
theta = linspace(0, 2*pi, 360);
radii = 500e3 : 500e3 : 3000e3;
for rVal = radii
    xCirc = rVal * cos(theta);
    yCirc = rVal * sin(theta);
    plot(ax, xCirc, yCirc, 'k:', 'HandleVisibility','off', 'Color', [0.5, 0.5, 0.5]);
end

%radial lines every 30 degrees
angles = 0 : 30 : 330;
rMax   = max(radii);
for aVal = angles
    xRad = [0, rMax * cosd(aVal)];
    yRad = [0, rMax * sind(aVal)];
    plot(ax, xRad, yRad, 'k:', 'HandleVisibility','off', 'Color', [0.5, 0.5, 0.5]);
end

%polar coordinate labels
labelOffset = 100e3;

%labels at each 30° interval
for aVal = angles
    xLabel = (rMax - labelOffset) * cosd(aVal);
    yLabel = (rMax - labelOffset) * sind(aVal);
    text(ax, xLabel, yLabel, sprintf('%d°', mod(90 - aVal, 360)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', 10, 'Color', 'k');
end

%500-km scale bar
scaleLen     = 500e3; 
offset_right = 100e3; 
offset_bottom= 100e3; 
sx = mesh.xmax - offset_right - scaleLen;  
sy = mesh.ymin + offset_bottom; 

%plot the scale bar
plot(ax, [sx, sx + scaleLen], [sy, sy], 'k-', 'LineWidth', 4, 'HandleVisibility','off'); 
text(ax, sx + scaleLen/2, sy + 10e3, '500 km', ... 
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...  
    'FontSize',12, ...
    'Color','k');

%execute basins and shelves function
overlayBasinsAndShelves(ax, ...
    "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Data\Input\Basins\Basins_Antarctica_v02.shp", ...
    "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Data\Input\Basins\IceShelf_Antarctica_v02.shp");

hold(ax, 'off');

%---plot ocean S---
S_ocean_t1 = ncread(filename, 'S_ocean', [1, depth_level, 1],  [Inf, 1, 1]);
S_ocean_t2 = ncread(filename, 'S_ocean', [1, depth_level, ti], [Inf, 1, 1]);
S_diff     = S_ocean_t2 - S_ocean_t1;

H = plot_mesh_data(mesh, S_diff);
title(H.Ax, ['Ocean salinity difference']);
%title(H.Ax, ['Ocean salinity difference (depth=', num2str(depth_level), ')']);
ylabel(H.Cbar, 'Salinity (psu)');

%position colorbar labels
H.Cbar.Label.Units = 'normalized';
posLabel = H.Cbar.Label.Position;
posLabel(1) = posLabel(1) - 1;
H.Cbar.Label.Position = posLabel;

%position colorbar
pos = get(H.Cbar, 'Position');
pos(1) = pos(1) - 0.005;
set(H.Cbar, 'Position', pos);

H.Cbar.FontSize = 12;
H.Cbar.Label.FontSize = 12;

%---overlay basins and shelves---
ax = H.Ax;
hold(ax, 'on');

%polar grid lines
theta = linspace(0, 2*pi, 360);
radii = 500e3 : 500e3 : 3000e3;
for rVal = radii
    xCirc = rVal * cos(theta);
    yCirc = rVal * sin(theta);
    plot(ax, xCirc, yCirc, 'k:', 'HandleVisibility','off', 'Color', [0.5, 0.5, 0.5]);
end

%radial lines every 30 degrees
angles = 0 : 30 : 330;
rMax   = max(radii);
for aVal = angles
    xRad = [0, rMax * cosd(aVal)];
    yRad = [0, rMax * sind(aVal)];
    plot(ax, xRad, yRad, 'k:', 'HandleVisibility','off', 'Color', [0.5, 0.5, 0.5]);
end

%polar coordinate labels
labelOffset = 100e3;

%labels at each 30° interval
for aVal = angles
    xLabel = (rMax - labelOffset) * cosd(aVal);
    yLabel = (rMax - labelOffset) * sind(aVal);
    text(ax, xLabel, yLabel, sprintf('%d°', mod(90 - aVal, 360)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', 10, 'Color', 'k');
end

%500-km scale bar
scaleLen     = 500e3; 
offset_right = 100e3; 
offset_bottom= 100e3; 
sx = mesh.xmax - offset_right - scaleLen;  
sy = mesh.ymin + offset_bottom; 

%plot the scale bar
plot(ax, [sx, sx + scaleLen], [sy, sy], 'k-', 'LineWidth', 4, 'HandleVisibility','off'); 
text(ax, sx + scaleLen/2, sy + 10e3, '500 km', ... 
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...  
    'FontSize',12, ...
    'Color','k');

%execute basins and shelves function
overlayBasinsAndShelves(ax, ...
    "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Data\Input\Basins\Basins_Antarctica_v02.shp", ...
    "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Data\Input\Basins\IceShelf_Antarctica_v02.shp");

hold(ax, 'off');