clc;
clear all;
close all;

%filename
filename = "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\results_linear_time_higher_resolution_enable_jourdain_2000_test_implicit\main_output_ANT_00001.nc";

%read mesh from file
mesh = read_mesh_from_file(filename);

%read time (time)
time = ncread(filename, 'time');
ti   = length(time);

%ice thickness (time, vi)
Hi_initial = ncread(filename, 'Hi', [1, 1], [Inf, 1]);
Hi_final   = ncread(filename, 'Hi', [1, ti], [Inf, 1]);
Hi_diff    = Hi_final - Hi_initial;

%depth level
depth_level = 1;

%read ocean temperature at depth_level for initial and final time steps
T_ocean_t1 = ncread(filename, 'T_ocean', [1, depth_level, 1], [Inf, 1, 1]); 
T_ocean_t2 = ncread(filename, 'T_ocean', [1, depth_level, ti], [Inf, 1, 1]); 

%temperature difference
T_diff = T_ocean_t2 - T_ocean_t1;

% Mask values:
% icefree_land                        = 1
% icefree_ocean                       = 2
% grounded_ice                        = 3
% floating_ice                        = 4
% groundingline_gr                    = 5
% groundingline_fl                    = 6
% calvingfront_gr                     = 7
% calvingfront_fl                     = 8
% margin                              = 9
% coastline                           = 10

%read masks (time, vi)
mask_initial = ncread(filename, 'mask', [1, 1], [Inf, 1]);
mask_final   = ncread(filename, 'mask', [1, ti], [Inf, 1]);

%mask values 3 to 10 represent ice-covered regions
mask_initial_ice = (mask_initial >= 3) & (mask_initial <= 10);
mask_final_ice   = (mask_final   >= 3) & (mask_final   <= 10);

%mask value counts
mask_values = 1:10;
counts_initial = arrayfun(@(x) sum(mask_initial == x), mask_values);
counts_final   = arrayfun(@(x) sum(mask_final   == x), mask_values);

%ice thickness difference, overlain with initial and final margins
f = figure('Position',[100,100,1400,900],'Color','w');
ax = axes('Parent',f);
hold(ax, 'on');

% --- Ice thickness --

%plot ice thickness difference as background
plot_submesh_data(ax, mesh, Hi_diff, true);
title(ax, 'Ice thickness difference between LGM and PI', 'FontSize', 22);

%plot_submesh_data(ax, mesh, T_diff);

%overlay initial margin (black)
x_init = mesh.V(mask_initial_ice, 1);
y_init = mesh.V(mask_initial_ice, 2);
if numel(x_init) > 2
    K_init = boundary(x_init, y_init, 1);
    plot(ax, x_init(K_init), y_init(K_init), 'k-', 'LineWidth', 1, 'DisplayName', 'PI');
end

%overlay final margin (red)
x_final = mesh.V(mask_final_ice, 1);
y_final = mesh.V(mask_final_ice, 2);
if numel(x_final) > 2
    K_final = boundary(x_final, y_final, 1);
    plot(ax, x_final(K_final), y_final(K_final), 'r-', 'LineWidth', 1, 'DisplayName', 'LGM');
end

%legend
legend(ax, 'show', 'Location', 'best', 'FontSize', 16);

%grid on, hold off
grid(ax, 'on');
hold(ax, 'off');

% --- Masks ---

%figure for the masks themselves, sharing same color scale
figMask = figure('Position',[150,150,1800,900],'Color','w');

%positions subplots and tables
subplotWidth = 0.4;
subplotHeight = 0.5;
tableHeight = 0.2;
margin = 0.02;

%subplot initial mask
axMask1 = axes('Position',[margin, 1 - subplotHeight - margin - tableHeight, subplotWidth, subplotHeight]);
plot_submesh_data_mask(axMask1, mesh, double(mask_initial));
title(axMask1, 'Initial mask (PI)', 'FontSize', 20);

%overlay boundary
hold(axMask1,'on');
if numel(x_init) > 2
   plot(axMask1, x_init(K_init), y_init(K_init), 'k-', 'LineWidth', 1);
end
hold(axMask1,'off');

%mask counts table below initial mask subplot
uitable('Parent', figMask, ...
        'Data', [mask_values; counts_initial]', ...
        'ColumnName', {'Mask Value', 'Count'}, ...
        'RowName', arrayfun(@mask_label, mask_values, 'UniformOutput', false), ...
        'Units', 'normalized', ...
        'Position', [margin, margin, subplotWidth, tableHeight], ...
        'FontSize', 12);

%subplot final mask
axMask2 = axes('Position',[2*margin + subplotWidth, 1 - subplotHeight - margin - tableHeight, subplotWidth, subplotHeight]);
plot_submesh_data_mask(axMask2, mesh, double(mask_final));
year = (ti-1) * 1000
title(axMask2, ['Final mask (time= ' num2str(year) ' year)'], 'FontSize', 20);

%overlay boundary
hold(axMask2,'on');
if numel(x_final) > 2
   plot(axMask2, x_final(K_final), y_final(K_final), 'k-', 'LineWidth', 1);
end
hold(axMask2,'off');

%mask counts table below final mask subplot
uitable('Parent', figMask, ...
        'Data', [mask_values; counts_final]', ...
        'ColumnName', {'Mask Value', 'Count'}, ...
        'RowName', arrayfun(@mask_label, mask_values, 'UniformOutput', false), ...
        'Units', 'normalized', ...
        'Position', [2*margin + subplotWidth, margin, subplotWidth, tableHeight], ...
        'FontSize', 12);

%function to map mask values to labels
function label = mask_label(value)
    switch value
        case 1
            label = 'icefree land';
        case 2
            label = 'icefree ocean';
        case 3
            label = 'grounded ice';
        case 4
            label = 'floating ice';
        case 5
            label = 'groundingline gr';
        case 6
            label = 'groundingline fl';
        case 7
            label = 'calvingfront gr';
        case 8
            label = 'calvingfront fl';
        case 9
            label = 'margin';
        case 10
            label = 'coastline';
        otherwise
            label = 'unknown';
    end
end

%function to mimic 'plot_mesh_data_a' on a given axis,
%with color scale suitable for mask values
function plot_submesh_data_mask(axHandle, mesh, dataVals)
    %"plot_mesh_data_a" patch logic for mask values
    patch('Parent', axHandle,...
          'Vertices', mesh.V(1:mesh.nV,:),...
          'Faces',    mesh.Tri(1:mesh.nTri,:),...
          'FaceColor','interp',...
          'FaceVertexCData', dataVals,...
          'EdgeColor','none');

    axis(axHandle, [mesh.xmin mesh.xmax mesh.ymin mesh.ymax]);
    set(axHandle, 'XTick', [], 'YTick', [], 'FontSize', 14);

    %color scale [1 10]
    caxis(axHandle, [1 10]);
    cb = colorbar(axHandle, 'Location', 'eastoutside');
    set(cb, 'FontSize', 14);

    %colormap jet with 10 colors
    colormap(axHandle, jet(10));
    daspect(axHandle, [1 1 1]);

    %colorbar ticks 1-10 with labels
    cb.Ticks = 1:10;
    cb.TickLabels = arrayfun(@mask_label, 1:10, 'UniformOutput', false);
end 

%function to plot ice thickness difference or similar data
function plot_submesh_data(axHandle, mesh, dataVals, addColorbar)
    if nargin < 4
        addColorbar = true; 
    end

    patch('Parent', axHandle,...
          'Vertices', mesh.V(1:mesh.nV,:),...
          'Faces',    mesh.Tri(1:mesh.nTri,:),...
          'FaceColor','interp',...
          'FaceVertexCData', dataVals,...
          'EdgeColor','none',...
          'HandleVisibility','off');

 axis(axHandle, [mesh.xmin mesh.xmax mesh.ymin mesh.ymax]);
    set(axHandle, 'XTick', [], 'YTick', [], 'FontSize', 14);

    caxis(axHandle, [min(dataVals(:)), max(dataVals(:))]);
    
    if addColorbar
        cb = colorbar(axHandle, 'Location', 'eastoutside');
        set(cb, 'FontSize', 14);
    end

    daspect(axHandle, [1 1 1]);
end 
