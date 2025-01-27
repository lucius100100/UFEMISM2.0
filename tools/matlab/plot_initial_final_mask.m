clc;
clear all;
close all;

%filename
filename = "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\results_test_realistic_ocean_WOA_realistic_climate\main_output_ANT_00001.nc";

%read mesh from file
mesh = read_mesh_from_file(filename);

%read time (time)
time = ncread(filename, 'time');
ti   = length(time);

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

%compute boundary for initial and final mask
x_init = mesh.V(mask_initial_ice, 1);
y_init = mesh.V(mask_initial_ice, 2);
if numel(x_init) > 2
    K_init = boundary(x_init, y_init, 1);
end

x_final = mesh.V(mask_final_ice, 1);
y_final = mesh.V(mask_final_ice, 2);
if numel(x_final) > 2
    K_final = boundary(x_final, y_final, 1);
end

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

% --- Mask counts table ---
figMaskTables = figure('Position',[200,200,800,400],'Color','w');

combinedData = [mask_values', counts_initial', counts_final'];
columnNames = {'Mask Value', 'Initial Count', 'Final Count'};
rowNames = arrayfun(@mask_label, mask_values, 'UniformOutput', false);

%combined uitable
uitable('Parent', figMaskTables, ...
        'Data', combinedData, ...
        'ColumnName', columnNames, ...
        'RowName', rowNames, ...
        'Units', 'normalized', ...
        'Position', [0.1, 0.1, 0.8, 0.8], ...
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

%function to mimic 'plot_mesh_data_a' with color scale
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
