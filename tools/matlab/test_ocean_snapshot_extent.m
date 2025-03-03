clear; close all; clc;

%file paths
meshFile = 'C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Test_ice_shelf_1\main_output_ANT_00001.nc';
dataDir = 'C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Test_ice_shelf\';
files = { 'ocean_T.nc', ...
           'ocean_matrix_T_initialise_0.nc', ...
          'ocean_matrix_T_initialise_1.nc', ...
    %'ocean_main_T_initialise.nc', ...
          %'ocean_main_T_run.nc', ...
          %'ocean_matrix_T_initialise_final.nc', ...
          %'ocean_matrix_T_run.nc', ...
          %'UFEMISM_main_ocean_T_initialise_A.nc', ...
          %'UFEMISM_main_ocean_T_initialise_B.nc', ...
          %'UFEMISM_main_ocean_T_initialise_run.nc', ...
          %'UFEMISM_main_ocean_T_run_A.nc', ...
          %'UFEMISM_program_ocean_T_initialise.nc', ...
          %'UFEMISM_program_ocean_T_run.nc', ... 
          %'ocean_netcdf_before_mesh_remapping.nc', ...
          %'ocean_netcdf_after_mesh_remapping.nc' };
          };

nFiles = length(files);

%load Mesh
mesh = read_mesh_from_file(meshFile);
% mesh.nV from the size of mesh.V.
sz = size(mesh.V);
if sz(1) < sz(2)
    mesh.nV = sz(2);
else
    mesh.nV = sz(1);
end
%number of triangles
mesh.nTri = size(mesh.Tri, 1);
%dummy edge count if missing
if ~isfield(mesh, 'nE')
    mesh.nE = 0;
end
%mesh boundaries
if sz(1) < sz(2)
    mesh.xmin = min(mesh.V(1,:));
    mesh.xmax = max(mesh.V(1,:));
    mesh.ymin = min(mesh.V(2,:));
    mesh.ymax = max(mesh.V(2,:));
else
    mesh.xmin = min(mesh.V(:,1));
    mesh.xmax = max(mesh.V(:,1));
    mesh.ymin = min(mesh.V(:,2));
    mesh.ymax = max(mesh.V(:,2));
end

%one big figure with subplots
nRowsLayout = 3;
nColsLayout = 5;
figure('Position',[50,50,2000,1000],'Color','w');

%loop over each ocean T file
for k = 1:nFiles
    subplot(nRowsLayout, nColsLayout, k);
    
    %full file path
    curFile = fullfile(dataDir, files{k});
    [~, varName, ~] = fileparts(curFile);
    
    %read in ocean T
    T = ncread(curFile, varName);
    
    if size(T,1) == mesh.nV
        d = T(:,1);  
    elseif size(T,2) == mesh.nV
        d = T(1,:);  
    else
        %error('Data dimensions in %s do not match mesh.nV', curFile);
    end
    
    %plot mesh and overlay data
    edgecolor = 'none';
    patch('Vertices', mesh.V, 'Faces', mesh.Tri, 'FaceColor', 'interp', ...
          'FaceVertexCData', d, 'EdgeColor', edgecolor);
    axis equal;
    xlim([mesh.xmin, mesh.xmax]);
    ylim([mesh.ymin, mesh.ymax]);
    title(varName, 'Interpreter', 'none', 'FontSize',8);
    colorbar;
end

%---top level t_an from two snapshots interpolated between---

%file paths
extraFile1 = 'C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Data\Data_Meike\PMIP4_GCM_data_Interpolated\Ens_AWI_INM_MIROC_MPI_PMIP4_Ocean_LGM.nc';
extraFile2 = 'C:\Users\luciu\Documents\Guided research\UFEMISM2.0\Data\Input\Ocean\WOA\woa18_decav_ts00_04_remapcon_r360x180_NaN.nc';

figure('Position',[100,100,1200,600],'Color','w');

%---PMIP4---
subplot(1,2,1);

t_an1 = ncread(extraFile1, 't_an');


if ndims(t_an1)==3
    topT1 = t_an1(:,:,1);
elseif ndims(t_an1)==4
    topT1 = t_an1(:,:,1,1);
else
    topT1 = t_an1;
end

%coordinates
lon1 = ncread(extraFile1, 'lon');
lat1 = ncread(extraFile1, 'lat');

if isvector(lon1) && isvector(lat1)

    [LON1, LAT1] = meshgrid(lon1, lat1);

    if ~isequal(size(topT1), size(LON1))

        if isequal(size(topT1), fliplr(size(LON1)))
            topT1 = topT1';
        else
            error('Dimensions of topT1 do not match the lon/lat grid for extraFile1.');
        end
    end
else

    LON1 = lon1;
    LAT1 = lat1;
    if ~isequal(size(topT1), size(LON1))
        error('Dimensions of topT1 do not match the lon/lat grid for extraFile1.');
    end
end

%plot
h1 = pcolor(LON1, LAT1, topT1);
set(h1, 'AlphaData', ~isnan(topT1), 'FaceAlpha', 'flat'); 
shading interp;
colorbar;
title('PMIP4 Ocean LGM (Top Level t\_an)');
xlabel('Longitude');
ylabel('Latitude');

%---WOA---
subplot(1,2,2);

t_an2 = ncread(extraFile2, 't_an');

if ndims(t_an2)==3
    topT2 = t_an2(:,:,1);
elseif ndims(t_an2)==4
    topT2 = t_an2(:,:,1,1);
else
    topT2 = t_an2;
end

lon2 = ncread(extraFile2, 'lon');
lat2 = ncread(extraFile2, 'lat');

if isvector(lon2) && isvector(lat2)
    [LON2, LAT2] = meshgrid(lon2, lat2);
    if ~isequal(size(topT2), size(LON2))
        if isequal(size(topT2), fliplr(size(LON2)))
            topT2 = topT2';
        else
            error('Dimensions of topT2 do not match the lon/lat grid for extraFile2.');
        end
    end
else
    LON2 = lon2;
    LAT2 = lat2;
    if ~isequal(size(topT2), size(LON2))
        error('Dimensions of topT2 do not match the lon/lat grid for extraFile2.');
    end
end

pcolor(LON2, LAT2, topT2);
shading interp;
colorbar;
title('WOA Ocean Data (Top Level t\_an)');
xlabel('Longitude');
ylabel('Latitude');

%---plotting of ocean fields from netcdf input module---

beforeFile = fullfile(dataDir, 'ocean_netcdf_before_mesh_remapping.nc');
afterFile  = fullfile(dataDir, 'ocean_netcdf_after_mesh_remapping.nc');

varBefore = 'ocean_netcdf_before_mesh_remapping';
varAfter  = 'ocean_netcdf_after_mesh_remapping';

%read data
dataBefore = ncread(beforeFile, varBefore);
lonBefore  = ncread(beforeFile, 'lon');
latBefore  = ncread(beforeFile, 'lat');

if isvector(lonBefore) && isvector(latBefore)
    [LONBefore, LATBefore] = meshgrid(lonBefore, latBefore);
else
    LONBefore = lonBefore;
    LATBefore = latBefore;
end

% Read data and coordinate variables from the "after" file
dataAfter = ncread(afterFile, varAfter);
lonAfter  = ncread(afterFile, 'lon');
latAfter  = ncread(afterFile, 'lat');

if isvector(lonAfter) && isvector(latAfter)
    [LONAfter, LATAfter] = meshgrid(lonAfter, latAfter);
else
    LONAfter = lonAfter;
    LATAfter = latAfter;
end

%plotting
figure('Position', [200, 200, 1200, 600], 'Color', 'w');

subplot(1,2,1);
pcolor(LONBefore, LATBefore, dataBefore);
shading interp;
colorbar;
title('Before Mesh Remapping','Interpreter','none');
xlabel('Longitude');
ylabel('Latitude');

subplot(1,2,2);
pcolor(LONAfter, LATAfter, dataAfter);
shading interp;
colorbar;
title('After Mesh Remapping','Interpreter','none');
xlabel('Longitude');
ylabel('Latitude');
