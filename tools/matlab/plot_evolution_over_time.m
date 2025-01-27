clc;
clear all;
close all;

%file path
ncFile = "C:\Users\luciu\Documents\Guided research\UFEMISM2.0\results_test_realistic_ocean_WOA_realistic_climate\main_output_ANT_00001.nc";

time = ncread(ncFile, 'time');  %[time]
z_ocean = ncread(ncFile, 'depth');

%number of time steps
numTimeSteps = length(time);

%ice thickness Hi
Hi = ncread(ncFile, 'Hi');  % [vi, time]
Hi = double(Hi');           % [time, vi]

%surface elevation Hs
Hs = ncread(ncFile, 'Hs');  % [vi, time]
Hs = double(Hs');           % [time, vi]

%rate of change dHi_dt
dHi_dt = ncread(ncFile, 'dHi_dt'); % [vi, time]
dHi_dt = double(dHi_dt');          % [time, vi]

%ocean temperature
T_ocean = ncread(ncFile, 'T_ocean');   % [time, k, vi]
T_ocean = double(T_ocean);             % [time, k, vi]

%ocean salinity
S_ocean = ncread(ncFile, 'S_ocean');  % [time, k, vi]
S_ocean = double(S_ocean);            % [time, k, vi]

%average over all vertices
avg_Hi = mean(Hi, 2);   % [time]
avg_Hs = mean(Hs, 2);   % [time]
avg_T = mean(T_ocean, 3);  % [time, k]
avg_T_overall = mean(avg_T, 2);  % [time, 1]
avg_S = mean(S_ocean, 3);  % [time, k]
avg_S_overall = mean(avg_S, 2);  % [time, 1]

%plotting

%ice thickness
figure('Name', 'Average ice thickness over time', 'NumberTitle', 'off');
plot(time, avg_Hi, 'b-', 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Average Ice Thickness (m)');
title('Temporal evolution of average ice thickness');
grid on;
set(gca, 'FontSize', 12);

%surface elevation
figure('Name', 'Average surface elevation over time', 'NumberTitle', 'off');
plot(time, avg_Hs, 'g-', 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Average surface elevation (w.r.t. PD sea level)');
title('Temporal evolution of average surface elevation');
grid on;
set(gca, 'FontSize', 12);

%rate of change of ice thickness
 if exist('dHi_dt', 'var')
    figure('Name', 'Rate of change of ice thickness over time', 'NumberTitle', 'off');
     avg_dHi_dt = mean(dHi_dt, 2);  
     plot(time, avg_dHi_dt, 'm-', 'LineWidth', 2);
     xlabel('Time (years)');
     ylabel('Average dHi/dt (m yr^{-1})');
     title('Temporal evolution of average rate of change of ice thickness');
     grid on;
     set(gca, 'FontSize', 12);
 end

%ocean T
figure('Name', 'Average ocean temperature profile over time', 'NumberTitle', 'off');
imagesc(time, z_ocean, avg_T.'); 
colorbar;
xlabel('Time (years)');
ylabel('Depth (m)');
title('Average ocean temperature profile over time');
set(gca, 'YDir', 'reverse'); 
set(gca, 'FontSize', 12);
axis tight;

%ocean S
figure('Name', 'Average ocean salinity profile over time', 'NumberTitle', 'off');
imagesc(time, z_ocean, avg_S.'); 
colorbar;
xlabel('Time (years)');
ylabel('Depth (m)');
title('Average ocean salinity profile over time');
set(gca, 'YDir', 'reverse'); 
set(gca, 'FontSize', 12);
axis tight;

