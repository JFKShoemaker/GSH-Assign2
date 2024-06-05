% makefile for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Model
% Load previous saved model

%model_name = 'Crust01_crust';
%load(model_name);

% Construct new model
inputModel 
make_topo
make_grav

%plot(Topo)
figure; % Create a new figure
imagesc(Topo); % Display the data as a heatmap

% Customize the colormap
bwr_colormap = [
    0.0, 0.0, 1.0;  % Blue
    1.0, 1.0, 1.0;  % White
    1.0, 0.0, 0.0   % Red
];

% Interpolate to get a smooth colormap with 256 colors
n_colors = 256;
bwr_colormap = interp1([1, round(n_colors / 2), n_colors], bwr_colormap, 1:n_colors);

% Use the colormap in your figure
colormap(bwr_colormap);

% Add a colorbar
colorbar;

% Add labels
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Topography'); % Replace with appropriate title

% Adjust axis properties if needed
axis equal; % Ensures the aspect ratio is equal

factor = 4;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels







%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [-180 180 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    225000.0; % height of computation above spheroid
SHbounds =  [0 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%
%gravfield = GSHS(grav_data,latLim,lonLim,SHbounds)


latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0 360 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0; % height of computation above spheroid
SHbounds =  [0 90;];

tic
%[V] = model_SH_analysis(Model);
%V(1,3) = 0;
%V(3,3) = 0;
[GF_generated] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
toc


%%%%%%%%%%%%%%%%%%%%%%% aglagj %%%%%%%%%%%%%%%%%%%%%%%



gdata = flip(sqrt(GF_generated.vec.R.^2 + GF_generated.vec.T.^2 + GF_generated.vec.L.^2));
figure; % Create a new figure
imagesc(gdata); % Display the data as a heatmap

% Customize the colormap
colormap(bwr_colormap);

% Add a colorbar
colorbar;

% Add labels
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Gravity Field'); % Replace with appropriate title

% Adjust axis properties if needed
axis equal; % Ensures the aspect ratio is equal


% Adjust tick labels
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels




Bouguer = 2*pi*6.67430E-11*2900*Topo;

Bouguer_resized = imresize(Bouguer,[180, 361]);
figure; % Create a new figure
imagesc(Bouguer_resized  ); % Display the data as a heatmap

% Customize the colormap
colormap(bwr_colormap);
% Add a colorbar
colorbar;

% Add labels
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Bouguer_resized'); % Replace with appropriate title

% Adjust axis properties if needed
axis equal; % Ensures the aspect ratio is equal

factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels

figure; % Create a new figure
imagesc(gdata-Bouguer_resized  ); % Display the data as a heatmap

% Customize the colormap
colormap(bwr_colormap);
% Add a colorbar
colorbar;

% Add labels
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('gravity-Bouguer'); % Replace with appropriate title

% Adjust axis properties if needed
axis equal; % Ensures the aspect ratio is equal

factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels








%%%%%%%%%%%%%%%%%%%%%%% aglagj %%%%%%%%%%%%%%%%%%%%%%%


%% Global Spherical Harmonic Analysis 

tic;
[V] = model_SH_analysis(Model);
toc

save(['Results/' Model.name '.mat'],'V')

%% Global Spherical Harmonic Synthesis

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
toc

%% Save data

DATE = datestr(now);
save(['Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '_' DATE '.mat'],'GF_generated')