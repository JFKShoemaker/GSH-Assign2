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
%inputModel 
make_topo

%plot(Topo)
figure; % Create a new figure
imagesc(Topo); % Display the data as a heatmap

% Customize the colormap
colormap('hot'); % Choose a colormap (e.g., 'hot', 'jet', 'parula', etc.)

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
save(['Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '_' DATE '.mat'],'data')