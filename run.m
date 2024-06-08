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
flexural


%%%%%%%%%%%%%% Topography from PDS %%%%%%%%%%%%%%%%%%%%%%%

figure; % Create a new figure
imagesc(Topo/1e3); % Display the data as a heatmap

bwr_colormap = [
    0.0, 0.0, 1.0;  % Blue
    1.0, 1.0, 1.0;  % White
    1.0, 0.0, 0.0   % Red
];
n_colors = 256;
bwr_colormap = interp1([1, round(n_colors / 2), n_colors], bwr_colormap, 1:n_colors);

colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'Topography (km)';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Topography'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 4;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels

%%%%%%%%%%%%%% gravity field from Observations (PDS) %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0.5 359.5 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0; % height of computation above spheroid
SHbounds =  [0 90;];

V(1,3) = 0;
V(3,3) = 0;
[GF_generated] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);

g_obs = 1e5* flip(sqrt(GF_generated.vec.X.^2 + GF_generated.vec.Y.^2 + GF_generated.vec.Z.^2)); %1e5 for converting into mGal

figure; % Create a new figure
imagesc(g_obs); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('g_{obs}'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels


Bouguer = 1e5*2*pi*6.67430E-11*2900*Topo;  % 1e5 for converting into mGal
Bouguer_resized = imresize(Bouguer,[180, 360]);

figure; % Create a new figure
imagesc(Bouguer_resized  ); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Bouguer resized'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels

figure; % Create a new figure
imagesc(g_obs-Bouguer_resized  ); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('g_{obs}-g_{Bouguer}'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels


% Finding D:

delta_g = 1e-5*(g_obs-Bouguer_resized);
delta_rho = 600;
Delta_h = delta_g / (2* pi *6.67430E-11* delta_rho);



%%%%%%%%%%%%%%%%%%%%%%% MODEL 1 %%%%%%%%%%%%%%%%%%%%%%%
tic
FC = 0.3;
D = 75000; %max(Delta_h(:))-min(Delta_h(:));
V_m1 = segment_2layer_model(imresize(Topo,[180, 360]), -ones(180, 360)*D, -200000, 2900, 3500, 25000, Model );
V_m1(1,3) = 0;
V_m1(4,3) = 0;

[GF_generated_M1] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_m1,Model);
gM1 = 1e5 * flip(sqrt(GF_generated_M1.vec.R.^2 + GF_generated_M1.vec.T.^2 + GF_generated_M1.vec.L.^2));

dg = g_obs - gM1;
dr1 = dg*FC;
dr2 = 2*dr1;
iter = 1;
eps = 1e3;
diff = dr1;
tab1 = g_obs - gM1;

while abs(1-max(dr1(:))/max(dr2(:))) > 0.01 & iter < 0
    dr2 = dr1;
    V_m1 = segment_2layer_model(imresize(Topo,[180, 360]), -ones(180, 360)*D-dr1, -200000, 2900, 3500, 25000, Model );
    V_m1(1,3) = 0;
    V_m1(3,3) = 0;
    
    [GF_generated_M1] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_m1,Model);
    gM1 = 1e5*flip(sqrt(GF_generated_M1.vec.X.^2 + GF_generated_M1.vec.Y.^2 + GF_generated_M1.vec.Z.^2));
    
    dg = g_obs - gM1;
    dr1 = dr1+dg*FC;
    diff = dr2-dr1;
    iter = iter+1;

    disp(iter);
    disp(['max ratio ', num2str(abs(1-max(dr1(:))/max(dr2(:)))), ' ------ ', 'min ratio ', num2str(abs(1-min(dr1(:))/min(dr2(:))))])
end
toc
tab2 = g_obs-gM1;

figure; % Create a new figure
imagesc(g_obs-gM1); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('g_{obs} - g_{M1}'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels

figure; % Create a new figure
imagesc(gM1); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('g_{M1}'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels

%%%%%%%%%%%%%%%%%%%%%%% M2 isostasy %%%%%%%%%%%%%%%%%%%%%%%

rho_c = 2900;
rho_m = 3500;

delta_rm2 = Topo*rho_c/(rho_m-rho_c);
Tc = Topo + ones(720, 1440)*D + delta_rm2;

figure; % Create a new figure
imagesc(Tc); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'Thickness (km)';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Tc'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels

V_M2 = segment_2layer_model(Topo, -Tc, -200000, 2900, 3500, 25000, Model );
V_M2(1,3) = 0;
V_M2(3,3) = 0;

[GF_generated_M2] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_M2,Model);
g_M2 = 1e5* flip(sqrt(GF_generated_M2.vec.R.^2 + GF_generated_M2.vec.T.^2 + GF_generated_M2.vec.L.^2)); %1e5 for converting into mGal

figure; % Create a new figure
imagesc(g_obs-g_M2); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title(['g_{obs}-g_{M2}']); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels


%%%%%%%%%%%%%%%%%%%%%%% M3 flexural %%%%%%%%%%%%%%%%%%%%%%%

Te = 400000;
rho_c = 2900;
rho_m = 3500;
g = 3.72;
R = 3396000;
D_M3 = (65e9*Te^3)/(12*(1-0.25*0.25));

function [phi] = flexural_inf (n, D, rho_c, rho_m, g, R)
    phi = ( 1+ (D/((rho_m-rho_c)*g))*((2*n+1)/(2*R))^4) ^ -1;
end

V = segment_2layer_model(Topo, -Tc, -200000, 2900, 3500, 25000, Model );

for i = 1:7381
    V(i,3)=V(i,3)*flexural_inf(V(i,1), D_M3, rho_c, rho_m, g, R);
    V(i,4)=V(i,4)*flexural_inf(V(i,1), D_M3, rho_c, rho_m, g, R);
end

V(1,3) = 0;
V(3,3) = 0;

[GF_generated_M3] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
g_M3 = 1e5* flip(sqrt(GF_generated_M3.vec.R.^2 + GF_generated_M3.vec.T.^2 + GF_generated_M3.vec.L.^2)); %1e5 for converting into mGal

figure; % Create a new figure
imagesc(g_M3); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('Flexural'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels


figure; % Create a new figure
imagesc(g_obs - g_M3); % Display the data as a heatmap
colormap(bwr_colormap);
hColorbar = colorbar;
hColorbar.Label.String = 'mGal';
xlabel('Longitude'); % Replace with appropriate label
ylabel('Latitude'); % Replace with appropriate label
title('g_{obs} - g_{M3}'); % Replace with appropriate title
axis equal; % Ensures the aspect ratio is equal
factor = 1;
xticks = get(gca, 'XTick'); % Get current x-axis tick values
xticklabels = xticks / factor; % Compute new tick labels
set(gca, 'XTickLabel', xticklabels); % Set new tick labels
yticks = get(gca, 'YTick'); % Get current y-axis tick values
yticklabels = yticks / factor; % Compute new tick labels
set(gca, 'YTickLabel', yticklabels); % Set new tick labels


%%%%%%%%%%%%%%%%%%%%%%% aglagj %%%%%%%%%%%%%%%%%%%%%%%



