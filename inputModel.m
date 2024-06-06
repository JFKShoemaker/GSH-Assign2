% This is an input file for the GSHA procedure
%
% It should contain a list of names and location of geometry boundaries followed by a
% list of names for density values
% 
% Single values for density files are allowed.

HOME = pwd;
Model = struct();
Model.number_of_layers = 2;
Model.name = 'Model_001';
% Additional variables

Model.GM =  42828.3748574E9;
Model.Re_analyse = 3396000;
Model.Re = 3396000;
Model.geoid = 'none';
Model.nmax = 120;     
Model.correct_depth = 0;


% First layer
Model.l1.bound = [HOME '/Data/crust1.bd1.gmt'];
Model.l1.dens = 2650;%[HOME '/Data/crust1.rho1.gmt'];
% Second layer
Model.l2.bound = -75000.*ones(size(gmt2matrix(importdata(Model.l1.bound))));
Model.l2.dens = 3750;
% Bottom boundary
Model.l3.bound = -300000;
save([HOME '/Data/' Model.name '.mat'],'Model')