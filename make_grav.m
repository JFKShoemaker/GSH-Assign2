
V = importdata("jgmro_120f_sha.tab");


%V = [1st column 2nd column 3rd column 4th
%column];
V = V([2:end],[1:4]);

V = [0 0 1 0; V];
V = sortrows(V,1);

grav_data = V