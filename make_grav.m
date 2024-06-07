
V = importdata("jgmro_120f_sha.tab");


%V = [1st column 2nd column 3rd column 4th
%column];
V = V([2:end],[1:4]);

V = [0 0 1 0; V];
V = sortrows(V,2);
V(1,3) = 0;
V(4,3) = 0;

grav_data = V;