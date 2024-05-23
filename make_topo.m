filename = 'megt90n000cb.img';
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
Topo = fread(f,[360*resolution Inf ],'int16')';
fclose(f );
