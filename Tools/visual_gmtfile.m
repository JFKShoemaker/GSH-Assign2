function [A,Lon,Lat] = visual_gmtfile(filename,zunit,file_type)
%
% This function visualizes a .gmt file. 
%
% ONLY for inspection purposes

if nargin == 1
    file_type = 'block';
    zunit = '-';
end

if nargin == 2
    file_type = 'block';    
end

d = load(filename);

if strcmp(file_type,'block')
    [A,Lon,Lat] = gmt2matrix(d);
elseif strcmp(file_type,'gauss')
    [A,Lon,Lat] = gmt2matrix_gauss(d);
else
    error('File type must be string: block or gauss')
end

B = Europe_centered(A);

lon = Lon(1,:);
lats = Lat(:,1);
lon = lon - 180;

load coast;

figure;
hold on
imagesc(lon,lats,(B));c=colorbar; 
plot(long,lat,'k');
xlim([min(lon) max(lon)])
ylim([min(lats) max(lats)])
hold off
title(['Visualization of ' filename])
xlabel('LON degree')
ylabel('LAT degree')
ylabel(c,zunit)

function [n,DV] = degreeVariance(gmt)

nmax = max(gmt(:,1));

A = zeros(nmax+1,nmax+1);
B = zeros(nmax+1,nmax+1);
k = 0;

for s=0:nmax
    
    i = s + 1;
    l = k + 1;
    k = l + nmax-s;
    
    A(i:end,i) = gmt(l:k,3);
    B(i:end,i) = gmt(l:k,4);
    
end

A = A.^2;
B = B.^2;

n = 0:nmax;
n = n';

DV = sum(A,2)+sum(B,2);

function [d] = matrix2gmt(A,Lon,Lat)
%
% This function converts the matrix into gmt based vectors
%
% Written by Bart Root, 26-07-2012
%
% input:
%           - A = value matrix
%           - Lon = longitude matrix in degree
%           - Lat = latitude matrix
%
%%%%%%%%%%%%%%% start of routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll = length(Lat(:,1));
ww = length(Lat(1,:));
axl = size(A,1)*size(A,2);
d = zeros(axl,3);

for i = 1:ll
   
    bb = 1 + (i-1)*ww;
    ee = ww + (i-1)*ww;
    
    d(bb:ee,3) = A(i,:);
    d(bb:ee,1) = Lon(i,:); 
    d(bb:ee,2) = Lat(i,:);

end

