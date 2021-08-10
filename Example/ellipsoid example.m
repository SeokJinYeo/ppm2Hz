# ppm2Hz usage example to calculate B0 map induced by magnetized ellipsoid

# making ellipsoid model
x=[1:80]-(1+80)/2;
y=[1:40]-(1+40)/2;
z=[1:20]-(1+20)/2;
[xx,yy,zz] = meshgrid(x,y,z);
S = (xx.^2)/1600 + (yy.^2)/400 + (zz.^2)/100 <=1;   % finally 'S' is input source model
dr = [0.001 0.001 0.001];  % voxel size[m]

# conventional case(n=l=z)
n=[0 0 1];
l=[0 0 1];
