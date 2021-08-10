# making ellipsoid model
x=[1:80]-(1+80)/2;
y=[1:40]-(1+40)/2;
z=[1:20]-(1+20)/2;
[xx,yy,zz] = meshgrid(x,y,z);
S = (xx.^2)/1600 + (yy.^2)/400 + (zz.^2)/100 <=1;
