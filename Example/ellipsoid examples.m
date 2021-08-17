% ppm2Hzs usage examples to calculate B0 map induced by magnetized ellipsoid

% making ellipsoid model
x=[1:80]-(1+80)/2;
y=[1:40]-(1+40)/2;
z=[1:20]-(1+20)/2;
[xx,yy,zz] = meshgrid(x,y,z);
chi = (xx.^2)/1600 + (yy.^2)/400 + (zz.^2)/100 <=1;   % finally 'chi' is input source model
dr = [0.001 0.001 0.001];  % voxel size[m]

% other parameters
r = [0 0 0];      
t = [40 80 20];
ff = 1; % fine-grain factor, could be larger in order to decrease ringing artifact
bf = 2; % zero-filling buffer, could be larger in order to reduce aliasing artifact
B0 = 3; % strength of applied b0 field(3T)

% conventional case(n=l=z)
n=[0 0 1];
l=[0 0 1];
b0ellip = ellipsoid_field3(80,40,20,n,l,chi,B0);  % analytical solution of b0 field inside the ellipsoid(Osborn JA. 1945), should download at example folder
ellipanalytical = b0ellip(20,40,10);  % analytical b0 value along x-axis
b01 = ppm2Hz(chi,dr,r,t,B0);  %produced by gSVC
b02 = ppm2Hz_KD(chi,dr,ff,bf,B0); %produced by KD
y1=squeeze(b01(20, :, 10)); y2=squeeze(b02(20, :, 10)); % calculated b0 value along x-axis of each method
x=[1:80]-(1+80)/2;
figure; plot(x,y1,'b'); 
hold on; plot(x,y2,'g');
axis([-40 40 -40 0]);
xlabel('x[mm]');
xticks([-40 -20 0 20 40]);
ylabel('B_{0}[Hz]'); 
yticks([-40 -30 -20 -10 0]);
hold on;
fplot(ellipanalytical,[-40,40],'r--');
set(gca,'FontSize',23);

% arbitrary directon case
n=[1/sqrt(3),1/sqrt(3),1/sqrt(3)];
l=[0,0,1];
b0ellip = ellipsoid_field3(80,40,20,n,l,chi,B0);  % analytical solution of b0 field inside the ellipsoid(Osborn JA. 1945), should download at example folder
ellipanalytical = b0ellip(20,40,10);  % analytical b0 value along x-axis
b01 = ppm2Hznl(chi,dr,r,t,n,l,B0);  %produced by gSVC
b02 = ppm2Hznl_KD(chi,dr,ff,bf,n,l,B0); %produced by KD
y1=squeeze(b01(20, :, 10)); y2=squeeze(b02(20, :, 10)); % calculated b0 value along x-axis of each method
x=[1:80]-(1+80)/2;
figure; plot(x,y1,'b'); 
hold on; plot(x,y2,'g');
axis([-40 40 -40 0]);
xlabel('x[mm]');
xticks([-40 -20 0 20 40]);
ylabel('B_{0}[Hz]'); 
yticks([-40 -30 -20 -10 0]);
hold on;
fplot(ellipanalytical,[-40,40],'r--');
set(gca,'FontSize',23);

% spatially varying applied field case
l=[0,0,1];
Bx = -yy.^2;
By = xx.^2;
Bz = zeros(40,80,20); % make spatially varying magnetic field
b01 = ppm2Hznl2(chi,dr,r,t,l,Bx,By,Bz);  %produced by gSVC
b02 = ppm2Hznl2_KD(chi,dr,ff,bf,l,Bx,By,Bz); %produced by KD
b03 = ppm2Hznl2_KD(chi,dr,2,bf,l,Bx,By,Bz); %produced by KD when ff=2
y1=squeeze(b01(20, :, 10)); y2=squeeze(b02(20, :, 10)); y3=squeeze(b03(20, :, 10)); % calculated b0 value along x-axis of each method
x=[1:80]-(1+80)/2;
figure; plot(x,y1,'b'); 
hold on; plot(x,y2,'g');
hold on; plot(x,y3,'r');
axis([-40 40 -100 100]);
xlabel('x[mm]');
xticks([-40 -20 0 20 40]);
ylabel('B_{0}[Hz]'); 
yticks([-100 -50 0 50 100]);
