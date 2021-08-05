    
function b0 = ppm2Hznl2(chi,dr,r,t,l,Bx,By,Bz)  %B field should have same grid size with chi field. 
[s2,s1,s3] = size(chi);                 % (2,1,3) correspond to (1st,2nd,3rd fast dim.)
dy = dr(1); dx = dr(2); dz = dr(3);     %[m]    (y,x,z) correspond to the same.
r2 = r(1);  r1 = r(2);  r3 = r(3);      %[m]
% dx = dr(1); dy = dr(2); dz = dr(3);     %[m]
% r1 = r(1);  r2 = r(2);  r3 = r(3);      %[m]
t2 = t(1);  t1 = t(2);  t3 = t(3);      % watch for 1 <-> 2
l1 = l(1);  l2 = l(2);  l3 = l(3);
chix = chi.*Bx;   % Mx
chiy = chi.*By;   % My
chiz = chi.*Bz;   % Mz

%% Grid size = S + T - 1    % Note negative sign on the r shift!
x = single([0:s1-1,-t1+1:-1]*dx - r1);            % first element is zero - displ.
y = single([0:s2-1,-t2+1:-1]*dy - r2);
z = single([0:s3-1,-t3+1:-1]*dz - r3);

%% Construct the core matrix (Can be a separate function)
[X,Y,Z] = meshgrid(x,y,z);
X1 = X-dx/2;
X2 = X+dx/2;
Y1 = Y-dy/2;
Y2 = Y+dy/2;
Z1 = Z-dz/2;
Z2 = Z+dz/2;
Kx = l1.*atan1(Z1,Y1,X1,Z2,Y2,X2)+...
    l2.*sinh1(X1,Y1,Z1,X2,Y2,Z2)+l3.*sinh1(X1,Z1,Y1,X2,Z2,Y2); 
Ky = l2.*atan1(X1,Z1,Y1,X2,Z2,Y2)+...
    l1.*sinh1(X1,Y1,Z1,X2,Y2,Z2)+...
    l3.*sinh1(Z1,Y1,X1,Z2,Y2,X2); % 2l
Kz = l3.*atan1(X1,Y1,Z1,X2,Y2,Z2)+...
    l1.*sinh1(X1,Z1,Y1,X2,Z2,Y2)+...
    l2.*sinh1(Z1,Y1,X1,Z2,Y2,X2); % 3l
% F = fftn(K)/4/pi;                 % [OLD], source-centered kernel
Fx = conj(fftn(Kx))/4/pi;             %  [NEW], target-centered kernel
Fy = conj(fftn(Ky))/4/pi; 
Fz = conj(fftn(Kz))/4/pi; 

factor = 42.578e6*1e-6;         % ppm to Hz
Fx = Fx*factor;                       % k-space dipola r field core
Fy = Fy*factor;  
Fz = Fz*factor;  
% assignin('base','K',K);disp('K > workplace');
% assignin('base','F',F);disp('F > workplace');

%% Align the source and target matrices
Sx = zeros(size(X)); Sy = zeros(size(X)); Sz = zeros(size(X));
Sx(1:s2,1:s1,1:s3) = chix; Sy(1:s2,1:s1,1:s3) = chiy; Sz(1:s2,1:s1,1:s3) = chiz;        % zero-filled for target-size expansion
T = ifftn(fftn(Sx).*Fx + fftn(Sy).*Fy + fftn(Sz).*Fz);
b0 = T(1:t2,1:t1,1:t3); 
end

function K=atan1(X1,Y1,Z1,X2,Y2,Z2)  %z-z
K = atan(X1.*Y1./Z1./sqrt(X1.^2+Y1.^2+Z1.^2))+...
    atan(X1.*Y2./Z2./sqrt(X1.^2+Y2.^2+Z2.^2))+...
    atan(X2.*Y1./Z2./sqrt(X2.^2+Y1.^2+Z2.^2))+...
    atan(X2.*Y2./Z1./sqrt(X2.^2+Y2.^2+Z1.^2))-...
    atan(X2.*Y2./Z2./sqrt(X2.^2+Y2.^2+Z2.^2))-...
    atan(X2.*Y1./Z1./sqrt(X2.^2+Y1.^2+Z1.^2))-...
    atan(X1.*Y2./Z1./sqrt(X1.^2+Y2.^2+Z1.^2))-...
    atan(X1.*Y1./Z2./sqrt(X1.^2+Y1.^2+Z2.^2));
X = (X1+X2)/2;
Y = (Y1+Y2)/2;
Z = (Z1+Z2)/2;
ind0 = (X==0 & Y==0 & Z==0);
K(ind0) = K(ind0) + 4*pi/3;         % delta function integral term
return
end

function K=sinh1(X1,Y1,Z1,X2,Y2,Z2)  %x-y
K = asinh(Z2./sqrt(X2.^2+Y2.^2))+asinh(Z2./sqrt(X1.^2+Y1.^2))+...
    asinh(Z1./sqrt(X1.^2+Y2.^2))+asinh(Z1./sqrt(X2.^2+Y1.^2))-...
    asinh(Z2./sqrt(X2.^2+Y1.^2))-asinh(Z2./sqrt(X1.^2+Y2.^2))-...
    asinh(Z1./sqrt(X2.^2+Y2.^2))-asinh(Z1./sqrt(X1.^2+Y1.^2));
return
end
