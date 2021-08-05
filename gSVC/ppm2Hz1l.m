
function [b0,S,T] = ppm2Hz1l(chi,dr,r,t,l)  % 'n' or 'l' is [ 1 0 0 ]
[s2,s1,s3] = size(chi);                 % (2,1,3) correspond to (1st,2nd,3rd fast dim.)
dy = dr(1); dx = dr(2); dz = dr(3);     %[m]    (y,x,z) correspond to the same.
r2 = r(1);  r1 = r(2);  r3 = r(3);      %[m]
% dx = dr(1); dy = dr(2); dz = dr(3);     %[m]
% r1 = r(1);  r2 = r(2);  r3 = r(3);      %[m]
t2 = t(1);  t1 = t(2);  t3 = t(3);      % watch for 1 <-> 2
l1 = l(1);  l2 = l(2);  l3 = l(3);

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
K = l1.*atan1(Z1,Y1,X1,Z2,Y2,X2)+...
    l2.*sinh1(X1,Y1,Z1,X2,Y2,Z2)+l3.*sinh1(X1,Z1,Y1,X2,Z2,Y2); 
% F = fftn(K)/4/pi;                 % [OLD], source-centered kernel
F = conj(fftn(K))/4/pi;             %  [NEW], target-centered kernel

factor = 42.578e6*3*1e-6;         % ppm to Hz
F = F*factor;                       % k-space dipola r field core
% assignin('base','K',K);disp('K > workplace');
% assignin('base','F',F);disp('F > workplace');

%% Align the source and target matrices
S = zeros(size(X));
S(1:s2,1:s1,1:s3) = chi;        % zero-filled for target-size expansion
T = ifftn(fftn(S).*F);
b0 = T(1:t2,1:t1,1:t3);      % carve out target-size block
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
