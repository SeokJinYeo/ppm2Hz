%% 2016-12-11 Changed index order: (x,y,z) or (y,x,z) --> 1st, 2nd, 3rd dimensions
%% 2016-02-01 Cleaned up and renamed from "genppm2Hz_SplusT_1.m" of 2015-01-14.

%% Input [revised 2016-12-11]
% chi: susceptibility source block [ppm], size(chi) gives (1st,2nd,3rd dim sizes) 
% grid size: dr = [1st dim size,2nd dim size,3rd dim size] [m]
% source-to-target displacement at (1,1,1): r = [in 1st dim,in 2nd dim,in 3rd dim] [m]
% target block's grid size: t = [1st dim,2nd dim,3rd dim size]
% b0 field direction: n = [nx,ny,nz] (nx^2=ny^2=nz^2=1) 
% k : strength of b0 (3T,7T) 

%% Output
% b0 [Hz] with size [t2,t1,t3].
% No masking in this version.

%% Assumption
% 1. If source and target blocks overlap, then the two blocks must be part
% of a single regular Cartesian grid (no subvoxel displacement)
% 2. If they don't overlap, I believe the subvoxel displacement is okay.
% 3. Applied field is 3 Tesla. This can be easily changed (line 67).

%% Example
%{
chi=zeros(8,8);chi(2:7,[2:3])=1; chi([6:7],2:7) = 1;
dr = [0.001 0.001 0.001];
r = [0 0.2 0.3];
t = [8,8,1];            % Y x X x Z size     
n = [0,0,1];
[b0,S,T] = nppm2Hz(chi,dr,r,t,n);
figure;imagesc(a);axis image;colorbar;title chi(ppm);
figure;imagesc(b0);axis image;colorbar;title b0(Hz);
%}

function b0 = nppm2Hz(chi,dr,r,t,n,B0)
[s2,s1,s3] = size(chi);                 % (2,1,3) correspond to (1st,2nd,3rd fast dim.)
dy = dr(1); dx = dr(2); dz = dr(3);     %[m]    (y,x,z) correspond to the same.
r2 = r(1);  r1 = r(2);  r3 = r(3);      %[m]
% dx = dr(1); dy = dr(2); dz = dr(3);     %[m]
% r1 = r(1);  r2 = r(2);  r3 = r(3);      %[m]
t2 = t(1);  t1 = t(2);  t3 = t(3);      % watch for 1 <-> 2
n1 = n(1);  n2 = n(2);  n3 = n(3);

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
K=n1^2*atan1(Z1,Y1,X1,Z2,Y2,X2)+n2^2*atan1(X1,Z1,Y1,X2,Z2,Y2)+...
  n3^2*atan1(X1,Y1,Z1,X2,Y2,Z2)+2*n1*n2*sinh1(X1,Y1,Z1,X2,Y2,Z2)+...
  2*n2*n3*sinh1(Z1,Y1,X1,Z2,Y2,X2)+2*n1*n3*sinh1(X1,Z1,Y1,X2,Z2,Y2);

% F = fftn(K)/4/pi;                 % [OLD], source-centered kernel
F = conj(fftn(K))/4/pi;             % [NEW], target-centered kernel

factor = 42.578e6*B0*1e-6;         % ppm to Hz
F = F*factor;                       % k-space dipolar field core
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

