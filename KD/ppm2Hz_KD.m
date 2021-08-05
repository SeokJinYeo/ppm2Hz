
%% Input:
%   chi: 3D grid [ppm]
%   dr: voxel size (dx1,dx2,dx3), [any unit] 
%   ff: fine-grain factor (integer)
%   bf: buffer factor (1st, 2nd, 3rd dimension, floating number >= 1)
%   n : (future extension) B0 direction [0,0,1]
%   B0: [T], defaults to 1.
%% Output:
%   b0 : field shift [Hz]
%% Example
%{
chi=zeros(8,8);chi(2:7,[2:3])=1; chi([6:7],2:7) = 1;
% chi = zeros(8,8);chi(4,4)=1;
dr = [0.001 0.001 0.001];
r = [0 0 0];
t = [8,8,1];            % Y x X x Z size     
b0 = ppm2Hz_KD(chi,dr,2,4);
figure;imagesc(chi);axis image;colorbar;title chi(ppm);
figure;imagesc(b0);axis image;colorbar;title b0(Hz);

[c0,S,T] = ppm2Hz(chi,dr,r,t);  %compare original ppm2Hz
figure;imagesc(c0);axis image;colorbar;title c0(Hz);
comp2maps(b0,c0,'RD','SVC');grid;

%% Compare with single dipole field
y1=[-1:6]*1e-3;x1=[-3:4]*1e-3;[X1,Y1]=meshgrid(x1,y1);
M=(4*pi*1e-7/4/pi)*(-1)./(eps+(X1.^2+Y1.^2).^(3/2))*1e-9*1e-6*3/(4*pi*1e-7)*42.578e6.*~logical(and(X1==0,Y1==0));
max(abs(b0(:)-M(:)))
%}

%% padarray: image processing toolbox

function b0 = ppm2Hz_KD(chi,dr,ff,bf,B0)
if nargin < 6
    B0 = 3;
 if nargin < 5
    n = [0 0 1];
  if nargin < 4
      bf = [1 1 1];
   if nargin < 3
       ff = 1;
   end
  end
 end
end
chi0 = chi;
[n1,n2,n3] = size(chi0);

if ff > 1
    M1 = ceil([1:n1*ff]/ff);
    M2 = ceil([1:n2*ff]/ff);
    M3 = ceil([1:n3*ff]/ff);
    chi = chi0(M1,M2,M3);       % fine-grained (over-sampled)
end

chi0 = chi; %fine-grained but not buffered

if max(bf) > 1
    if length(bf)==1    % isotropic buffer
        bf = [1 1 1]*bf;      % if chi is 2D, still bf = [bf, bf, bf]
    end
    chi = padarray(chi0,round([size(chi0,1),size(chi0,2),size(chi0,3)].*(bf-1)),'post');  
                % may make odd-numbered size
end

% chi: fine-grained and buffered
[N1,N2,N3] = size(chi);
%% The kernel in k space (flavor of Bouwman's "demonstration.m")
FOV = [N1,N2,N3].*dr;           % extended FOV
kx_sq = ( ifftshift( [1:N2]-(floor(N2/2)+1) )/FOV(2) ).^2;  % omit 2*pi in k space
ky_sq = ( ifftshift( [1:N1]-(floor(N1/2)+1) )/FOV(1) ).^2;
kz_sq = ( ifftshift( [1:N3]-(floor(N3/2)+1) )/FOV(3) ).^2;
[KX_SQ,KY_SQ,KZ_SQ] = meshgrid(kx_sq,ky_sq,kz_sq);
kernel = 1/3 - KZ_SQ./(KX_SQ + KY_SQ + KZ_SQ);
kernel(1,1,1) = 0;
dField = ifftn(fftn(chi).*kernel);  % [T] for B0=1T, chi in [1].
temp = dField*B0*10^(-6)*42.576e6;  % [Hz] for B0[T], chi[ppm].    

%% Cut out unbuffered size
temp2 = temp(1:size(chi0,1),1:size(chi0,2),1:size(chi0,3));
%% Downsampling by averaging

if ff>1
     b0 = C(temp2,ff);
else
     b0 = temp2;
end

function DX = C(X,ff)
ker = ones(ff,ff,ff)/ff^3;
tem = convn(X,ker);                         % Matlab's n-Dimensional convolution function
DX = tem(ff:ff:end,ff:ff:end,ff:ff:end);    % undersampling from the full matrix
end     



end