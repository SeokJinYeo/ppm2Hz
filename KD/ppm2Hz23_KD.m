
function b0 = ppm2Hz23_KD(chi,dr,ff,bf,B0) % 'n' and 'l' are [0 1 0] and [0 0 1], changable
if nargin < 5
    B0 = 3;
   if nargin < 4
       bf = [1 1 1];
    if nargin < 3
        ff = 1;
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
%% The kernel in k  space (flavor of Bouwman's "demonstration.m")
FOV = [N1,N2,N3].*dr;           % extended FOV
kx_sq = ( ifftshift( [1:N2]-(floor(N2/2)+1) )/FOV(2) );  % omit 2*pi in k space
ky_sq = ( ifftshift( [1:N1]-(floor(N1/2)+1) )/FOV(1) );
kz_sq = ( ifftshift( [1:N3]-(floor(N3/2)+1) )/FOV(3) );
[KX_SQ,KY_SQ,KZ_SQ] = meshgrid(kx_sq,ky_sq,kz_sq);
kernel = - (KZ_SQ.*KY_SQ)./(KX_SQ.^2 + KY_SQ.^2 + KZ_SQ.^2);
if mod(N1,2) == 0;
    kernel(N1/2+1,:,:) = 0;
end
if mod(N2,2) == 0;
    kernel(:,N2/2+1,:) = 0;
end
if mod(N3,2) == 0;
    kernel(:,:,N3/2+1) = 0;
end

kernel(1,1,1) = 0; 
dField = ifftn(fftn(chi).*kernel);  % [T] for B0=1T, chi in [1].
temp = dField*B0*10^(-6)*42.576e6;  % [Hz] for B0[T], chi[ppm].  
%% Cut out unbuffered size
temp2 = temp(1:size(chi0,1),1:size(chi0,2),1:size(chi0,3));
%% Downsampling by averaging
if ff>1
     b0 = D(temp2,ff);
else
     b0 = temp2;
end

function DX = C(X,ff)
ker = ones(ff,ff,ff)/ff^3;
tem = convn(X,ker);                         % Matlab's n-Dimensional convolution function
DX = tem(ff:ff:end,ff:ff:end,ff:ff:end);    % undersampling from the full matrix
end    

end