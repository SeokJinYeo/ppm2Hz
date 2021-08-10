 function b0 = ellipsoid_field3(a,b,c,n,l,chi,B) % a>=b>=c 
if a == b
    m = a/c;
    L = ((m^2)*(m^2-1)^(-0.5)*asin((m^2-1)^(0.5)/m)-1)*4*(pi)/(2*(m^2-1));
    M=L; 
    N=(m^2)/(m^2-1)*(1-1/(m^2-1)^0.5*asin((m^2-1)^0.5/m))*4*pi; 
elseif b == c
    m = a/c;   
    L = (m/(2*(m^2-1)^0.5)*log((m+(m^2-1)^0.5)/(m-(m^2-1)^0.5))-1)/(m^2-1)*4*pi;
    M = m*(m-1/(2*(m^2-1)^0.5)*log((m+(m^2-1)^0.5)/(m-(m^2-1)^0.5)))/(2*(m^2-1))*4*pi;
    N = M;
else  
v=acos(c/a);
u=acos(b/a);
k=sin(u)/sin(v);
alp=asin(k); 
syms x;
F=ellipticF(v,k^2);
E=ellipticE(v,k^2);
L=cos(u)*cos(v)*(F-E)/(sin(v))^(3)/(sin(alp)^(2))*4*pi;
M=cos(u)*cos(v)*(E-(cos(alp))^(2)*F-(sin(alp))^(2)*sin(v)*cos(v)/cos(u))/(sin(v))^(3)/(sin(alp))^(2)/(cos(alp)^(2))*4*pi;
N=cos(u)*cos(v)*(sin(v)*cos(u)/cos(v)-E)/(sin(v))^(3)/(cos(alp)^(2))*4*pi;
end
bx=42.576*B*n(1)*(1-L/4/pi-2/3)*chi;
by=42.576*B*n(2)*(1-M/4/pi-2/3)*chi;
bz=42.576*B*n(3)*(1-N/4/pi-2/3)*chi;
b0 = bx*l(1) + by*l(2) + bz*l(3);
end