function [x,V] =  CI_R(Hr,Hd,G,u,order)

H=Hr*G+Hd;
PHI = pi/order;
arg_u=pi/order+u*pi/2;  % psk mod angle
rotM=diag(exp(-1j.*arg_u));
N = size(Hd,2);
cvx_begin quiet
    variable x_raw(N) complex
    maximize (min(real(rotM*H*x_raw)*tan(PHI)-abs(imag(rotM*H*x_raw))))
    subject to 
    norm(x_raw)<=1;
cvx_end
x = x_raw;
V= eye(size(Hr,2));