function [xc,theta,iter_v,iter_x]= A2_DAO3(Hr,Hd,G,u,k_max,sigma,x_raw,L)
[K,P] = size(Hr);     % Hr should be a matrix whose size is K*P
[~,M] = size(G);      % G  should be a matrix whose size is P*M
PHI=pi/L;
arg_u=pi/L+u*pi/2;    % psk mod angle
rotM=diag(exp(-1j.*arg_u));
theta = diag(ones(P,1));
% cvx_optval_i(1)=0;
xc=x_raw;
theta = diag(randn(P,1)+1j*randn(P,1));
for i=1:k_max
%% new
H=Hr*theta*G+Hd;
vc = diag(theta);
vc = opt_v(Hr,vc,G,Hd,xc,u,L,1);
theta = diag(vc);
iter_v(i) = min(real(rotM*H*xc)*tan(PHI)-abs(imag(rotM*H*xc))-1)*cos(PHI);


%%
H=Hr*theta*G+Hd;
xc = opt_x(Hr,vc,G,Hd,xc,u,L,1);
iter_x=min(real(rotM*H*xc)*tan(PHI)-abs(imag(rotM*H*xc))-1);

cvx_optval_i(i+1)=iter_x;
% disp(iter_x)
if abs((cvx_optval_i(i+1)-cvx_optval_i(i)) < sigma)
    break
end

end

end