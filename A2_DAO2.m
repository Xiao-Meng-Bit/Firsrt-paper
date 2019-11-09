function [x,theta]= A2_DAO2(Hr,Hd,G,u,k_max,sigma,Pt,x_raw,L)
[K,P] = size(Hr);     % Hr should be a matrix whose size is K*P
[~,M] = size(G);      % G  should be a matrix whose size is P*M
PHI=pi/L;
arg_u=pi/L+u*pi/2;    % psk mod angle
rotM=diag(exp(-1j.*arg_u));
cvx_optval_i(1)=0;
x=x_raw;
theta = diag(randn(P,1)+1j*randn(P,1));
for i=1:k_max
%% new


vc = diag(theta);
vc = opt_v(Hr,vc,G,Hd,x,u,L,1);
theta = diag(vc);

%%
theta=theta/abs(theta);
H=Hr*theta*G+Hd;
% iter_v(i)=min(real(rotM*H*x)*tan(PHI)-abs(imag(rotM*H*x))-1)*cos(PHI);
cvx_begin quiet
    variable x(M) complex
    maximize (min(real(rotM*H*x)*tan(PHI)-abs(imag(rotM*H*x))))
    subject to 
    norm(x)<=Pt;
cvx_end


% iter_x(i)=min(real(rotM*H*x)*tan(PHI)-abs(imag(rotM*H*x))-1)*cos(PHI);
disp(cvx_optval);
cvx_optval_i(i+1)=cvx_optval;
if abs((cvx_optval_i(i+1)-cvx_optval_i(i)) < sigma)
    break
end

end

end