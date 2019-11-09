function [x,theta,cput,cvx_optval_i]= A2_DAO(Hr,Hd,G,u,k_max,sigma,x_raw,L)
[K,P] = size(Hr);     % Hr should be a matrix whose size is K*P
[~,M] = size(G);      % G  should be a matrix whose size is P*M
PHI=pi/L;
arg_u=pi/L+u*pi/2;    % psk mod angle
rotM=diag(exp(-1j.*arg_u));
theta = diag(ones(P,1));
cvx_optval_i(1)=0;
x=x_raw;
for i=1:k_max
    
    cvx_begin quiet
    variable theta(P,P) diagonal complex
    maximize (min(real(rotM*(Hr*theta*G+Hd)*x)*tan(PHI)...
        -abs( imag(rotM*(Hr*theta*G+Hd)*x))...
        ))
    subject to
    abs(theta)<=1;
    cvx_end
    
    cput1(i)=cvx_cputime;
    theta=theta/abs(theta);
    H=Hr*theta*G+Hd;
    cvx_begin quiet
    variable x(M) complex
    maximize (min(real(rotM*H*x)*tan(PHI)-abs(imag(rotM*H*x))))
    subject to
    norm(x)<=1;
    cvx_end
    
    
    cvx_optval_i(i+1)=cvx_optval;
%     disp(cvx_optval)
    if abs((cvx_optval_i(i+1)-cvx_optval_i(i)) < sigma)
        break
    end
    
end
end