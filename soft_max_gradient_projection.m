function [x,V]=soft_max_gradient_projection(varargin)
% initialize
if isempty(varargin)
    Nt = 8;
    P  = 256;
    M  = 4;
    x  = randn(Nt,1) + 1j * randn(Nt,1);
    Hr = randn(M,P)  + 1j * randn(M,P) ;
    V  = diag( randn(P,1) + 1j * randn(P,1) );
    G  = randn(P,Nt) + 1j * randn(P,Nt);
    Hd = randn(M,Nt) + 1j * randn(M,Nt);
    s  = randi([0 3] , M ,1);
    um = 1;
    L  = 4;
    epm= pskmod(s,4,pi/L);
    PHI= pi/L;
    iter_max = 100;
else
    Hr = varargin{1};
    V  = varargin{2};
    G  = varargin{3};
    Hd = varargin{4};
    x  = varargin{5};
    s  = varargin{6};
    L  = varargin{7};
    um = varargin{8};
    iter_max = varargin{9};
    epm= pskmod(s,4,pi/4);
    PHI = pi/L;
    P = size(Hr,2);
    M = size(Hr,1);
    Nt = size(Hd,2);
end

% matrices transformation
v= diag(V);
    z = [v ; x];
    z = [real(z);imag(z)];
    Qpre = eye(P+Nt);
    Q(1:Nt,:) = Qpre (P+1:P+Nt,:);
    Q(Nt+1:P+Nt,:) = Qpre(1:P,:);
    A=zeros(P+Nt,P+Nt);
    for m=1:M
        Hrex= diag(Hr(m,:));
        A_pre = Hrex*G/epm(m);
        A(1:P,end-Nt+1:end) =A_pre;
        B = [Hd(m,:) zeros(1,P)]*Q/epm(m);
        Ar=real(A);Ai=imag(A);
        Br=real(B);Bi=imag(B);
        para.Kr{m} = [Ar -Ai; -Ai -Ar];para.Ki{m}=[Ai  Ar;  Ar -Ai];
        para.Fr{m} = [Br -Bi]         ;para.Fi{m}=[Bi  Br];
        para.Kp{m} =  para.Kr{m}*tan(PHI) + para.Ki{m} ;
        para.Kn{m} =  para.Kr{m}*tan(PHI) - para.Ki{m} ;
        para.KKp{m} = para.Kp{m} + (para.Kp{m})';
        para.KKn{m} = para.Kn{m} + (para.Kn{m})';
    end
    
    
    [Margin,partial_M_z]=gradient_cal_quick(para,z,um,M,L);
    part1               =exp(-Margin);
    soft_max_gradient   = 1/sum(part1)*(-partial_M_z'*exp(-Margin));
    for i=1:iter_max        
        [x,v] = armijo_rule(soft_max_gradient/norm(soft_max_gradient),x,v,Hr,G,Hd,s,L,um,epm);
        x = x/norm(x);
        v = v./abs(v);
        z = [real(v) ;real(x) ;imag(v); imag(x)];
        z = [real(v) ;real(x) ;imag(v); imag(x)];
        [Margin,partial_M_z]= gradient_cal_quick(para,z,um,M,L);
        part1               = exp(-Margin);
        soft_max_new(i)     = log(sum(exp(-Margin)));
        soft_max_gradient   = 1/sum(part1)*(-partial_M_z'*exp(-Margin));
    end
    V = diag(v);