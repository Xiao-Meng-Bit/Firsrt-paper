function xc = opt_x(varargin)


if isempty(varargin)
    Nt = 8;
    P  = 256;
    M  = 4;
    x  = [randn(Nt,1) ;  randn(Nt,1)];
    Hr = randn(M,P)  + 1j * randn(M,P) ;
    v  = [randn(P,1) ;  randn(P,1)];
    G  = randn(P,Nt) + 1j * randn(P,Nt);
    Hd = randn(M,Nt) + 1j * randn(M,Nt);
    s  = randi([0 3] , M ,1);
    um = 1;
    L  = 4;
    epm= pskmod(s,4,pi/L);
    PHI= pi/L;
    iter_max = 100;
    xc = randn(Nt,1)+1j*randn(Nt,1);
    vc = randn(P,1)+1j*randn(P,1);
else
    Hr = varargin{1};
    vc  = varargin{2};
    G  = varargin{3};
    Hd = varargin{4};
    xc  = varargin{5};
    s  = varargin{6};
    L  = varargin{7};
    um = varargin{8};
    epm= pskmod(s,4,pi/4);
    PHI = pi/L;
    P = size(Hr,2);
    M = size(Hr,1);
    Nt = size(Hd,2);
end
Qpre = eye(P+Nt);
Q(1:Nt,:) = Qpre (P+1:P+Nt,:);
Q(Nt+1:P+Nt,:) = Qpre(1:P,:);
% A=zeros(P+Nt,P+Nt);
x = [real(xc) ; imag(xc)];
v = [real(vc) ; imag(vc)];
para.v = v;
para.epm= epm;
para.s=s;
para.L=L;
para.um=um;
for m=1:M
    Hrex= diag(Hr(m,:));
    A = Hrex*G/epm(m);
    %     A(1:P,end-Nt+1:end) =A_pre;
    B = [Hd(m,:)]/epm(m);
    Ar=real(A);Ai=imag(A);
    Br=real(B);Bi=imag(B);
    para.Kr{m} = [Ar -Ai; -Ai -Ar];para.Ki{m}=[Ai  Ar;  Ar -Ai];
    para.Fr{m} = [Br -Bi]         ;para.Fi{m}=[Bi  Br];
    para.Kp{m} =  para.Kr{m}*tan(PHI) + para.Ki{m} ;
    para.Kn{m} =  para.Kr{m}*tan(PHI) - para.Ki{m} ;
    para.vKn{m} = v'*para.Kn{m};
    para.vKp{m} = v'*para.Kp{m};
    %     para.KKp{m} = para.Kp{m} + (para.Kp{m})';
    %     para.KKn{m} = para.Kn{m} + (para.Kn{m})';
end
vc = v(1:end/2) + 1j * v(end/2+1:end);
xc = x(1:end/2) + 1j * x(end/2+1:end);
para.Hr = Hr;
para.Hd = Hd;
para.G  = G;
para.vr_ind = 1:length(v)+1:length(v)^2;
para.vc_ind = 1:length(vc)+1:(length(vc))^2;
tm=((Hr*diag(vc)*G*xc+Hd*xc)./epm-1);
MM=(real(tm)*tan(PHI)-abs(imag(tm)))*cos(PHI);
for i=1:200
    [Margin,pMpx]=gdx(para,x,v,um,M,L,epm);
    obj(i)= log(sum(exp(-Margin)));
    part1               =exp(-Margin);
    gd = 1/sum(part1)*(-pMpx'*exp(-Margin));
    x = -gd/norm(gd)*0.1 +x;
    xc = x(1:end/2) + 1j * x(end/2+1:end);
    xc= xc/norm(xc);
    x = [real(xc) ; imag(xc)];

end
tm=((Hr*diag(vc)*G*xc+Hd*xc)./epm-1);
MM=(real(tm)*tan(PHI)-abs(imag(tm)))*cos(PHI);
