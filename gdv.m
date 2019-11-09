function [Margin,pMpv]=gdv(para,x,v,um,M,L,epm)
%Margin = ([vr vi] * [Kp or Kn] * [xr;xi] + [Fp or Fn]*x - um)*cos(PHI)
%partial Margin/ partial v = transpose([Kp or Kn]*[xr;xi]) = transpose(Kpx or Knx)
Margin = zeros(M,1);
pMpv = zeros(M,length(v));
VC=zeros(length(v)/2);
xc = x(1:end/2) + 1j*x(end/2+1:end);
vc = v(1:end/2) + 1j*v(end/2+1:end);
VC(para.vc_ind) =vc;
PHI = pi/L;
tm_imag = imag(((para.Hr*VC*para.G)*xc + para.Hd*xc)./epm);
for m=1:M
%     tm_imag = v'*para.Ki{m}*x + v'*para.Fi{m}*x;
    if tm_imag(m) > 0 
        Margin(m) = (v' * para.Knx{m} + (para.Fr{m}*tan(PHI)-para.Fi{m})*x - um)*cos(pi/L);
        pMpv(m,:) = transpose(para.Knx{m});
    else
        Margin(m) = (v' * para.Kpx{m} + (para.Fr{m}*tan(PHI)+para.Fi{m})*x - um)*cos(pi/L);
        pMpv(m,:) = transpose(para.Kpx{m});
    end
end