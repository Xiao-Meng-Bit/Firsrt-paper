function [Margin,pMpx]=gdx(para,x,v,um,M,L,epm)
%Margin = ([vr vi] * [Kp or Kn] * [xr;xi] + [Fp or Fn]*x - um)*cos(PHI)
%partial Margin/ partial v = transpose([Kp or Kn]*[xr;xi]) = transpose(Kpx or Knx)
Margin = zeros(M,1);
pMpx = zeros(M,length(x));
VC=zeros(length(v)/2);
xc = x(1:end/2) + 1j*x(end/2+1:end);
vc = v(1:end/2) + 1j*v(end/2+1:end);
VC(para.vc_ind) =vc;
PHI = pi/L;
tm_imag = imag(((para.Hr*VC*para.G)*xc + para.Hd*xc)./epm);
for m=1:M
    if tm_imag(m) > 0 
        Margin(m) = (para.vKn{m}*x + (para.Fr{m}*tan(PHI)-para.Fi{m})*x - um)*cos(pi/L);
        pMpx(m,:) = para.vKn{m}+(para.Fr{m}*tan(PHI)-para.Fi{m});
    else
        Margin(m) = (para.vKp{m}*x + (para.Fr{m}*tan(PHI)+para.Fi{m})*x - um)*cos(pi/L);
        pMpx(m,:) = para.vKp{m}+(para.Fr{m}*tan(PHI)+para.Fi{m});
    end
end