function [x_dx,v_dv,i] = armijo_rule(gradient,x,v,Hr,G,Hd,s,L,um,epm)
t = 1;
alpha = 0.2;
beta = 0.9;
iter_max = 100;
dlt= -gradient;
fx = Margin_cal(Hr,G,Hd,x,v,s,L,um,epm);
flag = 0;
P = size(Hr,2);
Nt = size(Hd,2);
for i=1:iter_max
   t_dlt = t*dlt ;
   x_dx = x + t_dlt(P+1:P+Nt) + 1j*t_dlt(Nt+2*P+1:end);
   v_dv = v + t_dlt(1:P) + 1j * t_dlt(Nt+P+1:Nt+P+P);
   fx_dltx = Margin_cal(Hr,G,Hd,x_dx,v_dv,s,L,um,epm);
   if soft_cal(fx_dltx) < soft_cal(fx) + alpha * t *gradient' * dlt 
      step_size = t * dlt;
      if soft_cal(fx_dltx) > soft_cal(fx)
          warning
      end
      flag = 1;
      break
   else
        t = beta *t;
   end
end
if flag == 0
%     step_size = t * dlt_x;
%     error('break')
end




function y = soft_cal (M)
y = log( sum( exp(-M) ) );
