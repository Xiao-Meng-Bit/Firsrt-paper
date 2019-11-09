function [Margin,partial_M_z]=gradient_cal_quick(varargin)

% gradient calculation function
% If there is no input the function will give a demo of how to calculate
% gradient of Margin with defalt parameter and random value,
% else the inputs must have the following format [Hr,V,G,Hd,x,s] . And the
% dimension of matrices must be accordant.
% The output Margin is a M*1 vector where M is the number of users.
% The output partial_M_z is the  $\frac{\partial M}{\partial \tlide z}$
% whose dimension is M*Nt


para = varargin{1};
z = varargin{2};
um = varargin{3};
M = varargin{4};
L = varargin{5};
% Margin_right = varargin{6};
% partial_right=varargin{7};
PHI  = pi/L;
Margin = zeros(M,1);
partial_M_z = zeros(M,length(z));
% KK = sparse(length(z));
for m=1:M
    tm_imag =  z' * para.Ki{m}*z + para.Fi{m}*z;
    if tm_imag > 0 % K = Kn KK = KKn
%     K = Kr{m}*tan(PHI) - sign(tm_imag)*Ki{m} ;
        F = para.Fr{m} *tan(PHI) - para.Fi{m};
        Margin(m) = (z'*para.Kn{m}*z + F*z )*cos(PHI)-um*sin(PHI);
%     KK = (K+K');
    partial_M_z(m,:) = z'*para.KKn{m}+F;
    else        % K = Kp   KK  = KKp
        F = para.Fr{m} *tan(PHI) + para.Fi{m};
        Margin(m) = (z'*para.Kp{m}*z + F*z )*cos(PHI)-um*sin(PHI);
        partial_M_z(m,:) = z'*para.KKp{m}+F;
    end
end

