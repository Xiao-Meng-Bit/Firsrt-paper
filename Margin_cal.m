function [Margin] = Margin_cal(varargin)

Hr = varargin{1};
G  = varargin{2};
Hd = varargin{3};
x  = varargin{4};
v  = varargin{5};
s  = varargin{6};
L  = varargin{7};
um = varargin{8};
epm= varargin{9};
PHI = pi/L;

% P = size(Hr,2);
% Nt= size(Hd,2);
% x = z(P+1:P+Nt) + 1j * z(P+Nt+P+1:end);
% v = z(1:P)      + 1j * z(P+Nt+1:P+Nt+P);
V=diag(v);
y = Hr*V*G*x+ Hd* x;
y_rot = y./epm;
%           real part                          imag part          PHI factor
Margin =( (real(y_rot)-abs(um))*tan(PHI)  -  abs(imag(y_rot)) ) * cos(PHI);
