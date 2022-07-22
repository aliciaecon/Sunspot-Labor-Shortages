%% ========================================================================
% 
%  CONSTRUCTS TRUNCATED NORMAL VARIABLES FOR DISUTILITY IN SURPLUS FUNCTION
% 
% =========================================================================
function [phi, phi_d]  = understaff(xx, mu, sigma, nbar) 
 syms x
 Phi   = (1/sigma)*normpdf((x-mu)/sigma)/(normcdf((nbar-mu)/sigma) - normcdf((-mu)/sigma)); 
 phi_d = vpa(subs(diff(Phi),x,xx));
 phi   = vpa(subs(Phi,x,xx));
end

