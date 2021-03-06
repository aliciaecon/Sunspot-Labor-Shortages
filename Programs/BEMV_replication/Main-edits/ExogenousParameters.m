%% ========================================================================
% 
%                       DEFINES PRE-SET PARAMETERS
% 
% =========================================================================

ExogParams.rho   = log(1+.05)/12 ; % Discount rate
ExogParams.beta  = 0.5 ;           % Matching function exponent on vacancies
ExogParams.MF    = (1-0.05)/22 ;   % Fix mass of firms
ExogParams.cbar  = 100 ;           % Scalar in vacancy posting cost
ExogParams.gamma = 3.45 ;          % Curvature of vacancy post cost
ExogParams.d     = 0.0002 ;        % Exogenous firm exit rate

ExogParams.chi         = 0 ;      % Exogenous disutility
ExogParams.nbar        = 1000 ;    % Upper bound of truncated normal 
ExogParams.mu_u        = 0;        % Mean of truncated normal 
ExogParams.sigma_u     = 10 ;      % St dev of truncated normal 

% plot the functions to make sure they make sense
xgrid=linspace(0,ExogParams.nbar/2,250);
for i=1:size(xgrid,2)
    [x(i), y(i)] = understaff(xgrid(i), ExogParams.mu_u, ExogParams.sigma_u, ExogParams.nbar);
end
plot(xgrid,ExogParams.chi*x)
hold on
plot(xgrid,ExogParams.chi*y)


