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