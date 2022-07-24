%% ========================================================================
% 
%                       CONSTRUCT EXOGENOUS MATRICES
% 
% =========================================================================
% For the evolution of productivity z, which is exogenous. Notice that
% upwinding dictates that for value function updating, when drift is
% positive forward difference is used, and when drift is negative backward
% difference is used. For distribution updating, upwinding dictates that
% when drift is positive backward difference is used, and when drift is
% negative forward difference is used.

function [ ZMatExog ] = MatricesExogenousZ( Params, ~, Dmatrices, ~ ) 

    % Load derivatives
% For the value function
DzBack          = Dmatrices.DzBack ;  % Backward approximation to derivative w.r.t. log z, for value function
Dzz             = Dmatrices.Dzz ;     % Second derivative w.r.t. log z, for value function
% For the the distribution
DzForwd         = Dmatrices.DzForwd ; % Forward approximation to derivative w.r.t. log z, for distribution
Dzzd            = Dmatrices.Dzzd ;    % Second derivative w.r.t. log z, for distribution

    % Construct matrices for the evolution of z
sigma   = Params.sigma ; % Productivity volatility
mu      = Params.mu ; % Drift of productivity
Zv      =  mu * DzBack  + sigma^2/2 * Dzz ;          % Matrix for value function updating
Zd      = -mu * DzForwd + sigma^2/2 * Dzzd ;         % Matrix for distribution updating

    % Build struct
ZMatExog.Zd = Zd ;
ZMatExog.Zv = Zv ;

end