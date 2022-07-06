%% ========================================================================
% 
%             DEFINES SETTINGS FOR NUMERICAL APPROXIMATIONS
% 
% =========================================================================

function Numerical = NumericalApproximations( Params, ExogParams, options )
%%  LOAD PARAMETERS
%__________________________________________________________________________

zeta    = Params.zeta ;     % Productivity distribution of entrants
sigma   = Params.sigma ;    % Productivity volatility
alpha   = Params.alpha ;    % Labor exponent in production function
b       = Params.b ;        % Unemployment benefit
rho     = ExogParams.rho ;  % Time discount rate
mu      = Params.mu ;       % Drift of productivity

%%  BASELINE SETTINGS
%__________________________________________________________________________

    % Basics
Numerical.maxtime       = 240 ; % Max time in seconds
Numerical.maxiter       = 200 ; % Max total number of iterations
Numerical.maxiter_d     = 200 ; % Max number of iterations for distribution
Numerical.maxiter_v     = 200 ; % Max number of iterations for value functions
Numerical.ExitDrift     = 1e2 ; % Drift at exit boundary ; better not to exceed 1e2 or 1e3
Numerical.LayoffDrift   = 1e2 ; % Drift at layoff boundary ; better not to exceed 1e2 or 1e3
Numerical.Nn            = 50 ;  % Gridpoints in n dimension; should be weakly larger than Nz
Numerical.Nz            = 50 ;  % Gridpoints in z dimension

    % Tolerances for convergence
Numerical.tol   = 1e-3; % Tolerance for overall convergence
Numerical.tol_d = 1e-4; % Tolerance for convergence of distributions
Numerical.tol_v = 1e-4; % Tolerance for convergence of value

    % Time step size
Numerical.Delta = 10 ; % Smoothing parameter in HJB -> *** NOTE: THIS CAN EITHER BE TOO LARGE NOR TOO SMALL (for algorithm to converge) ***

    % Discretionary settings
factor      = 100 ;  % Scales maximum Zgrid point, and hence also maximum Ngrid point
Zmin        = 0.1 ;  % Lower bound on z grid
Nmax_uppcap = 1e8 ;  % Upper cap on maximum Ngrid point
Nmin_uppcap = 0.01 ; % Upper cap on minimum Ngrid point
Nmin_lowcap = 0.1 ;  % Lower cap on minimum Ngrid point

    % Previous calculations
EntryBound      = exp( -log(0.001)/zeta ) ;
kappa           = -2 * mu / sigma^2 ;
InvariantBound  = min( exp( -log(0.001)/kappa ) , 10000 ) ;
Zmax            = factor * max(EntryBound,InvariantBound) ; % From the a priori invariant productivity distribution and entry
Nmax            = ( alpha * Zmax / b )^(1/(1-alpha)) ;      % From maximum static unconstrained size
Nmin            = alpha / (rho*10^3 + (1-alpha)*b) * (1-alpha)^(1/(1-alpha)) ;

    % Grid limits (now all in logs)
Numerical.zmin  = log(Zmin) ; % Lower bound on z grid
Numerical.zmax  = log(Zmax) ; % Upper bound on z grid
Numerical.nmin  = log(max(Nmin_lowcap,min(Nmin,Nmin_uppcap))); % Lower bound on n grid
Numerical.nmax  = log(min(Nmax,Nmax_uppcap)); % Upper bound on n grid

    % Changes for hopenhayn limit
if options.HopenhaynNumerical == 1
    Numerical.maxtime   = 120;
    Numerical.Delta     = 1 ;
end

%%  TRANSITION SETTINGS
%__________________________________________________________________________

if options.Transition == 1
    
        % Transition-specific
    Numerical.tolTransition     = 5 * 1e-2 ;
    Numerical.maxiterTransition = 40 ; 
    
        % Basics
    Numerical.tol          	= 1e-3; % Tolerance for overall convergence
    Numerical.tol_d       	= 1e-4; % Tolerance for convergence of distributions
    Numerical.tol_v        	= 1e-4; % Tolerance for convergence of value
    Numerical.maxiter      	= 200 ; % Max total number of iterations
    Numerical.maxiter_d    	= 200 ; % Max number of iterations for distribution
    Numerical.maxiter_v   	= 200 ; % Max number of iterations for value functions
    Numerical.ExitDrift     = 1e2 ; % Drift at exit boundary ; better not to exceed 1e2 or 1e3
    Numerical.LayoffDrift   = 1e2 ; % Drift at layoff boundary ; better not to exceed 1e2 or 1e3
    Numerical.Nn           	= 50  ; % Gridpoints in n dimension; should be weakly larger than Nz
    Numerical.Nz           	= 50  ; % Gridpoints in z dimension
    
        % Time step size
    Numerical.Delta = 1 ; % Smoothing parameter in HJB -> *** NOTE: THIS CAN NEITHER BE TOO LARGE NOR TOO SMALL (for algorithm to converge) ***
    
        % Previous calculations
    factor          = 100 ;
    EntryBound      = exp( - log(0.001) / Params.zeta ) ;
    kappa           = - 2 * Params.mu / Params.sigma^2 ;
    InvariantBound  = min( exp( - log(0.001) / kappa ) , 10000 ) ;
    Zmax            = factor * max(EntryBound,InvariantBound) ;
    Nmax            = ( Params.alpha * Zmax / Params.b )^(1/(1-Params.alpha)) ;
    Nmin            = Params.alpha / (ExogParams.rho * 10^3 + (1-Params.alpha)*Params.b) ;
    FactorNmin      = (1-Params.alpha)^( 1 / (1-Params.alpha) ) ;
    factor          = 1 ;

        % Grid limits
    Numerical.zmin = log(0.1);                                                      % Lower bound on z grid
    Numerical.zmax = log( factor * Zmax ) ;                                         % Upper bound on z grid, based on invariant distribution of productivity process and entry of productivities
    Numerical.nmin = log( max( 1e-1 , factor * min( FactorNmin * Nmin , 0.01 ) ) ); % Lower bound on n grid
    Numerical.nmax = log( min( factor * Nmax , factor * 1e8 ) );                   	% Upper bound on n grid, based on z n^(1-alpha) = b
end

end