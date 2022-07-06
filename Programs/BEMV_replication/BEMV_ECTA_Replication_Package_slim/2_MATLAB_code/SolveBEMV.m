%% ========================================================================
% 
%       SOLVES THE SS PROBLEM FOR A PARTICULAR CHOICE OF PARAMETERS
% 
% =========================================================================

%   (1) File InitialGuess.m constructs an initial guess for the surplus
%       function and distribution;
%   
%   (2) File IndividualBehavior.m takes as given aggregate outcomes--a
%       worker finding rate, q, a share of unemployed job
%       seekers, phi, and distribution of workers Gn--and solves the
%       coalition's problem to obtain separation, exit and vacancy
%       policies. It does so by iterating to convergence on the surplus
%       function using an implicit method;
%
%   (3) File AggregateBehavior.m updates the aggregate state given optimal
%       coalition behavior by first updating the distribution of workers 
%       ranked by marginal surplus (file Cdf_Gn.m); then the distribution 
%       of vacancies ranked by marginal surplus (file Cdf_Gv.m); then the 
%       distribution of firms g over productivity and size (file 
%       Distribution.m), and finally the finding rates (done by file
%       FindingRates.m). It iterates on this until these aggregate outcomes 
%       have converged for given coalition behavior in (2);
%
%   (4) We subsequently return to (2) given the updated aggregate
%       outcomes from (3), and go back and forth between (2)-(3) until both
%       the surplus function S, the hazard rates q, p and phi, the
%       distribution of workers Gn, the distribution of vacancies Gv, and the
%       distribution of firms g have converged.

function [error_Ind, error_Agg, error, SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat,...
    v_znmat, g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, NumGrids, Dmatrices, Numerical, it]...
    = SolveBEMV( Params, ExogParams, options )
%%  NUMERICAL OBJECTS
%__________________________________________________________________________

Numerical   = NumericalApproximations( Params, ExogParams, options );           % get settings for numerical approximations     -> (zmin,zmax,nmin,nmax)
NumGrids    = Grids( Numerical, Params ) ;                                      % define grids and useful objects from them     -> (Nz,Nn,z,n,dz_z,dn_n,dzdn_zn,dzdn_znmat,...)
Dmatrices   = MatricesDerivatives( NumGrids );                                  % numerical approximation to derivatives        -> (Dn,Dz,Dzz) [forward and backward versions]
Zmatrices   = MatricesExogenousZ( Params, NumGrids, Dmatrices , options ) ;   % construct exogenous (productivity) matrices   -> (Zv,Zd)

%%  LOAD SPECIFIC OBJECTS
%__________________________________________________________________________

alpha       = Params.alpha ;       % Labor exponent in production function
z           = NumGrids.z;          % z grid in log space    (1 x Nz)
n           = NumGrids.n;          % n grid in log space    (1 x Nn)
n_znmat     = NumGrids.n_znmat;    % Replica of the n grid
dzdn_znmat  = NumGrids.dzdn_znmat; % Measure (dz*dn in a matrix of size Nz * Nn)
tol         = Numerical.tol;       % Tolerance for overall convergence
maxiter     = Numerical.maxiter;   % Max total number of iterations
tol_v       = Numerical.tol_v;     % Tolerance for convergence of value
maxiter_v   = Numerical.maxiter_v; % Max number of iterations for value
tol_d       = Numerical.tol_d;     % Tolerance for overall distributions
maxiter_d   = Numerical.maxiter_d; % Max number of iterations for distribution
Delta       = Numerical.Delta;     % Time step size
maxtime     = Numerical.maxtime;   % Max time in seconds

%%  INITIAL GUESS
%__________________________________________________________________________

PerPeriod.y_znmat  = kron( exp(z)' , exp(alpha*n) ) ;               % Flow output in levels, recall z and n are in logs
PerPeriod.yn_znmat = alpha * kron( exp(z)', exp((alpha-1)*n) ) ;    % Marginal product

    % Initialize distribution and surplus functions
if options.Hopenhayn == 0
    [ g_znmat, u, p, q, phi, Gn_znmat, ~, Gv_znmat, S_znmat, Sn_znmat, v_znmat, exitrate ] ...
        = InitialGuess( ExogParams, NumGrids, Params, Zmatrices, PerPeriod ) ;
elseif options.Hopenhayn == 1
    g_znmat     = ExogParams.InitialGuess.g_znmat ;
    u           = ExogParams.InitialGuess.u ;
    p           = ExogParams.InitialGuess.p ;
    q           = ExogParams.InitialGuess.q ;
    phi         = ExogParams.InitialGuess.phi ;
    Gn_znmat    = ExogParams.InitialGuess.Gn_znmat ;
    Gv_znmat    = ExogParams.InitialGuess.Gv_znmat ;
    S_znmat     = ExogParams.InitialGuess.S_znmat ;
    Sn_znmat    = ExogParams.InitialGuess.Sn_znmat ;
    v_znmat     = ExogParams.InitialGuess.v_znmat ;
    exitrate    = ExogParams.InitialGuess.exitrate ;
else
    warning('options.Hopenhayn must be either 0 or 1.')
end

%%  MAIN LOOP
%__________________________________________________________________________

    % Setup
error    = 1; % distance to convergence for entire algorithm
it       = 1; % iteration for entire algorithm
timeBEMV = tic ; % solution time

while error > tol && it < maxiter % In the paper, each of these iterations is "t"
    %%  STEP I - SOLVE FOR INDIVIDUAL BEHAVIOR GIVEN AGGREGATES (INCLUDING DISTRIBUTIONS)
    %______________________________________________________________________
    
        % Setup
    [ Gn_znmat, Ghatn_znmat, inds, ~ ] = Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ; % CDF weighted by employment
    S_old_znmat     = S_znmat ;                                 % Surplus function from previous main iteration
    GhatnSort_zn    = Ghatn_znmat(inds) ;                       % CDF sorted by marginal surplus
    SnGrid_zn       = Sn_znmat(inds) .* exp(-n_znmat(inds)) ;   % Set up a sorted grid for level Sn (need to divide by n because Sn is differentiated wrt logs) 
    % Note: Grid needed because the iterations over the individual problem
    % change Sn, and for convergence we need a fixed grid
    error_v = 1 ; % distance to convergence for value (surplus) function
    it_v 	= 1 ; % iteration for value function loop
    
        % Loop until individual behavior (value function problem) converges
    while error_v > tol_v && it_v < maxiter_v
        
            % Iterate once over surplus function (in the paper, iteration "tau")
        [ error_v , S_znmat , Sn_znmat , SUntrunc_znmat , SnUntrunc_znmat , v_znmat ] = ...
            IndividualBehavior( S_znmat, Sn_znmat, q, phi, GhatnSort_zn, SnGrid_zn, Delta, NumGrids,...
                                    ExogParams, Params, Zmatrices, PerPeriod, Dmatrices, options ) ;
        % Note: Sn numerically differentiates S w.r.t. log(n)
        
        it_v = it_v + 1 ;
    end
    
        % Compute Individual error, relative to previous main iteration
    error_Ind = max(max( abs(S_znmat - S_old_znmat) )) / sum(sum(S_old_znmat));
     
    %% STEP II - SOLVE FOR AGGREGATES GIVEN INDIVIDUAL BEHAVIOR
    %______________________________________________________________________
    
        % Setup
    error_d         = 1;
    it_d            = 1;
    g_old_znmat     = g_znmat;
    Gv_old_znmat    = Gv_znmat;
    Gn_old_znmat    = Gn_znmat;
    p_old           = p;
    phi_old         = phi;
    q_old           = q;
    u_old           = u;
    
        % Loop until distributions converge
    while error_d > tol_d && it_d < maxiter_d
        
            % Iterate once over aggregates (in the paper, iteration "tau")
        [ error_g, error_f, error_Gn, error_Gv, q, phi, p, u, Gn_znmat, ~,...
            Gv_znmat, g_znmat, exitrate, L, pi0Trunc_zn, ~ ] ...
            = AggregateBehavior( q, phi, p, Gn_znmat, Gv_znmat, g_znmat, SUntrunc_znmat,...
                                   Sn_znmat, SnUntrunc_znmat, v_znmat, exitrate, Numerical, NumGrids,...
                                   Dmatrices, Zmatrices, Params, ExogParams, options ) ;
        
            % Get error and update iteration "tau"
        error_d = max([ error_g error_f error_Gn error_Gv ]) ;
        it_d = it_d + 1 ;
    end
    
        % Compute errors relative to previous main iteration
    error_q     = abs( q   - q_old )   / q_old ;
    error_phi   = abs( phi - phi_old ) / phi_old ;
    error_p     = abs( p   - p_old )   / p_old ;
    error_Gn    = max(max( abs( Gn_znmat                          - Gn_old_znmat                          ) )) ;
    error_Gv    = max(max( abs( Gv_znmat                          - Gv_old_znmat                          ) )) ;
    error_g1    = max(max( abs( g_znmat                           - g_old_znmat                           ) )) ;
    error_g2    = max(max( abs( g_znmat.*dzdn_znmat               - g_old_znmat.*dzdn_znmat               ) )) ;
    error_g3    = max(max( abs( g_znmat.*exp(n_znmat).*dzdn_znmat - g_old_znmat.*exp(n_znmat).*dzdn_znmat ) )) ;
    
        % Aggregate error
    error_Agg = max([ error_q error_phi error_p error_Gn error_Gv error_g1 error_g2 error_g3 ]) ;
    
    %%  CLOSE MAIN ITERATION "t"
    %______________________________________________________________________
    
        % Error of main iteration
    error = max( [ error_Agg , error_Ind ] );
    if options.Hopenhayn == 2
        error = max( [ error_q , error_phi , error_p ] );
    end
    it = it + 1 ;
    
        % Choose dampening parameter
    damp = 0.25 ; % Baseline dampening
    if it > 50
        damp = 0.75 ; % Increased dampening in case algorithm is not converging
    end
    
        % Update aggregates
    g_znmat = (1-damp) * g_znmat + damp * g_old_znmat ;
    q       = (1-damp) * q       + damp * q_old ;
    p       = (1-damp) * p       + damp * p_old ;
    phi     = (1-damp) * phi     + damp * phi_old ;
    u       = (1-damp) * u       + damp * u_old ;
    
        % Time from the start of the algorithm
    time = toc(timeBEMV) ;
    if time > maxtime
        error = - 1 ;
    end
end

end