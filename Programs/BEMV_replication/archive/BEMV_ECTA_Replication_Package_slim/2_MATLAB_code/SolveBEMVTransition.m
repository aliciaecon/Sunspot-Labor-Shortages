%% ========================================================================
% 
%   SOLVES THE TRANSITION PROBLEM FOR A PARTICULAR CHOICE OF PARAMETERS
% 
% =========================================================================

function [error, TransitionEq, TransitionEqOld] = SolveBEMVTransition...
    (Params, ExogParams, options, TransitionNum, TransitionEq, TransitionShock, SmoothPath)
%%  NUMERICAL OBJECTS
%__________________________________________________________________________

Numerical   = NumericalApproximations( Params, ExogParams, options );           % get settings for numerical approximations     -> (zmin,zmax,nmin,nmax)
NumGrids    = Grids( Numerical, Params ) ;                                      % define grids and useful objects from them     -> (Nz,Nn,z,n,dz_z,dn_n,dzdn_zn,dzdn_znmat,...)
Dmatrices   = MatricesDerivatives( NumGrids );                                  % numerical approximation to derivatives        -> (Dn,Dz,Dzz) [forward and backward versions]
Zmatrices   = MatricesExogenousZ( Params, NumGrids, Dmatrices , options ) ;     % construct exogenous (productivity) matrices   -> (Zv,Zd)

%%  LOAD SPECIFIC OBJECTS
%__________________________________________________________________________

ExitDrift   = Numerical.ExitDrift; 	% Drift at exit boundary
LayoffDrift = Numerical.LayoffDrift;% Drift at layoff boundary
d           = ExogParams.d;         % Exogenous firm exit rate
alpha       = Params.alpha ;        % Labor exponent in production function
delta       = Params.delta;         % Exogenous separation rate
Nz          = NumGrids.Nz;          % Number of points in z grid
Nn          = NumGrids.Nn;          % Number of points in n grid
z           = NumGrids.z;           % z grid in log space    (1 x Nz)
n           = NumGrids.n;           % n grid in log space    (1 x Nn)
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn in a matrix of size Nz * Nn)
in0         = NumGrids.in0;         % Index of n gridpoint closest to n_0
Nt          = TransitionNum.Nt;     % Number of periods in transition

%%  SETUP
%__________________________________________________________________________

    % Initialize distribution and surplus functions
PerPeriod.y_znmat  = kron( exp(z)' , exp(alpha*n) ) ;               % Flow output in levels, recall z and n are in logs
PerPeriod.yn_znmat = alpha * kron( exp(z)', exp((alpha-1)*n) ) ;    % Marginal product

    % Entry distribution
pi0_z         	= zeros( Nz , Nn ) ; 
pi0_z(:,in0)   	= NumGrids.pi0_z ;
pi0_z           = pi0_z / sum(pi0_z(:) .* dzdn_znmat(:)) ;

    % Loop setup
ExcessEntryValue    = zeros( Nt , 1 ) ;
EntryValue         	= zeros( Nt , 1 );
error              	= 1;
it                	= 1;

%%  MAIN LOOP
%__________________________________________________________________________

while error > Numerical.tolTransition && it < Numerical.maxiterTransition 

TransitionEqOld = TransitionEq ; % Save old time-dependent variables
t               = Nt - 1 ; % Backward pass to update the value functions

while t > 0
    %%  BACKWARDS RECURSION FOR SURPLUS FUNCTION
    %______________________________________________________________________
    
        % Load aggregates from last iteration and period t+1
    g_zn = TransitionEq.g(:,t+1) ;
    g_znmat = reshape( g_zn , Nz , Nn );
 
        % Load values from current iteration and period t+1
    S_znmat  = reshape( TransitionEq.S(:,t+1)  , Nz , Nn ) ;
    Sn_znmat = reshape( TransitionEq.Sn(:,t+1) , Nz , Nn ) ;
    
        % Compute Gn as a function of any marginal surplus
    % use last iteration's distribution h, at period t+1
    % use the CURRENT marginal surplus at period t
    [ ~ , Ghatn_znmat , inds , ~ ]  = Cdf_Gn( Sn_znmat , g_znmat , NumGrids ) ;
    GG                              = Ghatn_znmat(inds) ;
    SnGrid                          = Sn_znmat(inds) .* exp( -n_znmat(inds) ) ;
    
        % Solve for optimal individual behavior given aggregates
    % Can use exactly the same function as before
    % Now simply do not loop until fixed point, just evaluate it once
    [ S1, Sn1, SUntrunc1, SnUntrunc1, v1 ] = IndividualBehaviorTransition...
        ( S_znmat, Sn_znmat, GG, NumGrids, ExogParams, Params, Zmatrices, options,...
        PerPeriod, Dmatrices, SnGrid, TransitionNum, TransitionEq, TransitionShock, t ) ;
    
        % Save results and update time
    TransitionEq.S(:,t)         = S1(:) ;
    TransitionEq.Sn(:,t)        = Sn1(:) ;
    TransitionEq.SUntrunc(:,t)  = SUntrunc1(:) ;
    TransitionEq.SnUntrunc(:,t) = SnUntrunc1(:) ;
    TransitionEq.v(:,t)         = v1(:) ;
    t                           = t-1 ;
end

    % Compute excess value of entry and update mass of entrants
for t = 1:Nt
    EntryValue(t)       = sum( TransitionEq.S(:,t) .* pi0_z(:) .* dzdn_znmat(:) ) ;
    ExcessEntryValue(t) = EntryValue(t) / Params.entrycost - 1 ;
end
TransitionEq.ME = ExogParams.M0 * EntryValue.^ExogParams.EntryElast...
    ./ ( Params.entrycost^ExogParams.EntryElast + EntryValue.^ExogParams.EntryElast ) ; % Use logit rule

    % Update aggregates for period 1
[ q1, phi1, p1, ~, V1 ] = FindingRatesTransition(         ...
                            TransitionEqOld.q(1)        , ...
                            TransitionEq.phi(1)         , ...
                            TransitionEq.g(:,1)         , ...
                            TransitionEq.v(:,1)         , ...
                            TransitionEq.Gn(:,1)        , ...
                            TransitionEq.SUntrunc(:,1)  , ...
                            TransitionEq.MF(1)          , ...
                            TransitionEq.ME(1)          , ...
                            NumGrids, Params, ExogParams ) ;
TransitionEq.q(1)	= q1 ;
TransitionEq.phi(1)	= phi1 ;
TransitionEq.p(1)	= p1 ;
TransitionEq.V(1) 	= V1 ;
t = 2 ;

while t <= Nt
    %%  FORWARD RECURSION FOR DISTRIBUTION
    %______________________________________________________________________
    
        % Step size
    Delta = TransitionNum.dT(t) ;

        % Load policies from current backward pass and period t
    Sn_zn           = TransitionEq.Sn(:,t) ;
    SUntrunc_zn     = TransitionEq.SUntrunc(:,t) ;
    SnUntrunc_zn    = TransitionEq.SnUntrunc(:,t) ;
    v_zn            = TransitionEq.v(:,t) ;
    
        % Load distribution and finding rates from current iteration and period t-1
    g_zn    = TransitionEq.g(:,t-1) ;
    MassF1  = TransitionEq.MF(t-1) ;
    Gn      = TransitionEq.Gn(:,t-1) ;
    Gv      = TransitionEq.Gv(:,t-1) ;
    q       = TransitionEq.q(t-1)  ;
    phi     = TransitionEq.phi(t-1)  ;
    p       = TransitionEq.p(t-1) ;
    
        % Load mass of entrants from previous iteration and period t   
    MassE = TransitionEq.ME(t) ;
    
        % Reshape
    g_znmat         = reshape(g_zn,Nz,Nn);
    Gn_znmat        = reshape(Gn,Nz,Nn);
    Gv_znmat        = reshape(Gv,Nz,Nn);
    SUntrunc_znmat  = reshape(SUntrunc_zn,Nz,Nn);
    SnUntrunc_znmat = reshape(SnUntrunc_zn,Nz,Nn);
    Sn_znmat        = reshape(Sn_zn,Nz,Nn);
    v_znmat         = reshape(v_zn,Nz,Nn);
   
        % Compute distribution at time t
    [ g1, exitrate, L, MassF1 ] = DistributionTransition( q, phi, p, Gn_znmat,...
        Gv_znmat, g_znmat, SUntrunc_znmat, SnUntrunc_znmat, v_znmat, MassE,...
        MassF1, Delta, Numerical, NumGrids, Dmatrices, Zmatrices, Params, ExogParams ) ;
    
    TransitionEq.g(:,t)     = g1(:)  ; 
    TransitionEq.MF(t)      = MassF1 ;                       
    TransitionEq.L(:,t-1)   = L(:) ;
    
        % Update distribution of employment
    [ Gnnew, ~, inds, inds_inv ] = Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ;
    TransitionEq.Gn(:,t)         = Gnnew(:) ;

        % Update the distribution of poaching firms
    Gvnew = Cdf_Gv( inds, inds_inv, v_znmat, g_znmat, SUntrunc_znmat, exitrate, q, Gn_znmat, phi, NumGrids) ;
    TransitionEq.Gv(:,t) = Gvnew(:) ;

        % Update finding rates
    [ q1, phi1, p1, u1, V1 ] = FindingRatesTransition(q, phi, g1, v_znmat, Gn_znmat,...
                                SUntrunc_znmat, TransitionEq.MF(t), TransitionEq.ME(t),...
                                NumGrids, Params, ExogParams ) ;
    TransitionEq.u(t)        = u1 ;
    TransitionEq.q(t)        = q1 ;
    TransitionEq.phi(t)      = phi1 ;
    TransitionEq.p(t)        = p1 ;
    TransitionEq.V(t)        = V1 ;
    
        % Alternative computation of labor flows
    sep     = ( delta + d ) * ones( Nz*Nn , 1 ) ;
    sep( (SnUntrunc_zn<0) & ( reshape(repmat(1:Nn,Nz,1),Nn*Nz,1) >= 2 ) ) ...
            = sep( (SnUntrunc_zn<0) & ( reshape(repmat(1:Nn,Nz,1),Nn*Nz,1) >= 2 ) ) + LayoffDrift ;
    sep     = sep + ExitDrift * (SnUntrunc_zn<0) ;
    G       = g_zn .* dzdn_znmat(:) ; 
    G       = G / sum(G(:)) ;
	EUflow  = MassF1 * ( sep(:) .* exp(n_znmat(:)) )' * G(:) ;
    um1     = TransitionEq.u(t-1) ;
    dt      = TransitionNum.dT(t) ;
    
        % Save in struct
    TransitionEq.uflow(t)   = um1 + EUflow*dt - TransitionEq.p(t)*um1*dt ;
    TransitionEq.EUrate(t)  = EUflow / (1-um1) ; 
    
        % Go into next period
    t = t + 1 ;
end

    % Last transition matrix
TransitionEq.L(:,Nt) = TransitionEq.L(:,Nt-1) ;

    % Calculate errors
dzdn        = repmat(                    dzdn_znmat(:) , [Nt , 1] ) ;
ndzdn       = repmat( exp(n_znmat(:)) .* dzdn_znmat(:) , [Nt , 1] ) ;
error_u     = max( abs( TransitionEq.u              - TransitionEqOld.u )     ./  TransitionEqOld.u ) ;
error_q     = max( abs( TransitionEq.q              - TransitionEqOld.q )     ./  TransitionEqOld.q ) ;
error_phi   = max( abs( TransitionEq.phi            - TransitionEqOld.phi ) ./  TransitionEqOld.phi ) ;
error_p     = max( abs( TransitionEq.p              - TransitionEqOld.p )    ./  TransitionEqOld.p ) ;
error_Gn    = max( abs( TransitionEq.Gn(:)          - TransitionEqOld.Gn(:) ) ) ;
error_Gv    = max( abs( TransitionEq.Gv(:)          - TransitionEqOld.Gv(:) ) ) ;
error_g1    = max( abs( TransitionEq.g(:)           - TransitionEqOld.g(:) ) ) ;
error_g2    = max( abs( TransitionEq.g(:) .* dzdn   - TransitionEqOld.g(:)    .* dzdn(:) ) ) ;
error_g3    = max( abs( TransitionEq.g(:) .* ndzdn  - TransitionEqOld.g(:)    .* ndzdn(:) ) ) ;
error_mf    = max( abs( TransitionEq.MF(:)          - TransitionEqOld.MF(:) ) ./ TransitionEqOld.MF(:) ) ;
error_me    = max( abs( TransitionEq.ME(:)          - TransitionEqOld.ME(:) ) ./ ( 1e-4 + TransitionEqOld.ME(:) ) ) ;
error       = max( [ error_u error_q error_phi error_p error_Gn error_Gv error_g1 error_g2 error_g3 error_mf error_me ] ) ;
disp([ 'it ' num2str(it) ' ; error = ' num2str(error) ])
disp([ '   error_u  = ' num2str(error_u)  ])
disp([ '   error_q  = ' num2str(error_q)  ])
disp([ '   error_ph = ' num2str(error_phi)])
disp([ '   error_p  = ' num2str(error_p)  ])
disp([ '   error_Gn = ' num2str(error_Gn) ])
disp([ '   error_Gv = ' num2str(error_Gv) ])
disp([ '   error_g1 = ' num2str(error_g1) ])
disp([ '   error_g2 = ' num2str(error_g2) ])
disp([ '   error_g3 = ' num2str(error_g3) ])
disp([ '   error_mf = ' num2str(error_mf) ])
disp([ '   error_me = ' num2str(error_me) ])
disp('')

    % Update whole path
smoothPath          = SmoothPath.s;
x0                  = 1 - smoothPath ;
x1                  = smoothPath ;
TransitionEq.u      = x0 * TransitionEq.u   + x1 * TransitionEqOld.u ;
TransitionEq.q      = x0 * TransitionEq.q 	+ x1 * TransitionEqOld.q ;
TransitionEq.phi    = x0 * TransitionEq.phi + x1 * TransitionEqOld.phi ;
TransitionEq.p      = x0 * TransitionEq.p   + x1 * TransitionEqOld.p ;
TransitionEq.Gn     = x0 * TransitionEq.Gn  + x1 * TransitionEqOld.Gn ;
TransitionEq.Gv     = x0 * TransitionEq.Gv  + x1 * TransitionEqOld.Gv ;
TransitionEq.g      = x0 * TransitionEq.g   + x1 * TransitionEqOld.g ;
TransitionEq.MF     = x0 * TransitionEq.MF  + x1 * TransitionEqOld.MF ;
TransitionEq.ME     = x0 * TransitionEq.ME  + x1 * TransitionEqOld.ME ;

    % Next iteration
it = it + 1 ;
    
end

TransitionEq.ExcessEntryValue = ExcessEntryValue ;

end



