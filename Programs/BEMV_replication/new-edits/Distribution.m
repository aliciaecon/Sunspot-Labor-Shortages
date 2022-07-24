%% ========================================================================
% 
%  UPDATE FIRM DISTRIBUTION GIVEN FINDING RATES AND WEIGHTED DISTRIBUTIONS
% 
% =========================================================================

function [ error, g_znmat, pi0Trunc_zn, exitrate, L, gprenorm_znmat ] = ...
    Distribution( q, phi, p, Gn_znmat, Gv_znmat, g_znmat, SUntrunc_znmat, SnUntrunc_znmat, v_znmat,...
    Numerical, NumGrids, Derivatives, MatExog, Params, ExogParams, options )
%%  LOAD NUMERICAL OBJECTS AND PARAMETERS
%__________________________________________________________________________

ExitDrift   = Numerical.ExitDrift; 	% Drift at exit boundary
LayoffDrift = Numerical.LayoffDrift;% Drift at layoff boundary
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid, in logs
Nz          = NumGrids.Nz;          % Number of points in z grid
Nn          = NumGrids.Nn;          % Number of points in n grid
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn as matrix, where d is log difference)
in0         = NumGrids.in0;         % Index of n gridpoint closest to n_0
pi0_z       = NumGrids.pi0_z;       % Distribution of entering firms in z grid
DnForwd     = Derivatives.DnForwd;  % Forward approximation to derivative w.r.t. log n, for distribution
DnBackd     = Derivatives.DnBackd;  % Backward approximation to derivative w.r.t. log n, for distribution
Zd          = MatExog.Zd;           % Matrix for distribution updating considering derivatives w.r.t. log z
delta       = Params.delta;         % Exogenous separation rate
xi          = Params.xi;            % Relative search efficiency of the employed
d           = ExogParams.d;         % Exogenous firm exit rate

%%  GET TRANSITION MATRIX AND ENTRY DISTRIBUTION
%__________________________________________________________________________

    % Drift of employment
ndotn_znmat                     = q * v_znmat ./ exp(n_znmat) .* ( phi + (1-phi)*Gn_znmat ) - delta - xi * p * (1-Gv_znmat); % drift of log labor, (dn/dt)/n
ndotn_znmat(SnUntrunc_znmat<0)  = ndotn_znmat(SnUntrunc_znmat<0) - LayoffDrift; % Layoffs when Sn is negative
ndotn_zn                    	= reshape( ndotn_znmat , Nz*Nn , 1 ); % Vector form
[ ndotnForw, ndotnBack, ~,~]    = FillSparse( ndotn_zn , NumGrids );
Nd                              = abs(ndotnBack') .* DnForwd + abs(ndotnForw') .* DnBackd; % Important to transpose and switch combination. Note that differentiation is w.r.t. log n

    % Add exit
iexit   = reshape( SUntrunc_znmat , Nz*Nn , 1 ) < 0;            % Indexes exiting firms
X       = ExitDrift * spdiags( iexit , 0 , Nz*Nn , Nz*Nn );     % Construct exit matrix X

    % Construct transition matrix excluding the exit state
L       = Nd + Zd - X - d*speye(Nz*Nn); % This is a slice of the true, larger transition matrix inclusive of the exit state

    % Construct distribution of entry
pi0_znmat                   = zeros( Nz , Nn ); 
pi0_znmat(:,in0)            = pi0_z;                                    % Entering distribution of productivity
pi0_znmat(SUntrunc_znmat<0) = 0;                                        % Truncate distribution: no firm enters with negative surplus
pi0Trunc_zn                 = reshape( pi0_znmat , Nz*Nn , 1 );         % Vector form
dzdn_zn                     = reshape( dzdn_znmat , Nz*Nn , 1);         % Vector form
pi0Trunc_zn                 = pi0Trunc_zn / ( pi0Trunc_zn' * dzdn_zn ); % ensure that this integrates to one

%%  GET NEW DISTRIBUTION CANDIDATE THROUGH KFE
%__________________________________________________________________________
gNew_zn         = - L \ ( 1 * pi0Trunc_zn) ; % Kolmogorov Forward Equation

    % Sanity checks
gprenorm_znmat  = reshape( gNew_zn , Nz , Nn ) ;        % Distribution pre-normalization
gNew_zn         = max( gNew_zn , 0 );                   % To ensure stability
gNew_zn         = gNew_zn / sum( gNew_zn .* dzdn_zn );  % Normalize to unit mass
gNew_znmat      = reshape( gNew_zn , Nz , Nn );         % Matrix form

%%  CALCULATE ERRORS AND UPDATE DISTRIBUTION
%__________________________________________________________________________

    % Assign errors as maximum absolute difference
errorg1 = max(max( abs(g_znmat                           - gNew_znmat                          ) )); % Error as level of pdf
errorg2 = max(max( abs(g_znmat.*dzdn_znmat               - gNew_znmat.*dzdn_znmat              ) )); % Error as measure
errorg3 = max(max( abs(g_znmat.*exp(n_znmat).*dzdn_znmat - gNew_znmat.*exp(n_znmat).*dzdn_znmat) )); % Error as measure weighted by employment
error   = max([ errorg1 errorg2 errorg3 ]);

    % Dampen the update
if options.HopenhaynNumerical == 0
    damp = 0.2; 
elseif options.HopenhaynNumerical == 1
    damp = 0.5;
end

    % Update distribution and get exit rate
g_znmat     = (1-damp) * gNew_znmat + damp * g_znmat;
g_znmat     = g_znmat / sum(sum( g_znmat .* dzdn_znmat)); % Ensure distribution integrates to 1
g_zn        = reshape( g_znmat , Nz*Nn , 1 ); % Vector form
exitrate    = sum(-L) * ( g_zn .* dzdn_zn ); % In a stationary equilibrium, exit = entry

end