%% ========================================================================
% 
%          UPDATES SURPLUS FUNCTIONS ONLY ONCE IN VFI
% 
% =========================================================================

function [ error, S_znmat, Sn_znmat, SUntrunc_znmat, SnUntrunc_znmat, v_znmat ] = ...
    IndividualBehavior( S_znmat, Sn_znmat, q, phi, GhatnSort_zn, SnGrid_zn, Delta, NumGrids,...
    ExogParams, Params, MatExog, PerPeriod, Derivatives, options )
%%  LOAD NUMERICAL OBJECTS AND PARAMETERS
%__________________________________________________________________________

Nz      = NumGrids.Nz;        % Number of points in z grid
Nn      = NumGrids.Nn;        % Number of points in n grid
n_znmat = NumGrids.n_znmat;   % Replica of the n grid in logs
n_zn    = NumGrids.n_zn;      % Replica of the n grid in logs in vector
cbar    = ExogParams.cbar;    % Scalar in vacancy posting cost
gamma   = ExogParams.gamma;   % Curvature of vacancy post cost
rho     = ExogParams.rho;     % Time discount rate
d       = ExogParams.d;       % Exogenous firm exit rate
delta   = Params.delta;       % Exogenous separation rate
b       = Params.b;           % Unemployment benefit
c_f     = Params.c_f;         % Fixed cost of operation
y_znmat = PerPeriod.y_znmat;  % Flow output
DnForw  = Derivatives.DnForw; % Forward approximation to derivative w.r.t. logn, for value function
DnBack  = Derivatives.DnBack; % Backward approximation to derivative w.r.t. logn, for value function
Zv      = MatExog.Zv ;        % Matrix for value function updating considering derivatives w.r.t. logz

%%  SETUP
%__________________________________________________________________________

    % Save old surplus functions and get vector forms
S0_znmat    = S_znmat ; 
Sn0_znmat   = Sn_znmat ; % Recall this differentiates S by logn
S_zn        = reshape( S_znmat , Nz*Nn , 1 ) ;
Sn_zn       = reshape( Sn_znmat , Nz*Nn , 1 ) ;

    % Interpolate Ghatn and recover unsorted integral
[ iBot, iTop, ~,~, frac ]   = Brackets2fast( SnGrid_zn, Sn_zn.*exp(-n_zn), Nz*Nn, 1, Nz*Nn ) ; % Locate indices on SnGrid such that Sn lies in between, and specify position (frac)
frac                        = min( frac , 1 ) ; % Sanity check
frac                        = max( frac , 0 ) ; % Sanity check
Ghatn_zn                    = (1-frac) .* GhatnSort_zn(iBot) + frac .* GhatnSort_zn(iTop) ; % Recover unsorted integral from sorted one
Ghatn_znmat                 = reshape( Ghatn_zn , Nz , Nn ) ; % Matrix form
Ghatn_znmat                 = min( Ghatn_znmat , Sn_znmat ) ; % Sanity check
Ghatn_znmat                 = max( Ghatn_znmat , 0 ) ;% Sanity check
% Note: must use Sn/n to get S differentiated by n instead of logn, so use
% .*exp(-n_zn) in code, as n grid is in logs

%%  GET OBJECTS FOR UPDATING SURPLUS FUNCTION
%__________________________________________________________________________

    % Get the hiring rate and vacancy policy
        % Note that as we are in logs, then Sn(z,n)*e^(-n) computes Sn(z,n) in levels
H_znmat = q * ( phi * (Sn_znmat.*exp(-n_znmat)) + (1-phi) * Ghatn_znmat ) ;     % Hires per posted vacancy
v_znmat = cbar^(-1/gamma) * H_znmat.^(1/gamma) .* exp(n_znmat) ;                % Vacancy policy function

    % Drift of labor (minus vacancy cost)
xi          = cbar^(-1/gamma) * gamma / (1+gamma);
mu_n0_znmat = xi * H_znmat.^((1+gamma)/gamma) ./ (Sn_znmat.*exp(-n_znmat)) - delta ; % Original formula
mu_n1_znmat = xi * H_znmat.^(1/gamma) .* exp(n_znmat) .* q .* (phi .* exp(-n_znmat)) - delta ; % Sn=0 is possible, mu_n1 avoids dividing by Sn for numerical purposes

    % Numerical adjustments
mu_n_znmat = mu_n0_znmat ;
if options.HopenhaynNumerical == 0
    mu_n_znmat(Sn_znmat <= 0) =  - 1000 ;
elseif options.HopenhaynNumerical == 1
    mu_n_znmat(Sn_znmat <= 0) =  - 1e4 ;
end
mu_n_znmat(0 < Sn_znmat & Sn_znmat < 1e-4) = mu_n1_znmat(0 < Sn_znmat & Sn_znmat < 1e-4) ;

% Note: DnForw has  values on the diagonal and above the diagonal. DnBack 
% has values on the diagonal and below the diagonal.
    % Define objects for updating in implicit scheme
mu_n_zn                     = reshape ( mu_n_znmat, Nz*Nn , 1 ) ;               % Vector form
[ mu_nForw, mu_nBack, ~,~ ] = FillSparse( mu_n_zn , NumGrids ) ;                % Reshape mu_n into sparse matrices
Nv                          = mu_nForw .* DnForw + mu_nBack .* DnBack ;         % Employment drift matrix. Note that differentiation is w.r.t. log n
B                           = (rho + d + 1/Delta) * speye(Nz*Nn) - Nv - Zv ;    % Key matrix in scheme

%%  UPDATE SURPLUS FUNCTION - IMPLICIT SCHEME
%__________________________________________________________________________

y_zn    = reshape ( y_znmat , Nz*Nn , 1 );                      % Vector form
n_zn    = reshape ( n_znmat , Nz*Nn , 1 );                      % Vector form
pi_zn   = y_zn - b*exp(n_zn) - c_f;                             % Vector form
S_zn    = B \ ( pi_zn + 1/Delta * S_zn ) ;                      % Update Surplus Function
Sn_zn   = 0.5 * ( DnBack * S_zn + DnForw * S_zn ) ;             % Update marginal surlplus: average of backward and forward differentiation w.r.t. log n

%%  CONSIDER LAYOFF AND EXIT
%__________________________________________________________________________

    % Reshape as matrix
S_znmat             = reshape( S_zn , Nz , Nn );
Sn_znmat            = reshape( Sn_zn , Nz , Nn );

    % Adjust S for endogenous separation
[ Smax_z , im_z ]   = max( S_znmat , [] , 2 ) ; % For each z, get n that maximizes surplus in grid
icap                = repmat( 1:Nn , Nz , 1 ) >=  repmat( im_z , 1 , Nn ) ; % Index regions where lower n would increase surplus
cap                 = repmat( Smax_z , 1 , Nn ) ; % Maximum S for each z
S_znmat(icap)       = cap(icap) ; % If lower n would increase surplus, adjust S to be the cap

    % Adjust S for exit
SUntrunc_znmat      = S_znmat ; % saved to calculate exit region later
S_znmat             = real( max(S_znmat,0) ) ; % 0 surplus for exiting firms

    % Adjust Sn for endogenous separation
SnUntrunc_znmat     = Sn_znmat ; % saved to calculate separation region later
Sn_znmat            = real( max(Sn_znmat,0) );

    % Adjust Sn for exit
Sn_znmat(SUntrunc_znmat < 0 ) = 0 ;

%%  ASSIGN NEW ERROR
%__________________________________________________________________________

    % Mean absolute difference between old and updated value, over (absolute) average
error_level = max(max(abs(S_znmat -S0_znmat) / mean(abs(S0_znmat(:))) )) ;
error_deriv = max(max(abs(Sn_znmat-Sn0_znmat) / mean(abs(Sn0_znmat(:))) )) ;
error       = max( error_level , error_deriv ) ;

end










