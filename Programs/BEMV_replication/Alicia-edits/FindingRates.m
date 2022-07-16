%% ========================================================================
% 
%    UPDATES FINDING RATES GIVEN SURPLUS FUNCTIONS AND DISTRIBUTIONS
% 
% =========================================================================

function [ error, qnew, phinew, pnew, u ] = ...
    FindingRates( q, phi, p, g_znmat, v_znmat, exitrate, Gn_znmat, SUntrunc_znmat, NumGrids, Params, ExogParams )

    % Load numerical objects and parameters
Nz          = NumGrids.Nz ;         % Number of points in z grid
Nn          = NumGrids.Nn ;         % Number of points in n grid
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid, in logs
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn in a matrix, where d is log difference)
in0         = NumGrids.in0;         % Index of n gridpoint closest to n_0
pi0_z       = NumGrids.pi0_z;       % Distribution of entering firms in z grid
xi          = Params.xi;            % Relative search efficiency of the employed
A           = Params.A;             % Matching function efficiency
MF          = ExogParams.MF;        % Fixed mass of firms
beta        = ExogParams.beta;      % Matching function exponent on vacancies

    % Update units of search efficiency and share phi
n_agg   = MF * sum(sum( g_znmat .* exp(n_znmat) .* dzdn_znmat )) ; % Aggregate employment, recall n is in logs
u       = min( max( 1-n_agg , 0 ), 1 ) ; % Unemployment (cap for algorithm stability)
Search  = abs( u + xi*(1-u) ) ; % Units of search efficiency in the labor market
phinew  = max( abs( u / Search ) , 0.00001 ) ; % Given new match, probability worker is unemployed (with cap for stability)

    % Distribution of entrants
pi0_znmat                   = zeros( Nz , Nn ) ; % Distribution of entering firms
pi0_znmat(:,in0)            = pi0_z ;
pi0_znmat(SUntrunc_znmat<=0)= 0 ; % Truncate by negative surplus
pi0Trunc_znmat              = pi0_znmat / sum(sum(pi0_znmat .* dzdn_znmat)) ; % Ensure integration to 1 

    % Aggregate vacancies
HirePerVacancy_znmat= q * ( phi + (1-phi) * Gn_znmat ) ; % Hiring rate
VacEnt_znmat        = exp(n_znmat) ./ HirePerVacancy_znmat .* pi0Trunc_znmat .* dzdn_znmat ; % Vacancies posted by entrants
V_agg               = MF * ( sum(sum( v_znmat .* g_znmat .* dzdn_znmat)) + exitrate * sum(sum(VacEnt_znmat)) ) ; % Aggregate vacancies, recall exitrate=entryrate

    % Update finding rates
pnew = abs( A * Search^(1-beta) * V_agg^beta / Search ) ; % Rate at which unemployed find vacancies
qnew = abs( A * Search^(1-beta) * V_agg^beta / V_agg ) ; % Rate at which vacancies find workers

    % Assign errors
error_p     = abs( p   - pnew );
error_q     = abs( q   - qnew );
error_phi   = abs( phi - phinew );
error       = max([ error_p error_q error_phi ]) ;

end