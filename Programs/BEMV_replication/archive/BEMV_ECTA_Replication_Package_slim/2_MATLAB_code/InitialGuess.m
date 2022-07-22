%% ========================================================================
% 
%  CONSTRUCTS AN INITIAL GUESS FOR THE SURPLUS FUNCTIONS AND DISTRIBUTIONS
% 
% =========================================================================

function [ g_znmat, u, p, q, phi, Gn_znmat, Ghatn_znmat, Gv_znmat, S_znmat, Sn_znmat, v_znmat, exitrate ]...
    = InitialGuess( ExogParams, NumGrids, Params, MatExog, PerPeriod )
%%  LOAD OBJECTS
%__________________________________________________________________________

z           = NumGrids.z;           % z grid in log space
n           = NumGrids.n;           % n grid in log space
Nz          = NumGrids.Nz ;         % number of points in z grid
Nn          = NumGrids.Nn ;         % number of points in n grid
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid
dn_znmat    = NumGrids.dn_znmat;    % Replica of the n step size
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn in a matrix of size Nz * Nn)
y_znmat     = PerPeriod.y_znmat;    % Flow output
yn_znmat    = PerPeriod.yn_znmat;   % Marginal product
b           = Params.b;             % Unemployment benefit
xi          = Params.xi;            % Relative search efficiency of the employed
rho         = ExogParams.rho;       % Time discount rate
MF          = ExogParams.MF;        % Fixed mass of firms
Zd          = MatExog.Zd ;          % Matrix for distribution updating

%%  UNWEIGHTED DISTRIBUTION g(z,n)
%__________________________________________________________________________

    % Fix one entry
ifix2 = 2;                   % Normalization index
ifixn = ifix2+(0:(Nn-1))*Nz; % Replicate normalization index
g0_zn = zeros(Nz*Nn,1);      % Initialize right-hand-side
Zd2   = Zd;                  % Initialize matrix
for i = 1:Nn % Correct matrix for non-degeneracy
    ifixloc         = ifixn(i);
    g0_zn(ifixloc)  = 0.1;
    Zd2(ifixloc,:)  = [ zeros(1,ifixloc-1) , 1 , zeros(1,Nz*Nn-ifixloc) ];
end
g_zn = Zd2 \ g0_zn; % Compute distribution

    % To construct an initial guess, multiply the productivity distribution
% by a lognormal distribution in firm size: a normal in log size, with 
% E[log n] = 2 and Var[log n] = 1.5 
mm0     = 0.1;
s0      = 1;
g1_znmat= 1 / sqrt(2*pi) / s0 .* exp( -(n_znmat-mm0).^2 / 2 / s0^2 );
g1_zn   = reshape( g1_znmat , Nz*Nn , 1 );
g_zn    = g_zn .* g1_zn;
g_znmat = reshape( g_zn , Nz , Nn );
g_znmat = g_znmat / sum(sum( g_znmat .* dzdn_znmat)); % Normalize to unit mass

%%  SURPLUS FUNCTIONS
%__________________________________________________________________________

    % Initialize surplus S and marginal surplus Sn
S_znmat     = (y_znmat - exp(n_znmat)*b) / rho; % Present value of constant future output minus total benefit of unemployment, exp(n) because n is in logs
[ ~ , iz1 ] = min(abs(z)) ;                     % index for gridpoint closest to z = 1
[ ~ , in1 ] = min(abs(n)) ;                     % index for gridpoint closest to n = 1
FixedCost   = S_znmat(iz1,in1) ;                % Set exit subsidy such that firms exit at (z,n) = (1,1)
S_znmat     = S_znmat - FixedCost ;             % Incorporate fixed cost in surplus calculation
Sn_znmat    = (yn_znmat - b) / rho ;            % Marginal surplus = Present value of constant marginal product minus un. benefit

    % Impose boundary conditions
S_znmat     = max( S_znmat , 0 ) ;              % No firm in exit region
Sn_znmat    = max( Sn_znmat , 0 ) ;             % No firm in layoff region

    % Integrate up marginal to avoid issues with endogenous separations
ISn_znmat   = cumsum( Sn_znmat .* dn_znmat , 2) ;           % Integral of marginal surplus
ISn_znmat   = [ zeros(Nz,1) , ISn_znmat(:,1:end-1) ] ;      % Shift the integral by one
baseline_z  = S_znmat(:,1) ;                                % First value
S_znmat     = repmat( baseline_z , 1 ,  Nn) + ISn_znmat ;   % Re-define as integral of marginal

%%  AGGREGATES AND WEIGHTED DISTRIBUTIONS
%__________________________________________________________________________

    % Initial guess for aggregates
u                   = 0.05 ;                    % Unemployment
Search              = u + xi*(1-u) ;            % Units of search efficiency
v_znmat             = exp(n_znmat) ;            % Vacancies by (z,n), equal to level employment as guess
v_znmat(S_znmat==0) = 0 ;                       % No vacancies posted by exiting firms
V_agg               = MF * sum(sum(v_znmat .* g_znmat .* dzdn_znmat)) + 0.0001 ; % Aggregate vacancies
p                   = 0.2 ;                     % Job finding probability (matches/searchers)
q                   = p / V_agg * Search ;      % Vacancy filling probability (matches/vacancies)
phi                 = max( u/Search , 0.05 ) ;  % Probability a matched vacancy finds an unemployed worker insted of employed

    % Initial guess for distributions
[ Gn_znmat, Ghatn_znmat, inds, inds_inv ] = Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ; % Employment-weighted marginal surplus distribution
exitrate = 0 ;
Gv_znmat = Cdf_Gv( inds, inds_inv, v_znmat, g_znmat, S_znmat, exitrate, q, Gn_znmat, phi, NumGrids ) ; % Vacancy-weighted marginal surplus distribution

end