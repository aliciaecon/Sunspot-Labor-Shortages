%% ========================================================================
% 
%            COMPUTES RELEVANT MOMENTS BASED ON THE SOLUTION
% 
% =========================================================================

function [ Moments, g_znamat, gdzdn_znmat, gndzdn_znmat, L, T1, T12, T1_12, sep_zn, pi0Trunc_zn,...
        TotalNetPoach_zn, gdzdn_znamat, gndzdn_znamat, EntryRateUW ] ...
    = ComputeMoments( SUntrunc_znmat, SnUntrunc_znmat, ~, ~, v_znmat, g_znmat,...
    Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, Params, ExogParams, ~, options )
%%  NUMERICAL OBJECTS
%__________________________________________________________________________

Numerical = NumericalApproximations( Params , ExogParams, options ); % get settings for numerical approximations
NumGrids  = Grids( Numerical , Params ); % define grids

%%  LOAD SPECIFIC OBJECTS
%__________________________________________________________________________

MF          = ExogParams.MF;        % Fixed mass of firms
d           = ExogParams.d;         % Exogenous firm exit rate
cbar        = ExogParams.cbar;      % Scalar in vacancy posting cost
gamma       = ExogParams.gamma;     % Curvature of vacancy post cost
delta       = Params.delta;         % Exogenous separation rate
xi          = Params.xi;            % Relative search efficiency of the employed
alpha       = Params.alpha ;        % Labor exponent in production function
b           = Params.b;             % Unemployment benefit
Nz          = NumGrids.Nz;          % Number of points in z grid
Nn          = NumGrids.Nn;          % Number of points in n grid
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid
z_znmat     = NumGrids.z_znmat;     % Replica of the z grid
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn in a matrix of size Nz * Nn)
ExitDrift   = Numerical.ExitDrift; 	% Drift at exit boundary
LayoffDrift = Numerical.LayoffDrift;% Drift at layoff boundary

%%  SETUP FOR MONTHLY RATES
%__________________________________________________________________________

    % Miscellaneous
v_znmat     = real(v_znmat) ;                           % Throw out small imaginary numbers
v_zn        = reshape( v_znmat , Nz*Nn , 1 );           % Vector form
gdzdn_znmat = g_znmat .* dzdn_znmat ;                   % Density of firms in grid space
gdzdn_zn    = reshape( gdzdn_znmat , Nz*Nn , 1 );       % Vector form
gndzdn_znmat= g_znmat .* exp( n_znmat ) .* dzdn_znmat ; % Employment-weighted density of firms in grid space
gndzdn_zn   = reshape( gndzdn_znmat , Nz*Nn , 1 );      % Vector form
n_zn        = reshape( n_znmat , Nz*Nn , 1 );           % Vector form
z_zn        = reshape( z_znmat , Nz*Nn , 1 );           % Vector form
SizeNewFirm = repmat( exp( n_zn )  , 1     , Nz*Nn ) ;  % Useful to calculate change in size
SizeOldFirm = repmat( exp( n_zn' ) , Nz*Nn , 1     ) ;  % Useful to calculate change in size
Gv_zn       = reshape( Gv_znmat , Nz*Nn , 1 );          % Vector form

    % Entrants
Gn_zn               = reshape( Gn_znmat , Nz*Nn , 1 );  % Vector form
HirePerVacancy_zn   = q * ( phi + (1-phi)*Gn_zn ) ;     % Hiring rate
VacancyPerHire_zn   = 1 ./ HirePerVacancy_zn ;
vE_zn               = exp(n_zn) .* VacancyPerHire_zn ;  % Vacancies posted by (z,n) entrants
dzdn_zn             = reshape( dzdn_znmat , Nz*Nn , 1 );% Vector form
pi0dzdn_zn          = pi0Trunc_zn .* dzdn_zn ;          % Density of entrants in grid space

    % Vacancies actually posted (resulting in hires)
dlogndt_znmat       = q * v_znmat ./ exp(n_znmat) .* ( phi + (1-phi)*Gn_znmat ) - delta - xi * p * (1-Gv_znmat) ; % dn/dt/n
vHires_znmat        = v_znmat ; % Vacancies actually posted must be slightly modified from policy function...
top                 = dlogndt_znmat > 0 & repmat(1:Nn,Nz,1) == Nn ; % ...if n is increasing and firm is at maximum n in grid...
vHires_znmat(top)   = ( delta + xi * p * (1-Gv_znmat(top)) ) .* exp(n_znmat(top)) ...
    / q ./ ( phi + (1-phi)*Gn_znmat(top) ) ; % ...then vacancies posted are such that n is constant
bottom              = dlogndt_znmat < 0 & repmat(1:Nn,Nz,1) == 1 ; % ...if n is decreasing and firm is at minimum n in grid...
vHires_znmat(bottom)= ( delta + xi * p * (1-Gv_znmat(bottom)) ) .* exp(n_znmat(bottom)) ...
    / q ./ ( phi + (1-phi)*Gn_znmat(bottom) ) ; % ...then vacancies posted are such that n is constant
vHires_zn           = reshape( vHires_znmat , Nz*Nn , 1 );% Vector form

    % Separations
SUntrunc_zn = reshape( SUntrunc_znmat , Nz*Nn , 1 );% Vector form
sep_zn         = ( delta + d ) * ones( Nz*Nn , 1 ) ; % Exogenous job and firm destruction
sep_zn( ( SnUntrunc_znmat < 0 ) & ( repmat(1:Nn,Nz,1) >= 2 ) ) = LayoffDrift ; % Endogenous separations
% Note: When size is at its lowest value, there can be no endogenous 
% separations since we are imposing a reflecting boundary or absorbing 
% state at the lowest size gridpoint
sep_zn = sep_zn + ExitDrift * ( SUntrunc_zn < 0 ) ; 
% Note: the ordering for the definition of the separation rate is consistent with the construction of L in Distribution.m    

%%  MONTHLY EXACT RATES
%__________________________________________________________________________

% (1) Moment 1, Unweighted Exit Rate: whatever mass is not flown into a new
% state must be exiting
ExitRateUW = sum(-L) * gdzdn_zn ;

% (2) Moment 2, Unweighted Entry Rate: entry = exit in SS
EntryRateUW = ExitRateUW ;

% (3) Moment 3, Weighted Exit Rate: weighted exit rate is mass that is not 
% flown into a new state weighted by employment in the initial state
ExitRateW = MF * sum(-L) * gndzdn_zn / (1-u) ;

% (4) Moment 4, Weighted Entry Rate: weighted entry rate is the unweighted 
% entry rate times the pre-assigned size of entrant firms
EntryRateW = EntryRateUW * MF * sum( exp(n_zn) .* pi0Trunc_zn .* dzdn_zn ) / (1-u) ;

% (5) Moment 5, Job Creation of Incumbent firms: average change in 
% employment of expanding firms, as a fraction of total employment
JobCreationIncumbents = MF ...
    * sum( ( L .* max( SizeNewFirm - SizeOldFirm , 0 ) ) * gdzdn_zn ) / (1-u) ;

% (6) Moment 6, Job Destruction of Incumbent firms: average change in 
% employment of contracting firms
JobDestructionIncumbents = MF ...
    * sum( ( L .* max( SizeOldFirm - SizeNewFirm , 0 ) ) * gdzdn_zn ) / (1-u) ; 

% (7) Moment 7, Total Job Creation: incumbents + entrants
JobCreation = JobCreationIncumbents + EntryRateW ;

% (8) Moment 8, Total Job Destruction: incumbents + exiting firms
JobDestruction = JobDestructionIncumbents + ExitRateW ;

% (9) Moment 9, Total Job Reallocation
JobReallocation = JobCreation + JobDestruction ;

% (10) Moment 10, Hires from Employment: hires from other firms is the flow 
% of vacancies posted times the rate at which they contact an employed 
% worker times the probability that the worker will accept the offer.
% Consider entrants as well as incumbents.
HiresE = MF * q * (1-phi) * ( ( v_zn .* Gn_zn )' * gdzdn_zn ...
    + EntryRateUW * ( vE_zn .* Gn_zn )' * pi0dzdn_zn ) / (1-u) ;

% (11) Moment 11, Separations to Employment: Probability of finding any
% other firm times probability of finding a firm with higher marginal
% surplus. In steady state, equal to hires from employment. 
SepE_zn = MF * xi * p * ( (1-Gv_zn) .* exp(n_zn) )' * gdzdn_zn / (1-u) ;

% (12) Moment 12, Hires from Unemployment: include "free hires" from 
% unemployment of entrant firms
HiresU = MF * q * phi * ( vHires_zn' * gdzdn_zn + EntryRateUW * vE_zn' * pi0dzdn_zn ) / (1-u) ;

% (10) and (12) Quarterly moments: 3*monthly
Moments.UEoverEquarterly = 3 * HiresU ;
Moments.UEoverUquarterly = 3 * HiresU * (1-u) / u ;
Moments.EEoverEquarterly = 3 * HiresE ;
Moments.EUoverEquarterly = 3 * HiresU ;

% (13) Moment 13, Separations to Unemployment
SepU_zn = MF * ( sep_zn .* exp(n_zn) )' * gdzdn_zn / (1-u) ;
% Discrepancy between SeparationsU and HiresU seems mostly numerical. As we 
% increase the number of gridpoints (particularly Nn), the discrepancy goes 
% to 0. It would simply require a large Nn (e.g. 500) to get something as 
% good as aggregate net poaching = 0
sep_zn = sep_zn + xi * p * (1-Gv_zn) ; % Include separations to employment

% (14) Moment 14, U-E Transition Rate
Moments.UE = HiresU * (1-u) / u ;

% (15) Moment 15, E-E Transition Rate
Moments.EE = HiresE ;

% (16) Moment 16, E-U Transition Rate
Moments.EU = SepU_zn ;

% (17) Moment 17, Unemployment Rate
Moments.URate = u ;

%%  CROSS SECTIONAL MOMENTS
%__________________________________________________________________________

% (18) Moment 18, Mean of log Productivity
MeanLogZ        = z_zn' * gdzdn_zn ;
Moments.MeanTFP = MeanLogZ ;

% (19) Moment 19, Standard Deviation of log Productivity
Moments.StdLogZ = sqrt( z_zn.^2' * gdzdn_zn - MeanLogZ.^2 ) ;

% (20) Moment 20, Mean of log Size
MeanLogSize         = n_zn' * gdzdn_zn ;
Moments.MeanLogSize = MeanLogSize ;

% (21) Moment 21, Standard Deviation of log Productivity
Moments.StdLogSize = sqrt( n_zn.^2' * gdzdn_zn - MeanLogSize.^2 ) ;

% (22) Moment 22, Output
Moments.MeanY = MF * sum( exp(z_zn) .* exp(n_zn).^alpha .* gdzdn_zn ) ;

% (23) Moment 23, Mean Value Added Per Worker and Unemployment Replacement Rate
MeanVAPW                = Moments.MeanY / (1-u) ;
Moments.MeanVAPW        = MeanVAPW ;
Moments.ReplacementRate = b / MeanVAPW ; % flow value of leisure relative to average output

% (24) Moment 24, Mean log Productivity of Entering Firms
Moments.MeanLogZAge0 = sum( z_zn .* pi0Trunc_zn .* dzdn_zn ) ;

% (25) Moment 25, Standard Deviation of log Productivity of Entering Firms
Moments.StdLogZAge0 = sqrt( sum( z_zn.^2 .* pi0Trunc_zn .* dzdn_zn ) - Moments.MeanLogZAge0.^2 ) ;

% (26) Moment 26, Employment Share of Firms with at least 500 workers
Moments.EmpShare500 = MF * gndzdn_zn' * (exp(n_zn)>=500) / (1-u) ;

% (27) Moment 27, Employment Share of Firms with at least 50 workers, but less than 500
Moments.EmpShare50 = MF * gndzdn_zn' * ( (exp(n_zn)<500) .* (exp(n_zn)>=50) ) / (1-u) ;

% (28) Moment 28, Employment Share of Firms with less than 50 workers
Moments.EmpShare0 = MF * gndzdn_zn' * (exp(n_zn)<50) / (1-u) ;

% (29) Moment 29, Share of Firms with at least 500 workers
Moments.FirmShare500 = gdzdn_zn' * (exp(n_zn)>=500) ;

% (30) Moment 30, Share of Firms with at least 50 workers, but less than 500
Moments.FirmShare50 = gdzdn_zn' * ( (exp(n_zn)<500) .* (exp(n_zn)>=50) )  ;

% (31) Moment 31, Share of Firms with less than 50 workers
Moments.FirmShare0 = gdzdn_zn' * (exp(n_zn)<50)  ;

% (32) Moment 32, Mean Growth Rate of Firms
Moments.MeanGrowthRate = sum( L .* (SizeNewFirm - SizeOldFirm) ...
    ./ (0.5*(SizeNewFirm + SizeOldFirm)) * gdzdn_zn ) ;

% (33) Moment 33, Std. Dev. of Growth Rate of Firms
Moments.StdGrowthRate = sqrt( sum( L .* ( (SizeNewFirm - SizeOldFirm) ...
    ./ (0.5*(SizeNewFirm + SizeOldFirm)) ).^2 * gdzdn_zn ) - Moments.MeanGrowthRate^2 ) ;

%%  ANNUAL RATES AND LIFE CYCLE MOMENTS
%__________________________________________________________________________
if options.Hopenhayn == 0
    
        % Transition Matrices, setup for annual rates
    T1      = fastExpm(L) ; % Monthly
    T12     = fastExpm(12*L) ; % Annual
    T1_12   = inv(L) * ( T12 - speye(Nz*Nn) ) ;
    
    % (34) Moment 34 Unweighted Exit Rate: take all firms inear t - their 
    % mass is one. Compute the mass of firms that remain one year later as 
    % the discrete time annual transition matrix times the initial 
    % distribution of firms. Exit is the former minus the latter.
    Moments.BDSExitRateUW = ( 1 - sum(T12) ) * gdzdn_zn  ;
    
    % (35) Moment 35 Unweighted Entry Rate: equal to exit in steady state
    Moments.BDSEntryRateUW = Moments.BDSExitRateUW ;
    
    % (36) Moment 36 Weighted Exit Rate: the fraction of firms that exit
    % given each starting point weighted by their initial employment
    Moments.BDSExitRateW = MF * ( 1 - sum(T12) ) * gndzdn_zn / (1-u) ;
    
    % (37) Moment 37 Weighted Entry Rate: the unweighted monthly entry rate
    % times the annual entry matrix
    Moments.BDSEntryRateW = MF / (1-u) ...
        * exp(n_zn)' * ( T1_12 * ( pi0Trunc_zn .* dzdn_zn * EntryRateUW ) ) ;
    
    % (38) Moment 38 Mean Job Creation of Incumbent Firms: average change in 
    % employment of expanding firms
    Moments.BDSJobCreationIncumbents = MF ...
        * sum( T12 .* ( max(SizeNewFirm - SizeOldFirm , 0) ) * gdzdn_zn ) / (1-u) ;
    
    % (39) Moment 39 Mean Job Destruction of Incumbent Firms: average change in 
    % employment of contracting firms
    Moments.BDSJobDestructionIncumbents = MF ...
        * sum( T12 .* ( max(SizeOldFirm - SizeNewFirm , 0) ) * gdzdn_zn ) / (1-u) ; 
    
    % (40) Moment 40, Total Job Creation: incumbents + entrants
    Moments.BDSJobCreation = Moments.BDSJobCreationIncumbents + Moments.BDSEntryRateW ;
    
    % (41) Moment 41, Total Job Destruction: incumbents + exiting firms
    Moments.BDSJobDestruction = Moments.BDSJobDestructionIncumbents + Moments.BDSExitRateW ;
    
    % (42) Moment 42, Total Job Reallocation
    Moments.BDSJobReallocation = Moments.BDSJobCreation + Moments.BDSJobDestruction ;
    
    % (43) Moment 43, Mean Growth Rate of Firms: the average
    % change of all firms that remain, rescaled to adjust for the fact that some
    % firms exit
    Moments.BDSMeanGrowthRate = sum( T12 .* ( (SizeNewFirm - SizeOldFirm) ...
        ./ (0.5*(SizeNewFirm + SizeOldFirm)) ) * gdzdn_zn ) / sum( T12 * gdzdn_zn ) ;
    
    % (44) Moment 44, Std. Dev. of Growth Rate of Firms:
    % the squared change of all firms that remain rescaled to adjust for the
    % fact that some firms exit
    Moments.BDSStdGrowthRate = sqrt( sum( T12 .* ( (SizeNewFirm - SizeOldFirm) ...
        ./ (0.5*(SizeNewFirm + SizeOldFirm)) ).^2 * gdzdn_zn ) / sum( T12 * gdzdn_zn ) - Moments.MeanGrowthRate^2 ) ;

        % Choices for Life Cycle Moments
    Na              = 22 ; % Number of points for age grid
    AgeGridLevels   = [ 0:6 , 11 , 16 , 21 ] + 1 ; % Compute moments on these ages
    SizeGridLevels  = [ 0 , 20 , 50 , 250 , 500 ] ; % Compute moments on these sizes
    
        % Transition matrices with age
    AA = spdiags( [ ones(Na*Nz*Nn,1)/(12) , -ones(Na*Nz*Nn,1)/(12) ] , [-Nz*Nn,0] , Na*Nz*Nn , Na*Nz*Nn ) ; % Matrix for evolution of Age
    AA( (Na-1)*Nz*Nn+1:end , (Na-1)*Nz*Nn+1:end ) = 0 ; % Oldest group stays
    LL = kron( speye(Na) , L ) ; % Expand size of L to consider age
    
        % Distribution with age
    g0_zn           = full( T1_12 * pi0Trunc_zn ) ; % Initial distribution
    g_zna           = - ( LL + AA ) \ ( EntryRateUW * [ g0_zn ; zeros((Na-1)*Nz*Nn,1) ] ) ; % Compute distribution with age
    g_znamat        = reshape( g_zna , Nz*Nn , Na ) ; % Matrix form for distribution
    g_znamat        = max( g_znamat , 0 ) ; % Sanity check
    gdzdn_znamat    = g_znamat .* repmat( dzdn_zn , 1 , Na ) ; % Distribution times measure
    g_znamat        = g_znamat / sum( gdzdn_znamat ,'all') ;
    gdzdn_znamat    = gdzdn_znamat / sum( gdzdn_znamat ,'all') ; % Normalize so density sums to one
    gndzdn_znamat   = g_znamat .* repmat( exp(n_zn).*dzdn_zn , 1 , Na ) ; % Employment weighted distribution
    gndzdn_znamat   = gndzdn_znamat / sum(gndzdn_znamat ,'all') ; % Normalize so density sums to one
    
        % Calculate Moments by Age and Size
    [ ~,~,~,~,~,~,~,~, JCRateAge, JDRateAge, ~,~, ExitRateAge, ExitRateUWAge, FirmShareAge, EmpShareAge ] = ...
        ComputeMomentsLifeCycleSize( NumGrids, T12, q, v_znmat, phi, Gn_znmat, gdzdn_znmat,...
        gndzdn_znmat, gdzdn_znamat, gndzdn_znamat, Na, AgeGridLevels, SizeGridLevels ) ;

        % Convert from quarterly to annual rates
    JCRateAge       = 4*JCRateAge ;
    JDRateAge       = 4*JDRateAge ;
    ExitRateAge     = 4*ExitRateAge ;
    ExitRateUWAge   = 4*ExitRateUWAge ;
    
    % (45)-(48) Save Moments 45-48
    for str = { 'ExitRateUW' , 'ExitRate' , 'JCRate' , 'JDRate' , 'EmpShare' , 'FirmShare' } 
        j = 1 ;
        for i = [ 0:6 , 11 , 16 , 21 ]
            eval([ 'Moments.BDS' num2str(cell2mat(str)) 'Age' num2str(i) ' = ' num2str(cell2mat(str)) 'Age(' num2str(j) ') ;' ])
            j = j+1 ;
        end
    end

    % (49) Moment 49, Average log productivity of firms less than one years
    % old, relative to incumbents
    meanzinc = z_zn' * gdzdn_zn ; % Mean log productivity of incumbents
    Moments.BDSMeanLogZYoungRelToInc = z_zn' * ( g_znamat(:,1) .* dzdn_zn ) ...
        / ( g_znamat(:,1)' * dzdn_zn ) - meanzinc ;

    % (50) Moment 50, Average log productivity of firms less than one years
    % old, relative to incumbents
    Moments.BDSMeanLogZOldRelToInc = z_zn' * ( sum(g_znamat(:,11:end),2) .* dzdn_zn ) ...
        / ( sum(g_znamat(:,11:end),2)' * dzdn_zn ) - meanzinc ;

    % (51) Moment 51, Standard Deviation of log productivity of firms 1 to 5 years old
    mean = z_zn' * ( sum(g_znamat(:,1:6),2) .* dzdn_zn ) ...
        / ( sum(g_znamat(:,1:6),2)' * dzdn_zn ) ;
    Moments.BDSStdLogZAge15 = sqrt( sum( z_zn.^2 .* ( sum(g_znamat(:,1:6),2) .* dzdn_zn ) ...
        / ( sum(g_znamat(:,1:6),2)' * dzdn_zn ) ) - mean.^2 ) ;
    
    % (52) Moment 52, Standard Deviation of log productivity of firms over 5 years old
    mean = z_zn' * ( sum(g_znamat(:,7:end),2) .* dzdn_zn ) ...
        / ( sum(g_znamat(:,7:end),2)' * dzdn_zn ) ;
    Moments.BDSStdLogZAge5p = sqrt( sum( z_zn.^2 .* ( sum(g_znamat(:,7:end),2) .* dzdn_zn ) ...
        / ( sum(g_znamat(:,7:end),2)' * dzdn_zn ) ) - mean.^2 ) ;
end

%%  VACANCY MOMENTS
%__________________________________________________________________________
        
% (53) Moment 53, Total vacancy costs paid in the economy as a fraction of output
Moments.VacCostToY = MF * ( cbar * ( v_zn ./ exp(n_zn) ).^(1+gamma) / (1+gamma) .* exp(n_zn) )'...
    * gdzdn_zn / Moments.MeanY ;

% (54) Moment 54, Vacancy Cost as a fraction of hires: must multiply by
% (1-u) because HiresE is calculated as fraction of (1-u)
Moments.VacCostPerHires = Moments.VacCostToY / ( (HiresE+HiresU) * (1-u) ) ;

% (55) Moment 55, Vacancy Yield: total hires over total vacancies
Moments.VacYield = ( (HiresE+HiresU) * (1-u) ) / ( MF * v_zn' * gdzdn_zn ) ;

% (56) Moment 56, Vacancy Rate: total vacancies over total employment
Moments.VacRate = MF * v_zn' * gdzdn_zn / (1-u) ;

    % Vacancy yield at the firm level
Wyield_zn   = q * ( phi + (1-phi) * Gn_zn ) .* gdzdn_zn ; % weighted
yield_zn    = q * ( phi + (1-phi) * Gn_zn ) ; % unweighted

    % Log employment growth
ndotn_znmat = q * v_znmat ./ exp(n_znmat) .* ( phi + (1-phi)*Gn_znmat ) - delta - xi * p * (1-Gv_znmat) ;
ndotn_zn    = reshape( ndotn_znmat , Nz*Nn , 1 ) ;

% (57)-(69) Moments 57-69, Mean Vacancy Yields by employment growth
Moments.VacYield0       = sum( Wyield_zn( abs( ndotn_zn ) < .01            ) ) ./ sum( gdzdn_znmat( abs( ndotn_zn ) < .01            ) ) ;
Moments.VacYield0_5     = sum( Wyield_zn( ndotn_zn >= .01 & ndotn_zn < .05 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .01 & ndotn_zn < .05 ) ) ;
Moments.VacYield5_10    = sum( Wyield_zn( ndotn_zn >= .05 & ndotn_zn < .10 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .05 & ndotn_zn < .10 ) ) ;
Moments.VacYield10_15   = sum( Wyield_zn( ndotn_zn >= .10 & ndotn_zn < .15 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .10 & ndotn_zn < .15 ) ) ;
Moments.VacYield15_20   = sum( Wyield_zn( ndotn_zn >= .15 & ndotn_zn < .20 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .15 & ndotn_zn < .20 ) ) ;
Moments.VacYield20_25   = sum( Wyield_zn( ndotn_zn >= .20 & ndotn_zn < .25 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .20 & ndotn_zn < .25 ) ) ;
Moments.VacYield25_30   = sum( Wyield_zn( ndotn_zn >= .25 & ndotn_zn < .30 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .25 & ndotn_zn < .30 ) ) ;
ndotn_zn                = -ndotn_zn ;
Moments.VacYieldm0_5    = sum( Wyield_zn( ndotn_zn >= .01 & ndotn_zn < .05 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .01 & ndotn_zn < .05 ) ) ;
Moments.VacYieldm5_10   = sum( Wyield_zn( ndotn_zn >= .05 & ndotn_zn < .10 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .05 & ndotn_zn < .10 ) ) ;
Moments.VacYieldm10_15  = sum( Wyield_zn( ndotn_zn >= .10 & ndotn_zn < .15 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .10 & ndotn_zn < .15 ) ) ;
Moments.VacYieldm15_20  = sum( Wyield_zn( ndotn_zn >= .15 & ndotn_zn < .20 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .15 & ndotn_zn < .20 ) ) ;
Moments.VacYieldm20_25  = sum( Wyield_zn( ndotn_zn >= .20 & ndotn_zn < .25 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .20 & ndotn_zn < .25 ) ) ;
Moments.VacYieldm25_30  = sum( Wyield_zn( ndotn_zn >= .25 & ndotn_zn < .30 ) ) ./ sum( gdzdn_znmat( ndotn_zn >= .25 & ndotn_zn < .30 ) ) ;

% (70-75) Moments 70-75, Vacancy Yields by size
Moments.VacYield0_10        = Hmean( yield_zn , ( exp(n_zn) <  10 )                          , gdzdn_znmat ) ;
Moments.VacYield10_50       = Hmean( yield_zn , ( exp(n_zn) >= 10 )   & ( exp(n_zn) < 50 )   , gdzdn_znmat ) ;
Moments.VacYield50_250      = Hmean( yield_zn , ( exp(n_zn) >= 50 )   & ( exp(n_zn) < 250 )  , gdzdn_znmat ) ;
Moments.VacYield250_1000    = Hmean( yield_zn , ( exp(n_zn) >= 250 )  & ( exp(n_zn) < 1000 ) , gdzdn_znmat ) ;
Moments.VacYield1000_5000   = Hmean( yield_zn , ( exp(n_zn) >= 1000 ) & ( exp(n_zn) < 5000 ) , gdzdn_znmat ) ;
Moments.VacYield5000_Inf    = Hmean( yield_zn , ( exp(n_zn) >= 5000 )                        , gdzdn_znmat ) ;

% (76-81) Moments 76-81, Vacancy Rates by size
VacRate_zn                  = MF * v_zn ./ exp(n_zn) ;
Moments.VacRate0_10         = Hmean( VacRate_zn , ( exp(n_zn) < 10 )                           , gdzdn_znmat ) ;
Moments.VacRate10_50        = Hmean( VacRate_zn , ( exp(n_zn) >= 10 )   & ( exp(n_zn) < 50 )   , gdzdn_znmat ) ;
Moments.VacRate50_250       = Hmean( VacRate_zn , ( exp(n_zn) >= 50 )   & ( exp(n_zn) < 250 )  , gdzdn_znmat ) ;
Moments.VacRate250_1000     = Hmean( VacRate_zn , ( exp(n_zn) >= 250 )  & ( exp(n_zn) < 1000 ) , gdzdn_znmat ) ;
Moments.VacRate1000_5000    = Hmean( VacRate_zn , ( exp(n_zn) >= 1000 ) & ( exp(n_zn) < 5000 ) , gdzdn_znmat ) ;
Moments.VacRate5000_Inf     = Hmean( VacRate_zn , ( exp(n_zn) >= 5000 )                        , gdzdn_znmat ) ;

% (82-87) Moments 82-87, Vacancy Shares by size
VacShare_zn                 = v_zn / ( v_zn' * gdzdn_zn ) ;
Moments.VacShare0_10        = Hmean( VacShare_zn , ( exp(n_zn) < 10 )                           , gdzdn_znmat ) ;
Moments.VacShare10_50       = Hmean( VacShare_zn , ( exp(n_zn) >= 10 )   & ( exp(n_zn) < 50 )   , gdzdn_znmat ) ;
Moments.VacShare50_250      = Hmean( VacShare_zn , ( exp(n_zn) >= 50 )   & ( exp(n_zn) < 250 )  , gdzdn_znmat ) ;
Moments.VacShare250_1000    = Hmean( VacShare_zn , ( exp(n_zn) >= 250 )  & ( exp(n_zn) < 1000 ) , gdzdn_znmat ) ;
Moments.VacShare1000_5000   = Hmean( VacShare_zn , ( exp(n_zn) >= 1000 ) & ( exp(n_zn) < 5000 ) , gdzdn_znmat ) ;
Moments.VacShare5000_Inf    = Hmean( VacShare_zn , ( exp(n_zn) >= 5000 )                        , gdzdn_znmat ) ;

%%  SIZE AND PRODUCTIVITY MOMENTS
%__________________________________________________________________________

% (88) Moment 88, Std. Dev. of Value Added Per Worker, firm-weighted
vapw_zn                 = z_zn - (1-alpha) * n_zn ; % in logs
Moments.StdLogVAPWfirm  = sqrt( vapw_zn.^2' * gdzdn_zn - ( vapw_zn' * gdzdn_zn )^2 ) ;

% (89) Moment 89, Std. Dev. of Value Added Per Worker, employment-weighted
Moments.StdLogVAPWemp = sqrt( vapw_zn.^2' * gndzdn_zn / sum(gndzdn_zn) ...
    - ( vapw_zn' * gndzdn_zn )^2  / sum(gndzdn_zn)^2 ) ;

% (90) Moment 90, Std. Dev. of Positive Growth Rate
if options.Hopenhayn == 0
    mean = sum( ( T12 .* (SizeNewFirm-SizeOldFirm>0) ...
        .* ( (SizeNewFirm-SizeOldFirm) ./ (0.5*(SizeNewFirm+SizeOldFirm)) ) ) * gdzdn_zn ) ...
        / sum( ( T12 .* (SizeNewFirm-SizeOldFirm>0) ) * gdzdn_zn ) ;
    Moments.BDSStdGrowthRatePos = sqrt( sum( ( T12 .* (SizeNewFirm-SizeOldFirm>0) ...
        .* ( (SizeNewFirm-SizeOldFirm) ./ (0.5*(SizeNewFirm + SizeOldFirm)) ).^2 ) * gdzdn_zn ) ...
        / sum( ( T12 .* (SizeNewFirm-SizeOldFirm>0) ) * gdzdn_zn ) - mean^2 ) ;
end

% (91) Moment 91, Correlation between Size and Productivity
meanz                   = z_zn' * gdzdn_zn ; 
Moments.CorrLogZLogSize = ( ( z_zn' - meanz ) .* ( n_zn' - MeanLogSize ) ) * gdzdn_zn / ( Moments.StdLogSize * Moments.StdLogZ ) ;

    % Turn gndzdn into density and recalculate moments
gnden_zn    = gndzdn_zn / sum(gndzdn_zn) ;
MeanLogZ  	= z_zn' * gnden_zn ;
StdLogZ  	= sqrt( z_zn.^2' * gnden_zn - MeanLogZ.^2 ) ;   % (19) Moment 19
MeanLogSize = n_zn' * gnden_zn ;                            % (20) Moment 20
StdLogSize  = sqrt( n_zn.^2' * gnden_zn - MeanLogSize.^2 ); % (21) Moment 21

% (92) Moment 92, Density-Weighted Correlation between Size and Productivity
Moments.CorrLogZLogSizeW = ( ( z_zn' - MeanLogZ ) .* ( n_zn' - MeanLogSize ) ) * gnden_zn / ( StdLogSize * StdLogZ ) ;

%%  NET POACHING MOMENTS
%__________________________________________________________________________
if options.Hopenhayn == 0

        % Hires and Separations
    SnUntrunc_zn    = reshape( SnUntrunc_znmat , Nz*Nn , 1 ); % Vector form
    SUntrunc_zn     = reshape( SUntrunc_znmat , Nz*Nn , 1 ); % Vector form
    HiresE_Inc_zn  	= q * (1-phi) * v_zn  .* Gn_zn ; % Hires from employment by incumbents
    HiresE_Ent_zn   = q * (1-phi) * vE_zn .* Gn_zn ; % Hires from employment by entrants
    SepE_zn         = xi * p * (1-Gv_zn) .* exp(n_zn) ; % Separations to employment
    
        % Separations to unemployment
    SepU_zn = delta                               * exp(n_zn) ... exogenous separations (layoffs)
                + d                               * exp(n_zn) ... exogenous exit
              	+ LayoffDrift * (SnUntrunc_zn<0) .* exp(n_zn) ... endogenous separations (layoffs)
              	+ ExitDrift   * (SUntrunc_zn<0)  .* exp(n_zn) ; % endogenous exit
    
        % Poaching results
    GrowthRate_zn       = ( HiresE_Inc_zn - SepU_zn - SepE_zn ) ./ exp(n_zn) ;
    NetPoachByInc_zn    = 3 * ( HiresE_Inc_zn - SepE_zn ) ; 
    NetPoachByEnt_zn    = 3 *                             EntryRateUW * HiresE_Ent_zn .* pi0dzdn_zn ./ max(g_znamat(:,1),1e-20) ;
    TotalNetPoach_zn    = 3 * ( HiresE_Inc_zn - SepE_zn + EntryRateUW * HiresE_Ent_zn .* pi0dzdn_zn ./ max(gdzdn_zn,1e-20) ) ;
    
    % (93) Moment 93, Aggregate Net Poaching as fraction of employment
    Moments.AggregateNetPoaching = ( TotalNetPoach_zn' * gdzdn_zn ) / ( exp(n_zn)' * gdzdn_zn ) ;

        % Expand net poaching by age
    NetPoaching_znamat          = repmat( reshape(NetPoachByInc_zn,Nz,Nn) , 1 , 1 , Na ) ;
    NetPoaching_znamat(:,:,1)   = NetPoaching_znamat(:,:,1) + reshape( NetPoachByEnt_zn , Nz , Nn ) ; % add poaching by entrants to age 0 group
    n_znamat                    = repmat( exp(n_znmat) , 1 , 1 , Na ) ; % Expand employment by age

        % Bins
    SizeGridLevels  = [0 , 20 , 50 , 500 ] ; % Levels of size over which to get moments
    AgeGridLevels   = [0 ,  2 ,  4 ,   6 , 11 ] ; % Levels of age over which to get moments
    ranks           = 0:1:9 ; % Ranks over which to get moments
    ProdGridRanks   = ranks ;
    SizeGridRanks   = ranks ;
    AgeGridRanks    = ranks ;
    VapwGridRanks   = ranks ;
    GrowthGridRanks = ranks ;
    scalerank       = 10 ;
    
        % Marginal surplus rank
    Snrank_znmat = 100 * Gn_znmat ;
   
    % (X) Moments, Net Poaching aggregated by Age bins in Levels
    NetPoachingAgeLevels = AggregateBin(NetPoaching_znamat,...
                                           repmat(1:Na,Nz*Nn,1),...
                                           n_znamat,...
                                           AgeGridLevels,...
                                           gdzdn_znamat,...
                                           'Levels','mean') ;
    for i=1:length(NetPoachingAgeLevels) % assign in loops
        str = [ 'NetPoachingAgeLevels' num2str(AgeGridLevels(i)) ] ;
        Moments.(str) = NetPoachingAgeLevels(i) ;
    end

    % (X) Moments, Net Poaching aggregated by Age Ranks
    NetPoachingAgeRanks = AggregateBin(NetPoaching_znamat,...
                                          repmat(1:Na,Nz*Nn,1),...
                                          n_znamat,...
                                          AgeGridRanks/scalerank,...
                                          gdzdn_znamat,...
                                          'Ranks','mean') ;
    for i=1:length(NetPoachingAgeRanks)
        str = [ 'NetPoachingAgeRanks' num2str(AgeGridRanks(i)) ] ;
        Moments.(str) = NetPoachingAgeRanks(i) ;
    end           
    
    % (X) Moments, Net Poaching aggregated by Size Levels
    NetPoachingSizeLevels = AggregateBin(TotalNetPoach_zn,...
                                           exp(n_znmat),...
                                           exp(n_znmat),...
                                           SizeGridLevels,...
                                           gdzdn_znmat,...
                                           'Levels','mean','WeightDistribution') ;
    for i=1:length(NetPoachingSizeLevels)
        str = [ 'NetPoachingSizeLevels' num2str(SizeGridLevels(i)) ] ;
        Moments.(str) = NetPoachingSizeLevels(i) ;
    end                                       
    
    % (X) Moments, Net Poaching aggregated by Size Ranks
    NetPoachingSizeRanks = AggregateBin(TotalNetPoach_zn,...
                                          exp(n_znmat),...
                                          exp(n_znmat),...
                                          SizeGridRanks/scalerank,...
                                          gdzdn_znmat,...
                                          'Ranks','mean','WeightDistribution') ;
    for i=1:length(NetPoachingSizeRanks)
        str = [ 'NetPoachingSizeRanks' num2str(SizeGridRanks(i)) ] ;
        Moments.(str) = NetPoachingSizeRanks(i) ;
    end                                    

    % (X) Moments, Net Poaching aggregated by Productivity Ranks
    NetPoachingProdRanks = AggregateBin(TotalNetPoach_zn,...
                                          exp(z_znmat),...
                                          exp(n_znmat),...
                                          ProdGridRanks/scalerank,...
                                          gdzdn_znmat,...
                                          'Ranks','mean','WeightDistribution') ;
    for i=1:length(NetPoachingProdRanks)
        str = [ 'NetPoachingProdRanks' num2str(ProdGridRanks(i)) ] ;
        Moments.(str) = NetPoachingProdRanks(i) ;
    end    

    % (X) Moments, Net Poaching aggregated by Value Added per Worker Ranks
    NetPoachingVapwRanks = AggregateBin(TotalNetPoach_zn,...
                                           vapw_zn,...
                                           exp(n_znmat),...
                                           VapwGridRanks/scalerank,...
                                           gdzdn_znmat,...
                                           'Ranks','mean','WeightDistribution') ;
    for i=1:length(NetPoachingVapwRanks)
        str = [ 'NetPoachingVapwRanks' num2str(VapwGridRanks(i)) ] ;
        Moments.(str) = NetPoachingVapwRanks(i) ;
    end

    % (X) Moments, Net Poaching aggregated by Employment Growth Ranks
    NetPoachingGrowthRanks = AggregateBin(TotalNetPoach_zn,...
                                          GrowthRate_zn,...
                                          exp(n_znmat),...
                                          GrowthGridRanks/scalerank,...
                                          gdzdn_znmat,...
                                          'Ranks','mean','WeightDistribution') ;
    for i=1:length(NetPoachingGrowthRanks)
        str = [ 'NetPoachingGrowthRanks' num2str(GrowthGridRanks(i)) ] ;
        Moments.(str) = NetPoachingGrowthRanks(i) ;
    end

    % (X) Moments, Mean Marginal Surplus aggregated by Size Ranks
    MeanSnSizeRanks = AggregateBin(Snrank_znmat,...
                                        exp(n_znmat),...
                                        ones(Nz,Nn),...
                                        SizeGridRanks/scalerank,...
                                        gnden_zn,...
                                        'Ranks','mean') ;
    for i=1:length(MeanSnSizeRanks)
        str = [ 'MeanSnSizeRanks' num2str(SizeGridRanks(i)) ] ;
        Moments.(str) = MeanSnSizeRanks(i) ;
    end  
    
    % (X) Moments, Mean Marginal Surplus aggregated by Size Levels
    MeanSnSizeLevels = AggregateBin(Snrank_znmat,...
                                        exp(n_znmat),...
                                        ones(Nz,Nn),...
                                        SizeGridLevels,...
                                        gnden_zn,...
                                        'Levels','mean') ;
    for i=1:length(MeanSnSizeLevels)
        str = [ 'MeanSnSizeLevels' num2str(SizeGridLevels(i)) ] ;
        Moments.(str) = MeanSnSizeLevels(i) ;
    end
    
    % (X) Moments, Mean Marginal Surplus aggregated by Age Ranks
    MeanSnAgeRanks = AggregateBin(repmat(Snrank_znmat(:),1,Na),...
                                      repmat(1:Na,Nz*Nn,1),...
                                      ones(Nz*Nn,Na),...
                                      AgeGridRanks/scalerank,...
                                      gdzdn_znamat,...
                                      'Ranks','mean') ;
    for i=1:length(MeanSnAgeRanks)
        str = [ 'MeanSnAgeRanks' num2str(AgeGridRanks(i)) ] ;
        Moments.(str) = MeanSnAgeRanks(i) ;
    end 
    
    % (X) Moments, Mean Marginal Surplus aggregated by Age Levels
    MeanSnAgeLevels = AggregateBin(repmat(Snrank_znmat(:),1,Na),...
                                      repmat(1:Na,Nz*Nn,1),...
                                      ones(Nz*Nn,Na),...
                                      AgeGridLevels,...
                                      gdzdn_znamat,...
                                      'Levels','mean') ;
    for i=1:length(MeanSnAgeLevels)
        str = [ 'MeanSnAgeLevels' num2str(AgeGridLevels(i)) ] ;
        Moments.(str) = MeanSnAgeLevels(i) ;
    end 

    % (X) Moments, Mean Marginal Surplus aggregated by Value Added per Worker Ranks
    MeanSnVapwRanks = AggregateBin(Snrank_znmat,...
                                   vapw_zn,...
                                   ones(Nz,Nn),...
                                   VapwGridRanks/scalerank,...
                                   gnden_zn,...
                                   'Ranks','mean') ;
    for i=1:length(MeanSnVapwRanks)
        str = [ 'MeanSnVapwRanks' num2str(VapwGridRanks(i)) ] ;
        Moments.(str) = MeanSnVapwRanks(i) ;
    end                            

    % (X) Moments, Mean Marginal Surplus aggregated by Growth Ranks
    MeanSnGrowthRanks = AggregateBin(Snrank_znmat,...
                                          GrowthRate_zn,...
                                          ones(Nz,Nn),...
                                          GrowthGridRanks/scalerank,...
                                          gnden_zn,...
                                          'Ranks','mean') ;
    for i=1:length(MeanSnGrowthRanks)
        str = [ 'MeanSnGrowthRanks' num2str(GrowthGridRanks(i)) ] ;
        Moments.(str) = MeanSnGrowthRanks(i) ;
    end
    
else 
    TotalNetPoach_zn = 0 ; 
end

%%  LUTTMER MOMENTS
%__________________________________________________________________________

    % Setup new age grid and distribution
Na      = 22 ; % Points in age grid
AA      = spdiags( [ ones(Na*Nz*Nn,1)/(12*10) , -ones(Na*Nz*Nn,1)/(12*10) ] ...
            , [-Nz*Nn,0] , Na*Nz*Nn , Na*Nz*Nn ) ; % Matrix for evolution of Age
AA( (Na-1)*Nz*Nn+1:end , (Na-1)*Nz*Nn+1:end ) = 0 ; % Oldest group stays
LL      = kron( speye(Na) , L ) ; % Expand size of L to consider age
g_zna   = - ( LL + AA ) \ ( EntryRateUW * [ pi0Trunc_zn ; zeros((Na-1)*Nz*Nn,1) ] ) ; % Compute distribution with age
g_znamat= reshape( g_zna , Nz*Nn , Na ) ; % Matrix form

% (X) Moments, Mean and Median Age by Size
for thres = [ 1000 , 10000 ]
    
        % Get Mean
    ind_zn          = ( exp(n_zn) >= thres ) ; % Index for employment over the threshold
    a_znamat        = repmat( 10*((1:Na)-1) , Nz*Nn , 1 ) ; % Age in a matrix
    a_zna           = reshape( a_znamat , Nz*Nn*Na , 1 ) ; % Age in a vector
    dist_znamat     = repmat( ind_zn , 1 , Na ) .* ( g_znamat .* repmat( dzdn_zn , 1 , Na ) ) ; % Distribution of firms with correct age
    dist_zna        = reshape( dist_znamat , Nz*Nn*Na , 1 ) ; % Distribution in a vector
    str             = [ 'MeanAgeSize' num2str(thres) ] ;
    Moments.(str)   = ( a_zna' * dist_zna ) / sum( dist_zna ) ; % Mean age for size bin
        
        % Get Median
    pdf_zna         = dist_zna / sum(dist_zna) ; % Normalize so density integrates to one
    [ ~ , ind ]     = sort( a_zna ) ; % Sort firms by age
    pdfSort_zna     = pdf_zna(ind) ; % Sorted density
    cdfSort_zna     = cumsum( pdfSort_zna ) ; % Integrate, get CDF
    [ ~ , mm ]      = min(abs( cdfSort_zna - 0.5 )) ; % Get index of point closest to median
    aSort_znamat    = a_znamat(ind) ; % Sort to use median index
    str             = [ 'MedianAgeSize' num2str(thres) ] ;
    Moments.(str)   = aSort_znamat( mm ) ; % Median age for size bin
end

% (X) Moment, Churn
Moments.Churn = 12*2*(Moments.EE + Moments.EU) - Moments.BDSJobReallocation ;

    % Setup for elasticity of employment growth to change in TFP
SizeNewFirm = repmat( n_zn  , 1 , Nz*Nn ) ; % Useful to calculate size change, recall n is in logs
SizeOldFirm = repmat( n_zn' , Nz*Nn , 1 ) ; % Useful to calculate size change, recall n is in logs
TFPNewFirm  = repmat( z_zn  , 1 , Nz*Nn ) ; % Useful to calculate productivity change, recall z is in logs
TFPOldFirm  = repmat( z_zn' , Nz*Nn , 1 ) ; % Useful to calculate productivity change, recall z is in logs
dlogn       = T12 .* ( SizeNewFirm - SizeOldFirm ) ; % Change in log size
dlogz       = T12 .* ( TFPNewFirm  - TFPOldFirm  ) ; % Change in log productivity
W           = repmat( gdzdn_zn' , Nz*Nn , 1 ) ; % Repeated distribution, weights in regression
beta        = lscov( [ ones(Nz^2*Nn^2,1) , dlogz(:) ] , dlogn(:) , W(:) ) ; % Weighted Least Squares of dlogz over dlogn

% (X) Moment, Change in log size over change in log productivity, weighted
% by firm
Moments.dLogSizedLogZFirm = beta(2) ;

% (X) Moment, Change in log size over change in log productivity, weighted
% by employment
W                           = repmat( gnden_zn' , Nz*Nn , 1 ) ;
beta                        = lscov( [ ones(Nz^2*Nn^2,1) , dlogz(:) ] , dlogn(:) , W(:) ) ;
Moments.dLogSizedLogZEmp    = beta(2) ;

end