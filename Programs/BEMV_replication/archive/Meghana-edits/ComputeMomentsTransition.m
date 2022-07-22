%% ========================================================================
% 
%           COMPUTES MOMENTS ALONG A TRANSITION EQUILIBRIUM
% 
% =========================================================================

function TransitionMoments = ComputeMomentsTransition( TransitionNum, TransitionEq,...
    ~, SsEq, NTmin, NTmax, ExogParams, Params, NumGrids, options )
%%  LOAD SPECIFIC OBJECTS
%__________________________________________________________________________

d           = ExogParams.d;         % Exogenous firm exit rate
alpha       = Params.alpha ;        % Labor exponent in production function
delta       = Params.delta;         % Exogenous separation rate
xi          = Params.xi;            % Relative search efficiency of the employed
Nz          = NumGrids.Nz;          % Number of points in z grid
Nn          = NumGrids.Nn;          % Number of points in n grid
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid
z_znmat     = NumGrids.z_znmat;     % Replica of the z grid
n_zn        = reshape(n_znmat,Nz*Nn,1);
z_zn        = reshape(z_znmat,Nz*Nn,1);
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn in a matrix of size Nz * Nn)
in0         = NumGrids.in0;         % Index of n gridpoint closest to n_0
pi0_z       = NumGrids.pi0_z;       % Distribution of entering firms in z grid
Nt          = TransitionNum.Nt;     % Number of periods in transition


%%  SETUP
%__________________________________________________________________________

    % Various useful objects
Numerical   = NumericalApproximations( Params, ExogParams, options );
ExitDrift   = Numerical.ExitDrift; 	% Drift at exit boundary
LayoffDrift = Numerical.LayoffDrift;% Drift at layoff boundary
SizeNewFirm = repmat( exp(n_zn)  , 1                  , length(n_zn) ) ;
SizeOldFirm = repmat( exp(n_zn') , length(n_zn) , 1                  ) ;
size1000    = ( exp(n_zn) >= 1000 ) ;
size500     = ( exp(n_zn) >= 500 ) ;
size50      = ( ( exp(n_zn) < 500 ) .* ( exp(n_zn) >= 50 ) ) ;
size0       = ( exp(n_zn) < 50 ) ;
size00      = ( exp(n_zn) < 10 ) ;
t           = floor( size(SsEq.g,2) / 2 ) ;
g_zn        = SsEq.g(:,t) ;
q           = SsEq.q(t) ;
phi         = SsEq.phi(t) ;
exitrate    = SsEq.ME(t) / SsEq.MF(t) ;
SUntrunc_zn = SsEq.SUntrunc(:,t) ;
Sn_zn       = TransitionEq.Sn(:,end) ;
cutlow      = 0.2 ;
cuthigh     = 0.8 ;

    % Reshape
g_znmat         = reshape(g_zn,Nz,Nn);
Sn_znmat        = reshape(Sn_zn,Nz,Nn);
SUntrunc_znmat  = reshape(SUntrunc_zn,Nz,Nn);

    % Marginal surplus rank
[ Gn_znmat ,~, inds, inds_inv ] = Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ;
SnRankFinal                     = Cdf_Gv( inds, inds_inv, ones(Nz,Nn), g_znmat, SUntrunc_znmat, exitrate, q, Gn_znmat, phi, NumGrids) ;


%%  LOOP OVER TRANSITION PERIODS
%__________________________________________________________________________
for t=NTmin:NTmax
disp(['Computing moments for time t = ' num2str(t) ' out of ' num2str(NTmax)])

    % Construct distribution of entry
SUntrunc_zn                 = TransitionEq.SUntrunc(:,t) ;
SUntrunc_znmat              = reshape(SUntrunc_zn,Nz,Nn);
pi0_znmat                   = zeros(Nz,Nn);
pi0_znmat(:,in0)            = pi0_z;                                    % Entering distribution of productivity
pi0_znmat(SUntrunc_znmat<0) = 0;                                        % Truncate distribution: no firm enters with negative surplus
pi0Trunc_zn                 = reshape( pi0_znmat , Nz*Nn , 1 );         % Vector form
dzdn_zn                     = reshape( dzdn_znmat , Nz*Nn , 1);         % Vector form
pi0Trunc_zn                 = pi0Trunc_zn / ( pi0Trunc_zn' * dzdn_zn ); % ensure that this integrates to one

    % Load equilibrium variables
q               = TransitionEq.q(t) ;
p               = TransitionEq.p(t) ;
phi             = TransitionEq.phi(t) ;
u               = TransitionEq.u(t) ;
MF              = TransitionEq.MF(t) ;
ME              = TransitionEq.ME(t) ;
v_zn            = TransitionEq.v(:,t) ;
Gn_zn           = TransitionEq.Gn(:,t) ;
Gv_zn           = TransitionEq.Gv(:,t) ;
g_zn            = TransitionEq.g(:,t) ;
gdzdn_zn        = g_zn .* dzdn_zn ;
gndzdn_zn       = g_zn .* exp(n_zn) .* dzdn_zn ;
EntrantRetention= q * ( phi + (1-phi) * Gn_zn ) ;
EntrantVacancy 	= exp(n_zn) ./ EntrantRetention .* pi0Trunc_zn .* dzdn_zn ;
gvdzdn_zn       = MF * g_zn .* v_zn .* dzdn_zn + ME .* EntrantVacancy ;
gvdzdnInc_zn    = MF * g_zn .* v_zn .* dzdn_zn ; 
Sn_zn           = TransitionEq.Sn(:,t) ;
sn_zn           = Sn_zn ./ exp(n_zn) ;
SnUntrunc_zn    = TransitionEq.SnUntrunc(:,t) ;

    % Save into moments struct
TransitionMoments.q(t)      = q ;
TransitionMoments.p(t)      = p ;
TransitionMoments.phi(t)    = phi ;
TransitionMoments.u(t)      = u ;
TransitionMoments.MF(t)     = MF ;
TransitionMoments.ME(t)     = ME ;

    % Weighted distributions
[ GSn, ~,~,~ ]  = Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ;
[ Gz,  ~,~,~ ]  = Cdf_Gn( z_znmat .* exp(n_znmat), g_znmat, NumGrids ) ;
[ Gva, ~,~,~ ]  = Cdf_Gn( ( exp( z_znmat + (alpha-1) * n_znmat ) ) .* exp(n_znmat), g_znmat, NumGrids ) ;

    % Shares
LowSn = GSn(:) <  cutlow ; 
MidSn = GSn(:) >= 0.33 & GSn(:) < 0.67 ;
TopSn = GSn(:) >= cuthigh ; 
Lowz  = Gz(:)  <  cutlow ; 
Midz  = Gz(:)  >= 0.33 & Gz(:)  < 0.67 ;
Topz  = Gz(:)  >= cuthigh ; 
Lowva = Gva(:) <  cutlow ; 
Midva = Gva(:) >= 0.33 & Gva(:) < 0.67 ;
Topva = Gva(:) >= cuthigh ; 
LOWva = Gva(:) <  0.6 ;
TOPva = Gva(:) >= 0.6 ;


%%  LABOR FLOWS
%__________________________________________________________________________

if t > 1
    L   = reshape( TransitionEq.L(:,t-1) , Nz*Nn , Nz*Nn ) ;
    um1 = TransitionMoments.u(t-1) ;
else
    L   = reshape( TransitionEq.L(:,Nt-1) , Nz*Nn , Nz*Nn ) ;
    um1 = TransitionEq.u(Nt) ;
end
TransitionMoments.V(t)              = sum( gvdzdn_zn ) ;
TransitionMoments.MF1(t)            = MF * sum ( gdzdn_zn .* ( exp(n_zn) >= 1 ) )  ;
TransitionMoments.HiresE(t)         = MF * q * (1-phi) * ( v_zn .* Gn_zn )' * gdzdn_zn ; 
TransitionMoments.SeparationsE(t)   = MF * xi * p * sum( (1-Gv_zn) .* gndzdn_zn ) ;
TransitionMoments.EntryRateUW(t)    = ME / MF ;
TransitionMoments.EntryRateW(t)     = ME * sum(sum( exp(n_znmat) .* pi0_znmat .* dzdn_znmat )) / (1-u) ;                                   
TransitionMoments.HiresU(t)         = (1-u) * TransitionMoments.EntryRateW(t) + MF * q * phi * v_zn' * gdzdn_zn ; 
TransitionMoments.Hires(t)          = TransitionMoments.HiresU(t) + TransitionMoments.HiresE(t) ;
TransitionMoments.SeparationsU(t)   = ( u - um1 ) / TransitionNum.dT(t) + TransitionMoments.HiresU(t) ;
TransitionMoments.EU(t)             = TransitionMoments.SeparationsU(t) / (1-um1)  ;
TransitionMoments.UE(t)             = TransitionMoments.HiresU(t) / u ;
TransitionMoments.Tightness(t)      = TransitionMoments.V(t) / ( u + xi*(1-u) ) ;
TransitionMoments.EEh(t)            = TransitionMoments.HiresE(t) / (1-u) ;
TransitionMoments.EEs(t)            = TransitionMoments.SeparationsE(t) / (1-u);
TransitionMoments.URate(t)          = u ;

%%  EE SHARE BY DESTINATION
%__________________________________________________________________________

    % Marginal surplus
TransitionMoments.EETopSn(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* TopSn )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EEMidSn(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* MidSn )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EELowSn(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* LowSn )' * gdzdn_zn / ( 1 - u ) ;

    % Value added
TransitionMoments.EETopva(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* Topva )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EEMidva(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* Midva )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EELowva(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* Lowva )' * gdzdn_zn / ( 1 - u ) ;

    % Productivity
TransitionMoments.EETopz(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* Topz )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EEMidz(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* Midz )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EELowz(t) = MF * q * (1-phi) * ( v_zn .* Gn_zn .* Lowz )' * gdzdn_zn / ( 1 - u ) ;

    % Size
TransitionMoments.EE500(t)  = MF * q * (1-phi) * ( v_zn .* Gn_zn .* size500 )' * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EE50(t) 	= MF * q * (1-phi) * ( v_zn .* Gn_zn .* size50 )'  * gdzdn_zn / ( 1 - u ) ;
TransitionMoments.EE0(t) 	= MF * q * (1-phi) * ( v_zn .* Gn_zn .* size0 )'   * gdzdn_zn / ( 1 - u ) ;


%%  CROSS-SECTIONAL FIRM CHARACTERISTICS
%__________________________________________________________________________

TransitionMoments.MeanLogZ(t)               = sum( z_zn .* gdzdn_zn );
TransitionMoments.MeanZ(t)                  = sum( exp(z_zn) .* gdzdn_zn );
TransitionMoments.StdLogZ(t)                = sqrt( sum( (z_zn).^2 .* gdzdn_zn ) - TransitionMoments.MeanLogZ(t).^2 ) ;
TransitionMoments.StdZ(t)                   = sqrt( sum( exp(z_zn).^2 .* gdzdn_zn ) - TransitionMoments.MeanZ(t).^2 ) ;
TransitionMoments.MeanLogSize(t)            = sum( n_zn .* gdzdn_zn ) ;
TransitionMoments.MeanSize(t)               = sum( exp(n_zn) .* gdzdn_zn ) ;
TransitionMoments.StdLogSize(t)             = sqrt( sum( n_zn.^2 .* gdzdn_zn ) - TransitionMoments.MeanLogSize(t).^2 ) ;
TransitionMoments.StdSize(t)                = sqrt( sum( exp(n_zn).^2 .* gdzdn_zn ) - TransitionMoments.MeanSize(t).^2 ) ;
TransitionMoments.MeanY(t)                  = MF * sum( exp(z_zn) .* exp(n_zn).^alpha .* gdzdn_zn ) ;
TransitionMoments.YPW(t)                    = MF * sum( exp(z_zn) .* exp(n_zn).^alpha .* gdzdn_zn ) / (1-u) ;
TransitionMoments.MeanVAPW(t)               = TransitionMoments.MeanY(t) / (1-u) ;
TransitionMoments.CorrLogZLogSize(t)        = ( ( z_zn' - TransitionMoments.MeanLogZ(t) ) .* ( n_zn' - TransitionMoments.MeanLogSize(t) ) ) * gdzdn_zn / ( TransitionMoments.StdLogSize(t) * TransitionMoments.StdLogZ(t) ) ;
TransitionMoments.CorrZSize(t)              = ( ( exp(z_zn') - TransitionMoments.MeanZ(t) ) .* ( exp(n_zn') - TransitionMoments.MeanSize(t) ) ) * gdzdn_zn / ( TransitionMoments.StdSize(t) * TransitionMoments.StdZ(t) ) ;
lvapw                                       = z_zn - (1-alpha) * n_zn ;
vapw                                        = exp(lvapw) ;
keep                                        = Sn_zn > 0 ;
MeanVAPWfirm                                = sum( vapw .* gdzdn_zn ) ;
MeanVAPWemp                                 = sum( vapw .* gndzdn_zn ) / sum( gndzdn_zn ) ;
MeanLogVAPWKeepemp                          = sum( lvapw(keep) .* gndzdn_zn(keep) ) / sum( gndzdn_zn(keep) ) ;
MeanSnfirm                                  = sum( Sn_zn .* gdzdn_zn ) ;
MeanSnemp                                   = sum( Sn_zn .* gndzdn_zn ) / sum( gndzdn_zn ) ;
MeanLogSnKeepemp                            = sum( log(Sn_zn(keep)) .* gndzdn_zn(keep) ) / sum( gndzdn_zn(keep) ) ;
TransitionMoments.StdLogVAPWfirm(t)         = sqrt( sum( lvapw.^2 .* gdzdn_zn ) - sum( lvapw .* gdzdn_zn )^2 ) ;
TransitionMoments.StdLogVAPWemp(t)          = sqrt( sum( lvapw.^2 .* gndzdn_zn ) / sum(gndzdn_zn) - sum( lvapw .* gndzdn_zn )^2  / sum(gndzdn_zn)^2 ) ;
TransitionMoments.StdLogVAPWKeepemp(t)      = sqrt( sum( lvapw(keep).^2 .* gndzdn_zn(keep) ) / sum(gndzdn_zn(keep)) - sum( lvapw(keep) .* gndzdn_zn(keep) )^2  / sum(gndzdn_zn(keep))^2 ) ;
TransitionMoments.StdVAPWfirm(t)            = sqrt( sum( vapw.^2 .* gdzdn_zn ) - sum( vapw .* gdzdn_zn )^2 ) ;
TransitionMoments.StdVAPWemp(t)             = sqrt( sum( vapw.^2 .* gndzdn_zn ) / sum(gndzdn_zn) - sum( vapw .* gndzdn_zn )^2  / sum(gndzdn_zn)^2 ) ;
StdSnfirm                                   = sqrt( sum( Sn_zn.^2 .* gdzdn_zn ) / sum(gdzdn_zn) - sum( Sn_zn .* gdzdn_zn )^2  / sum(gdzdn_zn)^2 ) ;
StdSnemp                                    = sqrt( sum( Sn_zn.^2 .* gndzdn_zn ) / sum(gndzdn_zn) - sum( Sn_zn .* gndzdn_zn )^2  / sum(gndzdn_zn)^2 ) ;
StdLogSnKeepemp                             = sqrt( sum( log(Sn_zn(keep)).^2 .* gndzdn_zn(keep) ) / sum(gndzdn_zn(keep)) - sum( log(Sn_zn(keep)) .* gndzdn_zn(keep) )^2  / sum(gndzdn_zn(keep))^2 ) ;
TransitionMoments.CorrSnVAPWfirm(t)         = ( ( Sn_zn' - MeanSnfirm ) .* ( vapw' - MeanVAPWfirm ) ) * gdzdn_zn / ( TransitionMoments.StdVAPWfirm(t) * StdSnfirm ) ;
TransitionMoments.CorrSnVAPWemp(t)          = ( ( Sn_zn' - MeanSnemp ) .* ( vapw' - MeanVAPWemp ) ) * gndzdn_zn / ( TransitionMoments.StdVAPWemp(t) * StdSnemp ) / sum( gndzdn_zn) ;
TransitionMoments.CorrLogSnLogVAPWemp(t)    = ( ( log(Sn_zn(keep))' - MeanLogSnKeepemp ) .* ( lvapw(keep)' - MeanLogVAPWKeepemp ) ) * gndzdn_zn(keep) / ( TransitionMoments.StdLogVAPWKeepemp(t) * StdLogSnKeepemp ) / sum( gndzdn_zn ) ;


%%  FIRM SHARES
%__________________________________________________________________________

TransitionMoments.FirmShare500(t)   = gdzdn_zn' * ( exp(n_zn) >= 500 )   ;
TransitionMoments.FirmShare50(t)	= gdzdn_zn' * ( ( exp(n_zn) < 500 ) .* ( exp(n_zn) >= 50 ) ) ;
TransitionMoments.FirmShare0(t) 	= gdzdn_zn' * ( exp(n_zn) < 50 )  ;


%%  EMPLOYMENT SHARES
%__________________________________________________________________________

    % By size
TransitionMoments.EmpShare500(t)  	= gndzdn_zn' * size500 / sum( gndzdn_zn ) ;
TransitionMoments.EmpShare50(t)   	= gndzdn_zn' * size50  / sum( gndzdn_zn ) ;
TransitionMoments.EmpShare0(t)      = gndzdn_zn' * size0   / sum( gndzdn_zn ) ;

    % By productivity
TransitionMoments.EmpShareLowz(t)  	= gndzdn_zn' * Lowz / sum( gndzdn_zn ) ;
TransitionMoments.EmpShareMidz(t)   = gndzdn_zn' * Midz / sum( gndzdn_zn ) ;
TransitionMoments.EmpShareTopz(t)   = gndzdn_zn' * Topz / sum( gndzdn_zn ) ;

    % By labor prod
TransitionMoments.EmpShareLowva(t)  = gndzdn_zn' * Lowva / sum( gndzdn_zn ) ;
TransitionMoments.EmpShareMidva(t)  = gndzdn_zn' * Midva / sum( gndzdn_zn ) ;
TransitionMoments.EmpShareTopva(t)  = gndzdn_zn' * Topva / sum( gndzdn_zn ) ;

    % By Sn
TransitionMoments.EmpShareLowSn(t)  = gndzdn_zn' * LowSn / sum( gndzdn_zn ) ;
TransitionMoments.EmpShareMidSn(t)  = gndzdn_zn' * MidSn / sum( gndzdn_zn ) ;
TransitionMoments.EmpShareTopSn(t)  = gndzdn_zn' * TopSn / sum( gndzdn_zn ) ;
                                                      

%%  VACANCY SHARES
%__________________________________________________________________________

    % By size
TransitionMoments.VacShare1000(t)       = gvdzdn_zn' * size1000  / sum( gvdzdn_zn ) ;
TransitionMoments.VacShare500(t)        = gvdzdn_zn' * size500   / sum( gvdzdn_zn ) ;
TransitionMoments.VacShare50(t)         = gvdzdn_zn' * size50    / sum( gvdzdn_zn ) ;
TransitionMoments.VacShare0(t)          = gvdzdn_zn' * size0     / sum( gvdzdn_zn ) ;
TransitionMoments.VacShare00(t)         = gvdzdn_zn' * size00    / sum( gvdzdn_zn ) ;
TransitionMoments.VacTot(t)             = sum( gvdzdn_zn ) ;
TransitionMoments.VacIncTot(t)          = sum( gvdzdnInc_zn ) ;

    % By productivity
TransitionMoments.VacShareLowz(t)       = gvdzdn_zn' * Lowz / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareMidz(t)       = gvdzdn_zn' * Midz / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareTopz(t)       = gvdzdn_zn' * Topz / sum( gvdzdn_zn ) ;

    % By labor prod
TransitionMoments.VacShareLowva(t)      = gvdzdn_zn' * Lowva / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareMidva(t)      = gvdzdn_zn' * Midva / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareTopva(t)      = gvdzdn_zn' * Topva / sum( gvdzdn_zn ) ;

    % By marginal surplus
TransitionMoments.VacShareLowSn(t)      = gvdzdn_zn' * LowSn / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareMidSn(t)      = gvdzdn_zn' * MidSn / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareTopSn(t)  	= gvdzdn_zn' * TopSn / sum( gvdzdn_zn ) ;
TransitionMoments.VacShareEmpLowSn(t) 	= 1 / MF * gvdzdn_zn' * LowSn / ( gndzdn_zn' * LowSn ) ;
TransitionMoments.VacShareEmpTopSn(t)   = 1 / MF * gvdzdn_zn' * TopSn / ( gndzdn_zn' * TopSn ) ;


%%  NET POACHING
%__________________________________________________________________________

    % Previous calculations
hiresU  = q * phi * v_zn ;
hiresE  = q * (1-phi) * v_zn .* Gn_zn ;
quitsE  = xi * p * (1-Gv_zn) .* exp(n_zn) ;
np      = hiresE - quitsE ;

    % By size
TransitionMoments.NetPoach1000(t) 	= sum( np .* gdzdn_zn .* size1000 ) / sum( gndzdn_zn .* size1000 ) ;
TransitionMoments.NetPoach500(t)    = sum( np .* gdzdn_zn .* size500  ) / sum( gndzdn_zn .* size500  ) ;
TransitionMoments.NetPoach50(t)     = sum( np .* gdzdn_zn .* size50   ) / sum( gndzdn_zn .* size50   ) ;
TransitionMoments.NetPoach0(t)      = sum( np .* gdzdn_zn .* size0    ) / sum( gndzdn_zn .* size0    ) ;
TransitionMoments.NetPoach00(t)     = sum( np .* gdzdn_zn .* size00   ) / sum( gndzdn_zn .* size00   ) ;

    % By productivity
TransitionMoments.NetPoachTopz(t)   = sum( np .* gdzdn_zn .* Topz ) / sum( gndzdn_zn .* Topz ) ;
TransitionMoments.NetPoachMidz(t)   = sum( np .* gdzdn_zn .* Midz ) / sum( gndzdn_zn .* Midz ) ;
TransitionMoments.NetPoachLowz(t)   = sum( np .* gdzdn_zn .* Lowz ) / sum( gndzdn_zn .* Lowz ) ;

    % By labor prod
TransitionMoments.NetPoachTopva(t)  = sum( np .* gdzdn_zn .* Topva ) / sum( gndzdn_zn .* Topva ) ;
TransitionMoments.NetPoachMidva(t)  = sum( np .* gdzdn_zn .* Midva ) / sum( gndzdn_zn .* Midva ) ;
TransitionMoments.NetPoachLowva(t)  = sum( np .* gdzdn_zn .* Lowva ) / sum( gndzdn_zn .* Lowva ) ;

    % By marginal surplus
TransitionMoments.NetPoachTopSn(t)  = sum( np .* gdzdn_zn .* TopSn ) / sum( gndzdn_zn .* TopSn ) ;
TransitionMoments.NetPoachMidSn(t)  = sum( np .* gdzdn_zn .* MidSn ) / sum( gndzdn_zn .* MidSn ) ;
TransitionMoments.NetPoachLowSn(t)  = sum( np .* gdzdn_zn .* LowSn ) / sum( gndzdn_zn .* LowSn ) ;

    % Other stuff by Sn
TransitionMoments.MeanGTopSn(t)     = sum( Gn_zn .* gndzdn_zn .* TopSn ) / sum( gndzdn_zn .* TopSn ) ;
TransitionMoments.MeanGMidSn(t)     = sum( Gn_zn .* gndzdn_zn .* MidSn ) / sum( gndzdn_zn .* MidSn ) ;
TransitionMoments.MeanGLowSn(t) 	= sum( Gn_zn .* gndzdn_zn .* LowSn ) / sum( gndzdn_zn .* LowSn ) ;
TransitionMoments.HireRateTopSn(t)  = sum( hiresE .* gdzdn_zn .* TopSn ) / sum( gndzdn_zn .* TopSn ) ;
TransitionMoments.HireRateLowSn(t)  = sum( hiresE .* gdzdn_zn .* LowSn ) / sum( gndzdn_zn .* LowSn ) ;
TransitionMoments.QuitRateTopSn(t)  = sum( quitsE .* gdzdn_zn .* TopSn ) / sum( gndzdn_zn .* TopSn ) ;
TransitionMoments.QuitRateLowSn(t)  = sum( quitsE .* gdzdn_zn .* LowSn ) / sum( gndzdn_zn .* LowSn ) ;


%%  HALTIWANGER MOMENTS
%__________________________________________________________________________

    % Previous calculations
sep     = ( delta + d ) * ones(Nz*Nn,1) ;
sep( ( reshape(SnUntrunc_zn,Nz,Nn) < 0 ) & ( repmat(1:Nn,Nz,1) >= 2 ) ) ...
        = LayoffDrift ;
sep     = sep + ExitDrift * (SUntrunc_zn<0) ;   
sepU    = sep .* exp(n_zn) ;

    % Hires from employment
TransitionMoments.HireRateE_TOPVA(t)   = sum( hiresE .* gdzdn_zn .* TOPva ) / sum( gndzdn_zn .* TOPva ) ;
TransitionMoments.HireRateE_LOWVA(t)   = sum( hiresE .* gdzdn_zn .* LOWva ) / sum( gndzdn_zn .* LOWva ) ;

    % Hires from nonemployment
TransitionMoments.HireRateU_TOPVA(t)   = sum( hiresU .* gdzdn_zn .* TOPva ) / sum( gndzdn_zn .* TOPva ) ;
TransitionMoments.HireRateU_LOWVA(t)   = sum( hiresU .* gdzdn_zn .* LOWva ) / sum( gndzdn_zn .* LOWva ) ;

    % Separations to employment
TransitionMoments.SepRateE_TOPVA(t)    = sum( quitsE .* gdzdn_zn .* TOPva ) / sum( gndzdn_zn .* TOPva ) ;
TransitionMoments.SepRateE_LOWVA(t)    = sum( quitsE .* gdzdn_zn .* LOWva ) / sum( gndzdn_zn .* LOWva ) ;

    % Separations to nonemployment
TransitionMoments.SepRateU_TOPVA(t)    = sum( sepU .* gdzdn_zn .* TOPva ) / sum( gndzdn_zn .* TOPva ) ;
TransitionMoments.SepRateU_LOWVA(t)    = sum( sepU .* gdzdn_zn .* LOWva ) / sum( gndzdn_zn .* LOWva ) ;


%%  CHURN-TYPE MEASURES
%__________________________________________________________________________

TransitionMoments.ExitRateUW(t)                  = sum(-L) * gdzdn_zn ; 
TransitionMoments.ExitRateW(t)                   = MF * sum(-L) * gndzdn_zn / (1-u) ;
TransitionMoments.JobCreationIncumbents(t)       = MF * sum( ( L .* ( max(SizeNewFirm-SizeOldFirm,0) ) ) * gdzdn_zn ) / (1-u) ;
TransitionMoments.JobDestructionIncumbents(t)    = MF * sum( ( L .* ( max(SizeOldFirm-SizeNewFirm,0) ) ) * gdzdn_zn ) / (1-u) ; 
TransitionMoments.JobCreation(t)                 = TransitionMoments.JobCreationIncumbents(t) + TransitionMoments.EntryRateW(t) ;
TransitionMoments.JobDestruction(t)              = TransitionMoments.JobDestructionIncumbents(t) + TransitionMoments.ExitRateW(t) ;
TransitionMoments.JobReallocation(t)             = TransitionMoments.JobCreation(t) + TransitionMoments.JobDestruction(t) ;
TransitionMoments.Hires(t)                       = TransitionMoments.HiresE(t) + TransitionMoments.HiresU(t) ;
TransitionMoments.Separations(t)                 = TransitionMoments.SeparationsE(t) + TransitionMoments.SeparationsU(t) ;
TransitionMoments.Churn(t)                       = (   TransitionMoments.Hires(t) / (1-u) ...
                                                     + TransitionMoments.Separations(t) / (1-u) ...
                                                     - TransitionMoments.JobCreation(t) ...
                                                     - TransitionMoments.JobDestruction(t) ...
                                                    ) / 2 ;

                                                
%%  HIRES BY SIZE
%__________________________________________________________________________

TransitionMoments.Hires1000(t) = ... 
    (    MF * q * ( ( phi + (1-phi) * v_zn .* Gn_zn ) .* size1000  )' * gdzdn_zn ...
      +  ME * sum(  exp( n_zn ) .* size1000 .* pi0Trunc_zn .* dzdn_zn ) ) ;   
TransitionMoments.Hires500(t) = ... 
    (    MF * q * ( ( phi + (1-phi) * v_zn .* Gn_zn ) .* size500  )' * gdzdn_zn ...
      +  ME * sum(  exp( n_zn ) .* size500 .* pi0Trunc_zn .* dzdn_zn ) ) ;
TransitionMoments.Hires50(t) = ... 
    (    MF * q * ( ( phi + (1-phi) * v_zn .* Gn_zn ) .* size50  )' * gdzdn_zn ...
      +  ME * sum(  exp( n_zn ) .* size50 .* pi0Trunc_zn .* dzdn_zn ) ) ;
TransitionMoments.Hires0(t) = ... 
     (    MF * q * ( ( phi + (1-phi) * v_zn .* Gn_zn ) .* size0  )' * gdzdn_zn ...
       +  ME * sum(  exp( n_zn ) .* size0 .* pi0Trunc_zn .* dzdn_zn ) ) ;
TransitionMoments.Hires00(t) = ... 
    (     MF * q * ( ( phi + (1-phi) * v_zn .* Gn_zn ) .* size00  )' * gdzdn_zn ...
       +  ME * sum(  exp( n_zn ) .* size00 .* pi0Trunc_zn .* dzdn_zn ) ) ;
TransitionMoments.VacYield(t) = ... 
    ( TransitionMoments.HiresE(t) + TransitionMoments.HiresU(t) ) ...
    / TransitionMoments.VacTot(t) ;
TransitionMoments.VacYield1000(t) = ... 
    TransitionMoments.Hires1000(t) ...
  / TransitionMoments.VacShare1000(t) ...
  / TransitionMoments.VacTot(t) ;
TransitionMoments.VacYield500(t) = ... 
    TransitionMoments.Hires500(t) ...
  / TransitionMoments.VacShare500(t) ...
  / TransitionMoments.VacTot(t) ;
TransitionMoments.VacYield50(t) = ... 
    TransitionMoments.Hires50(t) ...
  / TransitionMoments.VacShare50(t) ...
  / TransitionMoments.VacTot(t) ;
TransitionMoments.VacYield0(t) = ... 
    TransitionMoments.Hires0(t) ...
  / TransitionMoments.VacShare0(t) ...
  * TransitionMoments.VacTot(t)  ;
TransitionMoments.VacYield00(t) = ... 
    TransitionMoments.Hires00(t) ...
  / TransitionMoments.VacShare00(t) ...
  / TransitionMoments.VacTot(t)  ;


%%  JOB LADDER
%__________________________________________________________________________

    % Rank-rank correlation
[ ~,~, inds, inds_inv ]             = Cdf_Gn( reshape(Sn_zn,Nz,Nn), reshape(g_zn,Nz,Nn), NumGrids ) ;
SnRank                              = Cdf_Gv( inds, inds_inv, ones(Nz,Nn), reshape(g_zn,Nz,Nn), SUntrunc_znmat, exitrate, q, reshape(Gn_zn,Nz,Nn), phi, NumGrids) ;
MeanSnRank                          = sum( SnRank(:) .* gdzdn_zn ) ;
StdSnRank                           = sqrt( sum( SnRank(:).^2 .* gdzdn_zn ) - MeanSnRank^2 ) ;
MeanSnRankFinal                     = sum( SnRankFinal(:) .* gdzdn_zn ) ;
StdSnRankFinal                      = sqrt( sum( SnRankFinal(:).^2 .* gdzdn_zn ) - MeanSnRankFinal^2 ) ;
TransitionMoments.RankRankCorr(t)   = ( ( SnRank(:)' - MeanSnRank ) .* ( SnRankFinal(:)' - MeanSnRankFinal ) ) * gdzdn_zn / ( StdSnRank * StdSnRankFinal ) ;

    % Percentiles of Sn distribution
[ ~, inds ]                 = sort( Sn_zn ) ;
[ ~, inds_inv ]             = sort( inds ,'ascend' ) ;
hS                          = g_zn(inds) .* dzdn_znmat(inds) ;
CumH                        = cumsum( hS ) ;
p10                         = CumH >= 0.05 & CumH < 0.15 ;
p50                         = CumH >= 0.45 & CumH < 0.55 ;
p90                         = CumH >= 0.85 & CumH < 0.95 ; 
P10                         = p10(inds_inv) ;
P50                         = p50(inds_inv) ;
P90                         = p90(inds_inv) ;
TransitionMoments.SnP10(t)  = sum( Sn_zn(P10) .* gdzdn_zn(P10) ) / sum( gdzdn_zn(P10) ) ;
TransitionMoments.SnP50(t)  = sum( Sn_zn(P50) .* gdzdn_zn(P50) ) / sum( gdzdn_zn(P50) ) ;
TransitionMoments.SnP90(t)  = sum( Sn_zn(P90) .* gdzdn_zn(P90) ) / sum( gdzdn_zn(P90) ) ;

    % Rank-rank correlation
[ ~ , indsVA ]                      = sort(lvapw) ;
[ ~ , indsinvVA ]                   = sort(indsVA) ;
hSortVA                             = gdzdn_zn(indsVA) ;
rankVA                              = cumsum( hSortVA ) ;
rankVA                              = rankVA(indsinvVA) ;
[ ~ , indsSN ]                      = sort(sn_zn) ;
[ ~ , indsinvSN ]                   = sort(indsSN) ;
hSortSN                             = gdzdn_zn(indsSN) ;
rankSN                              = cumsum( hSortSN ) ;
rankSN                              = rankSN(indsinvSN) ;
mVA                                 = rankVA' * gdzdn_zn ;
mSN                                 = rankSN' * gdzdn_zn ;
varVA                               = ( ( rankVA - mVA ).^2 )' * gdzdn_zn ;
varSN                               = ( ( rankSN - mSN ).^2 )' * gdzdn_zn ;
corrVASN                            = ( ( (rankVA - mVA ) .* ( rankSN - mSN ) )' * gdzdn_zn ) / sqrt( varVA * varSN ) ;
TransitionMoments.RankRankVASN(t)   = corrVASN ;


end
end