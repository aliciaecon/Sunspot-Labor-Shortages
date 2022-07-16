%% ========================================================================
% 
%          UPDATES AGGREGATE STATES GIVEN SURPLUS FUNCTIONS
% 
% =========================================================================

function [ error_g, error_f, error_Gn, error_Gv, ...
            q, phi, p, u, Gn_znmat, Ghatn_znmat, Gv_znmat, g_znmat, exitrate, ...
            L, pi0Trunc_zn, gprenorm_znmat ] = ...
        AggregateBehavior( q, phi, p, Gn_znmat, Gv_znmat, g_znmat,...
                            SUntrunc_znmat, Sn_znmat, SnUntrunc_znmat, v_znmat, exitrate,...
                            Numerical, NumGrids, Derivatives, MatExog, Params, ExogParams, options )

    % STEP IIa - Update the employment-weighted distribution
[ GnNew_znmat, Ghatn_znmat, inds, inds_inv ]= Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ;
error_Gn                                    = max(max( abs(GnNew_znmat - Gn_znmat) )) ;
Gn_znmat                                    = GnNew_znmat ;

    % STEP IIb - Update the vacancy-weighted distribution
GvNew_znmat = Cdf_Gv( inds, inds_inv, v_znmat, g_znmat, SUntrunc_znmat, exitrate, q, Gn_znmat, phi, NumGrids ) ;
error_Gv    = max(max( abs(GvNew_znmat - Gv_znmat) )) ;
Gv_znmat    = GvNew_znmat ;

    % STEP IIc - Solve for stationary distribution of firms
[ error_g , g_znmat , pi0Trunc_zn , exitrate , L , gprenorm_znmat ] = ...
    Distribution( q, phi, p, Gn_znmat, Gv_znmat, g_znmat, SUntrunc_znmat, SnUntrunc_znmat, v_znmat,...
                    Numerical, NumGrids, Derivatives, MatExog, Params, ExogParams, options ) ;
exitrate = min(max( exitrate , 0.0001 ), 0.5 ) ; % caps for the exit rate

    %  STEP IId - Update finding rates (Uses updated: g_znmat, Gn_znmat)
[ error_f, q, phi, p, u ] = ...
    FindingRates( q, phi, p, g_znmat, v_znmat, exitrate, Gn_znmat, SUntrunc_znmat, NumGrids, Params, ExogParams ) ;

end