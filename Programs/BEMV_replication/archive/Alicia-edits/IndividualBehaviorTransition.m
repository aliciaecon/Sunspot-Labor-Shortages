%% ========================================================================
% 
%          UPDATES SURPLUS FUNCTIONS ONLY ONCE IN VFI
% 
% =========================================================================

function [ S, Sn, SUntrunc, SnUntrunc, v ] = IndividualBehaviorTransition...
        ( S_znmat, Sn_znmat, GG, NumGrids, ExogParams, Params, Zmatrices, options,...
        PerPeriod, Dmatrices, SnGrid, TransitionNum, TransitionEq, TransitionShock, t )

    % Get aggregates to compute flow payoffs and numerical scheme
DeltaT  = TransitionNum.dT(t) ;
GGT     = reshape( GG , NumGrids.Nz , NumGrids.Nn ) ;
qT      = TransitionEq.q(t+1) ;
phiT    = TransitionEq.phi(t+1) ;

    % Get shock to discount factor
ExogParams.rho = TransitionShock.rho(t) ;

    % Perform VFI once
[ ~, S, Sn, SUntrunc, SnUntrunc, v ] = IndividualBehavior( S_znmat,...
    Sn_znmat, qT, phiT, GGT, SnGrid, DeltaT, NumGrids, ExogParams,...
    Params, Zmatrices, PerPeriod, Dmatrices, options );

end


