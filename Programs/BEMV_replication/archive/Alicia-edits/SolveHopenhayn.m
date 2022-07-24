%% ========================================================================
% 
%               SOLVES THE PROBLEM OF COMPARATIVE STATICS
% 
% =========================================================================

function [ SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat,...
    g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, ExogParams ] ...
    = SolveHopenhayn( Params, entrycost, ExogParams, options, gridpoint, dir )
       
    % Loop over the mass of active firms, build grid
NumPoints = 100 ;
if strcmp(dir,'up') == 1
    if gridpoint <= 1.5
        lowerscale = 1 ;
        upperscale = 1.25 ;
    elseif gridpoint > 1.5
        lowerscale = 0.98 ;
        upperscale = 1.25 ;
    end
    MFgrid = linspace( lowerscale*ExogParams.MF , upperscale*ExogParams.MF , NumPoints ) ;
elseif strcmp(dir,'down') == 1
    if gridpoint > 0.6
        lowerscale = 0.75 ;
        upperscale = 1 ;
    elseif gridpoint <= 0.6
        lowerscale = 0.5 ;
        upperscale = 0.95 ;
    end
    MFgrid = fliplr(linspace( lowerscale*ExogParams.MF , upperscale*ExogParams.MF , NumPoints )) ;
end    
error_hop = nan(length(MFgrid),1) ;

    % Setup
disp('Inside Solve Hopenhayn')
disp(['*** Initial MF = ' num2str(ExogParams.MF,3)])
i     = 1 ;
error = 1 ;

    % Main Loop
while i <= length(MFgrid) && error > 0
    
        % Solve steady state for value of grid
    ExogParams.MF = MFgrid(i) ;
    [ ~,~,~,~,~, S_znmat, Sn_znmat, v_znmat, g_znmat, Gn_znmat, Gv_znmat, q,...
        phi, p, u, L, ~, NumGrids] = SolveBEMV( Params, ExogParams, options ) ;
     
        % Compute error
    pi0_znmat                   = zeros( NumGrids.Nz , NumGrids.Nn ) ; 
    pi0_znmat(:,NumGrids.in0)   = NumGrids.pi0_z ;
    valueofentry                = S_znmat(:)' * ( pi0_znmat(:) .* NumGrids.dzdn_znmat(:) )  ;
    error_hop(i)                = ( valueofentry - entrycost ) / entrycost ;
    error                       = error_hop(i) ;
    
        % Display progress
    if strcmp(dir,'up')==1
        if i == 1 && error < 0
            fprintf('negative initial error=%6.4f while trying to move up \n',error )
            MFgrid = 0.975 * MFgrid ;
            error  = 1 ;
        else
            i = i + 1 ;
            run AssignHopenhaynInitialGuess.m
        end
    elseif strcmp(dir,'down')==1
        if i == 1 && error > 0
            fprintf('positive initial error=%6.4f while trying to move down \n',error )
            MFgrid = 1.025 * MFgrid ;
            error  = - 1 ;
        else
            i = i + 1 ;
            run AssignHopenhaynInitialGuess.m
        end
    end
    fprintf('Chi/Chi0=%6.2f, i=%2.0f, error=%6.4f, MF=%6.4f \n',gridpoint,i,error,MFgrid(i))
    if strcmp(dir,'down')==1
        error = - error ; % to work with the while loop based on error>0
    end
end

    % Close loop
i = i - 1 ;
MF = MFgrid(i-1) + abs( ( 0 - error_hop(i-1) ) / ( error_hop(i) - error_hop(i-1) ) ) * ( MFgrid(i) - MFgrid(i-1) );
ExogParams.MF = MF ;
[ ~,~,~, SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat,...
    g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, ~,~,~,~,~ ] ...
    = SolveBEMV( Params, ExogParams, options ) ;

end





