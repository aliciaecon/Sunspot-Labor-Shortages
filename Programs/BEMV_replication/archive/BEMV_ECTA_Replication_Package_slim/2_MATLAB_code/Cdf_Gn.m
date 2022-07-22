%% ========================================================================
% 
%    CONSTRUCTS EMPLOYMENT-WEIGHTED DISTRIBUTION OF MARGINAL SURPLUSES
% 
% =========================================================================

function [ Gn_znmat, Ghatn_znmat, inds, inds_inv ] = Cdf_Gn( Sn_znmat, g_znmat, NumGrids )

% INPUTS:
% g_znmat   -> Unweighted distribution of firms over (z,n)
% Sn_znmat  -> Marginal surplus of firms over (z,n)

    % Load numerical objects
Nz          = NumGrids.Nz ;         % number of points in z grid
Nn          = NumGrids.Nn ;         % number of points in n grid
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid, in logs
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn as matrix, where d is log difference)

    % Sort marginal surpluses
    % Ranks depend on S(n,z) in levels, not Stilde(ntilde,ztilde) which is
    % in logs. Therefore transform back into levels using line from
    % computational appendix: S_n(n,z)n = Stilde_ntilde(ntilde,ztilde)
Sn_znmat            = Sn_znmat .* exp(-n_znmat) ;       % Transform marginal surplus: difference of S w.r.t. n instead of logn
Sn_zn               = reshape( Sn_znmat , Nz*Nn , 1 ) ; % Vector form
[ SnSort_zn , inds ]= sort( real(Sn_zn) , 'ascend' ) ;  % Sort
[ ~ , inds_inv ]    = sort( inds , 'ascend' ) ;         % Unsort mapping

    % Compute GnSort, employment-weighted CDF of marginal surpluses, in sorted space
gdzdn_znmat     = g_znmat .* dzdn_znmat ;                   % Integrate against g(z,n)dz*dn
ngdzdn_znmat    = exp(n_znmat) .* gdzdn_znmat;              % Employment times pdf
ngdzdn_zn       = reshape( ngdzdn_znmat , Nz*Nn , 1 ) ;     % Vector
GnSort_zn       = cumsum( ngdzdn_zn(inds) ) ;               % Get cdf by integrating employment weighted pdf, sorted by marginal surplus
GnSort_zn       = [ 0 ; GnSort_zn(1:end-1) ] ;              % Hiring from those strictly below
GnSort_zn       = GnSort_zn / nanmax(GnSort_zn) ;           % Normalize (equivalent to dividing by total employment)

    % Compute Ghatn, integral of Gn, in sorted space
dSN_zn       = [ SnSort_zn(2) - SnSort_zn(1) ; SnSort_zn(2:Nz*Nn) - SnSort_zn(1:(Nz*Nn-1)) ] ; % This has to be constructed winding downwards. Recall also Sn is difference w.r.t. n instead of logn
GhatnSort_zn = cumsum( GnSort_zn .* dSN_zn ); % Integral of Gn from 0 to Sn

    % Unsort and reshape as matrix
Gn_zn       = GnSort_zn(inds_inv);              % Gn_zn is unsorted CDF
Ghatn_zn    = GhatnSort_zn(inds_inv);           % Ghatn_zn is unsorted integral
    % Ghatn_zn is the object that appears in firm (z,n)'s Bellman equation
    % Gives probability of hiring conditional on a vacancy meeting an
    % employed worker at some other firm
Gn_znmat    = reshape( Gn_zn , Nz , Nn );
Ghatn_znmat = reshape( Ghatn_zn , Nz , Nn );

end