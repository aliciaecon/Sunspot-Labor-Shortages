%% ========================================================================
% 
%     CONSTRUCTS VACANCY-WEIGHTED DISTRIBUTION OF MARGINAL SURPLUSES
% 
% =========================================================================

function Gv_znmat = Cdf_Gv( inds, inds_inv, v_znmat, g_znmat, SUntrunc_znmat,...
    exitrate, q, Gn_znmat, phi, NumGrids)

    % Load numerical objects
Nz          = NumGrids.Nz ;         % number of points in z grid
Nn          = NumGrids.Nn ;         % number of points in n grid
n_znmat     = NumGrids.n_znmat;     % Replica of the n grid, in logs
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn as matrix, where d is log difference)
in0         = NumGrids.in0;         % index of n gridpoint closest to n_0
pi0_z       = NumGrids.pi0_z;       % Distribution of entering firms in z grid

    % Distribution of entrants
pi0_znmat                   = zeros( Nz , Nn ) ; % Distribution of entering firms
pi0_znmat(:,in0)            = pi0_z ;
pi0_znmat(SUntrunc_znmat<=0)= 0 ; % Truncate by negative surplus
pi0Trunc_znmat              = pi0_znmat / sum(sum( pi0_znmat .* dzdn_znmat)) ; % Ensure integration to 1 

    % Entrant vacancies
HirePerVacancy_znmat = q * ( phi + (1-phi) * Gn_znmat ) ;
VacancyPerHire_znmat = 1 ./ HirePerVacancy_znmat ;  
vgdzdn_znmat_Ent     = exitrate * exp(n_znmat) .* VacancyPerHire_znmat .* pi0Trunc_znmat .* dzdn_znmat ; % Vacancies of entrants, recall n is in logs

    % Incumbent vacancies
gdzdn_znmat         = g_znmat .* dzdn_znmat ;                   % Integrate against g(z,n)dz*dn
vgdzdn_znmat_Inc    = v_znmat .* gdzdn_znmat;                   

    % Total vacancies
vgdzdn_znmat        = vgdzdn_znmat_Inc + vgdzdn_znmat_Ent ;     % Incumbents + Entrants
vgdzdn_zn           = reshape( vgdzdn_znmat , Nz*Nn , 1 );      % Vector form to calculate integral

    % Compute Gv in space sorted by marginal surplus
GvSort_zn = cumsum(vgdzdn_zn(inds));

    % Sanity checks and normalize
mF = max(GvSort_zn);
if mF > 0
    GvSort_zn = GvSort_zn / mF;
else
    GvSort_zn = zeros( Nz*Nn , 1 );
end
GvSort_zn = GvSort_zn / nanmax(GvSort_zn) ; % Normalization, equivalent to dividing by total vacancies

    % Unsort and reshape as matrix
Gv_zn       = GvSort_zn(inds_inv) ; % Gv1 is unsorted CDF
Gv_znmat    = reshape( Gv_zn , Nz , Nn );

end