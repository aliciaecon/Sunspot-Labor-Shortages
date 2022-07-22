%% ========================================================================
% 
%                 CONSTRUCT MATRICES TO COMPUTE DERIVATIVES
% 
% =========================================================================
% Entries are simply 1/dn and 1/dz, where dn and dz are log differences

function [ Derivatives ] = MatricesDerivatives( NumGrids )
%%  LOAD NUMERICAL OBJECTS
%__________________________________________________________________________

Nz      = NumGrids.Nz ;     % Number of points in z grid
Nn      = NumGrids.Nn ;     % Number of points in n grid
dz_z    = NumGrids.dz_z ;   % Log step size in z grid
dn_n    = NumGrids.dn_n ;   % Log step size in n grid

%%  DERIVATIVES W.R.T. LOG N
%__________________________________________________________________________

    % Backward approximation for value computation
Diag        = kron( 1./dn_n' , ones(Nz,1) ) ;
Diag(1:Nz)  = 0 ; 
SubDiag     = kron( -1./dn_n' , ones(Nz,1) ) ; 
SubDiag     = [ SubDiag(Nz+1:end) ; zeros(Nz,1) ] ;
DnBack      = spdiags( [SubDiag,Diag] , [-Nz,0] , Nz*Nn , Nz*Nn ) ;

    % Forward approximation for value computation
Diag                = kron( -1./dn_n' , ones(Nz,1) ) ;
Diag(end-Nz+1:end)  = 0 ;
SupDiag             = kron( 1./dn_n' , ones(Nz,1) ) ; 
SupDiag             = [ zeros(Nz,1) ; SupDiag(1:(Nn-1)*Nz) ] ;
DnForw              = spdiags( [Diag,SupDiag] , [0,Nz] , Nz*Nn , Nz*Nn ) ;

    % Approximations for distribution (absorbing state at the bottom)
DnBackd =  DnForw' ;
DnForwd = -DnBack' ;

%%  DERIVATIVES W.R.T. LOG Z
%__________________________________________________________________________

    % Backward approximation for value function
Diag    = 1./dz_z' ;
Diag(1) = -1/dz_z(1) ; % approximate backward derivative with forward derivative at bottom 
SubDiag = -1./dz_z' ;
SubDiag = [ SubDiag(2:end) ; zeros(1,1) ] ;
SupDiag = [ zeros(1,1) ; 1./dz_z(1) ; zeros(Nz-2,1) ] ; % just for forward derivative at bottom 
base    = spdiags( [SubDiag,Diag,SupDiag] , [-1,0,1] , Nz , Nz ) ;
DzBack  = kron( speye(Nn) , base ) ;

    % Forward approximation for value function
Diag        = -1./dz_z' ;
Diag(end)   = 1/dz_z(Nz) ; % approximate forward derivative with backward derivative at top
SupDiag     = 1./dz_z' ;
SupDiag     = [ zeros(1,1) ; SupDiag(1:Nz-1) ] ;
SubDiag     = [ zeros(Nz-2,1) ; -1./dz_z(Nz) ; zeros(1,1) ] ; % just for backward derivative at top
base        = spdiags( [SubDiag,Diag,SupDiag] , [-1,0,1] , Nz , Nz ) ;
DzForw      = kron( speye(Nn) , base ) ;

    % Forward approximation for distribution when drift is negative (mass falls below diagonal)
SupDiag = 1./dz_z' ; % boundary condition at bottom
Diag    = -1./dz_z' ; 
SubDiag = [ 1/dz_z(1) ; zeros(Nz-1,1) ] ;
base    = spdiags( [SubDiag,Diag,SupDiag] , [-1,0,1] , Nz , Nz ) ;
DzForwd = kron( speye(Nn) , base ) ;

    % Backward approximation for distribution (mass falls ABOVE diagonal)
SubDiag = 1./dz_z' ; % boundary condition at top
Diag    = 1./dz_z' ;
SupDiag = [ zeros(Nz-1,1) ; 1/dz_z(end) ] ;
base    = spdiags( [SubDiag,-Diag,SupDiag] , [-1,0,1] , Nz , Nz ) ;
DzBackd = kron( speye(Nn) , base ) ;

    % Second derivative for value computation
Diag            = -2./dz_z'.^2 ;
SupDiag         = 1./dz_z'.^2 ; 
SupDiag         = [ 0 ; SupDiag(1:end-1) ] ;
SubDiag         = 1./dz_z'.^2 ; 
SubDiag         = [ SubDiag(2:end) ; 0 ] ;
base            = spdiags( [SubDiag,Diag,SupDiag] , [-1,0,1] , Nz , Nz );
base(1,1)       = -1/dz_z(1)^2;
base(1,2)       = 1/dz_z(1)^2;
base(Nz,Nz)     = -1/dz_z(Nz)^2;
base(Nz,Nz-1)   = 1/dz_z(Nz)^2;
Dzz             = kron(speye(Nn),sparse(base)) ;

    % Second derivative for distribution
Dzzd            = Dzz' ;

%%  BUILD STRUCT
%__________________________________________________________________________

Derivatives.DnBack  = DnBack ;
Derivatives.DnForw  = DnForw ;
Derivatives.DnBackd = DnBackd ;
Derivatives.DnForwd = DnForwd ;

Derivatives.DzBack  = DzBack ;
Derivatives.DzForw  = DzForw ;
Derivatives.DzBackd = DzBackd ;
Derivatives.DzForwd = DzForwd ;

Derivatives.Dzz     = Dzz ;
Derivatives.Dzzd    = Dzzd ;

end