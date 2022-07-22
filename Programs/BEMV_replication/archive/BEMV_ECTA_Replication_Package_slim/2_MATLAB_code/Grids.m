%% ========================================================================
% 
%                            DEFINES GRIDS
% 
% =========================================================================
% Note that the grids are in logs. Therefore, dz is log difference between
% different points in z grid for example.

function [ NumGrids ] = Grids( Numerical, Params )

    % Load numerical objects and parameters
Nz      = Numerical.Nz ;    % Number of z gridpoints
Nn      = Numerical.Nn ;    % Number of n gridpoints
zmin    = Numerical.zmin ;  % Lowest z gridpoint
zmax    = Numerical.zmax ;  % Highest z gridpoint
nmin    = Numerical.nmin ;  % Lowest z gridpoint
nmax    = Numerical.nmax ;  % Highest z gridpoint
n_0     = Params.n_0 ;      % Initial employment of entering firms
zeta    = Params.zeta ;     % Shape of productivity distribution of entering firms

    % Grids in logs
z = linspace(zmin,zmax,Nz); % Define z grid in log space
n = linspace(nmin,nmax,Nn); % Define n grid in log space

    % Useful objects from grids
    % _zn     -> Vector = ([z1,n1],[z2,n1],...,[zNz,n1],[z1,n2],...,[zNz,n2],...[z1,nNn],...,[zNz,nNn])
    % X_znmat -> Matrix = reshape(X_zn,Nz,Nn)
    % X_zn              = reshape(X_znmat,Nz*Nn,1)
dz_z        = [z(2) - z(1),z(2:Nz) - z(1:(Nz-1))];  % Step size for z
dn_n        = [n(2) - n(1),n(2:Nn) - n(1:(Nn-1))];  % Step size for n
dzdn_znmat  = kron(dz_z',dn_n);                     % Measure as a (z,n) matrix
n_znmat     = kron(ones(Nz,1),n);                   % Replica of the n grid as a (z,n) matrix
z_znmat     = kron(z',ones(1,Nn));                  % Replica of the z grid as a (z,n) matrix
dn_znmat    = kron(ones(Nz,1),dn_n);                % Replica of the n step size as a (z,n) matrix
dz_znmat    = kron(dz_z',ones(1,Nn));               % Replica of the z step size as a (z,n) matrix
dzdn_zn     = reshape( dzdn_znmat , Nz*Nn , 1 );    % Measure as a (z,n) vector
n_zn        = reshape( n_znmat , Nz*Nn , 1 );       % Replica of the n grid as a (z,n) vector
z_zn        = reshape( z_znmat , Nz*Nn , 1 );       % Replica of the z grid as a (z,n) vector
dn_zn       = reshape( dn_znmat , Nz*Nn , 1 );      % Replica of the n step size as a (z,n) vector
dz_zn       = reshape( dz_znmat , Nz*Nn , 1 );      % Replica of the z step size as a (z,n) vector

    % Characterize firm entry in grids
[ ~, in0 ]  = min( abs( n - log(n_0) ) ) ;   % index of n gridpoint closest to n_0
pi0_z       = exppdf( z-z(1) , 1/zeta ) ; % Distribution of entering firms in z grid (note: does not sum to 1)

    % Finer grids (for some computations)
Nz2         = 10*Nz ;                                    % Define how much finer
Nn2         = 10*Nn ;                                    % Define how much finer
z2          = linspace(zmin,zmax,Nz2);                   % Define z grid in log space
n2          = linspace(nmin,nmax,Nn2);                   % Define n grid in log space
dz2_z       = [z2(2) - z2(1),z2(2:Nz2) - z2(1:(Nz2-1))]; % Step size for z
dn2_n       = [n2(2) - n2(1),n2(2:Nn2) - n2(1:(Nn2-1))]; % Step size for n
dzdn2_znmat = kron(dz2_z',dn2_n);                        % Measure
n2_znmat    = kron(ones(Nz2,1),n2);                      % Replica of the n grid
z2_znmat    = kron(z2',ones(1,Nn2));                     % Replica of the z grid
dn2_znmat   = kron(ones(Nz2,1),dn2_n);                   % Replica of the n step size
dz2_znmat   = kron(dz2_z',ones(1,Nn2));                  % Replica of the z step size

    % Build struct
NumGrids.Nz         = Nz ;
NumGrids.Nn         = Nn ;
NumGrids.z          = z ;
NumGrids.n          = n ;
NumGrids.dz_z       = dz_z ;
NumGrids.dn_n       = dn_n ;
NumGrids.dzdn_znmat = dzdn_znmat ;
NumGrids.n_znmat    = n_znmat ;
NumGrids.z_znmat    = z_znmat ;
NumGrids.dn_znmat   = dn_znmat ;
NumGrids.dz_znmat   = dz_znmat ;
NumGrids.dzdn_zn    = dzdn_zn ;
NumGrids.n_zn       = n_zn ;
NumGrids.z_zn       = z_zn ;
NumGrids.dn_zn      = dn_zn ;
NumGrids.dz_zn      = dz_zn ;
NumGrids.in0        = in0 ;
NumGrids.pi0_z      = pi0_z ;
NumGrids.Nz2        = Nz2 ;
NumGrids.Nn2        = Nn2 ;
NumGrids.z2         = z2 ;
NumGrids.n2         = n2 ;
NumGrids.dz2_z      = dz2_z ;
NumGrids.dn2_n      = dn2_n ;
NumGrids.dzdn2_znmat= dzdn2_znmat ;
NumGrids.n2_znmat   = n2_znmat ;
NumGrids.z2_znmat   = z2_znmat ;
NumGrids.dn2_znmat  = dn2_znmat ;
NumGrids.dz2_znmat  = dz2_znmat ;

end