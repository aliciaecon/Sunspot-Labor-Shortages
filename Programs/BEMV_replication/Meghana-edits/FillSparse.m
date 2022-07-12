%% ========================================================================
% 
%                 DEFINES CONVENIENT SPARSE MATRICES
% 
% =========================================================================

% Fills in a sparse matrix with three bands in the same places as
% repmat(s), and consistent with upwinding scheme. The input s is a
% (Nz*Nn)x1 vector which forms the diagonal of the sparse matrices sForw
% and sBack, which contain the positive and negative values of s
% respectively.

function [ sForw , sBack , Forw , Back ] = FillSparse( s , NumGrids )

    % Load Nz and Nn
Nz = NumGrids.Nz; % number of points in z grid
Nn = NumGrids.Nn; % number of points in n grid
    
    % Positive and negative indicator
sPos = ( s >= 0 ) ;
sNeg = ( s <  0 ) ;

    % Positive and negative part of the s vector
sPlus =  s .* sPos ;
sMinus = s .* sNeg ;

    % Upward- and downward-shifted versions of the previous vectors to use spdiags
sUpPos     = [ sPos(Nz+1:end)   ; zeros(Nz,1)      ] ;
sDownPos   = [ zeros(Nz,1)      ; sPos(1:end-Nz)   ] ;

sUpNeg     = [ sNeg(Nz+1:end)   ; zeros(Nz,1)      ] ;
sDownNeg   = [ zeros(Nz,1)      ; sNeg(1:end-Nz)   ] ;

sUpPlus    = [ sPlus(Nz+1:end)  ; zeros(Nz,1)      ] ;
sDownPlus  = [ zeros(Nz,1)      ; sPlus(1:end-Nz)  ] ;

sUpMinus   = [ sMinus(Nz+1:end) ; zeros(Nz,1)      ] ;
sDownMinus = [ zeros(Nz,1)      ; sMinus(1:end-Nz) ] ;

    % Compute multiplier of forward and backward approximation matrices   
sForw = spdiags( [sUpPlus,sPlus,sDownPlus]    , [-Nz,0,Nz] , Nz*Nn , Nz*Nn ) ;
sBack = spdiags( [sUpMinus,sMinus,sDownMinus] , [-Nz,0,Nz] , Nz*Nn , Nz*Nn ) ;
Forw  = spdiags( [sUpPos,sPos,sDownPos]       , [-Nz,0,Nz] , Nz*Nn , Nz*Nn ) ;
Back  = spdiags( [sUpNeg,sNeg,sDownNeg]       , [-Nz,0,Nz] , Nz*Nn , Nz*Nn ) ;

end