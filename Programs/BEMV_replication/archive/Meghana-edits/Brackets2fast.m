function [iBot,iTop,iBotStack,iTopStack,frac] = Brackets2fast(BracketingPoints,BracketedPoints,N1,N2,M1)

% N1 N2 is the size of BracketingPoints, and M1 is the size of BracketedPoints
% The bracketing(basis) is done along the first dimension of BracketingPoints,
% for each index of the second dimension on BracketingPoints
% the BracketingPoints and BrackedPoints are both columns

BracketingSquare = reshape(BracketingPoints,N1,N2);

% change start here

%% sort each row of bracheting points, keep track of the mapping
% the build-in matrix for matlab use quick sort
[Bracketingsort,I] = sort(BracketingSquare,'ascend') ;
Bracketingsort = reshape(Bracketingsort, N1*N2,1);

% Define scale to shift all values, and vectorize
scale0 = 2 * ( 10 + max(max(BracketingPoints(:)),max(BracketedPoints(:))) - min(min(BracketingPoints(:)),min(BracketedPoints(:)))) ;
%scale0 = 100 ;
scale = scale0*kron(1:N2,ones(1,N1))';

Bracketing = Bracketingsort + scale;
scale2 = scale0*kron(1:N2,ones(1,M1))';
Bracketed = scale2 + repmat(BracketedPoints,N2, 1);

%find the index with scale
iBotscale = discretize(Bracketed, Bracketing);
iBotscale(Bracketed < Bracketing(1)) = 1;
iBotscale(Bracketed > Bracketing(end)) = length(Bracketing)-1;

iBotsquare = reshape(iBotscale,M1,N2) - repmat(linspace(0,N2-1,N2)*N1,M1,1);

% Correct ends of grid; deal with the boundary here
iBotsquare(iBotsquare == 0) = 1;
iBotsquare(iBotsquare == N1) = N1-1;

iTopsquare = iBotsquare +1;


iBotSort = reshape(iBotsquare, M1*N2,1);
iTopSort = reshape(iTopsquare, M1*N2,1);

%%unsort ibot, iTop, using the indices,
%baseindex = ( kron((1:N2)',ones(M1,1))- 1)*M1;
%typo here, baseindex should be N1

baseindex = ( kron((1:N2)',ones(M1,1))- 1)*N1;



iBotSortStack = iBotSort + baseindex;
iTopSortStack = iTopSort+ baseindex;


iBot = I(iBotSortStack);
iTop = I(iTopSortStack);

iBotStack = iBot + baseindex;
iTopStack = iTop + baseindex;

SolBot = BracketingPoints(iBotStack);
SolTop = BracketingPoints(iTopStack);

frac = ( repmat(BracketedPoints,N2,1) - SolBot )./(SolTop-SolBot) ;
frac(isnan(frac)) = 0.5;





end





