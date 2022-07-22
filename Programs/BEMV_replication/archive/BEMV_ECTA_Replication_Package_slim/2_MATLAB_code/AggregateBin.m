function [ BinnedYX , VecX ] = AggregateBin(y,x,w,bins,distrib,RankOption,option,WeightDistribOption)

% x is the grid (e.g. age grid)
% y is the variable we bin (e.g. density over age)
% w are weights for the denominator in case y is extensive and we want to make the output
%       be adimensional (e.g. equal weights in the case of age bins)
% bins are the bins we're grouping the variable into (e.g. agebins)
% distrib is a distribution to weight observations y

if nargin == 7
    Wdistrib = distrib ;
elseif nargin == 8
    if strcmp(WeightDistribOption,'WeightDistribution')
        Wdistrib = distrib .* w ;
        Wdistrib = Wdistrib / sum(Wdistrib(:)) ;
    else
        error('Invalid Weight Distribution argument')
    end
end

if strcmp(RankOption,'Ranks')
    [~ , Sort] = sort(x(:)) ;
    [~ , Unsort ] = sort(Sort) ;
    xRankSort = cumsum( Wdistrib(Sort) ) ;
    xRank = xRankSort(Unsort) ;
    x = xRank;
elseif strcmp(RankOption,'Levels')
else
    error('Invalid rank option');
end

VecX = ones(size(x))*length(bins) ;
for i=1:( max(size(bins)) - 1 )
    arg = x < bins(i+1) & x >= bins(i) ;
    VecX(arg) = i ;
end

if strcmp(option,'sum')
    BinnedYX = accumarray( VecX(:) , y(:) .* distrib(:) ) ;
elseif strcmp(option,'mean')
    BinnedYX = accumarray( VecX(:) , y(:) .* distrib(:) ) ./ accumarray( VecX(:) , w(:) .* distrib(:) ) ;
else
    error('AggregateBin option not well defined')
end
BinnedYX = [ BinnedYX ; nan(length(bins)-length(BinnedYX),1) ] ;

end

