function [ m ] = Hmean( X , selection , H )
% Computes mean of X weighted by H

m = sum( X(selection) .* H(selection) ) / sum( H(selection) ) ;

end

