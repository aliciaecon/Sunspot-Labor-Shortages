    % Re-compute relevant aggregates
[ Gn_znmat, Ghatn_znmat, inds, inds_inv ] = Cdf_Gn( Sn_znmat, g_znmat, NumGrids ) ;
H        = g_znmat(:) .* NumGrids.dzdn_znmat(:) ;
H        = H / sum(H(:)) ;
exitrate = sum(-L) * H(:) ; 

    % Assign initial guess in Hopenhayn
ExogParams.InitialGuess.g_znmat     = g_znmat ;
ExogParams.InitialGuess.u           = u ;
ExogParams.InitialGuess.p           = p ;
ExogParams.InitialGuess.q           = q ;
ExogParams.InitialGuess.phi         = phi ;
ExogParams.InitialGuess.Gn_znmat    = Gn_znmat ;
ExogParams.InitialGuess.Ghatn_znmat = Ghatn_znmat ;
ExogParams.InitialGuess.Gv_znmat    = Gv_znmat ;
ExogParams.InitialGuess.S_znmat     = S_znmat ;
ExogParams.InitialGuess.Sn_znmat    = Sn_znmat ;
ExogParams.InitialGuess.v_znmat     = v_znmat ;
ExogParams.InitialGuess.exitrate    = exitrate ;