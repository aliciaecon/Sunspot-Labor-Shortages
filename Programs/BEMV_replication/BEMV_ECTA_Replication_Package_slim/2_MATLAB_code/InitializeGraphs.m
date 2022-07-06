function TransitionMoments = InitializeGraphs(TransitionMoments,SsMoments)

TransitionMoments.u(1)      = SsMoments.u(end) ;  
TransitionMoments.MeanY(1)  = SsMoments.MeanY(end) ;
TransitionMoments.YPW(1)    = SsMoments.YPW(end) ;
TransitionMoments.ME(1)     = SsMoments.ME(end) ;
TransitionMoments.JobCreation(1) = SsMoments.JobCreation(end) ;  
TransitionMoments.JobDestruction(1) = SsMoments.JobDestruction(end) ;  
TransitionMoments.V(1)      = SsMoments.V(end) ;
TransitionMoments.EU(1)     = TransitionMoments.EU(end) ;
TransitionMoments.UE(1)     = SsMoments.UE(end) ;
TransitionMoments.EEh(1)          = SsMoments.EEh(end) ;
TransitionMoments.VacYield(1)     = SsMoments.VacYield(end) ;
TransitionMoments.VacYield1000(1) = SsMoments.VacYield1000(end) ;
TransitionMoments.VacYield00(1)   = SsMoments.VacYield00(end) ;
TransitionMoments.VacShareTopSn(1) = SsMoments.VacShareTopSn(end) ;
TransitionMoments.VacShareLowSn(1) = SsMoments.VacShareLowSn(end) ;
TransitionMoments.VacShareTopva(1) = SsMoments.VacShareTopva(end) ;
TransitionMoments.VacShareLowva(1) = SsMoments.VacShareLowva(end) ;
TransitionMoments.VacShare1000(1) = SsMoments.VacShare1000(end) ;
TransitionMoments.VacShare00(1) = SsMoments.VacShare00(end) ;
TransitionMoments.NetPoachTopSn(1) = SsMoments.NetPoachTopSn(end) ;
TransitionMoments.NetPoachLowSn(1) = SsMoments.NetPoachLowSn(end) ;
TransitionMoments.NetPoachTopva(1) = SsMoments.NetPoachTopva(end) ;
TransitionMoments.NetPoachLowva(1) = SsMoments.NetPoachLowva(end) ;
TransitionMoments.NetPoach1000(1) = SsMoments.NetPoach1000(end) ;
TransitionMoments.NetPoach00(1) = SsMoments.NetPoach00(end) ;
end