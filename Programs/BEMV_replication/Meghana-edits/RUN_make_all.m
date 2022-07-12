%% Replication: Firm and Worker Dynamics in a Frictional Labor Market
%__________________________________________________________________________
% FILE: Running this file produces all MATLAB based figures and tables in 
% Bilal, Engbom, Mongey and Violante (2021)
% INPUTS: Contained in directories:
    % - Input_data
    % - Input_estimation
% OUTPUTS: Can start off with the following *empty* directories which will
% then be populated by the code:
    % - Created_mat_files
    % - Created_figure_files
    % - Created_table_files
% These are then loaded into the tex document: 
    % - ../1_Draft/
    
    % Housekeeping
clear
close all
clc

    % Choose which parts of replication to execute
execute.SS          = 1; % Solve Steady State
execute.CompStatics = 1; % Solve Comparative Statics
execute.Transition  = 1; % Solve Transition
execute.Figures     = 1; % Plot all Figures and Tables

%% A. SOLVE MODEL STEADY STATE
%__________________________________________________________________________
if execute.SS
disp('Solving model steady state...')
tic

% Parameters
run ExogenousParameters.m % sets values for ExogParams struct
load Input_mat_files/EstimatedParams.mat % load estimated parameters, Params struct
Params.zbar = 0; % normalization

% Options
options.Hopenhayn           = 0;
options.HopenhaynNumerical  = 0;
options.Transition          = 0;

% Solve Steady State
[~,~,~, SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat, ...
    g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, NumGrids, Derivatives]...
    = SolveBEMV( Params, ExogParams, options ) ;

% Compute Moments
disp('...solved steady state, now computing moments...')
[ Moments, g_znamat, gdzdn_znmat, gndzdn_znmat, ~,~, T12, ~, sep_zn,...
    ~, TotalNetPoach_zn, gdzdn_znamat, gndzdn_znamat, ~ ] = ...
    ComputeMoments( SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat, g_znmat,...
    Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, Params, ExogParams, Derivatives, options ) ;

% Age and size settings for more moments
Na              = size(g_znamat,2);
AgeGridLevels   = [ 0 , 2 , 4 , 6 , 11 ] + 1 ;
agenames        = {'0-1','2-3','4-5','6-10','11+'};
xA              = 1:1:5;
SizeGridLevels  = [ 0 , 20 , 50 , 250 , 500 ] ;
sizenames       = {'0-19','20-49','50-249','250-499','500+'};
xS              = 1:1:5;

% Compute life cycle and size moments
[ JCRateSize_model, JDRateSize_model, HiresRateSize_model, SepsRateSize_model, ~,...
    ExitRateUWSize_model, FirmShareSize_model, EmpShareSize_model, JCRateAge_model, JDRateAge_model,...
    HiresRateAge_model, SepsRateAge_model, ~, ExitRateUWAge_model, FirmShareAge_model, EmpShareAge_model ] = ...
    ComputeMomentsLifeCycleSize( NumGrids, T12, q, v_znmat, phi, Gn_znmat, gdzdn_znmat,...
    gndzdn_znmat, gdzdn_znamat, gndzdn_znamat, Na, AgeGridLevels, SizeGridLevels ) ;

    % Save results
save Created_mat_files/SS_results.mat


disp('... done!')
toc
end

%{
%% B. SOLVE COMPARATIVE STATICS FOR FIGURE 11
%__________________________________________________________________________
if execute.CompStatics
disp('Solving model comparative statics...')
tic

    % Parameters
run ExogenousParameters.m % sets values for ExogParams struct
load Input_mat_files/EstimatedParams.mat % load estimated parameters, Params struct 

    % Options
options.Hopenhayn           = 0;
options.HopenhaynNumerical  = 1;
options.Transition          = 0;

    % Numerical objects
Numerical   = NumericalApproximations( Params, ExogParams, options );           % get settings for numerical approximations     -> (zmin,zmax,nmin,nmax)
NumGrids    = Grids( Numerical, Params ) ;                                      % define grids and useful objects from them     -> (Nz,Nn,z,n,dz_z,dn_n,dzdn_zn,dzdn_znmat,...)

    % Grid for matching efficiency
grid        = 1:0.05:5 ; %1:0.05:5
griddown    = fliplr(0.2:0.05:1) ; %0.2:0.05:1
ngrid       = length(grid) ;
ngriddown   = length(griddown) ;
EstimateA   = Params.A ;

    % Load useful objects
Nz          = NumGrids.Nz;        % Number of points in z grid
Nn          = NumGrids.Nn;        % Number of points in n grid
n_znmat     = NumGrids.n_znmat;   % Replica of the n grid in logs
z_znmat     = NumGrids.z_znmat;   % Replica of the z grid in logs
dzdn_znmat  = NumGrids.dzdn_znmat;  % Measure (dz*dn as matrix, where d is log difference)
in0         = NumGrids.in0;         % index of n gridpoint closest to n_0
pi0_z       = NumGrids.pi0_z;       % Distribution of entering firms in z grid

    % Loop for each direction
for direction = { 'up' , 'down' } % 'up' and/or 'down'

        % Setup
    dir = direction{1} ;
    if strcmp(dir,'down') 
        grid  = griddown ;
        ngrid = ngriddown ;
    end

        % Initialize
    MOM             = zeros( ngrid , 6 ) ;
    MOMex           = zeros( ngrid , 6 ) ;
    Distrib         = zeros( Nz * Nz , ngrid ) ;
    Surplus         = zeros( Nz * Nz , ngrid ) ;
    MarginalSurplus = zeros( Nz * Nz , ngrid ) ;

        % Loop on each gridpoint conditional on direction
    for i=1:ngrid
        Params.A = grid(i) * EstimateA ; % set matching efficiency
        
        if i == 1  % Original steady state
            options.Hopenhayn = 0 ;
            run ExogenousParameters.m % sets values for ExogParams struct
            [ ~,~,~, SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat,...
                v_znmat, g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, NumGrids, Dmatrices, Numerical]...
                = SolveBEMV( Params, ExogParams, options );
            pi0_znmat        = zeros(Nz,Nn) ; 
            pi0_znmat(:,in0) = pi0_z ;
            entrycost        = S_znmat(:)' * ( pi0_znmat(:) .* dzdn_znmat(:) ) ; % get value of entry
            run AssignHopenhaynInitialGuess.m
            options.Hopenhayn = 1 ;

        else    % Comparative statics
            [ SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat, ...
                g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, ExogParams ] ...
                = SolveHopenhayn( Params, entrycost, ExogParams, options, grid(i), dir ) ;
            run AssignHopenhaynInitialGuess.m

        end

            % Compute moments function
        options.Hopenhayn = 0 ;
        [ Moments, ~, H, Hn, ~,~,~,~,~,~,~,~,~ ] = ...
            ComputeMoments( SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat,...
            v_znmat, g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u,  L, pi0Trunc_zn,...
            Params, ExogParams, Dmatrices, options ) ;
        options.Hopenhayn = 1 ;
        MOM( i , 1:6 ) = [ Moments.URate, Moments.StdLogVAPWfirm, Moments.CorrLogZLogSize, Moments.MeanY, Moments.VacCostToY, ExogParams.MF ] ;
        
            % Compute additional moments
        Distrib(:,i)        = g_znmat(:) ;
        Surplus(:,i)        = S_znmat(:) ;
        MarginalSurplus(:,i)= Sn_znmat(:) ;
        SN                  = Sn_znmat ./ exp(n_znmat) ;
        arg                 = SN > 0.00001 ;
        MeanLogSn           = sum(sum( log(SN(arg(:))) .* H(arg(:)) ) ) / sum(H(arg(:))) ;
        StdLogSn            = sqrt( sum(sum( ( log(SN(arg(:))) - MeanLogSn ).^2 .* H(arg(:)) )) / sum(H(arg(:))) ) ;
        MeanLogSnW          = sum(sum( log(SN(arg(:))) .* Hn(arg(:)) ) ) / sum(Hn(arg(:))) ;
        StdLogSnW           = sqrt( sum(sum( ( log(SN(arg(:))) - MeanLogSnW ).^2 .* Hn(arg(:)) )) / sum(Hn(arg(:))) ) ;
        MeanLogSize         = sum(sum( n_znmat(:) .* Hn(:) ) ) ;
        MeanLogProd         = sum(sum( z_znmat(:) .* Hn(:) ) ) ;
        StdLogSize          = sqrt( sum(sum( ( n_znmat(:) - MeanLogSize ).^2 .* Hn(:) )) ) ;
        StdLogProd          = sqrt( sum(sum( ( z_znmat(:) - MeanLogProd ).^2 .* Hn(:) )) ) ;
        CorrLogZLogSize     = sum(sum( ( ( z_znmat(:) - MeanLogProd ) .* ( n_znmat(:) - MeanLogSize ) ) .* Hn(:) ) ) / ( StdLogSize * StdLogProd ) ;
        vapw                = z_znmat - (1-Params.alpha) * n_znmat ;
        StdLogVAPW          = sqrt( sum(sum( vapw(:).^2 .* Hn(:) ) ) - sum(sum( vapw(:) .* Hn(:) ))^2 ) ;
        
            % Save moments
        MOM( i , 7  ) = StdLogVAPW ;
        MOM( i , 8  ) = CorrLogZLogSize ;
        MOM( i , 9  ) = MeanLogSize ;
        MOM( i , 10 ) = StdLogSize ;
        MOM( i , 11 ) = StdLogSn ;
        MOM( i , 12 ) = StdLogSnW ; 

            % Print results for this gridpoint
        fprintf('\n')
        fprintf('\n')
        fprintf('Unemployment rate = %6.4f \n',MOM(i,1))
        fprintf('St.d. of MPL   = %6.4f \n',MOM(i,2))
        fprintf('St.d. of Sn    = %6.4f \n',MOM(i,11))
        fprintf('Corr size-TFP  = %6.4f \n',MOM(i,3))
        fprintf('Gross output   = %6.4f \n',MOM(i,4))
        fprintf('Vacancy cost   = %6.4f \n',MOM(i,5))
        fprintf('\n')
        fprintf('\n')
    end

        % Loop for additional moments
    for i=1:ngrid
        g_zn            = Distrib(:,i) ;
        g_znmat         = reshape(g_zn,Nz,Nn) ;
        select          = 1:(Nn-1) ;
        Hs              = g_znmat .* dzdn_znmat ;
        Hs              = Hs(:,select) / sum(sum(Hs(:,select))) ;
        MeanLogSize     = sum(sum( n_znmat(:,select) .* Hs ) ) ;
        MeanLogProd     = sum(sum( z_znmat(:,select) .* Hs ) ) ;
        StdLogSize      = sqrt( sum(sum( ( n_znmat(:,select) - MeanLogSize ).^2 .* Hs )) ) ;
        StdLogProd      = sqrt( sum(sum( ( z_znmat(:,select) - MeanLogProd ).^2 .* Hs )) ) ;
        CorrLogZLogSize = sum(sum( ( ( z_znmat(:,select) - MeanLogProd ) .* ( n_znmat(:,select) - MeanLogSize ) ) .* Hs ) ) / ( StdLogSize * StdLogProd ) ;
        vapw            = z_znmat - (1-Params.alpha) * n_znmat ;
        StdLogVAPW      = sqrt( sum(sum( vapw(:,select).^2 .* Hs ) ) - sum(sum( vapw(:,select) .* Hs ))^2 ) ;
        MOMex( i , 1 )  = MOM(i,1) ;
        MOMex( i , 2 )  = StdLogVAPW ;
        MOMex( i , 3 )  = CorrLogZLogSize ;
        MOMex( i , 4 )  = MOM(i,4) ;
        MOMex( i , 5 )  = MOM(i,5) ;
        MOMex( i , 6 )  = MOM(i,6) ;
    end

        % Save results for this direction
     save(['Created_mat_files/Hopenhayn_' dir '.mat'],'grid','MOM','MOMex','Distrib','Surplus','MarginalSurplus');
end

    % Merge both directions
load Created_mat_files/Hopenhayn_down
griddown            = fliplr(grid) ;
ng                  = length(griddown) ;
MOMdown             = flipud(MOM(1:ng,:)) ; 
MOMexdown           = flipud(MOMex(1:ng,:)) ;
Distribdown         = flipud(Distrib(:,1:ng)) ;
Surplusdown         = flipud(Surplus(:,1:ng)) ;
MarginalSurplusdown = flipud(MarginalSurplus(:,1:ng)) ;
load Created_mat_files/Hopenhayn_up
gridA               = [ griddown(1:end-1) grid ] ;
MOM                 = [ MOMdown(1:end-1,:) ; MOM ] ;
MOMex               = [ MOMexdown(1:end-1,:) ; MOMex ] ;
Distrib             = [ Distribdown(:,1:end-1) , Distrib ] ;
Surplus             = [ Surplusdown(:,1:end-1) , Surplus ] ;
MarginalSurplus     = [ MarginalSurplusdown(:,1:end-1) , MarginalSurplus ] ;

    % Last Moments
nn = size(MOM,1);
for i = 1:nn
    h               = Distrib(:,i) ;
    h               = reshape(h,Nz,Nn) ;
    Sn              = MarginalSurplus(:,i) ;
    Sn              = reshape(Sn, Nz,Nn) ;
    Hn              = h .* exp(n_znmat) ;
    Hn              = Hn / sum(Hn(:)) ;
    SN              = Sn ./ exp(n_znmat) ;
    arg             = SN > 0.00001 ;
    select          = 1:(Nn-1) ;
    Hs              = h .* dzdn_znmat ;
    Hs              = Hs(:,select) ; 
    args            = SN(:,select) > 0.00001 ;
    Hs              = Hs(args) / sum(Hs(args(:))) ;
    lSNs            = log(SN(:,select)) ; 
    lSNs            = lSNs(args) ;
    n               = exp(n_znmat(:,select)) ;
    Hns             = Hs .* n(args) / sum(Hs .* n(args)) ;
    MeanLogSnUWex   = sum(sum( lSNs .* Hs ) );
    StdLogSnUWex    = sqrt( sum(sum( ( lSNs - MeanLogSnUWex ).^2 .* Hs )  ) );
    MeanLogSnWex    = sum(sum( lSNs .* Hns ) );
    StdLogSnWex     = sqrt( sum(sum( ( lSNs - MeanLogSnWex ).^2 .* Hns ) ) ) ;
    MOMex( i , 7 )  = StdLogSnUWex ; 
    MOMex( i , 8 )  = StdLogSnWex ; 
end

    % Save results
save('Created_mat_files/Hopenhayn','gridA','MOM','MOMex','Distrib');


disp('... done!')
toc
end
%% C. SOLVE TRANSITION DYNAMICS FOR FIGURE 13
%__________________________________________________________________________
if execute.Transition
disp('Solving model transition dynamics...')
tic

    % Parameters
run ExogenousParameters.m % sets values for ExogParams struct
load Input_mat_files/EstimatedParams.mat % load estimated parameters, Params struct
Params.zbar = 0; % normalization

    % Options
options.Hopenhayn           = 0;
options.HopenhaynNumerical  = 0;
options.Transition          = 1;

    % New parameters
ExogParams.EntryElast = 0.05 ;
ExogParams.M0 = 1 ; % normalization of mass of potential entrants

    % Solve Steady State
[~,~,~, SUntrunc_znmat, SnUntrunc_znmat, S_znmat, Sn_znmat, v_znmat, ...
    g_znmat, Gn_znmat, Gv_znmat, q, phi, p, u, L, pi0Trunc_zn, NumGrids, Derivatives, Numerical]...
    = SolveBEMV( Params, ExogParams, options ) ;
ExitRate = sum(-L) * ( g_znmat(:) .* NumGrids.dzdn_znmat(:) ) ;
disp('...solved initial steady state, now computing transition equilibrium without shocks...')

    % Numerical settings for transition
TransitionNum.TFinal    = 12 * 50 ; % number of months
TransitionNum.T         =  (0:(Numerical.Delta):TransitionNum.TFinal)' ;
TransitionNum.dT        = [ TransitionNum.T(2:end) - TransitionNum.T(1:end-1) ; ...
                            TransitionNum.T(end) - TransitionNum.T(end-1) ] ;
TransitionNum.Nt        = length(TransitionNum.T) ;
Nt                      = TransitionNum.Nt;
SmoothPath.s            = 0.5 ;
SmoothPath.ME           = 1 ;

    % Initial guess for transition equilibrium:
TransitionEq.S          = repmat( S_znmat(:) , [1,Nt] ) ;
TransitionEq.Sn         = repmat( Sn_znmat(:) , [1,Nt] ) ;
TransitionEq.SUntrunc   = repmat( SUntrunc_znmat(:) , [1,Nt] ) ;
TransitionEq.SnUntrunc  = repmat( SnUntrunc_znmat(:) , [1,Nt] ) ;
TransitionEq.v          = repmat( v_znmat(:) , [1,Nt] ) ;
TransitionEq.u          = u * ones(Nt,1) ;
TransitionEq.q          = q * ones(Nt,1) ;
TransitionEq.phi        = phi * ones(Nt,1) ;
TransitionEq.p          = p * ones(Nt,1) ;
TransitionEq.Gv         = repmat( Gv_znmat(:) , [1,Nt] ) ;
TransitionEq.Gn         = repmat( Gn_znmat(:) , [1,Nt] ) ;
TransitionEq.g          = repmat( g_znmat(:) , [1,Nt] ) ;
TransitionEq.L          = repmat( L(:) , [1,Nt] ) ;
TransitionEq.MF         = ExogParams.MF * ones(Nt,1) ;
TransitionEq.ME         = ExitRate * ExogParams.MF * ones(Nt,1) ;
TransitionEq.theta      = (TransitionEq.p / Params.A).^(1/1-ExogParams.beta) ;

    % Get entry cost
pi0_z                   = zeros( NumGrids.Nz , NumGrids.Nn ) ; 
pi0_z(:,NumGrids.in0)   = NumGrids.pi0_z ;
pi0_z                   = pi0_z / sum(pi0_z(:) .* NumGrids.dzdn_znmat(:)) ;
ME0                     = ExitRate * ExogParams.MF ;
S0                      = sum( S_znmat(:) .* pi0_z(:) .* NumGrids.dzdn_znmat(:) ) ;
x0                      = ME0 / ExogParams.M0 ;
Params.entrycost        = ( (1-x0) / x0 ) ^ ( 1 / ExogParams.EntryElast ) * S0 ; % invert logit share

    % Initialize transition equilibrium: solve without shock
TransitionShock.rho     = ExogParams.rho * ones(Nt,1) ;
[~, TransitionEq ,~]    = SolveBEMVTransition( Params, ExogParams, options, ...
   TransitionNum, TransitionEq, TransitionShock, SmoothPath ) ;
disp('...solved transition equilibrium without shocks, now with shocks...')

    % Extract steady state
i       = floor( 0.6 * Nt ) ;
names   = fieldnames(TransitionEq) ;
for j=1:numel(names)
    name = names{j} ;
    x    = TransitionEq.(name) ;
    if size(x,2) > 1
        TransitionEq.(name) = repmat( x(:,i) , [1,Nt] ) ;
    else
        TransitionEq.(name) = repmat( x(i) , [1,Nt] )' ;
    end
end

    % Solve transition equilibrium with shock
SsEq                = TransitionEq ;
MultiplierRho       = 150 ;
HalfLife            = 2 * 12 ;
TransitionShock.rho = TransitionShock.rho ...
    + MultiplierRho * exp( -TransitionNum.T / HalfLife ) .* TransitionShock.rho ;
[~, TransitionEq ,~] = SolveBEMVTransition( Params, ExogParams, options, ...
   TransitionNum, TransitionEq, TransitionShock, SmoothPath ) ;
disp('...solved transition equilibrium with shocks, now computing transition moments...')

    % Compute moments in the transition
TransitionMoments   = ComputeMomentsTransition( TransitionNum, TransitionEq,...
                        TransitionShock, SsEq, 1, Nt, ExogParams, Params, NumGrids, options ) ;
                    
    % Save results
save Created_mat_files/Transition.mat TransitionMoments TransitionNum


disp('... done!')
toc
end
%% D. FIGURES AND TABLES
%__________________________________________________________________________
if execute.Figures 
disp('Producing figures and tables...')
    
    % Load results
load Created_mat_files/SS_results

    % Select which figures and tables to produce
produce.table1  = 1;
produce.fig6    = 1;
produce.fig7    = 1;
produce.fig8    = 1;
produce.fig9    = 1;
produce.fig10   = 1;
produce.fig11   = 1;
produce.fig12   = 1;
produce.fig13   = 1;
produce.fig14   = 1;
produce.figE1   = 0; % Needs to load a very large mat file, available upon request to the authors
produce.figE2   = 1;

    % Produce figures and tables
run FiguresTables


disp('... done!')
end
%}