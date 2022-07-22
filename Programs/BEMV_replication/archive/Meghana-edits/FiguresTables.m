%% ========================================================================
% 
%                       PRODUCES FIGURES AND TABLES
% 
% =========================================================================

%%  SETUP
%__________________________________________________________________________

    % Parameter names
ParameterName.mu        = '$\mu$';
ParameterName.b         = '$b$';
ParameterName.sigma     = '$\sigma$';
ParameterName.xi        = '$\xi$';
ParameterName.delta     = '$\delta$';
ParameterName.alpha     = '$\alpha$';
ParameterName.A         = '$A$';
ParameterName.zeta      = '$\zeta$';
ParameterName.n_0       = '$n_0$';
ParameterName.c_f       = '$c_f$';
ParameterTitle.mu       = 'Drift of productivity' ;
ParameterTitle.b        = 'Flow value of leisure';
ParameterTitle.sigma    = 'St.d of productivity shocks';
ParameterTitle.xi       = 'Relative search efficiency of employed';
ParameterTitle.delta    = 'Exogenous separation rate';
ParameterTitle.alpha    = 'Curvature of production';
ParameterTitle.A        = 'Matching efficiency';
ParameterTitle.zeta     = 'Shape of entry distribution';
ParameterTitle.n_0      = 'Size of entrants';
ParameterTitle.c_f      = 'Fixed cost of operation';

	% Target names for calibration
ChoiceOfTargets.mu      = 'BDSExitRateUW';
ChoiceOfTargets.b		= 'BDSJobDestructionIncumbents';
ChoiceOfTargets.sigma	= 'BDSStdGrowthRate';
ChoiceOfTargets.xi		= 'EEoverEquarterly';
ChoiceOfTargets.delta	= 'EUoverEquarterly';
ChoiceOfTargets.alpha	= 'EmpShare500';
ChoiceOfTargets.A		= 'URate';
ChoiceOfTargets.zeta	= 'BDSJCRateAge1';

	% Import targets for calibration
run Input_data/Targets_for_calibration.m;

    % BDS Data 
BDSsize = xlsread('Input_data/BDS_J2J_by_size.xlsx') ;
BDSage  = xlsread('Input_data/BDS_J2J_by_age.xlsx') ;

%%  TABLE 1: ESTIMATED PARAMETERS AND TARGETED MOMENTS
%__________________________________________________________________________
if produce.table1

fid = fopen('Created_table_files/Table1_ParametersMoments.tex','w');
fprintf(fid,'\\begin{tabular}{l l c l c c}\n');
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'Parameter &  & Value & Moment & Data & Model \\\\ \n');
fprintf(fid,'\\hline \\\\[-4mm]');
fprintf(fid,'\\multicolumn{6}{c}{A. Externally set/normalized parameters} \\\\ \n');
fprintf(fid,'\\hline \\\\[-4mm]');
fprintf(fid,'%s & %s & %3.3f & %s \\\\ \n','$\rho$',                    'Discount rate',                            ExogParams.rho, '5\% annual real interest rate') ;
fprintf(fid,'%s & %s & %3.1f & %s \\\\ \n','$\beta$',                   'Elasticity of matches w.r.t. vacancies',   ExogParams.beta,'Petrongolo and Pissarides (2001)') ;
fprintf(fid,'%s & %s & %3.0f & %s \\\\ \n','$c_f$',                     'Fixed cost of operation',                  Params.c_f,     'Normalization') ;
fprintf(fid,'%s & %s & %3.0f & %s \\\\ \n','$\overline{c}/(1+\gamma)$', 'Scalar in the cost of vacancies',          ExogParams.cbar,'Normalization') ;
fprintf(fid,'%s & %s & %3.0f & %s \\\\ \n','$n_0$',                     'Size of entrants',                         Params.n_0,     'Normalization') ;
fprintf(fid,'\\hline \\\\[-4mm]');
fprintf(fid,'\\multicolumn{6}{c}{B. Estimated offline} \\\\ \n');
fprintf(fid,'\\hline \\\\[-4mm]');
fprintf(fid,'%s & %s & %3.3f & %s & %2.3f & %2.3f \\\\ \n','$\mathtt{m}$',  'Number of active firms',               ExogParams.MF,      'Average firm size (BDS)',              23.340, (1-Moments.URate)/ExogParams.MF) ;
fprintf(fid,'%s & %s & %3.3f & %s & %2.3f & %2.3f \\\\ \n','$d$',           'Exogenous exit rate',                  ExogParams.d,       'Exit rate, 1000--2499 empl. firms',    0.002,  12*ExogParams.d) ;
fprintf(fid,'%s & %s & %3.3f & %s & %2.3f & %2.3f \\\\ \n','$\gamma$',      'Curvature of vacancy cost function',   ExogParams.gamma,   'Vacancy filling rate vs. hiring rate', 3.45,   ExogParams.gamma) ;
fprintf(fid,'\\hline \\\\[-4mm]');
fprintf(fid,'\\multicolumn{6}{c}{C. Internally Estimated} \\\\ \n');
fprintf(fid,'\\hline \\\\[-4mm]');
ParamsStr = { 'mu' , 'sigma' , 'alpha' , 'zeta' , 'A' , 'xi' , 'delta' , 'b' } ;
ChoiceOfTargetsTitle.mu     = 'Exit rate (annual)' ;
ChoiceOfTargetsTitle.sigma  = 'St.d. of log empl. growth (annual)' ;
ChoiceOfTargetsTitle.alpha  = 'Empl. share of 500+ firms' ;
ChoiceOfTargetsTitle.zeta   = 'JC rate, age 1 firms (annual)' ;
ChoiceOfTargetsTitle.A      = 'Nonemployment rate' ;
ChoiceOfTargetsTitle.xi     = 'EE rate (quarterly)' ;
ChoiceOfTargetsTitle.delta  = 'EN rate (quarterly)' ;
ChoiceOfTargetsTitle.b      = 'JD rate of incumbents (annual)' ;
for params = ParamsStr
    fprintf(fid,'%s & %s & %3.3f & %s & %2.3f & %2.3f \\\\ \n',ParameterName.(params{1}), ...
                                                               ParameterTitle.(params{1}), ...
                                                               Params.(params{1}), ...
                                                               ChoiceOfTargetsTitle.(params{1}), ...
                                                               Targets.(ChoiceOfTargets.(params{1})),...
                                                               Moments.(ChoiceOfTargets.(params{1}))) ;
end
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);

end
%%  FIGURE 6: VACANCY RATES BY FIRM GROSS HIRING RATE
%__________________________________________________________________________
if produce.fig6

    % Setup
hires       = q * v_znmat .* ( phi + (1-phi)*Gn_znmat ) ;
seps        = reshape(sep_zn,NumGrids.Nz,NumGrids.Nn) .* exp(NumGrids.n_znmat) ;
GrowthRate  = ( hires - seps ) ./ exp(NumGrids.n_znmat) ;
GrowthBins  = [ -0.30:0.01:-0.01 0 0.01:0.01:0.31 ]' ;  
xxlim       = 1.5 ;

    % Model
VacByGrowth_model   = AggregateBin( v_znmat ,...
                                  GrowthRate,...
                                  ones(size(v_znmat)),...
                                  GrowthBins,...
                                  gdzdn_znmat,...
                                  'Levels','mean') ;
VacByGrowth_model   = VacByGrowth_model(1:end-1) ; % to remove all firms that grow more than 30%
SizeByGrowth_model  = AggregateBin( exp(NumGrids.n_znmat) ,...
                                  GrowthRate,...
                                  ones(size(NumGrids.n_znmat)),...
                                  GrowthBins,...
                                  gdzdn_znmat,...
                                  'Levels','mean') ;  
SizeByGrowth_model  = SizeByGrowth_model(1:end-1) ;
HiresByGrowth_model = AggregateBin( hires ,...
                                  GrowthRate,...
                                  ones(size(hires)),...
                                  GrowthBins,...
                                  gdzdn_znmat,...
                                  'Levels','mean') ;                                     
HiresByGrowth_model = HiresByGrowth_model(1:end-1) ;
HireRate_model      = HiresByGrowth_model ./ SizeByGrowth_model ;
VacRate_model       = VacByGrowth_model   ./ ( VacByGrowth_model + SizeByGrowth_model )  ;                     
FillingRate_model   = HiresByGrowth_model ./ VacByGrowth_model ;

    % Find normalization point
[~,i0]              = min( GrowthBins(1:end-1).^2 + 10000 * isnan( HiresByGrowth_model ) ) ;
zzz_model           = (GrowthBins(1:end-1)>=0.025) & (GrowthBins(1:end-1)<=0.095) ;

    % Data
data                = csvread('Input_data/Mongey_Violante_BLS_JOLTS_Microdata_Output.csv',1,0);
HireRate_data       = data(:,2);
VacRate_data        = data(:,3);
FillingRate_data    = data(:,4);
zzz_data            = (HireRate_data(1:end-1)>=2.5) & (HireRate_data(1:end-1)<=9.5) ;

    % Rescale both model and data to point that has zero growth
for str1 = { 'HireRate' , 'VacRate' , 'FillingRate' } 
    for str2 = { 'model' , 'data' } 
        eval([ 'dlog' cell2mat(str1) '_' cell2mat(str2) ' = log(' cell2mat(str1) '_' cell2mat(str2) '/'  cell2mat(str1) '_' cell2mat(str2) '(i0)) ;'])
        
            % Regressions with all observations
        eval([ 'beta_' cell2mat(str1) '_' cell2mat(str2) '= regress(dlog' cell2mat(str1) '_' cell2mat(str2) ',dlogHireRate_' cell2mat(str2) ') ;'])
        eval([ 'fit_'  cell2mat(str1) '_' cell2mat(str2) '= beta_' cell2mat(str1) '_' cell2mat(str2) ' * dlogHireRate_' cell2mat(str2) ' ;'])
    
            % Regressions with selected observations
        eval([ 'y = dlog' cell2mat(str1) '_' cell2mat(str2) ';' ])
        eval([ 'x = dlogHireRate_' cell2mat(str2) ';' ])
        y = y(x >= 0 & x <= xxlim ) ;
        x = x(x >= 0 & x <= xxlim ) ;
        eval([ 'select_beta_' cell2mat(str1) '_' cell2mat(str2) '= regress(y,x) ;']);
        eval([ 'select_fit_'  cell2mat(str1) '_' cell2mat(str2) '= select_beta_' cell2mat(str1) '_' cell2mat(str2) ' * dlogHireRate_' cell2mat(str2) ' ;']);
    end 
end

    % Figure setup
Figure6 = figure(6);
set(Figure6,'Pos',[423.6667 443 1.0887e+03 415.3333]);
font = 16;
set(0,'defaultlinemarkersize',10);
x = linspace(0,xxlim,10) ; % for evenly spaced fit line

    % Vacancy Rate plot
subplot(1,2,1);
plot(dlogHireRate_model,dlogVacRate_model,'bo','linewidth',1); grid on; hold on;
plot(dlogHireRate_data ,dlogVacRate_data ,'rx','linewidth',1); grid on; hold on;
plot(dlogHireRate_model,select_fit_VacRate_model,'b-','linewidth',2); hold on;
plot(x ,select_beta_VacRate_data * x ,'r--','linewidth',2); hold on;
plot([0,100],[0,100],'k-','linewidth',1);
title('A. Vacancy rate: $v/(v+n)$','interpreter','latex'); 
xlabel('Relative log gross hiring rate','interpreter','latex');
ylabel('Relative log vacancy rate','interpreter','latex');
l = legend(sprintf('Model ($\\beta$ = %3.2f)',beta_VacRate_model),...
           sprintf('Data $\\:\\:\\,$($\\beta$ = %3.2f)',beta_VacRate_data));
set(l,'interpreter','latex','location','NorthWest');
set(gca,'fontsize',font,'ticklabelinterpreter','latex');
set(gca,'YTick',[0,1,2,3]);
ylim([0,xxlim]);
xlim([0,xxlim]);

    % Vacancy Yield plot
subplot(1,2,2);
plot(dlogHireRate_model,dlogFillingRate_model,'bo','linewidth',1); grid on; hold on;
plot(dlogHireRate_data,dlogFillingRate_data,'rx','linewidth',1); grid on; hold on;
plot(dlogHireRate_model,select_fit_FillingRate_model,'b-','linewidth',2);hold on;
plot(x ,select_beta_FillingRate_data * x,'r--','linewidth',2);hold on;
plot([0,100],[0,100],'k-','linewidth',1);
title('B. Vacancy yield: $h/v$','interpreter','latex');
xlabel('Relative log gross hiring rate','interpreter','latex');
ylabel('Relative log vacancy yield','interpreter','latex');
l = legend(sprintf('Model ($\\beta$ = %3.2f)',select_beta_FillingRate_model),...
           sprintf('Data $\\:\\:\\,$($\\beta$ = %3.2f)',select_beta_FillingRate_data));
set(l,'interpreter','latex','location','NorthWest');
set(gca,'fontsize',font,'ticklabelinterpreter','latex');
set(gca,'YTick',[0,1,2,3]);
ylim([0,xxlim]);
xlim([0,xxlim]); 

    % Save
print '-depsc' 'Created_figure_files/Fig6_VacancyRatesYields_eps'
print '-dpng'  'Created_figure_files/Fig6_VacancyRatesYields_png'

end
%%  FIGURE 7: DISTRIBUTION OF FIRMS AND EMPLOYMENT BY FIRM AGE AND SIZE IN DATA AND MODEL
%__________________________________________________________________________
if produce.fig7
    
    % Rename data variables
FdataA  = BDSage(:,2);
EdataA  = BDSage(:,4);
FdataS  = BDSsize(:,1);
EdataS  = BDSsize(:,3);

    % Rename model variables
FmodelA = FirmShareAge_model;
EmodelA = EmpShareAge_model;
FmodelS = FirmShareSize_model;
EmodelS = EmpShareSize_model/sum(EmpShareSize_model);

    % Percent
FmodelA = 100*FmodelA;
FdataA  = 100*FdataA;
EmodelA = 100*EmodelA;
EdataA  = 100*EdataA;
FmodelS = 100*FmodelS;
FdataS  = 100*FdataS;
EmodelS = 100*EmodelS;
EdataS  = 100*EdataS;

    % Figure setup
Figure7 = figure(7);
set(Figure7,'Pos',[786 918 1214 420]);
set(0,'defaultlinelinewidth',2);
font = 14;

    % Distribution of firms by size
subplot(1,2,1);
plot(xA,FmodelS,'o-','color',rgb('Red'),'markersize',10);hold on; grid on;
plot(xA,EmodelS,'o-','color',rgb('Blue'),'markersize',10);hold on; grid on;
plot(xA,FdataS,'x:','color',rgb('Red'),'markersize',12);hold on; grid on;
plot(xA,EdataS,'x:','color',rgb('Blue'),'markersize',12);hold on; grid on;
title('A. By firm size','interpreter','latex');
xlabel('Size group','interpreter','latex');
ylabel('Share (percent)','interpreter','latex');
set(gca,'Xtick',xS,'Xticklabel',sizenames);xlim([0.5,5.5]);
l = legend('Share of firms (red)','Share of employment (blue)');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.0f');
ylev = 55;
text(2,ylev,'Model (Circles, solid)','interpreter','latex','fontsize',font);
text(2,ylev-10,'Data (Crosses, dashed)','interpreter','latex','fontsize',font);

    % Distribution of firms by age
subplot(1,2,2);
plot(xA,FmodelA,'o-','color',rgb('Red'),'markersize',10);hold on; grid on;
plot(xA,EmodelA,'o-','color',rgb('Blue'),'markersize',10);hold on; grid on;
plot(xA,FdataA,'x:','color',rgb('Red'),'markersize',12);hold on; grid on;
plot(xA,EdataA,'x:','color',rgb('Blue'),'markersize',12);hold on; grid on;
title('B. By firm age','interpreter','latex');
xlabel('Age group','interpreter','latex');
ylabel('Share (percent)','interpreter','latex');
set(gca,'Xtick',xA,'Xticklabel',agenames);xlim([0.5,5.5]);
l = legend('Share of firms (red)','Share of employment (blue)','location','northwest');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.0f');

    % Save
print '-depsc' 'Created_figure_files/Fig7_DistributionSizeAge_eps'
print '-dpng'  'Created_figure_files/Fig7_DistributionSizeAge_png'

end
%%  FIGURE 8: JOB, WORKER AND FIRM REALLOCATION BY SIZE AND AGE
%__________________________________________________________________________
if produce.fig8

    % Rename data variables
JCdataA  = BDSage(:,5);
JDdataA  = BDSage(:,6);
HdataA   = BDSage(:,7);
SdataA   = BDSage(:,8);
ERdataA  = BDSage(:,13);
JCdataS  = BDSsize(:,5);
JDdataS  = BDSsize(:,6);
HdataS   = BDSsize(:,7);
SdataS   = BDSsize(:,8);
ERdataS  = BDSsize(:,13);

    % Rename model variables
JCmodelA = JCRateAge_model;
JDmodelA = JDRateAge_model;
HmodelA  = HiresRateAge_model;
SmodelA  = SepsRateAge_model;
ERmodelA = ExitRateUWAge_model ;
JCmodelS = JCRateSize_model;
JDmodelS = JDRateSize_model;
HmodelS  = HiresRateSize_model;
SmodelS  = SepsRateSize_model;
ERmodelS = ExitRateUWSize_model ;

    % Percent
JCmodelA    = 100*JCmodelA;
JCdataA     = 100*JCdataA;
JDmodelA    = 100*JDmodelA;
JDdataA     = 100*JDdataA;
HmodelA     = 100*HmodelA;
HdataA      = 100*HdataA;
SmodelA     = 100*SmodelA;
SdataA      = 100*SdataA;
ERmodelA    = 100*ERmodelA;
ERdataA     = 100*ERdataA;
JCmodelS    = 100*JCmodelS;
JCdataS     = 100*JCdataS;
JDmodelS    = 100*JDmodelS;
JDdataS     = 100*JDdataS;
HmodelS     = 100*HmodelS;
HdataS      = 100*HdataS;
SmodelS     = 100*SmodelS;
SdataS      = 100*SdataS;
ERmodelS    = 100*ERmodelS;
ERdataS     = 100*ERdataS;

    % Figure 8A setup
Figure8A = figure(8);
set(Figure8A,'Pos',[1866 828 1463 420]);
set(0,'defaultlinelinewidth',2);
font = 14;

    % Job reallocation by age
subplot(1,3,1);
plot(xA,JCmodelA,'o-','color',rgb('Green'));hold on; grid on;
plot(xA,JDmodelA,'o-','color',rgb('Red'));hold on; grid on;
plot(xA,JCdataA,'x-.','color',rgb('Green'),'markersize',10);hold on; grid on;
plot(xA,JDdataA,'x-.','color',rgb('Red'),'markersize',10);hold on; grid on;
title('(i) Job reallocation','interpreter','latex');
xlabel('Age group','interpreter','latex');
ylabel('Rate (percent)','interpreter','latex');
set(gca,'Xtick',xA,'Xticklabel',agenames);xlim([0.5,5.5]);
l = legend('Job creation rate (green)','Job destruction (red)');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.0f');
ylim([0,25]);
ylev = 13.5;
text(2,ylev,'Model (Circles, solid)','interpreter','latex','fontsize',font);
text(2,ylev-2,'Data (Crosses, dashed)','interpreter','latex','fontsize',font);

    % Worker reallocation by age
subplot(1,3,2);
plot(xA,HmodelA,'o-','color',rgb('Green'));hold on; grid on;
plot(xA,SmodelA,'o-','color',rgb('Red'));hold on; grid on;
plot(xA,HdataA,'x-.','color',rgb('Green'),'markersize',10);hold on; grid on;
plot(xA,SdataA,'x-.','color',rgb('Red'),'markersize',10);hold on; grid on;
title('(ii) Worker reallocation','interpreter','latex');
xlabel('Age group','interpreter','latex');
ylabel('Rate (percent)','interpreter','latex');
set(gca,'Xtick',xA,'Xticklabel',agenames);xlim([0.5,5.5]);
l = legend('Hire rate (green)','Separation rate (red)');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.0f');
ylim([0,25]);

    % Firm reallocation by age
subplot(1,3,3);
plot(xA,ERmodelA,'o-','color',rgb('Black'));hold on; grid on;
plot(xA,ERdataA,'x-.','color',rgb('Black'),'markersize',10);hold on; grid on;
title('(iii) Firm reallocation','interpreter','latex');
xlabel('Age group','interpreter','latex');
ylabel('Rate (percent)','interpreter','latex');
set(gca,'Xtick',xA,'Xticklabel',agenames);xlim([0.5,5.5]);
l = legend('Exit rate');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.1f');
ylim([0,3.5]);

    % Save
print '-depsc' 'Created_figure_files/Fig8A_ReallocationAge_eps'
print '-dpng'  'Created_figure_files/Fig8A_ReallocationAge_png'

    % Figure 8B setup
Figure8B = figure(88);
set(Figure8B,'Pos',[1867 315 1463 420]);
set(0,'defaultlinelinewidth',2);
font = 14;

    % Job reallocation by size
subplot(1,3,1);
plot(xS,JCmodelS,'o-','color',rgb('Green'));hold on; grid on;
plot(xS,JDmodelS,'o-','color',rgb('Red'));hold on; grid on;
plot(xS,JCdataS,'x-.','color',rgb('Green'),'markersize',10);hold on; grid on;
plot(xS,JDdataS,'x-.','color',rgb('Red'),'markersize',10);hold on; grid on;
title('(i) Job reallocation','interpreter','latex');
xlabel('Size group','interpreter','latex');
ylabel('Rate (percent)','interpreter','latex');
set(gca,'Xtick',xS,'Xticklabel',sizenames);xlim([0.5,5.5]);
l = legend('Job creation rate (green)','Job destruction (red)');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.0f');
ylim([0,25]);
ylev = 13.5;
text(2,ylev,'Model (Circles, solid)','interpreter','latex','fontsize',font);
text(2,ylev-2,'Data (Crosses, dashed)','interpreter','latex','fontsize',font);

    % Worker reallocation by size
subplot(1,3,2);
plot(xS,HmodelS,'o-','color',rgb('Green'));hold on; grid on;
plot(xS,SmodelS,'o-','color',rgb('Red'));hold on; grid on;
plot(xS,HdataS,'x-.','color',rgb('Green'),'markersize',10);hold on; grid on;
plot(xS,SdataS,'x-.','color',rgb('Red'),'markersize',10);hold on; grid on;
title('(ii) Worker reallocation','interpreter','latex');
xlabel('Size group','interpreter','latex');
ylabel('Rate (percent)','interpreter','latex');
set(gca,'Xtick',xS,'Xticklabel',sizenames);xlim([0.5,5.5]);
l = legend('Hire rate (green)','Separation rate (red)');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.0f');
ylim([0,25]);

    % Firm reallocation by size
subplot(1,3,3);
plot(xS,ERmodelS,'o-','color',rgb('Black'));hold on; grid on;
plot(xS,ERdataS,'x-.','color',rgb('Black'),'markersize',10);hold on; grid on;
title('(iii) Firm reallocation','interpreter','latex');
xlabel('Size group','interpreter','latex');
ylabel('Rate (percent)','interpreter','latex');
set(gca,'Xtick',xS,'Xticklabel',sizenames);xlim([0.5,5.5]);
l = legend('Exit rate');
set(l,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.1f');
ylim([0,3.5]);

    % Save
print '-depsc' 'Created_figure_files/Fig8B_ReallocationSize_eps'
print '-dpng'  'Created_figure_files/Fig8B_ReallocationSize_png'

end
%%  FIGURE 9: NET POACHING AND MARGINAL SURPLUS DISTRIBUTION
%__________________________________________________________________________
if produce.fig9

    % Define Sn bins by H quantiles
SnPerc              = [1,5,10,20:20:80,90,95,99]/100;
SnBins              = nan(length(SnPerc),1);
[SnSort,IndSort]    = sort( Sn_znmat(:) ./ exp(NumGrids.n_znmat(:)) ) ;
[~,IndUnsort]       = sort(IndSort ) ;
hSort               = g_znmat(IndSort) ;
hcdfSort            = cumsum( g_znmat(IndSort(:)).*NumGrids.dzdn_znmat(IndSort(:)) ) ;
for i=1:length(SnPerc)
    [~,j]       = min( (hcdfSort-SnPerc(i)).^2 );
    SnBins(i)   = SnSort(j);
end
SnBins      = [SnBins ; 100 ] ; 
LogSnBins   = log(SnBins) ;
LogSnBins   = [ - 100 ; LogSnBins ] ; % add a low point due to the structure of the AggregateBin function

    % Automatically adjust x-axis
xmin = max(floor(LogSnBins(2)),-1) ;
xmax = floor(LogSnBins(end-1)) + 1 ;
xvec = [xmin ; LogSnBins(2:(end-1)) ; xmax ] ;

    % Labels for x-axis
LogSnBinsGraph      = xmin:(0.5):xmax ;
LogSnBinsLabels     = {} ;
LogSnBinsLabels{1}  = [ '$\leq' num2str(xmin) '$'] ;
for i=2:(length(LogSnBinsGraph)-1)
    LogSnBinsLabels{i} = num2str(LogSnBinsGraph(i));   
end
LogSnBinsLabels{i+1} =[ '$\geq' num2str(xmax) '$'] ;

    % Rule out very low values of Sn
arg = Sn_znmat ./ exp(NumGrids.n_znmat) > 0.0001 ;

    % Bins for figure 9A
hBySn           = AggregateBin( gdzdn_znmat(arg) ,...
                      log( Sn_znmat(arg) ./ exp(NumGrids.n_znmat(arg)) ) ,...
                      ones(size(g_znmat(arg))),...
                      LogSnBins,...
                      ones(size(g_znmat(arg))),...
                      'Levels','sum') ;
hCDFBySn        = cumsum(hBySn) ;
NetPoachBySn    = AggregateBin( TotalNetPoach_zn(arg) ,...
                      log( Sn_znmat(arg) ./ exp(NumGrids.n_znmat(arg)) ) ,...
                      exp(NumGrids.n_znmat(arg)),...
                      LogSnBins,...
                      gdzdn_znmat(arg),...
                      'Levels','mean') ; 
 
    % Bins for figure 9B
pVacRate        = v_znmat(IndSort) ./ exp(NumGrids.n_znmat(IndSort)) ;
pHireRate       = Gn_znmat(IndSort) ;
pQuitRate       = 1 - Gn_znmat(IndSort) ;
PVacRate        =   cumsum( ( pVacRate(2:end) - pVacRate(1:(end-1)) )   .* pHireRate(2:end) ) ;
PVacRate(end+1) = PVacRate(end) ;
PHireRate       =   cumsum( ( pHireRate(2:end) - pHireRate(1:(end-1)) ) .* pVacRate(2:end) ) ;
PHireRate(end+1)= PHireRate(end) ;
PQuitRate       = - cumsum( pQuitRate(2:end) - pQuitRate(1:(end-1)) ) ;
PQuitRate(end+1)= PQuitRate(end) ;
PVacRate        = 3 * q * ( 1 - phi ) * PVacRate(IndUnsort) ;
PHireRate       = 3 * q * ( 1 - phi ) * PHireRate(IndUnsort) ;
PQuitRate       = 3 * Params.xi * p * PQuitRate(IndUnsort) ;
Ones            = ones(NumGrids.Nz,NumGrids.Nn) ;
PVacRateBySn    = AggregateBin( PVacRate(arg) ,...
                      log( Sn_znmat(arg) ./ exp(NumGrids.n_znmat(arg)) ) ,...
                      Ones(arg) ,...
                      LogSnBins,...
                      gdzdn_znmat(arg),...
                      'Levels','mean') ; 
PHireRateBySn   = AggregateBin( PHireRate(arg) ,...
                      log( Sn_znmat(arg) ./ exp(NumGrids.n_znmat(arg)) ) ,...
                      Ones(arg) ,...
                      LogSnBins,...
                      gdzdn_znmat(arg),...
                      'Levels','mean') ;
PQuitRateBySn   = AggregateBin( PQuitRate(arg) ,...
                      log( Sn_znmat(arg) ./ exp(NumGrids.n_znmat(arg)) ) ,...
                      Ones(arg) ,...
                      LogSnBins,...
                      gdzdn_znmat(arg),...
                      'Levels','mean') ; 
PTotalBySn      = PVacRateBySn + PHireRateBySn + PQuitRateBySn ;

    % Figure setup
Fig9 = figure(9);
set(Fig9,'Pos',[245.6667 364.3333 1410 390.6667]);
set(0,'defaultlinelinewidth',3);
font = 16;

    % Figure 9A plot
xvec_grid       = linspace( min(xvec(abs(xvec)<100)) , max(xvec(abs(xvec)<100)) , length(xvec) ) ;
hCDFBySn_sm     = interp1( xvec(abs(xvec)<1000)+1e-09*linspace(0,1,length(xvec(abs(xvec)<1000)))' , hCDFBySn(abs(xvec)<1000) , xvec_grid ) ;
hCDFBySn_sm     = smoothdata( hCDFBySn_sm , 'movmean' , 3 ) ;
NetPoachBySn_sm = interp1( xvec(abs(xvec)<1000)+1e-09*linspace(0,1,length(xvec(abs(xvec)<1000)))' , NetPoachBySn(abs(xvec)<1000) , xvec_grid ) ;
NetPoachBySn_sm = smoothdata( NetPoachBySn_sm , 'movmean' , 3 ) ;
subplot(1,2,1);
plot(xvec_grid,hCDFBySn_sm,'b-'); grid on; hold on;
xlabel('Log marginal surplus - $\log{S_n}$','interpreter','latex');
yyaxis left
set(gca,'YTick',0:0.2:1);
yyaxis right
plot(xvec_grid,NetPoachBySn_sm,'r--'); grid on; hold on;
ylim([-0.2,0.6]);
set(gca,'YTick',-0.2:0.1:0.6);
xlim([xmin,xmax]);
ax = gca;
ax.YAxis(1).Color = [0 0 1];
ax.YAxis(2).Color = [1 0 0];
l = legend('$H$ [Left]','Net poaching [Right]');
set(l,'interpreter','latex','location','NorthWest');
title('A. Net poaching rate and marginal surplus','interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
set(gca,'ticklabelinterpreter','latex','fontsize',font);
set(gca,'XTick',LogSnBinsGraph,'XTickLabel',LogSnBinsLabels);

    % Figure 9B plot
PVacRateBySn_sm = interp1( xvec(abs(xvec)<1000)+1e-09*linspace(0,1,length(xvec(abs(xvec)<1000)))' , PVacRateBySn(abs(xvec)<1000) , xvec_grid ) ;
PVacRateBySn_sm(isnan(PVacRateBySn_sm)) = 0 ;
PVacRateBySn_sm = smoothdata( PVacRateBySn_sm , 'movmean' , 3 ) ;
PHireRateBySn_sm = interp1( xvec(abs(xvec)<1000)+1e-09*linspace(0,1,length(xvec(abs(xvec)<1000)))' , PHireRateBySn(abs(xvec)<1000) , xvec_grid ) ;
PHireRateBySn_sm(isnan(PHireRateBySn_sm)) = 0 ;
PHireRateBySn_sm = smoothdata( PHireRateBySn_sm , 'movmean' , 3 ) ;
PQuitRateBySn_sm = interp1( xvec(abs(xvec)<1000)+1e-09*linspace(0,1,length(xvec(abs(xvec)<1000)))' , PQuitRateBySn(abs(xvec)<1000) , xvec_grid ) ;
PQuitRateBySn_sm(isnan(PQuitRateBySn_sm)) = 0 ;
PQuitRateBySn_sm = PQuitRateBySn_sm - 3*Params.xi*p ;
PQuitRateBySn_sm = smoothdata( PQuitRateBySn_sm , 'movmean' , 3 ) ;
subplot(1,2,2);
plot(xvec_grid,PVacRateBySn_sm,'r-');hold on;grid on;
plot(xvec_grid,PHireRateBySn_sm,'b:');hold on;grid on;
plot(xvec_grid,PQuitRateBySn_sm,'k--');hold on;grid on;
ylimdown = 1.1 * min(min(PHireRateBySn,PQuitRateBySn_sm')) ;
ylimup = 1.1 * max(max(PHireRateBySn,PQuitRateBySn)) ;
ax = gca;
ax.YAxis(1).Color = [0 0 0];
l = legend('$v(z,n)/n$','$H_n$','$H_v$ ');
set(l','interpreter','latex','location','NorthWest');
title('B. Components of the net poaching rate','interpreter','latex');
xlabel('Log marginal surplus - $\log{S_n}$','interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.3f');
xlim([xmin,xmax]);
ylim([ylimdown,2*ylimup])
set(gca,'XTick',LogSnBinsGraph,'XTickLabel',LogSnBinsLabels);
set(gca,'ticklabelinterpreter','latex','fontsize',font);
hold off

    % Save
print '-depsc' 'Created_figure_files/Fig9_PoachingSurplusDistribution_eps'
print '-dpng'  'Created_figure_files/Fig9_PoachingSurplusDistribution_png'

end
%%  FIGURE 10: NET POACHING RATES BY SIZE, AGE, LABOR PRODUCTIVITY AND NET EMPLOYMENT GROWTH RATE
%__________________________________________________________________________
if produce.fig10

    % Rename data
NetPoachingAge_data     = BDSage(:,9) ;
NetPoachingSize_data    = BDSsize(:,9) ;

    % Reconstruct model variables
momentnames = fieldnames(Moments);
for var = { 'NetPoachingSizeRanks' , 'MeanSnSizeRanks' , 'MeanSnSizeLevels' , 'NetPoachingAgeRanks' , ...
            'MeanSnAgeRanks' , 'MeanSnAgeLevels' , 'NetPoachingSizeLevels' , 'NetPoachingAgeLevels' , ...
            'MeanSnVapwRanks' , 'NetPoachingVapwRanks' , 'NetPoachingGrowthRanks' , ...
            'MeanSnGrowthRanks' }
    temp = [] ;
    for i=1:length(momentnames)
        if contains(momentnames{i},var{1})
            temp = [ temp , Moments.(momentnames{i}) ];     
        end
    end
    eval([ var{1} '= temp ;'])
end

    % Figure setup
Fig10 = figure(10);
set(Fig10,'Pos',[69 353 1352 643]);
set(0,'defaultlinelinewidth',3);
font = 14;
label_size = {'0-19','20-49','50-499','500+'};
label_age  = {'0-1','2-3','4-5','6-10','11+'};
ylim2 = [-.1,.1];

    % A. Net poaching by size
xvec_Size    = [1:1:numel(NetPoachingSizeLevels)]';
s3 = subplot(2,2,1);
i = ~isnan(NetPoachingSizeLevels);
NetPoachingSize_data_new = NetPoachingSize_data;
NetPoachingSize_data_new(3) = ( NetPoachingSize_data(3) + NetPoachingSize_data(4) ) / 2 ; % average 
NetPoachingSize_data_new(4) = [] ;
plot(xvec_Size,NetPoachingSize_data_new,'rx--'); hold on;
plot(xvec_Size(i),NetPoachingSizeLevels(i),'bo-'); hold on; grid on;
plot(xvec_Size,0*xvec_Size,'k-','linewidth',1);
xlabel('Size group (employees)','interpreter','latex','fontsize',font);
ylabel('Net poaching rate (quarterly)','interpreter','latex','fontsize',font);
title('A. Net poaching by size - Model / Data','interpreter','latex','fontsize',font);
l = legend('Data (Census J2J)','Model');
set(l,'interpreter','latex','location','NorthEast','fontsize',font);
set(gca,'XTick',xvec_Size,'XTickLabel',label_size);
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ylim([-0.01,0.07]);
ytickformat('%3.2f')
hAx = gca;             % handle to current axes
hAx.YAxis.Exponent=0;  % don't use exponent

    % B. Net poaching by age
xvec_Age    = [1:1:numel(NetPoachingAgeLevels)]';
s4 = subplot(2,2,2);
i = ~isnan(NetPoachingAgeLevels);
plot(xvec_Age,NetPoachingAge_data,'rx--'); hold on;
plot(xvec_Age(i),NetPoachingAgeLevels(i),'bo-'); hold on; grid on;
plot(xvec_Age,0*xvec_Age,'k-','linewidth',1);
xlabel('Age group (years)','interpreter','latex','fontsize',font);
ylabel('Net poaching rate (quarterly)','interpreter','latex','fontsize',font);
title('B. Net poaching by age - Model / Data','interpreter','latex','fontsize',font);
l = legend('Data (Census J2J)','Model');
set(l,'interpreter','latex','location','NorthEast','fontsize',font);
set(gca,'XTick',xvec_Age,'XTickLabel',label_age);
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ylim([-0.01,0.07]);
ytickformat('%3.2f')
hAx = gca;             % handle to current axes
hAx.YAxis.Exponent=0;  % don't use exponent

    % C. Net poaching by labor productivity rank
notnan = ~isnan(MeanSnVapwRanks) ;
NetPoachingVapwRanks_new = NetPoachingVapwRanks;
NetPoachingVapwRanks_new(6) = ( NetPoachingVapwRanks(5) + NetPoachingVapwRanks(7) ) / 2 ;
xvec = [5:10:95]';
ylim1 = [20,80] ;
Ytick1 = 0:10:100;
s1  = subplot(2,2,3);
i = ~isnan(NetPoachingVapwRanks_new);
plot(xvec(i),NetPoachingVapwRanks_new(i),'b-');hold on;
plot([-1,101],[0,0],'k-','linewidth',1);
xlim([0,100]);
grid on;
ax = gca;
ax.YAxis(1).Color = [0 0 0];
ylim(ylim2);
set(gca,'XTick',0:10:100);
xlabel('Firm rank in labor prod. distribution - $(y_i/n_i)$','interpreter','latex');
ylabel('Net poaching rate (quarterly)','interpreter','latex','fontsize',font);
title('C. Net poaching by rank in labor prod. distribution','interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.2f');

    % D. Net poaching by employment growth rank
xvec = [5:10:95]';
xvec = xvec(notnan) ;
s2 = subplot(2,2,4);
i = ~isnan(NetPoachingGrowthRanks);
plot(xvec(i),NetPoachingGrowthRanks(i),'b-');hold on;grid on;
plot([-1,101],[0,0],'k-','linewidth',1);
xlim([0,100]);
ylim(ylim2);
set(gca,'XTick',0:10:100);
xlabel('Firm rank in growth distribution - $\Delta\log n_i$','interpreter','latex');
ylabel('Net poaching rate (quarterly)','interpreter','latex','fontsize',font);
title('D. Net poaching by rank in growth distribution','interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',font);
ytickformat('%3.2f');
hold off

    % Save
print '-depsc' 'Created_figure_files/Fig10_PoachingByCharacteristics_eps'
print '-dpng'  'Created_figure_files/Fig10_PoachingByCharacteristics_png'

end
%%  FIGURE 11: FRICTIONLESS LIMIT: THE EFFECT OF INCREASING MATCH EFFICIENCY
%__________________________________________________________________________
if produce.fig11

    % Load previously computed mat file
load('Created_mat_files/Hopenhayn.mat');

    % Select points in grid to be displayed
iselect = gridA > 0.5 ;
gridselect = gridA(iselect);
MOMselect = MOM(iselect,:) ;
gridselect = log(gridselect) ; % choose logs for x axis
i0 = find(gridselect==0) ;

    % Figure setup
Fig11 = figure(11);
set(Fig11,'Pos',[0,0,1000,900]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
msize = 100 ;

    % A. Non-employment rate
subplot(2,2,1);
AX = plot( gridselect , MOMselect(:,1) ) ; hold on;
plot(0,MOMselect(i0,1),'color',[50 0 200]/255,'marker','.','markersize',msize)
set(gca,'fontsize',20);
ax = gca ;
ax.XGrid = 'on';
ax.YGrid = 'on';
title(sprintf('%s','A. Non-employment rate'),'interpreter','latex','fontsize',25); hold on;
xlabel(sprintf('%s','Log matching efficiency (rel. to benchmark)'),'interpreter','latex','fontsize',20); hold on;
set(gca,'GridLineStyle','--')
ax.GridAlpha = 1;
ax.LineWidth = 2 ;
set(gca,'GridColor',[ 0 0 0 ]) 
axis([ min(gridselect) , max(gridselect) , 0 , 1.05 * max(MOMselect( : , 1 )) ] )
AX(1).Color = [50 0 200]/255;
AX(1).LineStyle = '-';
AX(1).LineWidth = 5 ;
AX(1).Marker = 'none' ;
AX(1).MarkerSize = 15 ;

    % B. Standard Deviation of Marginal Productivity of Labor
subplot(2,2,2);
AX = plot( gridselect , MOMselect( : , 2 ) ) ;
hold on
plot(0,MOMselect(i0,2),'color',[50 0 200]/255,'marker','.','markersize',msize) ;
hold on
set(gca,'fontsize',20);
ax = gca ;
ax.XGrid = 'on';
ax.YGrid = 'on';
title(sprintf('%s','B. St.d. of MPL'),'interpreter','latex','fontsize',25); hold on;
xlabel(sprintf('%s','Log matching efficiency (rel. to benchmark)'),'interpreter','latex','fontsize',20); hold on;    
set(gca,'GridLineStyle','--')
ax.GridAlpha = 1;
ax.LineWidth = 2 ;
set(gca,'GridColor',[ 0 0 0 ]) 
axis([ min(gridselect) , max(gridselect) , 0 , 1.05 * max(MOMselect(:,2)) ] )
AX(1).Color = [50 0 200]/255;
AX(1).LineStyle = '-';
AX(1).LineWidth = 5 ;
AX(1).Marker = 'none' ;
AX(1).MarkerSize = 15 ;

    % C. Size-Productivity Correlation
subplot(2,2,3);
x = gridselect' ; 
y = MOMselect( : , 3 ) ;
b = regress(log(1-y),[ones(size(x,1),1) , x , x.^2  ]) ;
ysmooth = 1 - exp( b(1) + b(2) * x + b(3) * x.^2 ) ;
AX = plot( gridselect , ysmooth ) ;
hold on
plot(0,ysmooth(i0),'color',[50 0 200]/255,'marker','.','markersize',msize) ;
hold on
set(gca,'fontsize',20);
ax = gca ;
ax.XGrid = 'on';
ax.YGrid = 'on';
title(sprintf('%s','C. Correlation size-productivity'),'interpreter','latex','fontsize',25); hold on;
xlabel(sprintf('%s','Log matching efficiency (rel. to benchmark)'),'interpreter','latex','fontsize',20); hold on;    
set(gca,'GridLineStyle','--')
ax.GridAlpha = 1;
ax.LineWidth = 2 ;
set(gca,'GridColor',[ 0 0 0 ]) 
axis([ min(gridselect) , max(gridselect) , 0.95 * min(MOMselect( : , 3 )) , 1 ] )
AX(1).Color = [50 0 200]/255;
AX(1).LineStyle = '-';
AX(1).LineWidth = 5 ;
AX(1).Marker = 'none' ;
AX(1).MarkerSize = 15 ;

    % D. Output and TFP
subplot(2,2,4);
dlogY = log( MOMselect( : , 4 ) - MOMselect( : , 5 ) ) - log( ( MOMselect( i0 , 4 ) - MOMselect( i0 , 5 ) ) ) ;
dlogTFP = dlogY - log( ( 1 - MOMselect( : , 1 ) ).^Params.alpha ./ ( 1 - MOMselect( i0 , 1 ) ).^Params.alpha ) ;
AX = plot( gridselect , dlogY ) ; hold on;
AX2 = plot( gridselect , dlogTFP ) ;hold on
plot(0,0,'color',[50 0 200]/255,'marker','.','markersize',msize) ;
hold on
set(gca,'fontsize',20);
ax = gca ;
ax.XGrid = 'on';
ax.YGrid = 'on';
title(sprintf('%s','D. Log net output rel. to baseline'),'interpreter','latex','fontsize',25); hold on;
xlabel(sprintf('%s','Log matching efficiency (rel. to benchmark)'),'interpreter','latex','fontsize',20); hold on;        
set(gca,'GridLineStyle','--')
ax.GridAlpha = 1;
ax.LineWidth = 2 ;
set(gca,'GridColor',[ 0 0 0 ]) 
axis([ min(gridselect) , max(gridselect) , -0.3 , 1.05 * max(dlogY) ] )
legend('Output','TFP','Location','SouthEast')
AX(1).Color = [50 0 200]/255;
AX(1).LineStyle = '-';
AX(1).LineWidth = 5 ;
AX(1).Marker = 'none' ;
AX(1).MarkerSize = 15 ;
AX2(1).Color = 'r' ; 
AX2(1).LineStyle = '--';
AX2(1).LineWidth = 5 ;
AX2(1).Marker = 'none' ;
AX2(1).MarkerSize = 15 ;

    % Save
print '-depsc' 'Created_figure_files/Fig11_FrictionlessLimit_eps'
print '-dpng'  'Created_figure_files/Fig11_FrictionlessLimit_png'

end
%%  FIGURE 12: COMPARING THE AVERAGE FIRM LIFECYCLE TO THE LIFECYCLE OF SUPERSTARS
%__________________________________________________________________________
if produce.fig12

    % Simulate firms forward in time for 400 years
Nz          = NumGrids.Nz ;
Nn          = NumGrids.Nn ;
NT          = 400 ;
zrepstack 	= NumGrids.z_znmat(:) ;
firmsMiddle = zeros( Nz*Nn , NT ) ;
firmsTop    = zeros( Nz*Nn , NT ) ;
firmsMean   = zeros( Nz*Nn , NT ) ;

    % All firms in entering cohort
firmsMean(:,1) = pi0Trunc_zn/sum(pi0Trunc_zn);

    % Top productivity firms in entering cohort
shEntryTrunc    = cumsum( pi0Trunc_zn .* NumGrids.dzdn_znmat(:) ) ;
[ ~ , ind ]  	= min(abs( shEntryTrunc - .99999 )) ;
firmsTop(ind,1) = 1 ;

    % Evolve cohorts forward
A = fastExpm( 12*L ) ;
for t=2:NT
    firmsMean(:,t)  = A*firmsMean(:,t-1) ;
    firmsTop(:,t)   = A*firmsTop(:,t-1) ;
end
percentile  = 99.5;
sizeMean    = zeros( NT , 1 ) ;     % Panel A
sizeTop     = zeros( NT , 1 ) ;     % Panel B
sizeAll     = zeros( NT , 1 ) ;     % Panel B
for t=1:NT
    % Cohort of firms in the tail of entrants
    distTop         = reshape( firmsTop(:,t) , Nz , Nn ) ;
    distTop         = distTop/sum(sum(distTop));
    marginalTop     = cumsum( sum( distTop , 1 ) ) ;
    marginalTop     = marginalTop ./ marginalTop(end) ;
    % Cohort of all entering firms
    distMean        = reshape( firmsMean(:,t) , Nz , Nn ) ;
    distMean        = distMean/sum(sum(distMean));
    sizeMean(t)     = sum(exp(NumGrids.n_znmat(:)).*distMean(:)); % Average size (panel A)
    marginalAll     = cumsum( sum( distMean , 1 ) ) ;
    marginalAll     = marginalAll ./ marginalAll(end) ;
    % Save out size of percentile firm                        % Panel B
    ind             = sum( marginalTop < percentile/100 ) ;
    sizeTop(t)      = exp(NumGrids.n( ind+1 ) ) ; 
    ind             = sum( marginalAll < percentile/100 ) ;
    sizeAll(t)      = exp(NumGrids.n( ind+1 ) ) ; 
end

    % Figure setup
Fig12 = figure(12);
set(Fig12,'Pos',[1828 821 1508 421]);
font = 14;
order = 50;

    % A. Average lifecycle
subplot(1,3,1);
semilogy(sizeMean,'k-');hold on;grid on;
gvec        = (sizeMean(2:end)./sizeMean(1:end-1)-1);
fsize2      = zeros(size(sizeMean));
fsize2(1)   = sizeMean(1);
lambda      = 1.73;
gvec2       = lambda*gvec;
for aa = 2:numel(fsize2)
   fsize2(aa) = (1+gvec2(aa-1))*fsize2(aa-1);
end
semilogy(fsize2,'r--');
plot([1,NT],[10^4,10^4],'k-','linewidth',1);
l = legend('Average. Grows at rate $\overline{g}$','Growth at $1.73\overline{g}$ gets to $10^4$ in 110 years');
set(l,'interpreter','latex','location','Northwest','fontsize',11);
title('A. Average lifecycle','interpreter','latex'); 
ylim([10^0,2*10^6]);
ylabel('Employees','interpreter','latex','fontsize',font);
xlabel('Age','interpreter','latex');
set(gca,'ytick',10.^[0,1,2,3,4,5,6,7],'ticklabelinterpreter','latex','fontsize',font);
set(gca,'xtick',0:20:140);
xlim([0,140]);

    % B. Superstar lifecycle (model)
subplot(1,3,2);
semilogy(smooth(sizeAll,order),':','linewidth',3,'color',rgb('Blue')); grid on; hold on;
semilogy(smooth(sizeTop,order),'--','linewidth',3,'color',rgb('Blue')); grid on; hold on;
xlabel('Age','interpreter','latex');
plot([1,NT],[10^4,10^4],'k-','linewidth',1);
l = legend('Cohort of all entering firms','Cohort entering at 99$^{th}$ percentile');
set(l,'interpreter','latex','location','North','fontsize',12);
title('B. Superstar lifecycle ($99.5^{th}$ percentile of cohort)','interpreter','latex'); 
ylabel('Employees','interpreter','latex','fontsize',font);
ylim([10^0,2*10^6]);
set(gca,'ytick',10.^[0,1,2,3,4,5,6,7],'ticklabelinterpreter','latex','fontsize',font);
set(gca,'xtick',0:20:140);
xlim([0,140]);

    % C. Superstar lifecycles (data)
subplot(1,3,3);
age = [1906:2005]';
cd Input_data/rdq028_Supplementary_Data % Download Supplementary data to Luttmer (2011) in https://academic.oup.com/restud/article/78/3/1042/1566375?searchresult=1#supplementary-data

    % Load superstar data 
load     Bestbuy.txt; tbb =     Bestbuy(:,1); xbb =     Bestbuy(:,2);
load        Dell.txt; tde =        Dell(:,1); xde =        Dell(:,2);
load       Fedex.txt; tfx =       Fedex(:,1); xfx =       Fedex(:,2);
load          GM.txt; tgm =          GM(:,1); xgm =          GM(:,2);
load   Homedepot.txt; thd =   Homedepot(:,1); xhd =   Homedepot(:,2);
load       Intel.txt; tin =       Intel(:,1); xin =       Intel(:,2);
load       Lowes.txt; tlo =       Lowes(:,1); xlo =       Lowes(:,2);
load   McDonalds.txt; tmc =   McDonalds(:,1); xmc =   McDonalds(:,2);
load   Microsoft.txt; tmi =   Microsoft(:,1); xmi =   Microsoft(:,2);
load        Nike.txt; tnk =        Nike(:,1); xnk =        Nike(:,2);
load      Oracle.txt; toc =      Oracle(:,1); xoc =      Oracle(:,2); 
load         Sun.txt; tsn =         Sun(:,1); xsn =         Sun(:,2);
load  Tysonfoods.txt; ttf =  Tysonfoods(:,1); xtf =  Tysonfoods(:,2);
load      Abbott.txt; tab =      Abbott(:,1); xab =      Abbott(:,2);
load    Allstate.txt; tal =    Allstate(:,1); xal =    Allstate(:,2);
load      Amazon.txt; taz =      Amazon(:,1); xaz =      Amazon(:,2);
load      Disney.txt; tdi =      Disney(:,1); xdi =      Disney(:,2);
load DowChemical.txt; tdo = DowChemical(:,1); xdo = DowChemical(:,2);
load       Cisco.txt; tcs =       Cisco(:,1); xcs =       Cisco(:,2);
load        Ford.txt; tfo =        Ford(:,1); xfo =        Ford(:,2);
load         IBM.txt; tib =         IBM(:,1); xib =         IBM(:,2);
load         MMM.txt; t3m =         MMM(:,1); x3m =         MMM(:,2);
load     WalMart.txt; twa =     WalMart(:,1); xwa =     WalMart(:,2);
load  ProcterGamble.txt; tpg =  ProcterGamble(:,1); xpg =  ProcterGamble(:,2);
load HewlettPackard.txt; thp = HewlettPackard(:,1); xhp = HewlettPackard(:,2);
cd ../..
set(0,'defaultlinelinewidth',2);
semilogy(tbb,(xbb),'-');  legmat = [       'Bestbuy        ']; hold on;
semilogy(tde,(xde),'--');  legmat = [legmat;'Dell           '];
semilogy(tfx,(xfx),'-');  legmat = [legmat;'Fedex          '];
semilogy(tgm,(xgm),'--');  legmat = [legmat;'GM             '];
semilogy(thd,(xhd),'-');  legmat = [legmat;'Home Depot     '];
semilogy(tin,(xin),'--');  legmat = [legmat;'Intel          '];
semilogy(tlo,(xlo),'-' ); legmat = [legmat;'Lowes          '];
semilogy(tmc,(xmc),'--');  legmat = [legmat;'McDonalds      '];
semilogy(tmi,(xmi),'-');  legmat = [legmat;'Microsoft      '];
semilogy(tnk,(xnk),'--');  legmat = [legmat;'Nike           '];
semilogy(toc,(xoc),'-');  legmat = [legmat;'Oracle         '];
semilogy(tsn,(xsn),'--');  legmat = [legmat;'Sun            '];
semilogy(ttf,(xtf),'-');  legmat = [legmat;'Tyson Foods    '];
semilogy(tab,(xab),'-');
semilogy(tal,(xal),'--' );
semilogy(taz,(xaz),'-' );
semilogy(tdi,(xdi),'--' );
semilogy(tdo,(xdo),'-');
semilogy(tcs,(xcs),'--' );
semilogy(tfo,(xfo),'-' );
semilogy(thp,(xhp),'--' );
semilogy(tib,(xib),'-' );
semilogy(t3m,(x3m),'--' );
semilogy(twa,(xwa),'-' );
semilogy(tpg,(xpg),'--');

    % Plot
set(0,'defaultlinelinewidth',3);
ylim([10^0,2*10^6]);
xlim([min(age),max(age)]);
grid on;
ylabel('Employees','interpreter','latex','fontsize',font);
xlabel('Year','interpreter','latex','fontsize',font);
set(gca,'fontsize',font,'ticklabelinterpreter','latex');
set(gca,'ytick',10.^[0,1,2,3,4,5,6]);
set(gca,'xtick',1910:20:2000);
title('C. Superstar lifecycles (Luttmer 2011, Fig. 1)','interpreter','latex','fontsize',font);

    % Save
print '-depsc' 'Created_figure_files/Fig12_FirmLifecycleSuperstars_eps'
print '-dpng'  'Created_figure_files/Fig12_FirmLifecycleSuperstars_png'

end
%%  FIGURE 13: ENTRY AND JOB-TO-JOB HIRING RATES OVER THE GREAT RECESSION: AGGREGATE AND CROSS-SECTION
%__________________________________________________________________________
if produce.fig13

    % Load time series data
dat = csvread('Input_data/j2j_entry_matlab.csv',1,0); %J2Jhirerate and entryrate from 2001 to 2014 % 
dat = dat(1:14, :);

    % Unscaled
entry_J2J_dates             = datetime(2001,1,1):calyears(1):datetime(2014, 1,1);
J2Jt                        = dat(:, 1);
Entryt                      = dat(:, 2);
J2Jt_ts                     = timeseries(J2Jt);
J2Jt_ts.TimeInfo.Units      = 'Years';
J2Jt_ts.TimeInfo.StartDate  = '01-01-2001';     % Set start date.
J2Jt_ts.TimeInfo.Format     = 'yyyy';       % Set format for display on x-axis.
Entryt_ts                   = timeseries(Entryt);
Entryt_ts.TimeInfo.Units    = 'Years';
Entryt_ts.TimeInfo.StartDate= '01-01-2001';     % Set start date.
Entryt_ts.TimeInfo.Format   = 'yyyy';       % Set format for display on x-axis.

    % Scaled
J2Jt_scaled                         = (J2Jt - mean(J2Jt)) ./ std(J2Jt);
Entryt_scaled                       = (Entryt - mean(Entryt)) ./ std(Entryt);
J2Jt_scaled_ts                      = timeseries(J2Jt_scaled);
J2Jt_scaled_ts.TimeInfo.Units       = 'Years';
J2Jt_scaled_ts.TimeInfo.StartDate   = '01-01-2001';     % Set start date.
J2Jt_scaled_ts.TimeInfo.Format      = 'yyyy';       % Set format for display on x-axis.
Entryt_scaled_ts                    = timeseries(Entryt_scaled);
Entryt_scaled_ts.TimeInfo.Units     = 'Years';
Entryt_scaled_ts.TimeInfo.StartDate = '01-01-2001';     % Set start date.
Entryt_scaled_ts.TimeInfo.Format    = 'yyyy';       % Set format for display on x-axis.
Zero_ts                             = timeseries(0*Entryt_scaled);
Zero_ts.TimeInfo.Units              = 'Years';
Zero_ts.TimeInfo.StartDate          = '01-01-2001';     % Set start date.
Zero_ts.TimeInfo.Format             = 'yyyy';       % Set format for display on x-axis.

    % Load cross-section data
dat             = csvread('Input_data/matlab_ratechange_j2j_newfirm.csv',1,0);
yhat            = csvread('Input_data/predicted_ratechange_j2j.csv',1,0);
dlog_j2jhire    = dat(:, 1);
dlog_newfirm    = dat(:, 2);

    % Figure setup
Fig13 = figure(13);
set(Fig13,'Pos',[240 297 1437 432]);
set(0,'defaultlinelinewidth',3);
font = 16;

    % A. Time Series
subplot(1,2,1);
yyaxis left
plot(entry_J2J_dates, J2Jt,'b-'); hold on; grid on;
ylabel('$EE$ hire rate','interpreter','latex')
ytickformat('%4.3f');
set(gca,'Ytick',0.03:0.005:0.06);
yyaxis right
plot(entry_J2J_dates, Entryt,'r--')
ylabel('Establishment entry rate','interpreter','latex')
ylim([0.08,0.14]);
ytickformat('%4.2f');
set(gca,'Ytick',0.08:0.01:0.14);
ax = gca;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
l = legend('$EE$ hire rate [LHS]', 'Estab. entry rate [RHS]');
set(l,'interpreter','latex','fontsize',font-2,'location','southwest');
title('A. Time series','interpreter','latex');
xlabel('Year','interpreter','latex');
xticks(datetime(2000:2:2014,1,1));
xtickformat('yyyy');
set(gca,'fontsize',font,'ticklabelinterpreter','latex');

    % B. Changes at Metro Level
subplot(1,2,2);
plot(dlog_newfirm, dlog_j2jhire,'bx','markersize',10);hold on;
xvec            = -0.75:0.01:0.35;
p1              = polyfit(dlog_newfirm,yhat(:,1),2);     % Unweighted
p2              = polyfit(dlog_newfirm,yhat(:,2),2);     % Weighted
yhat_unweighted = polyval(polyfit(dlog_newfirm,yhat(:,1),2),xvec);
yhat_weighted   = polyval(polyfit(dlog_newfirm,yhat(:,2),2),xvec);
plot(xvec,yhat_unweighted,'g--','color',rgb('Black'));hold on;
plot(xvec,yhat_weighted,'g:','color',rgb('Red'));hold on;grid on;
l = legend('Metro area level data', ...
    sprintf('Fit (unweighted): $\\:\\:\\quad\\beta=%4.2f$',p1(2)), ...
    sprintf('Fit (emp. weighted): $\\beta=%4.2f$',p2(2)));
set(l,'interpreter','latex','fontsize',font,'location','south');
xlim([-0.08,0]);
ylim([-0.06,0]);
ylabel('Change in $EE$ hire rate','interpreter','latex');
xlabel('Change in establishment entry rate','interpreter','latex');
title('B. Changes at metro level: 2006 to 2009','interpreter','latex');
set(gca,'fontsize',font,'ticklabelinterpreter','latex');

    % Save
print '-depsc' 'Created_figure_files/Fig13_HiringRatesGreatRecession_eps'
print '-dpng'  'Created_figure_files/Fig13_HiringRatesGreatRecession_png'

end
%%  FIGURE 14: RESPONSE OF THE ECONOMY FOLLOWING A DISCOUNT RATE SHOCK
%__________________________________________________________________________
if produce.fig14

    % Load transition results
load Created_mat_files/Transition

    % Figure setup
figure(14)
nx = 3 ; ny = 3 ;
set(gcf,'Pos',[1500,1200,1100,600]);
color0 = [50 0 200]/255 ;
style0 = '-' ;
line1 = 'b-' ;
line2 = 'r--' ;
YLIM = 1 ; %0 or 1
legendgroup = 'quintile' ;
DisplayYears = 10; % Years displayed in the IRF
TransitionMoments = InitializeGraphs(TransitionMoments,TransitionMoments) ;
[ ~ , iYM ] = min( ( TransitionNum.T - 12 * DisplayYears ).^2  ) ;
years = 1:iYM;
years2 = 2:iYM;
maxyear = TransitionNum.T(iYM) / 12 ;
XTick = 0:2:20;
xticklabel = {'0','2','4','6','8','10','12','14','16','18','20'} ;

    % A. Unemployment rate
subplot(nx,ny,1)
plot(TransitionNum.T(years)/12,0.103 / TransitionMoments.u(1) * TransitionMoments.u(years),'Color',color0,'LineStyle',style0,'LineWidth',2); hold on;
title('A. Unemployment rate','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Level','interpreter','latex') 
xlim([0,maxyear])
if YLIM == 1
    ylim( [ 0.1 , 0.17 ] )
end
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % B. Number of entering firms
subplot(nx,ny,2)
plot(TransitionNum.T(years)/12,TransitionMoments.ME(years) / TransitionMoments.ME(1) - 1,'Color',color0,'LineStyle',style0,'LineWidth',2) ; hold on;
title('B. Number of entering firms','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Relative change','interpreter','latex') 
xlim([0,maxyear])
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % C. EU rate
subplot(nx,ny,3)
plot(TransitionNum.T(years)/12,TransitionMoments.EU(years),'Color',color0,'LineStyle',style0 ,'LineWidth',2); hold on;
title('C. EU rate','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Level','interpreter','latex') 
xlim([0,maxyear])
if YLIM == 1
    ylim( [ 0.015 , 0.04 ] )
end
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % D. Aggregate vacancies
subplot(nx,ny,4)
plot(TransitionNum.T(years)/12,TransitionMoments.V(years) / TransitionMoments.V(1) - 1,'Color',color0,'LineStyle',style0 ,'LineWidth',2); hold on;
title('D. Aggregate vacancies','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Relative change','interpreter','latex') 
xlim([0,maxyear])
if YLIM == 1
    ylim( [ -0.6 , 0 ] )
end
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % E. UE rate
subplot(nx,ny,5)
plot(TransitionNum.T(years)/12,TransitionMoments.UE(years),'Color',color0,'LineStyle',style0 ,'LineWidth',2); hold on;
title('E. UE rate','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Level','interpreter','latex') 
xlim([0,maxyear])
if YLIM == 1
    ylim( [ 0.1 , 0.17 ] )
end
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % F. Vacancy yield
subplot(nx,ny,6)
plot(TransitionNum.T(years)/12, ...
     TransitionMoments.VacYield(years) / TransitionMoments.VacYield(1) - 1,'Color',color0,'LineStyle',style0,'LineWidth',2); hold on;
plot(TransitionNum.T(years)/12, ...
     TransitionMoments.VacYield1000(years) / TransitionMoments.VacYield1000(1) - 1,'b-.','LineWidth',2); hold on;
plot(TransitionNum.T(years)/12, ...
     TransitionMoments.VacYield00(years) / TransitionMoments.VacYield00(1) - 1,line2 ,'LineWidth',2); hold on;
title('F. Vacancy yield','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Relative change','interpreter','latex') 
l = legend('Aggregate','$n>1,000$','$n<10$','Location','NorthEast') ;
set(l,'Interpreter','latex')
legend boxoff
xlim([0,maxyear])
if YLIM == 1
    ylim( [ 0 , 1.2 ] )
end
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % G. EE rate
subplot(nx,ny,7)
plot(TransitionNum.T(years)/12, ...
     TransitionMoments.EEh(years),'Color',color0,'LineStyle',style0 ,'LineWidth',2); hold on;
title('G. EE rate','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Level','interpreter','latex') 
xlim([0,maxyear])
ylim([0.008,0.014])
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % H. Net Poaching by Sn
subplot(nx,ny,8)
HPPARAMETER = 1600 ;
YSCALE = 0.7 ; 
TransitionMoments.NetPoachTopSn(2) = ( TransitionMoments.NetPoachTopSn(3) + TransitionMoments.NetPoachTopSn(1) ) / 2 ; % Correct numerical issues
TransitionMoments.NetPoachLowSn(2) = ( TransitionMoments.NetPoachLowSn(3) + TransitionMoments.NetPoachLowSn(1) ) / 2 ;
y1 = [ TransitionMoments.NetPoachTopSn(1) ; hpfilter(TransitionMoments.NetPoachTopSn(2:end),HPPARAMETER) ] ; % Smooth out numerical wiggles from input
y2 = [ TransitionMoments.NetPoachLowSn(1) ;  hpfilter(TransitionMoments.NetPoachLowSn(2:end),HPPARAMETER) ] ;
plot(TransitionNum.T(years)/12, ...
     ( y1(years) - TransitionMoments.NetPoachTopSn(1) ) ...
     / TransitionMoments.EEh(1),line1 ,'LineWidth',2); hold on;
plot(TransitionNum.T(years)/12, ...
     ( y2(years) - TransitionMoments.NetPoachLowSn(1) ) ...
     / TransitionMoments.EEh(1),line2 ,'LineWidth',2); hold on;
title('H. Net poaching by $S_n$','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Change, percent of EE','interpreter','latex') 
l = legend(['Top ' legendgroup],['Bottom ' legendgroup],'Location','SouthEast') ;
set(l,'Interpreter','latex')
legend boxoff
xlim([0,maxyear])
if YLIM == 1
    ylim( [ -YSCALE , YSCALE ] )
end
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % I. Output and TFP
subplot(nx,ny,9)
dlogY = log( TransitionMoments.MeanY(years) / TransitionMoments.MeanY(1) ) ;
dlogTFP = dlogY - log( (1-TransitionMoments.u(years)).^Params.alpha ./  (1-TransitionMoments.u(1)).^Params.alpha ) ;
plot(TransitionNum.T(years)/12,dlogY,'Color',color0,'LineStyle',style0 ,'LineWidth',2); hold on;
plot(TransitionNum.T(years)/12,dlogTFP,'b-.','LineWidth',2); hold on;
title('I. Output and TFP','interpreter','latex')
xlabel('Years','interpreter','latex') 
ylabel('Log change','interpreter','latex') 
xlim([0,maxyear])
if YLIM == 1
    ylim( [ -0.07 , 0.01 ] )
end
l = legend('Output','TFP','Location','SouthEast') ;
set(l,'Interpreter','latex')
legend boxoff
grid on
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTick',XTick,'xticklabel',xticklabel);
hold off

    % Save
print '-depsc' 'Created_figure_files/Fig14_DiscountRateIRF_eps'
print '-dpng'  'Created_figure_files/Fig14_DiscountRateIRF_png'

end
%%  FIGURE E.1: MINIMUM DISTANCE AS A FUNCTION OF EACH PARAMETER
%__________________________________________________________________________
if produce.figE1

    % Load estimation results
Params      = rmfield(Params,'c_f') ;
name        = fieldnames(Params)' ;
name(9)     = [] ;
load('Input_mat_files/EstimationResults.mat') ; % Heavy file available upon request to the authors
Params.xi   = Params.phi;
Params.A    = Params.chi;
Params      = rmfield(Params,{'phi','chi'});

    % Calculate Target-Moment Distance
DIST = zeros( 1 , length(Moments.EU) ) ;
for i=struct2cell(ChoiceOfTargets)'
    DIST = DIST + ((Moments.(i{1})-Targets.(i{1}))./Targets.(i{1})).^2 ;
end

    % Bounds
X.mu    = [ -.002   , -.0004] ;
X.b     = [ 0.5     , 1.6   ] ;
X.sigma = [ 0.008   , .04   ] ;
X.xi    = [ 0.1     , 0.2   ] ;
X.delta = [ 0.01    , 0.03  ] ;
X.alpha = [ .7      , .9    ] ;
X.A     = [ 0.1     , 0.25  ] ;
X.zeta  = [ 6       , 15    ] ;

    % Figure setup
FigE1 = figure(101);
set(FigE1,'Pos',[69 353 1352 643]);
k=1 ;

    % Plot
for i = name(1:end-1)
    xplot = Params.(i{1}) ;
    xgrid = linspace(X.(i{1})(1),X.(i{1})(2),6) ;
    mindist = zeros(length(xgrid)-1,1) ;
    for j = 2:length(xgrid)
        mindist(j-1) = nanmin( DIST( xplot >= xgrid(j-1) & xplot < xgrid(j) ) ) ;
    end
    mindist = mindist - min( mindist );
    xgrid = ( xgrid(2:end) + xgrid(1:end-1) ) / 2 ;
    yplot = mindist;
    subplot(2,4,k);
    AX = plot( xgrid , yplot ) ;
    hold on
    [ ~ , mmm ] = min( yplot.^2 ) ;
    plot([ xgrid(mmm),xgrid(mmm)],[0,.5*min(max(DIST),10)],'LineStyle','-','Color',[230 0 0]/255,'linewidth',4);
    set(gca,'fontsize',10);
    AX(1).Color = [0 0 0]/255;
    AX(1).LineStyle = '-.';
    AX(1).LineWidth = 4;
    AX(1).Marker = 'square';
    AX(1).MarkerSize = 10;
    AX(1).MarkerFaceColor = [0 0 0]/255;
    ax = gca ;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    set(ax,'XLim',[min(xgrid) max(xgrid)])
    xlabel(sprintf('%s',char(ParameterName.(i{1}))),'interpreter','latex','fontsize',18); hold on;
    ylabel('Minimum distance','interpreter','latex','fontsize',18);hold on;
    ylim([ 0 , .05 ]);
    set(gca,'GridLineStyle','--')
    ax.GridAlpha = 1;
    ax.LineWidth = 2 ;
    set(gca,'GridColor',[ 0 0 0 ])
    ax = ancestor(AX(1), 'axes') ;
    ax.XAxis.Exponent = 0 ;
    drawnow
    k = k + 1 ;
end
hold off

    % Save
print '-depsc' 'Created_figure_files/FigE1_MinimumDistance_eps'
print '-dpng'  'Created_figure_files/FigE1_MinimumDistance_png'

end
%%  FIGURE E.2: EACH TARGETED MOMENT AGAINST EACH PARAMETER
%__________________________________________________________________________
if produce.figE2
    
    % Save original Params struct
Params = rmfield(Params,'n_0') ;
Params = rmfield(Params,'c_f') ;
Params_orig = Params ;

    % Load Jacobian
load('Input_mat_files/EstimationJacobian.mat') ;
Jacobian.xi = Jacobian.phi;
Jacobian.A = Jacobian.chi;
Jacobian = rmfield(Jacobian,{'phi','chi'});

    % Bounds
X.mu    = [ -.002   , -.0004] ;
X.b     = [ 0.5     , 1.6   ] ;
X.sigma = [ 0.008   , .04   ] ;
X.xi    = [ 0.1     , 0.2   ] ;
X.delta = [ 0.01    , 0.03  ] ;
X.alpha = [ .7      , .9    ] ;
X.A     = [ 0.1     , 0.25  ] ;
X.zeta  = [ 6       , 15    ] ;

    % Target names
TargetName.BDSExitRateUW =                  'Exit rate' ;
TargetName.BDSJobDestructionIncumbents =    'JD rate of inc.' ;
TargetName.BDSStdGrowthRate =               'St.d. of growth rate' ;
TargetName.EEoverEquarterly =               'EE rate' ;
TargetName.EUoverEquarterly =               'EU rate' ;
TargetName.EmpShare500 =                    'Empl. share 500+' ;
TargetName.URate =                          'Unempl. rate' ;
TargetName.BDSJCRateAge1 =                  'JC rate, age 1';

    % Figure Setup
FigE2 = figure(102);
set(FigE2,'Pos',[100,500,1200,600]);
k = 1 ;
color2 = [230 0 30]/255 ;

    % Plot
for i=fieldnames(Jacobian)'
    xgrid = linspace(X.(i{1})(1),X.(i{1})(2),6) ;
    xgrid = ( xgrid(2:end) + xgrid(1:end-1) ) / 2 ;
    Moments = Jacobian.(i{1}).Moments ;
    xaxis = Jacobian.(i{1}).value ;
    subplot(2,4,k);
    for m=struct2cell(ChoiceOfTargets)'
        xaxis(imag(Moments.(m{1}))>0) = NaN ;
    end
    x1 = X.(i{1})(1) ;
    x2 = X.(i{1})(2) ;
    yaxis = Moments.(ChoiceOfTargets.(i{1})) ;
    if strcmp(i{1},'mu')
        yaxis(xaxis>-0.0008) = NaN ;
    end
    yaxis( imag( xaxis ) ~= 0 ) = NaN ;
    xaxis( imag( xaxis ) ~= 0 ) = NaN ;
    yaxis = yaxis( xaxis >= x1 & xaxis <= x2 ) ;
    xaxis = xaxis( xaxis >= x1 & xaxis <= x2 ) ;
    yaxis = smoothdata( yaxis , 'movmean' , 3 ) ;
    [~,m] = min(abs(xaxis-Params.(i{1}))) ;
    mom = yaxis(m) ;
    AX = plot(xaxis,yaxis); hold on;
    xlim([min(xgrid) max(xgrid)])
    AX(1).Color = [0 0 0]/255;
    AX(1).LineStyle = '-.';
    AX(1).LineWidth = 4;
    AX(1).Marker = 'square';
    AX(1).MarkerSize = 10;
    AX(1).MarkerFaceColor = [0 0 0]/255;
    ax = gca ;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    set(ax,'XLim',[min(xaxis) max(xaxis)])
    xlabel(sprintf('%s',char(ParameterName.(i{1}))),'interpreter','latex','fontsize',18); hold on;
    ylabel(sprintf('%s',TargetName.(ChoiceOfTargets.(i{1}))),'interpreter','latex','fontsize',18);
    set(gca,'GridLineStyle','--')
    ax.GridAlpha = 1;
    ax.LineWidth = 2 ;
    set(gca,'GridColor',[ 0 0 0 ])
    plot(exp([log(.25),log(4)]).*Params_orig.(i{1}),[mom,mom],'LineStyle','--','Color',color2,'linewidth',3);
    ax = ancestor(AX(1), 'axes') ;
    ax.XAxis.Exponent = 0 ;
    k = k+1 ;
end

    % Save
print '-depsc' 'Created_figure_files/FigE2_TargetMomentsParameters_eps'
print '-dpng'  'Created_figure_files/FigE2_TargetMomentsParameters_png'

end