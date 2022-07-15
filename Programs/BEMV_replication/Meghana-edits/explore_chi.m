%% Compare SS moments for different values of chi

%{
to-do:
- check for different equilbria (different intial distribution guesses)
- turn off z shocks
%}

%% Housekeeping
set(0, 'DefaultLineLineWidth', 2);

clear
close all
clc

vals = cell(7,1);
% Save results
for chi = 0:size(vals,1)-1
    vals{chi+1,1} = ['chi' num2str(chi)];
    file = ['Created_mat_files/SS_results_chi_' num2str(chi) '.mat'];
    data.(vals{chi+1}) = load(file);
end

% easy to view
chi0=data.chi0.Moments;
chi1=data.chi1.Moments;
chi2=data.chi2.Moments;
chi3=data.chi3.Moments;
chi4=data.chi4.Moments;
chi5=data.chi5.Moments;
chi6=data.chi6.Moments;

%% make some preliminary figures

lsize      = zeros(size(vals,1),1); % mean log size
std_lsize  = zeros(size(vals,1),1); % std log size
elasticity = zeros(size(vals,1),1); % dlogemp/dlogz - emp weighted
tfp        = zeros(size(vals,1),1); % mean tfp
gr         = zeros(size(vals,1),1); % mean growth rate
urate      = zeros(size(vals,1),1); % unemp rate

for chi = 0:size(vals,1)-1
    lsize(chi+1)      = data.(vals{chi+1}).Moments.MeanLogSize;
    elasticity(chi+1) = data.(vals{chi+1}).Moments.dLogSizedLogZEmp;
    tfp(chi+1)        = data.(vals{chi+1}).Moments.MeanTFP;
    gr(chi+1)         = data.(vals{chi+1}).Moments.MeanGrowthRate;
    urate(chi+1)      = data.(vals{chi+1}).Moments.URate;
    std_lsize(chi+1)  = data.(vals{chi+1}).Moments.StdLogSize;
end

figure
subplot(2,2,1);
plot(0:6, lsize)
title('mean log size')

subplot(2,2,2);
plot(0:6, elasticity)
title('dlog size / dlog z')

subplot(2,2,3);
plot(0:6, tfp)
title('mean tfp')

subplot(2,2,4);
plot(0:6, gr)
title('mean growth rate')

figure

subplot(2,1,1);
plot(0:6,urate)
title('unemployment rate')

subplot(2,1,2);
plot(0:6, std_lsize)
title('std log size')




