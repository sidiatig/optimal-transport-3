%%14-param ML data points
% %all, some data evaluated with old cost function
% t_ML=[2 1 0.5 1 1 0.5 122e-3 117.5e-3 104.6e-3 100.5e-3 88.5e-3 81.8e-3 ... 
%     73.6e-3 54e-3 66.4e-3 60.9e-3 54.2e-3 49.2e-3 46.2e-3 43.8e-3 41.5e-3 ...
%     39.9e-3 0.25];
% osc_ML=[2.2 1.4 2.4 1.7 1 3.7 7.5 6.7 8.8 5.7 54 6.5 8.5 NaN 16 3.9 8 9.5 30.7 19.7 78.3 1580 1.52];

% %outlier removed and MEAN of same data
% t_ML=[2 1 0.5 122e-3 117.5e-3 104.6e-3 100.5e-3 81.8e-3 ... 
%     73.6e-3 54e-3 66.4e-3 60.9e-3 54.2e-3 49.2e-3 46.2e-3 43.8e-3 41.5e-3 ...
%     39.9e-3 0.25];
% osc_ML=[2.2 1.37 3.05 7.5 6.7 8.8 5.7 6.5 8.5 NaN 16 3.9 8 9.5 30.7 19.7 78.3 1580 1.52];

% all data with updated cost
t_ML=[2 1 0.5 1 1 0.5 122e-3 117.5e-3 104.6e-3 100.5e-3 88.5e-3 81.8e-3 ... 
    73.6e-3 54e-3 66.4e-3 60.9e-3 54.2e-3 49.2e-3 46.2e-3 43.8e-3 41.5e-3 ...
    39.9e-3 0.25];
    
osc_ML=[2.2 1.5 2.4 1.7 1.0 3.7 7.5 6.7 8.8 5.7 54 6.5 8.5 NaN 16 3.9 8 9.5 30.7 19.7 78.3 1580 1.52];

%%exponential benchmarks
%Old data
% t_exp=[2 1 0.5 0.25 0.125 0.06 0.122];
% osc_exp=[1.2 1.2 2.3 1.6 NaN NaN 9.8];
%New data
t_exp=[2 1 0.5 0.25 0.122 117.5e-3];
osc_exp=[0.6 1.4 1.5 1.7 24.1 781];

hfig_ml_compare=figure();
hold on;

size=50;

scatter(t_ML,osc_ML,size,'bo','filled','MarkerEdgeColor','k');
scatter(t_exp,osc_exp,size,'r^','filled','MarkerEdgeColor','k');

grid on;

% axis tight;
% xlim([min([t_ML,t_exp]),1.1*max([t_ML,t_exp])]);  % nonzero min X
xlim([0,1.1*max([t_ML,t_exp])]);    %set min X to 0
ylim([0,1.1*max([osc_ML,osc_exp])]);    %set min Y to 0

xlabel('Shunt duration (s)','fontsize',13);
% ylabel('RMS oscillation (mm)','fontsize',13);
ylabel('Best cost found','fontsize',13);
legend({'ML-14','Best exp'},'fontsize',12);
box on;

saveas(hfig_ml_compare,'oscillation_comparison.png');


%create loglog
ax=gca;
ax.XScale='log';
ax.YScale='log';
saveas(hfig_ml_compare,'oscillation_comparison_loglog.png');