%%Summarise ML run
%30/03/2017 DK Shin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO
%   * save only vars specific to analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%configs
%experiment specific
dir_data='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\ML_Optimal_Transport\ml_14_param\run7';    % data dir
dir_log=[dir_data,'\log'];          % dir for logs (currently output_write mat files)
dir_output=[dir_data,'\output'];    % dir of output

fopts.filepath=[dir_data,'\d'];     %filepath to data
file_id=1:1462;  % file ids to analyse

% NEW COST CONFIG
fopts.xlim=[-35e-3, 25e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
fopts.ylim=[-15e-3, 20e-3];

fopts.dt=10e-3;
fopts.t_0=2.6075;
fopts.n_pulse=110;

%%%pass/fail/penalty
%pulse picking success rate
fopts.win_capture_rate=0.5;     %pulse/window capture rate for pass/fail (0 to turn off)
%pkpk amplitude penalty
fopts.pkpk_penalty=true;
%low number penalty
fopts.num_win_penalty=15;       %minimum counts per window (avgd) to pass number penalty
%width penalty
fopts.width_sat=[6e-3,3e-3,4e-3];
fopts.penalty_width_sat=0.005;      %set to [] to turn penalty OFF

% OLD COST CONFIG
% fopts.xlim=[-35e-3, 20e-3];
% fopts.ylim=[-10e-3, 18e-3];
% 
% fopts.dt=10e-3;
% fopts.t_0=2.6073;
% fopts.n_pulse=110;
% 
% %%%pass/fail/penalty
% %pulse picking success rate
% fopts.win_capture_rate=0.8;      %pulse/window capture rate for pass/fail
% %pkpk amplitude penalty
% fopts.pkpk_penalty=false;
% %low number penalty
% fopts.num_win_penalty=[];       %minimum counts per window (avgd) to pass number penalty
% %width penalty
% fopts.width_sat=[8e-3,5e-3,6e-3];
% fopts.penalty_width_sat=[];      %set to [] to turn penalty OFF

%output
fopts.log=false;
fopts.graphics=false;

%SHUNT BOUNDARIES
% 6 Hz trap
param_boundary=[3.4,0.14,0,0.87];   % shunt boundaries [Q_i,Q_f,Sh_i,Sh_f]

%%main
n_shots=length(file_id);

cost=zeros(n_shots,1);
cost_unc=zeros(n_shots,1);
bad=cell(n_shots,1);
output=cell(n_shots,1);
config=cell(n_shots,1);

parfor ii=1:n_shots
    [cost(ii),cost_unc(ii),bad{ii},output{ii},config{ii}]=cost_calculator('slosh',file_id(ii),fopts);
    progress_message=sprintf('%d/%d Complete',ii,n_shots);
    disp(progress_message);
end

%build time vector
t=fopts.t_0+fopts.dt*(0:fopts.n_pulse-1);

%% Get shunt params from logs
if isdir(dir_log)
    dir_orig=cd(dir_log);
    log_files=dir('output_write*.mat');
    
    idx_params=zeros(length(log_files),1);
    shunt_params=cell(length(log_files),1);
    
    for ii=1:length(log_files)
        [idx_params(ii),shunt_params{ii}]=get_params_from_log(log_files(ii).name);
    end
    
    cd(dir_orig);
else
    warning('Parameter log directory \"%s\" cannot be found. Do not interpret params given by this script.',dir_log);
    idx_params=NaN;
    shunt_params=NaN;
end

%% Summary
%%ML summary
%all costs
cost_not_bad=100;       % cost for "bad" run
hfig_ml_not_bad=figure();
hold on;
scatter(file_id(cost<cost_not_bad),cost(cost<cost_not_bad),10,'bo',...
    'filled');   % plot all runs that aren't "bad"
% scatter(file_id(cost<cost_max_disp_all),cost(cost<cost_max_disp_all),10,'go',...
%     'filled');   % plot "bad" runs
box on;
axis tight;
titlestr='ML: all shots except bad';
title(titlestr);
xlabel('Iteration');
ylabel('Cost function (m)');

%low costs / good shots only
cost_good=30e-3;
hfig_ml_good=figure();
hold on;
scatter(file_id(cost<cost_good),cost(cost<cost_good),10,'o',...
    'filled');   % don't plot bad costs
box on;
axis tight;
titlestr=sprintf('ML: cost $<$ %.3g m',cost_good);
title(titlestr);
xlabel('Iteration');
ylabel('Cost function (m)');

%%report on BEST file
[best_cost,best_id]=min(cost);
best_osc_std=output{best_id}.osc_std;
best_osc=output{best_id}.osc;

hfig_best=figure();
hold on;
for ii=1:3
    plot(t,best_osc(:,ii)-mean(best_osc(:,ii)),...
        'Linewidth',1.5);
end
box on;
axis tight;
titlestr=sprintf('Best 14-ML [%d]: cost=%0.3g',best_id,best_cost);
title(titlestr);
xlabel('Time (s)');
ylabel('Position (m)');
legend({'X','Y','Z'});

%plot the profile
hfig_best_profile=figure();
best_profile=shunt_params{idx_params==file_id(best_id)};
plot_shunt_profile(best_profile,param_boundary,hfig_best_profile);

box on;
axis tight;

titlestr=sprintf('Best 14-ML [%d]: cost=%0.3g',best_id,best_cost);
title(titlestr);
xlabel('Profile segment');
ylabel('DAQ voltage (V)');

%%report on exponential ramp - 1st run
if file_id(best_id)~=1
    exp_cost=cost(1);
    exp_osc_std=output{1}.osc_std;
    exp_osc=output{1}.osc;
    
    hfig_exp=figure();
    hold on;
    for ii=1:3
            plot(t,exp_osc(:,ii)-mean(exp_osc(:,ii)),...
        'Linewidth',1.2);
    end
    box on;
    titlestr=sprintf('First exponential profile: cost=%0.3g',exp_cost);
    title(titlestr);
    xlabel('Time (s)');
    ylabel('Position (m)');
    axis tight;
    legend({'X','Y','Z'});
    
    %plot the profile
    hfig_exp_profile=figure();
    exp_profile=shunt_params{idx_params==1};
    plot_shunt_profile(exp_profile,param_boundary,hfig_exp_profile);
    
    titlestr=sprintf('First exponential profile: cost=%0.3g',exp_cost);
    title(titlestr);
    xlabel('Profile segment');
    ylabel('DAQ voltage (V)');
end

%%Breathing mode
%best only
best_width=output{best_id}.width;

hfig_best_breathing=figure();
hold on;
for ii=1:3
    plot(t,best_width(:,ii),...
        'LineWidth',1.5);
end
box on;
titlestr=sprintf('Best 14-ML [%d]: cost=%0.3g',best_id,best_cost);
title(titlestr);
xlabel('Time (s)');
ylabel('PAL rms width (m)');
axis tight;
legend({'X','Y','Z'});

%%Breathing mode and oscillation optimisation
%evaluate breathing mode
bec_width_rms=zeros(n_shots,3);
for ii=1:n_shots
    if ~isempty(output{ii})
        %breathing mode magnitude as rms of pulse width (std)
        bec_width_rms(ii,:)=std(output{ii}.width,1);
    else
        bec_width_rms(ii,:)=nan;
    end
end
hfig_breathingmode=figure();
hold on;
%original cost
cost_max_disp_better=10e-3;
scatter(file_id(cost<cost_max_disp_better),cost(cost<cost_max_disp_better),10,'o',...
    'filled');   % don't plot bad costs
%breathing mode magnitude
for jj=1:3
    scatter(1:n_shots,bec_width_rms(:,jj),10,'^','filled');   %breathing mode in each axis
end
box on;
titlestr=sprintf('Cost function and breathing mode');
title(titlestr);
xlabel('Iteration');
ylabel('Cost function (m)');
legend({'original cost','$\Delta\sigma_X$','$\Delta\sigma_Y$','$\Delta\sigma_Z$'});

%% Save output
%output directory
if ~isdir(dir_output)
    mkdir(dir_output);
end

%%Figures PNG
saveas(hfig_ml_not_bad,[dir_output,'\','ml_not_bad.png']);       % ML progress - all costs except bad
saveas(hfig_ml_good,[dir_output,'\','ml_good.png']);  % ML progress - low costs
saveas(hfig_best,[dir_output,'\','best_oscillation.png']);        % best shot - oscillations
saveas(hfig_best_profile,[dir_output,'\','best_profile.png']);    % best shot - shunt profile
saveas(hfig_exp,[dir_output,'\','exp_oscillation.png']);          % exp benchmark - oscillations
saveas(hfig_exp_profile,[dir_output,'\','exp_profile.png']);      % exp benchmark - shunt profile
saveas(hfig_best_breathing,[dir_output,'\','best_breathing.png']);  	% best shot - rms width oscillation
saveas(hfig_breathingmode,[dir_output,'\','breathing_mode.png']); % breathing mode

%%Figures .fig
saveas(hfig_ml_not_bad,[dir_output,'\','ml_not_bad.fig']);       % ML progress - all costs except bad
saveas(hfig_ml_good,[dir_output,'\','ml_good.fig']);  % ML progress - low costs
saveas(hfig_best,[dir_output,'\','best_oscillation.fig']);        % best shot - oscillations
saveas(hfig_best_profile,[dir_output,'\','best_profile.fig']);    % best shot - shunt profile
saveas(hfig_exp,[dir_output,'\','exp_oscillation.fig']);          % exp benchmark - oscillations
saveas(hfig_exp_profile,[dir_output,'\','exp_profile.fig']);      % exp benchmark - shunt profile
saveas(hfig_best_breathing,[dir_output,'\','best_breathing.fig']);  	% best shot - rms width oscillation
saveas(hfig_breathingmode,[dir_output,'\','breathing_mode.fig']); % breathing mode

%%Data
save([dir_output,'\','ml_summary_',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat']);  %save all vars