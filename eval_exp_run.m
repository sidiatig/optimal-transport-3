%%Summarise exponential profile trap relaxation
%06/04/2017 DK Shin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Configs
%experiment specific
dir_data='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\ml_shunt_14param\exp_shunt\0_5s';    % data dir
dir_log=[dir_data,'\log'];          % dir for logs (currently output_write mat files)
dir_output=[dir_data,'\output'];    % dir of output

fopts.filepath=[dir_data,'\d'];     %filepath to data
file_id=1:1600;  % file ids to analyse

fopts.xlim=[-35e-3, 20e-3];
fopts.ylim=[-10e-3, 18e-3];

fopts.t_0=1.1074;

fopts.dt=10e-3;
fopts.n_pulse=280;

%general flags
% fopts.min_count=1000;
fopts.num_in_win=30;        % minimum counts per window (avgd) to pass number penalty
fopts.min_window_pass=10;   % minimum good windows to pass for oscillation analysis

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

for ii=1:n_shots
    [cost(ii),cost_unc(ii),bad{ii},output{ii},config{ii}]=cost_calculator('slosh',file_id(ii),fopts);
    progress_message=sprintf('%d/%d Complete',ii,n_shots);
    disp(progress_message);
end

%build time vector
t=fopts.t_0+fopts.dt*(0:fopts.n_pulse-1);

%% Get EXP shunt params from log
%TODO - for more shots than param number
%1. use directly the mat file
%!!!assumes Labview IterNum is synched with recorded shot ID!!!
PARAMLOG=load([dir_log,'/param_data.mat']);

%2. use the txt log - written by interface per shot

%parse data
param_values=PARAMLOG.param_values;
n_par_set=size(param_values,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO - we suffer here since I forgot to save ntaus with the params
Vquad_0=3.4;        %initial set point
Vquad_inf=0.14;     %asymptotic values of exponential ramp (boundary condition)
Vshunt_0=0;
Vshunt_inf=0.87;

ntau_array=zeros(n_par_set,2);   %ntau to calculate
for ii=1:n_par_set
    %get profile
    Vquad_this=param_values(ii,1:7);
    Vshunt_this=param_values(ii,8:end);
    
    %%get number of exp constants decayed
    %quad
    logv_temp=log(Vquad_this-Vquad_inf);
    ntau_array(ii,1)=log(Vquad_0-Vquad_inf)-logv_temp(end);
    
    %shunt
    logv_temp=log(Vshunt_inf-Vshunt_this);
    ntau_array(ii,2)=log(Vshunt_inf-Vshunt_0)-logv_temp(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Summary
%%EXP Summary
%all costs
cost_max_disp_all=100;       % cost for "bad" run
id_except_bad=file_id(cost<cost_max_disp_all);
hfig_exp_all=figure();
hold on;
scatter3(ntau_array(id_except_bad,1),ntau_array(id_except_bad,2),cost(id_except_bad),10,'bo',...
    'filled');   % plot all runs that aren't "bad"
% scatter(file_id(cost<cost_max_disp_all),cost(cost<cost_max_disp_all),10,'go',...
%     'filled');   % plot "bad" runs
box on;
titlestr='EXP: all shots except bad';
title(titlestr);
xlabel('NTAU1');
ylabel('NTAU2');
zlabel('Cost (m)');
view(3);

%low costs / good shots only
cost_max_disp_good=100e-3;
id_good=file_id(cost<cost_max_disp_good);
hfig_exp_lowcost=figure();
hold on;
scatter3(ntau_array(id_good,1),ntau_array(id_good,2),cost(id_good),10,'bo',...
    'filled');   % plot all runs that aren't "bad"
box on;
titlestr=sprintf('EXP - cost < %.3g m',cost_max_disp_good);
title(titlestr);
xlabel('NTAU1');
ylabel('NTAU2');
zlabel('Cost (m)');
view(3);

%%report on BEST file
[best_cost,best_id]=min(cost);
best_osc_std=output{best_id}.osc_std;
best_osc=output{best_id}.osc;

best_ntau=ntau_array(best_id,:);

hfig_best=figure();
hold on;
for ii=1:3
    plot(t,best_osc(:,ii)-mean(best_osc(:,ii)),...
        'Linewidth',1.5);
end
box on;
titlestr=sprintf('Best EXP shunt: ntau=[%0.3g,%0.3g]',best_ntau(1),best_ntau(2));
title(titlestr);
xlabel('Time (s)');
ylabel('Position (m)');
axis tight;
legend({'X','Y','Z'});

%plot the best EXP profile
hfig_best_profile=figure();
best_profile=param_values(best_id,:);
plot_shunt_profile(best_profile,param_boundary,hfig_best_profile);

titlestr=sprintf('Best EXP shunt: ntau=[%0.3g,%0.3g]',best_ntau(1),best_ntau(2));
title(titlestr);
xlabel('Profile segment');
ylabel('DAQ voltage (V)');
 
% %% Breathing mode
% %evaluate breathing mode
% bec_width_rms=zeros(n_shots,3);
% for ii=1:n_shots
%     if ~isempty(output{ii})
%         %breathing mode magnitude as rms of pulse width (std)
%         bec_width_rms(ii,:)=std(output{ii}.width,1);
%     else
%         bec_width_rms(ii,:)=nan;
%     end
% end
% hfig_breathingmode=figure();
% hold on;
% %original cost
% cost_max_disp_better=10e-3;
% scatter(file_id(cost<cost_max_disp_better),cost(cost<cost_max_disp_better),10,'o',...
%     'filled');   % don't plot bad costs
% %breathing mode magnitude
% for jj=1:3
%     scatter(1:n_shots,bec_width_rms(:,jj),10,'^','filled');   %breathing mode in each axis
% end
% box on;
% titlestr=sprintf('Cost function and breathing mode');
% title(titlestr);
% xlabel('Iteration');
% ylabel('Cost function (m)');
% legend({'original cost','\Delta\sigma_X','\Delta\sigma_Y','\Delta\sigma_Z'});
% 
% 
%% Save output
%output directory
if ~isdir(dir_output)
    mkdir(dir_output);
end

%%Figures
saveas(hfig_exp_all,[dir_output,'\','ml_progress_all.png']);       % ML progress - all costs except bad
saveas(hfig_exp_lowcost,[dir_output,'\','ml_progress_good.png']);  % ML progress - low costs
saveas(hfig_best,[dir_output,'\','best_oscillation.png']);        % best shot - oscillations
saveas(hfig_best_profile,[dir_output,'\','best_profile.png']);    % best shot - shunt profile

% saveas(hfig_breathingmode,[dir_output,'\','breathing_mode.png']); % breathing mode

%%Data
save([dir_output,'\','ml_summary_',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat']);  %save all vars