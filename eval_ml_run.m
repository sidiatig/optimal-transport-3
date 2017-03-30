%%Summarise ML run
%30/03/2017 DK Shin

%%configs
%experiment specific
dir_data='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\ml_shunt_14param\run4';    % data dir
dir_log=[dir_data,'\log'];      % dir for logs

fopts.filepath=[dir_data,'\d'];     %filepath to data
file_id=1:1885;  % file ids to analyse

fopts.xlim=[-40e-3, 35e-3];
fopts.ylim=[-15e-3, 20e-3];

fopts.t_0=1.9977;

fopts.dt=5e-3;
fopts.n_pulse=250;

%general flags
fopts.min_count=3000;

fopts.log=false;
fopts.graphics=false;

%SHUNT BOUNDARIES
param_boundary=[3.4,0.25,0,0.75];   % shunt boundaries [Q_i,Q_f,Sh_i,Sh_f]

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
t=fopts.t_0+fopts.dt*[0:fopts.n_pulse-1];

%% Get shunt params from logs
dir_orig=cd(dir_log);
log_files=dir('output_write*.mat');

idx_params=zeros(length(log_files),1);
shunt_params=cell(length(log_files),1);

for ii=1:length(log_files)
    [idx_params(ii),shunt_params{ii}]=get_params_from_log(log_files(ii).name);
end

cd(dir_orig);

%% Summary
%%ML summary
hfig_ml=figure();
hold on;
scatter(file_id(cost<100),cost(cost<100),10,'o',...
    'filled');   % don't plot bad costs
box on;
titlestr='Machine learning for optimal trap shunting';
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
titlestr=sprintf('Best shunt: shot %d',file_id(best_id));
title(titlestr);
xlabel('Time (s)');
ylabel('Position (m)');
axis tight;
legend({'X','Y','Z'});

%plot the profile
hfig_best_profile=figure();
best_profile=shunt_params{idx_params==file_id(best_id)};
plot_shunt_profile(best_profile,param_boundary,hfig_best_profile);

titlestr=sprintf('Best shunt: shot %d',file_id(best_id));
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
    titlestr=sprintf('Exponential shunt');
    title(titlestr);
    xlabel('Time (s)');
    ylabel('Position (m)');
    axis tight;
    legend({'X','Y','Z'});
end

%plot the profile
hfig_exp_profile=figure();
exp_profile=shunt_params{idx_params==1};
plot_shunt_profile(exp_profile,param_boundary,hfig_exp_profile);

titlestr=sprintf('Exponential shunt');
title(titlestr);
xlabel('Profile segment');
ylabel('DAQ voltage (V)');