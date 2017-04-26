%%Summarise exponential profile trap relaxation
%06/04/2017 DK Shin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Configs
%experiment specific
% data dir
dir_data='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\optimal_transport\2d_exp_search\0_25s_20x20';
dir_log=[dir_data,'\log'];          % dir for logs (currently output_write mat files)
dir_output=[dir_data,'\output'];    % dir of output

fopts.filepath=[dir_data,'\d'];     %filepath to data
file_id=1:800;  % file ids to analyse

fopts.xlim=[-35e-3, 20e-3];
fopts.ylim=[-10e-3, 18e-3];

fopts.dt=10e-3;
fopts.t_0=0.9974;
fopts.n_pulse=250;

%%%pass/fail/penalty
%pulse picking success rate
fopts.win_capture_rate=0.5;     %pulse/window capture rate for pass/fail
%pkpk amplitude penalty
fopts.pkpk_penalty=false;
%low number penalty
fopts.num_win_penalty=[];       %minimum counts per window (avgd) to pass number penalty
%width penalty
fopts.width_sat=[8e-3,5e-3,6e-3];
fopts.penalty_width_sat=[];      %set to [] to turn penalty OFF

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
% for ii=1:n_shots
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

if isfield(PARAMLOG,'ntau_perm')
ntau=PARAMLOG.ntau_perm;
else
    % we suffer here since I forgot to save ntaus with the params for first
    % run
    Vquad_0=param_boundary(1);          %initial set point
    Vquad_inf=param_boundary(2);        %asymptotic values of exponential ramp (boundary condition)
    Vshunt_0=param_boundary(3);
    Vshunt_inf=param_boundary(4);
    
    ntau=zeros(n_par_set,2);   %ntau to calculate
    for ii=1:n_par_set
        %get profile
        Vquad_this=param_values(ii,1:7);
        Vshunt_this=param_values(ii,8:end);
        
        %%get number of exp constants decayed
        %quad
        logv_temp=log(Vquad_this-Vquad_inf);
        ntau(ii,1)=log(Vquad_0-Vquad_inf)-logv_temp(end);
        
        %shunt
        logv_temp=log(Vshunt_inf-Vshunt_this);
        ntau(ii,2)=log(Vshunt_inf-Vshunt_0)-logv_temp(end);
    end
end
%% Summary
n_interp=40;    %for create grid data for cost lanscape/mesh plotting

%%EXP Summary
%all costs
cost_not_bad=100;       % cost for "bad" run
id_except_bad=file_id(cost<cost_not_bad);
[NTAU1q,NTAU2q,COST]=array2mesh(ntau(id_except_bad,1),ntau(id_except_bad,2),cost(id_except_bad),n_interp);

hfig_exp_not_bad=figure();
hold on;
im_not_bad=surf(NTAU1q,NTAU2q,COST,...
    'CDataMapping','scaled','LineStyle','none');
shading interp;
colormap(hfig_exp_not_bad,viridis);
scatter3(ntau(id_except_bad,1),ntau(id_except_bad,2),cost(id_except_bad),10,'o','filled');
colorbar;
box on;
axis tight;
titlestr='EXP: all shots except bad';
title(titlestr);
xlabel('$T/\tau_{Q}$');
ylabel('$T/\tau_{S}$');
zlabel('Cost (m)');
view(2);

%low costs / good shots only
cost_good=30e-3;
id_good=file_id(cost<cost_good);
[NTAU1q_good,NTAU2q_good,COST_good]=array2mesh(ntau(id_good,1),ntau(id_good,2),cost(id_good),n_interp);

hfig_exp_good=figure();
hold on;
im_good=surf(NTAU1q_good,NTAU2q_good,COST_good,...
    'CDataMapping','scaled','LineStyle','none');
shading interp;
colormap(hfig_exp_good,viridis);
scatter3(ntau(id_good,1),ntau(id_good,2),cost(id_good),10,'o','filled');
colorbar;
box on;
axis tight;
titlestr=sprintf('EXP: cost$<$%0.3g m',cost_good);
title(titlestr);
xlabel('$T/\tau_{Q}$');
ylabel('$T/\tau_{S}$');
zlabel('Cost (m)');
view(2);

%very good shots
cost_vgood=10e-3;
id_vgood=file_id(cost<cost_vgood);
[NTAU1q_vgood,NTAU2q_vgood,COST_vgood]=array2mesh(ntau(id_vgood,1),ntau(id_vgood,2),cost(id_vgood),n_interp);

hfig_exp_vgood=figure();
hold on;
im_vgood=surf(NTAU1q_vgood,NTAU2q_vgood,COST_vgood,...
    'CDataMapping','scaled','LineStyle','none');
shading interp;
colormap(hfig_exp_vgood,viridis);
scatter3(ntau(id_vgood,1),ntau(id_vgood,2),cost(id_vgood),10,'o','filled');
colorbar;
box on;
axis tight;
titlestr=sprintf('EXP: cost$<$%0.3g m',cost_vgood);
title(titlestr);
xlabel('$T/\tau_{Q}$');
ylabel('$T/\tau_{S}$');
zlabel('Cost (m)');
view(2);

%%report on BEST file
[best_cost,best_id]=min(cost);
best_osc_std=output{best_id}.osc_std;
best_osc=output{best_id}.osc;

best_ntau=ntau(best_id,:);

hfig_best=figure();
hold on;
for ii=1:3
    plot(t,best_osc(:,ii)-mean(best_osc(:,ii),'omitnan'),...
        'Linewidth',1.5);
end
box on;
titlestr=sprintf('Best EXP [%d]: $T/\\tau$=[%0.3g,%0.3g]: cost=%0.3g',best_id,best_ntau(1),best_ntau(2),best_cost);
title(titlestr);
xlabel('Time (s)');
ylabel('Position (m)');
axis tight;
legend({'X','Y','Z'});

%plot the best EXP profile
hfig_best_profile=figure();
best_profile=param_values(best_id,:);
plot_shunt_profile(best_profile,param_boundary,hfig_best_profile);
box on;
titlestr=sprintf('Best EXP [%d]: $T/\\tau$=[%0.3g,%0.3g]: cost=%0.3g',best_id,best_ntau(1),best_ntau(2),best_cost);
title(titlestr);
xlabel('Profile segment');
ylabel('DAQ voltage (V)');

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
titlestr=sprintf('Best EXP [%d]: $T/\\tau$=[%0.3g,%0.3g]: cost=%0.3g',best_id,best_ntau(1),best_ntau(2),best_cost);
title(titlestr);
xlabel('Time (s)');
ylabel('PAL rms width (m)');
axis tight;
legend({'X','Y','Z'});

%% Save output
%output directory
if ~isdir(dir_output)
    mkdir(dir_output);
end

%%Figures PNG
saveas(hfig_exp_not_bad,[dir_output,'\','ml_progress_not_bad.png']);       % ML progress - all costs except bad
saveas(hfig_exp_good,[dir_output,'\','ml_progress_good.png']);  % ML progress - low costs
saveas(hfig_exp_vgood,[dir_output,'\','ml_progress_vgood.png']);  % ML progress - very low cost
saveas(hfig_best,[dir_output,'\','best_oscillation.png']);        % best shot - oscillations
saveas(hfig_best_profile,[dir_output,'\','best_profile.png']);    % best shot - shunt profile
saveas(hfig_best_breathing,[dir_output,'\','best_breathing.png']);  	% best shot - rms width oscillation

%%Figures .fig
saveas(hfig_exp_not_bad,[dir_output,'\','ml_progress_not_bad.fig']);       % ML progress - all costs except bad
saveas(hfig_exp_good,[dir_output,'\','ml_progress_good.fig']);  % ML progress - low costs
saveas(hfig_exp_vgood,[dir_output,'\','ml_progress_vgood.fig']);  % ML progress - very low cost
saveas(hfig_best,[dir_output,'\','best_oscillation.fig']);        % best shot - oscillations
saveas(hfig_best_profile,[dir_output,'\','best_profile.fig']);    % best shot - shunt profile
saveas(hfig_best_breathing,[dir_output,'\','best_breathing.fig']);  	% best shot - rms width oscillation

%%Data
save([dir_output,'\','ml_summary_',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat']);  %save all vars