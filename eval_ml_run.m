%%Summarise ML run
%30/03/2017 DK Shin

%%configs
%experiment specific
fopts.filepath='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d';  %filepath to data
file_id=1:1000;  % file ids to analyse

fopts.xlim=[-40e-3, 35e-3];
fopts.ylim=[-15e-3, 20e-3];

fopts.t_0=1.9977;

fopts.dt=5e-3;
fopts.n_pulse=250;

%general flags
fopts.min_count=3000;

fopts.log=false;
fopts.graphics=false;


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