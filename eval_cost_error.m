%Script to evaluate cost variability

dir_data='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';    % data dir
dir_output=[dir_data,'\output'];    % dir of output

fopts.filepath=[dir_data,'\d'];     %filepath to data
file_id=1:7;  % file ids to analyse

fopts.xlim=[-35e-3, 35e-3];
fopts.ylim=[-35e-3, 35e-3];

fopts.t_0=2.6075;

fopts.dt=10e-3;
fopts.n_pulse=110;

%general flags
% fopts.min_count=1000;
fopts.num_in_win=30;        % minimum counts per window (avgd) to pass number penalty
fopts.min_window_pass=10;   % minimum good windows to pass for oscillation analysis

fopts.log=false;
fopts.graphics=false;

%% MAIN
num_files=length(file_id);
cost=zeros(1,num_files);
for i=1:num_files
    [cost_total,unc_total,bad,output,config_out]=cost_calculator('slosh',file_id(i),fopts);
    cost(i)=cost_total;
end

%statistics
cost_avg=mean(cost);
cost_std=std(cost);

disp(sprintf('cost=%0.3g ± %0.3g\n',cost_avg,cost_std));

%plot
figure();
plot(file_id,cost,'o');
xlabel('shot #');
ylabel('cost');

