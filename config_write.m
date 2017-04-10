%write out configuration for m-loop
%save('C:\Users\BEC Machine\Documents\labview control code\mloop_files\configvar.mat','max_wo_better','max_num','bounds','first','trust','mloop','param_select')

%input var from labview
%mloop double, set to 1 for mloop control
%max_num
%max_num=100

%input var from MloopUserConfig
%lower_bounds
%upper_bounds
%first  1 x n array of first values to use
%param_select 1 x n array of logicals(as doubles) if this variable should be used
%trust_region double, trust region between 0 and 1
%max_wo_better


dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file
f_log=fopen(path_log,'a');  % append to log-file

fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' config_write   : started. \n']);


path_user_config='..\MloopUserConfig.txt';		% TODO - assumes currently working dir (dir_this) but should work fine
user_config_fp=fopen(path_user_config,'r');  % read config file

%run all lines in Config files on MATLAB
while true
  this_line = fgetl(user_config_fp);
  if ~ischar(this_line)
      break; 
  end  %end of file
  eval(this_line);
end
fclose(user_config_fp);

fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' config_write   : Params read from MloopUserConfig. \n']);

if ~isequal(size(lower_bounds,2),size(upper_bounds,2),size(first,2))
    fprintf(f_log,[' config_write   : Param Dimentions do not match. \n']);
end

save(['.\logs\config_write',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat'])

%delete the param files so dont read in anything that is there from prev runs
delete('C:\Users\BEC Machine\Dropbox\labview control code\mloop_files\exp_input.txt');

file_name = 'C:\Users\BEC Machine\Dropbox\labview control code\mloop_files\exp_config.txt';
file = fopen(file_name,'w');

%which parameters to control
selecter = find(param_select);

%number of parameters to control
N = sum(param_select);
fprintf(file,'num_params = %d\n',N);

%vectors of the bounds for all parameters
lower_bounds=lower_bounds(logical(param_select));
lower_bounds_str=strjoin(arrayfun(@(x) num2str(x),lower_bounds,'UniformOutput',false),',');
fprintf(file,strcat('min_boundary = [',lower_bounds_str,']\n'));

upper_bounds=upper_bounds(logical(param_select));
upper_bounds_str=strjoin(arrayfun(@(x) num2str(x),upper_bounds,'UniformOutput',false),',');
fprintf(file,strcat('max_boundary = [',upper_bounds_str,']\n'));

%starting point for each parameter
first=first(logical(param_select));
first_str=strjoin(arrayfun(@(x) num2str(x),first,'UniformOutput',false),',');
param_select_str=strjoin(arrayfun(@(x) num2str(x),param_select,'UniformOutput',false),',');
fprintf(file,strcat('first_params = [',first_str,']\n'));

if ceil(trust_region) ~=1 || trust_region == 1
    trust_region = 0.45;
end
fprintf(file,'trust_region = %f\n',trust_region);

%Halting conditions
%have these ones optional
if max_num>0
    fprintf(file,'max_num_runs = %d\n',ceil(max_num)); 
end
if max_wo_better>0
    fprintf(file,'max_num_runs_without_better_params = %d\n',ceil(max_wo_better));
end
%target_cost = 0.01;                        %optimization halts when a cost below this target is found

fprintf(file,'cost_has_noise = True\n');

%no_delay = 'True';
%Timing options
if no_delay
    fprintf(file,'no_delay = True\n');
else
    fprintf(file,'no_delay = False\n');
end

%File format options
fprintf(file,'interface_file_type = ''txt''\n');
fprintf(file,'controller_archive_file_type = ''txt''\n');
fprintf(file,'learner_archive_file_type = ''txt''\n');

fprintf(file,'predict_global_minima_at_end  = True\n');
fprintf(file,'predict_local_minima_at_end = False\n');

fprintf(file,'visualizations = True');

fclose(file);



fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' config_write   : finished. \n']);
fclose(f_log);

%% noise_level = 0.1                      #initial noise level
%% update_hyperparameters = True          #whether noise level and lengths scales are updated
%% default_bad_cost = 10                  #default cost for bad run
%% default_bad_uncertainty = 1            #default uncertainty for bad run
