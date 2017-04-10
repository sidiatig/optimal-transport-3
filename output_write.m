%write out the result of the experiment for m-loop

%% USER CONFIGS
dir_data='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';    % default data dir
fopts.filepath=[dir_data,'\d'];     %filepath to data

fopts.xlim=[-35e-3, 20e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
fopts.ylim=[-10e-3, 18e-3];

fopts.t_0=1.1074;

fopts.dt=10e-3;
fopts.n_pulse=280;

%general flags
% fopts.min_count=1000;
fopts.num_in_win=30;        % minimum counts per window (avgd) to pass number penalty
fopts.min_window_pass=10;   % minimum good windows to pass for oscillation analysis

fopts.log=true;
fopts.graphics=true;
%%END USER CONFIGS

%% MAIN
%passed from labviews
%count
%i
%max_wo_better
%end_num

%returns
%bcost
%end_check
%count
dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file

%check which file number in dld output to read, this case strucutre is to
%avoid the problem when the iteration number is set back to zero on the
%last run
if i==0
    startfile = end_num;
else
    startfile = i;
end

%calcualte cost and uncertainty
%in future 'slosh' should be replaced with cost choice, so that the cost
%function can be easily selected by the user
% [cost unc bad] = cost_calculator('slosh',startfile);
[cost, unc, bad] = cost_calculator('slosh',startfile,fopts);

%writes out the file
file_name = 'C:\Users\BEC Machine\Dropbox\labview control code\mloop_files\exp_output.txt';
file = fopen(file_name,'w');

fprintf(file,'cost = %f\n',cost);
fprintf(file,'uncer = %f\n',unc);
fprintf(file,'bad = %s',bad);

fclose(file);

end_check = 0;

save(['.\logs\output_write',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat'])

f_log=fopen(path_log,'a');  % append to log-file
fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' output_write   : finished. \n']);
fclose(f_log);