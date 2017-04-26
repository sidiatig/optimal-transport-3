%write out the result of the experiment for m-loop

%% USER CONFIGS
dir_data='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';    % default data dir
fopts.filepath=[dir_data,'\d'];     %filepath to data


%thes options should not be defined in the script but rather a seprate config file
fopts.xlim=[-35e-3, 25e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
fopts.ylim=[-15e-3, 20e-3];

fopts.dt=5e-3;
fopts.t_0=2.06775;
fopts.n_pulse=300;

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
[cost, unc, bad] = cost_calculator('slosh',startfile,fopts);

%writes out the file
file_name = 'C:\Users\BEC Machine\Dropbox\labview control code\mloop_files\exp_output.txt';
file = fopen(file_name,'w');

fprintf(file,'cost = %f\n',cost);
%fprintf(file,'uncer = %f\n',unc); %turn off output uncert
fprintf(file,'bad = %s',bad);

fclose(file);

end_check = 0;

save(['.\logs\output_write',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat'])

f_log=fopen(path_log,'a');  % append to log-file
fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' output_write   : finished. \n']);
fclose(f_log);