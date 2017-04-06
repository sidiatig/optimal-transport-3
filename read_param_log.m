% Read param set from param_log.txt generated from param scanner code in
% interfacev5.m

function param_array = read_param_log(filepath,n_params)
%check file exists
if exist(filepath,'file')~=2
    error('File %s does not exist.\nEnter a valid path to log file.',filepath);
end

%read file
fid=fopen(filepath,'r');

%text search pattern
pattern_form=strcat('[',repmat('%f',[1,n_params]),']');  %how the param_log writer is formatted
data_cell=textscan(fid,pattern_form);

fclose(fid);  %close file

%create param array
param_array=cat(2,data_cell{:});