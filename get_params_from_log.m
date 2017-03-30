% extract shot number and params from output_write log
function [iter_num, params]=get_params_from_log(filename)
    S=load(filename,'i','line');
    
    % get params
    iter_num=S.i;
    
    params_str=strsplit(S.line(15:end-2),',');
    params=str2double(params_str);
end