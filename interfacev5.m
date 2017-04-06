%bryce:what needs to be in the path for this to run properly
%explain briefly what this program does and what the user will need to edit

%to do
%could read in number of itterations from settings file in order to be a
%bit intelegent

%NEEDS

%safety feature if path can't be found

%passes from labview
%mloop,boolean purpose: do you want to interface with mloop or use matlab to scan over variables
%file,string purpose: tells program where to look for setting files
%i
%Start user settings

%read from config file
% i=1; %DEBUG DO NOT LEAVE IN ---------------------------------------------------------

path_user_config='..\MloopUserConfig.txt';  % this is the user config file regardless of Mloop
user_config_fp=fopen(path_user_config,'r');  % read config file

while true
  this_line = fgetl(user_config_fp);
  if ~ischar(this_line); break; end  %end of file
  eval(this_line);
end
fclose(user_config_fp);
%input var from MloopUserConfig
%lower_bounds
%upper_bounds
%first  1 x n array of first values to use
%param_select 1 x n array of logicals(as doubles) if this variable should be used
%trust_region double, trust region between 0 and 1
%max_wo_better

%mloop=false;
%param_select,array purpose: tells the program which parameters to use, out of the list of variables in the paths array

save(['.\logs\interfacev5',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat'])

%List of all paths of control variables
paths = {
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',96},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',95},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',97},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',96},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',98},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',97},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',99},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',98},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',100},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',99},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',101},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',100},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',102},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',101},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',65},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',64},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',66},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',65},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',67},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',66},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',68},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',67},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',69},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',68},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',70},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',69},'Final value, Amplitude (exp)',''}},...
    {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',71},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',70},'Final value, Amplitude (exp)',''}}
    };      % 14-param shunting profile


% Previous settings
% paths = {{{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',96},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',95},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',97},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',96},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',98},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',97},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',99},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',98},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',100},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',99},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',101},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',100},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',102},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',101},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',103},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',102},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',104},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',103},'Final value, Amplitude (exp)',''}},...%shunt
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',65},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',64},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',66},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',65},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',67},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',66},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',68},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',67},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',69},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',68},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',70},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',69},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',71},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',70},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',72},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',71},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',73},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',72},'Final value, Amplitude (exp)',''}}
%         };
    
% paths = {
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',98},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',97},'Final value, Amplitude (exp)',''}},...
%          {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',99},'Initial Value ',''},...
%          {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',98},'Final value, Amplitude (exp)',''}},...
%         };  
%     
%     
    
%         {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',6},'Initial Value ',''},...
%         {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',2},'Initial Value ',''},...
%         {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array 6070,2 (Ch0,1)',{'<Cluster>',2},'Final value, Amplitude (exp)',''},...
%         {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config 6070E Device3 (Ch 0,1) ',{'<Cluster>',2},'Dwell time (ms,digital waveform)',''},...
%         {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',2},'Tau (ms, exp)',''},...
%         {{'<Cluster>',1},{'<Cluster>',4},'Digital Waveform Array (Channels 0-15 on 6289M) ',{'<Cluster>',2},'Pulse Stop',''},...
%         {{'<pCluster>',1},{'<Cluster>',4},'Digital Waveform Array (Channels 0-17 on 6733)',{'<Cluster>',2},'Pulse Start',''}};


%%%%NOTE: QUAD FIRST THEN SHUNT!!!! %%%%
%use a path of names to change variables, with the value you want to change
%the variable too at the end of the path
%eg path = {path to variable, value to change too};
%can also conjoin variables so they are set to the same parameter, simply
%put their paths into a cell array in athe paths master list
quad_lims = [0.14,3.4];
%shunt_lims = [0,0.75];
shunt_lims = [0,1];
param_limits=[repmat(quad_lims,7,1);repmat(shunt_lims,7,1)];    % NOTE: hard-coded parameter limits

%%Independent exponential ramp scan for various time constants (2D)
if i==1 && ~mloop && control
    %shunt boundaries (start and end points for exponential ramp design)
    Vquad_boundary=[3.4,0.14];
    Vshunt_boundary=[0,0.87];
    
    ntau_lim=[2,10];        %param lims - number of exponential constants
    n_search=40;            %interp points in lim to search for param
    n_shot_avg=1;          %number of shots to take per param set
    n_interp=7;             %14[7/7] param shunt
    ntau=linspace(ntau_lim(1),ntau_lim(2),n_search);    %single-dim param in linear space
    
    %create all permutations of param_values for experiment
    ntau_perm=transpose(combvec(ntau,ntau));
    ntau_perm=repmat(ntau_perm,[n_shot_avg,1]);         %averaging shots
    ntau_perm=ntau_perm(random(size(ntau_perm,1)),:);   %randomise the ntau
    
    %build param_values
    n_total_perm=size(ntau_perm,1);
    param_values=zeros(n_total_perm,n_interp*2);    %n_total_perm X 14 (7/7 ramp)
    for ii=1:n_total_perm
        param_values(ii,:)=[exp_ramp(Vquad_boundary(1),Vquad_boundary(2),n_interp,ntau_perm(ii,1)),...
            exp_ramp(Vshunt_boundary(1),Vshunt_boundary(2),n_interp,ntau_perm(ii,2))];
    end
%     param_values=repmat(param_values,[n_shot_avg,1]);   %averaging shots
%     param_values=param_values(randperm(size(param_values,1)),:);    %randomise param sets in run order
    
    %save params set and configs to disk
    %param_data is crucial since it is called every iteration by a new instance of interface script
    %but the permutation and ordering should remain fixed for current analysis
    filename_param_log=['param_data_',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat'];
    save(filename_param_log,'param_values','ntau_perm',...
        'ntau_lim','n_search','n_shot_avg','n_interp','Vquad_boundary','Vshunt_boundary');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETER SCANNER HINTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Parameter scan
% %for param scan only create an array of the values for each iteration

% %Combining variables
% %option one 
% %product
%a=[1 2]
%b=[5 6]
%transpose(combvec(a,b))
%     1     5
%     2     5
%     1     6
%     2     6
% %which maps through all posible combinations

% %Option two list
% %steps though both in sequenc
% %requires both lists to be same length !
%transpose([a; b])
%     1     5
%     2     6

%val1=linspace(0.6,3.4,30);
%val2=linspace(0.6,3.4,30);
%param_values=transpose([val1; val2]);
%
%param_values=param_values(randperm(size(param_values,1)),:);

% END PARAMETER SCANNER HINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END user settings

dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file
f_log=fopen(path_log,'a');  % append to log-file
fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : started with itteration=%u mloop=%u file=%s param_select=[%s] .\n'],i,mloop,file,num2str(param_select));
fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : param_limits size %u  size(paths) %u .\n'],size(param_limits,1),size(paths,2));
fclose(f_log);

%% read in the xml settings file
new_path = file;
file_id = fopen(file,'r');
j=1;
line = fgetl(file_id);
A{j} = line;
while ischar(line)
    j = j+1;
    line = fgetl(file_id);
    A{j} = line;
end
fclose(file_id);

% delete(file);   % DEBUG - delete xml file

%% make adjustments
if mloop
    %Check if M-LOOP has written out input file
    check=0;
    file = 'C:\Users\BEC Machine\Dropbox\labview control code\mloop_files\exp_input.txt';
    while check ~= 2
        check = exist(file,'file');
        pause(0.5);
    end
    %read in the mloop parameter file and then delete it
    file_id = fopen(file,'r');
    line = fgetl(file_id);
    fclose(file_id);
    param = str2num(line(14:end-1));
	
% 	copyfile(file,'C:\Users\BEC Machine\Dropbox\debug\exp_input_cp.txt');	% DEBUG
	
    delete(file);
    
    f_log=fopen(path_log,'a');  % append to log-file
    fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : Read in params from mloop=[%s] \n'],num2str(param));
    fclose(f_log);
    
    %which parameters to control
    selecter = find(param_select); 
    
    %selects the path for the different parameters
    paths = paths(selecter);

    %select the hard limits of those parameters
    param_limits = param_limits(selecter,:);
    
else
    %set params to current param values
    if control ~=0
        %load the mat file with predefine parameter sets
        data_param=load('param_data.mat');
        param_values=data_param.param_values;   %get the parameter set array
        num_param_set=size(param_values,1);     %number of parameter sets
        
        try
            idx_param=mod(i-1,num_param_set)+1;     %param idx to run - loops
            param=param_values(idx_param,:);
            file_name = 'param_log.txt';
            if i==1
                file_pointer = fopen(file_name,'w');
            else
                file_pointer = fopen(file_name,'a');
            end
            fprintf(file_pointer,[mat2str(param),'\n']);
            fclose(file_pointer);
        catch
            param=param_values(end,:);
        end
    end

    %if we do not wish for anything to be change the path file is set as
    %empty
    if control==0
        paths = {};
    end
end

conjoined_indx = [];

for idx=1:numel(paths)
    
    %check if param values are within the hard coded safety range
    if param(idx)>param_limits(idx,2)
        param(idx) = param_limits(idx,2);
    elseif param(idx)<param_limits(idx,1)
        param(idx) = param_limits(idx,1);
    end
    
    %place the variable inputs into the paths to be changed
    if iscell(paths{idx}{end})
        % variables are bunched (i.e. conjoined)
        conjoined_indx = [conjoined_indx, idx];
        for j = 1:numel(paths{idx})
            paths{idx}{j}{end} = num2str(param(idx));
        end
    else
        % solitary variable
        paths{idx}{end} = num2str(param(idx));
    end
end  

%split up conjoined variables
if numel(conjoined_indx)>0
    temp = paths;
    paths = {};
    r=1;
    for j = 1:numel(temp)
        if any(j==conjoined_indx)
            for k = 1:numel(temp{j})
                paths{r} = temp{j}{k};
                r=r+1;
            end
        else
            paths{r} = temp{j};
            r=r+1;
        end
    end
end

read_paths={};

%if we have multiple paths we can steup a looping structure as follows
%paths = {all the different paths}
%for j = 1:length(paths)
%then replace path with paths{j}
if i == 1
    exact_line =zeros(1,numel(paths));
    comment_line = zeros(1,numel(read_paths));
end
for k = 1:numel(paths)
    %bryce: what happens when the path string cant be found? is there an
    %error message?
    
    %select the path we want to change
    path = paths{k};
    
    %we only wish to find the exact path once, as it is computationally
    %expensive
    if i==1
        cline = 1;
        for j = 1:length(path)-1
            count = 0;
            toggle = 1; 
            test = 0;
            if iscell(path{j})
                name = path{j}{1};
                num = path{j}{2};
                end_name = strrep(name,'<','</');
            elseif ischar(path{j})
                name = path{j};
                num = 1;
                end_name = '';
            end
            while count<num
                cline = cline + 1;
                if ~isempty(strfind(A{1,cline},name))
                    if test == 0 && count>0
                        toggle = 1 - toggle;
                    end
                    if toggle
                        count = count + 1;
                        toggle = 1 - toggle;
                    end
                    test = test + 1;
                elseif ~isempty(strfind(A{1,cline},end_name))
                    test = test - 1;
                end
            end
        end
        exact_line(1,k) = cline+1;
    end

    %set the value to what we want it to be
    start_ps=strfind(A{exact_line(1,k)},'>');
    end_ps= strfind(A{exact_line(1,k)},'<');
    strt = A{exact_line(1,k)}(1:start_ps(1));
    fin = A{exact_line(1,k)}(end_ps(2):end);
    A{exact_line(1,k)} = strcat(strt,path{end},fin);
end

if i==1
    for k = 1:numel(read_paths)
        %bryce: what happens when the path string cant be found? is there an
        %error message?
        %select the path we want to change
        path = paths{k};
        %we only wish to find the exact path once, as it is computationally
        %expensive

        cline = 1;
        for j = 1:length(path)-1
            count = 0;
            toggle = 1; 
            test = 0;
            if iscell(path{j})
                name = path{j}{1};
                num = path{j}{2};
                end_name = strrep(name,'<','</');
            elseif ischar(path{j})
                name = path{j};
                num = 1;
                end_name = '';
            end
            while count<num
                cline = cline + 1;
                if ~isempty(strfind(A{1,cline},name))
                    if test == 0 && count>0
                        toggle = 1 - toggle;
                    end
                    if toggle
                        count = count + 1;
                        toggle = 1 - toggle;
                    end
                    test = test + 1;
                elseif ~isempty(strfind(A{1,cline},end_name))
                    test = test - 1;
                end
            end
        end
        comment_line(1,k) = cline+1;

        %or use the exact line is already known
        % exact_line = {};

        %print out what the comment says
        j=0;
        end_of_comment = 1;
        while 0
            %display this line somehow
            %A{comment_line(1,k)+j}
            if isempty(strfind(A{comment_line(1,k)+j},'</Val>'))
                j = j+1;
            else
                end_of_comment = 0; 
            end
        end
    end
end
%% write out new file
file_id = fopen(new_path,'w');
for j = 1:numel(A)
    if A{j+1} == -1
        fprintf(file_id,'%s', A{j});
        break
    else
        fprintf(file_id,'%s\n', A{j});
    end
end
fclose(file_id);

%%clean up
clear('A')

f_log=fopen(path_log,'a');  % append to log-file
fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : finished. \n']);
fclose(f_log);

% copyfile(new_path,'C:\Users\BEC Machine\Dropbox\debug\new_path_cp');    % DEBUG