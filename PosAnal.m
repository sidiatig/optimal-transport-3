function[out,err,frames]= PosAnal(filepath,filenumstart,files,windows,plots,rot_angle,isverbose)
%this function imports data,windows it and returns the avg pos data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% * out: NxM array: [t_window_cen, xavg, yavg, tavg, xsd, ysd, tsd,
%                   numcounts]
% * err: boolean flag triggered True when any error is encountered in
% script detrimental to callers such as "cost_calculator.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TO BE IMPROVED/FIXED

err=false;      % initilialise err flag
dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file

%Start User Options
%acceptable counts per file to still include as window
numcut=10;
min_count_low_file=1000;    % count threshold for lowcountfile categorisation
%End User Options

%this measures the position of a pulsed atom laser in x y t and plots a
%picture

xavg=NaN;
yavg=NaN;
tavg=NaN;
xsd=NaN;
ysd=NaN;
tsd=NaN;
numcounts=NaN;
t_window_cen_int=NaN;
% out=[t_window_cen_int xavg yavg tavg xsd ysd tsd numcounts];
    
lowcountfiles=0;
missingfiles=0;
filestotxy=0;
out=[];     %this is a dynamic array but it prevents having to "clean" it at the end

if isverbose
    disp('importing data')
    parfor_progress(files);
    tic
end

three_channel_multifile_output=cell(files);

for n=1:files
    current_file_str = num2str(filenumstart+n-1);
    
    % Check if TXY-file needs to be generated
    if ~fileExists([filepath,'_txy_forc',current_file_str,'.txt']) &&...
            fileExists([filepath,current_file_str,'.txt'])
    	%disp(['converting file # ',num2str(filenumstart+n-1)])
        filestotxy=filestotxy+1;
        % Generate TXY-file from raw DLD file
    	dld_raw_to_txy(filepath,filenumstart+n-1,filenumstart+n-1);
    end
    
    % If TXY-file is still missing
    if ~fileExists([filepath,'_txy_forc',current_file_str,'.txt'])
        missingfiles=missingfiles+1;
        %disp(['file "' filepath,current_file_str,'.txt" or ',filepath,'_txy_forc',current_file_str,'.txt not here will not output data'])
    else
        % TXY-file exists
        three_channel_output=importdata([filepath,'_txy_forc',current_file_str,'.txt']);
    
        if size(three_channel_output,1)<min_count_low_file
            %warning('selected shot has less than %d counts',min_count_low_file);
            lowcountfiles=lowcountfiles+1;
        else
            three_channel_output_rot=zeros(size(three_channel_output));
            sin_theta = sin(rot_angle);
            cos_theta = cos(rot_angle);
            three_channel_output_rot(:,1) = three_channel_output(:,1);
            three_channel_output_rot(:,2) = three_channel_output(:,2)*cos_theta...
                - three_channel_output(:,3)*sin_theta;
            three_channel_output_rot(:,3) = three_channel_output(:,2)*sin_theta...
                + three_channel_output(:,3)*cos_theta;
            three_channel_multifile_output{n} = three_channel_output_rot; 
        end
    end
    if isverbose
        parfor_progress;
    end
end

if isverbose
    parfor_progress(0);
    disp('finished importing data')
    disp([int2str(lowcountfiles),' files with low counts'])
    disp([int2str(missingfiles),' files missing'])
    disp([int2str(filestotxy),' files converted to txy'])
    toc
end

three_channel_multifile_output= vertcat(three_channel_multifile_output{:});

tic;    %time the masking

%here i premask spatialy and temporaly if the spatial mask is constant
%this gives a small speedup for the loop(inc this premasking)
%there is also premasking if window is not constant which looks at the
%min/max and winsows in there
premask=0;
if isverbose
    disp(['total counts :',int2str(size(three_channel_multifile_output,1))])
end
if size(three_channel_multifile_output,1)>=numcut
    % premask only files with counts at least numcut
    if all(windows(:,3)==windows(1,3))&&...
            all(windows(:,4)==windows(1,4))&&...
            all(windows(:,5)==windows(1,5))&&all(windows(:,6)==windows(1,6))
        %spatial and temporal mask
        premask=1;
        premaskval=three_channel_multifile_output(:,2)>windows(1,3) &...
            three_channel_multifile_output(:,2)<windows(1,4);
        premaskval=premaskval &  three_channel_multifile_output(:,3)>windows(1,5) &...
            three_channel_multifile_output(:,3)<windows(1,6);
        premaskval=premaskval & three_channel_multifile_output(:,1)>min(windows(:,1)) &...
            three_channel_multifile_output(:,1)<max(windows(:,2));
        three_channel_multifile_output=three_channel_multifile_output(premaskval,:);
    else
        premask=0;
        %here i mask only by the min and max in the 3d(xyt) volume
        %are therese two cases doing the same thing?
        premaskval=three_channel_multifile_output(:,2)>min(windows(:,3)) &...
            three_channel_multifile_output(:,2)<max(windows(:,4));
        premaskval=premaskval &  three_channel_multifile_output(:,3)>min(windows(:,5)) &...
            three_channel_multifile_output(:,3)<max(windows(:,6));
        premaskval=premaskval & three_channel_multifile_output(:,1)>min(windows(:,1)) &...
            three_channel_multifile_output(:,1)<max(windows(:,2));
        three_channel_multifile_output=three_channel_multifile_output(premaskval,:);
    end
end
if isverbose
    disp(['total count # after spatial window',...
        num2str(size(three_channel_multifile_output,1))])
    disp('windowing data')
    parfor_progress(size(windows,1));
end

lowcountwindows=0;
t_window_cen=zeros(1,size(windows,1));
%counts=cell(size(windows,1),1);
for p=1:size(windows,1)
    tmin=windows(p,1);
    tmax=windows(p,2);
    t_window_cen(p)=(tmax+tmin)/2;
    xmin=windows(p,3);
    xmax=windows(p,4);
    ymin=windows(p,5);
    ymax=windows(p,6);
    
    if size(three_channel_multifile_output,1)>=numcut
        % window only files with counts at least numcut
        mask=three_channel_multifile_output(:,1)>tmin &...
            three_channel_multifile_output(:,1)<tmax;
        if premask==0 %doe this ever get called?
            mask=mask &  three_channel_multifile_output(:,2)>xmin &...
                three_channel_multifile_output(:,2)<xmax;
            mask=mask &  three_channel_multifile_output(:,3)>ymin &...
                three_channel_multifile_output(:,3)<ymax;
        end
        
        three_channel_output_rot_mask=three_channel_multifile_output(mask,:);
        
        numcounts=size(three_channel_output_rot_mask,1);    % number of counts captured in window
        if numcounts<numcut
            if isverbose
                warning('selected window has less than %d counts',numcut);
            end
            %plots=0;
            lowcountwindows=lowcountwindows+1;
            
            %all results about this window is set to NaN (except numcounts)
            xavg=NaN;
            yavg=NaN;
            tavg=NaN;
            xsd=NaN;
            ysd=NaN;
            tsd=NaN;
        else
            tavg=mean(three_channel_output_rot_mask(:,1));
            xavg=mean(three_channel_output_rot_mask(:,2));
            yavg=mean(three_channel_output_rot_mask(:,3));
            
            tsd=std(three_channel_output_rot_mask(:,1));
            xsd=std(three_channel_output_rot_mask(:,2));
            ysd=std(three_channel_output_rot_mask(:,3));
%             numcounts=size(three_channel_output_rot_mask(:,3),1);
            
%             out=[out;[t_window_cen(p) xavg yavg tavg xsd ysd tsd numcounts]];
        end
        
        %build output vector
        out=[out;[t_window_cen(p) xavg yavg tavg xsd ysd tsd numcounts]];
        
    else
        
        %plots=0;
    end
    
    %Plot results
    if plots
        x_bin_centers_temp = linspace(xmin,xmax, 200);
        y_bin_centers_temp = linspace(ymin,ymax, 200);
        [counts{p}, bin_centers] = hist3(three_channel_output_rot_mask(:,[2,3]),{x_bin_centers_temp, y_bin_centers_temp});
        x_bin_centers(p,:) = bin_centers{1};
        y_bin_centers(p,:) = bin_centers{2};
    end
    if isverbose
        parfor_progress;
    end
end

if isverbose
    toc
    parfor_progress(0);
    disp([int2str(size(windows,1)),' windows done. '])
    disp([int2str(lowcountwindows),' windows with less than 10 counts'])
end

three_channel_multifile_output=[];
if plots
    if isverbose
        disp('building plot frames')
        parfor_progress(size(windows,1));
    end
    tic
    
    hfig = figure(1);
    %set(gca,'Units','pixels');
    %imagesc(x_bin_centers(1,:),y_bin_centers(1,:), counts{1});%,[0 20]);
    %set(gca, 'xlimmode','manual','ylimmode','manual','zlimmode','manual','climmode','manual','alimmode','manual');
    %set(hfig,'visible','off');
    %framepar.resolution = [200,200];
    for n=1:size(windows,1)
        imagesc(x_bin_centers(n,:),y_bin_centers(n,:), counts{n});%,[0 20]);
        %colormap(gray)
        colormap(hot)
        %colorbar;
        xlabel('y(m)');
        ylabel('x(m)');
        title(['Time=',num2str(t_window_cen(n),'%.4f')]);
        frames(n)=getframe(hfig);
        %clf(hfig);
        if isverbose
            parfor_progress;
        end
    end
    
    if isverbose
        parfor_progress(0);
        toc
        disp('plot frames finished')
    end
    
end

%Clean up the out matrix by removing empty
%keeplist=[]; %this defines the rows that will be kept (eg have non zero #)
%for p=1:size(out,1)
%    if out(p,8)>=numcut
%        keeplist=[keeplist p];
%    end
%end
%only output non zero slices
%not sure if this is any better than just having a dynamic array
%out= out(keeplist,:)

if isverbose
    disp('done with PosAnal')
end

%% Error check
%%% Needs to be improved if there are going to be lots of ways to throw error
% Check for empty array - e.g. when all files had low counts
if isempty(out)
    err=true; %what does this do?
    % log error message
    % % BUG
    % %     f_log=fopen(path_log,'a');  % append to log-file
    % %     fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' PosAnal        : Error low counts. ',...
    % %         'path=%s, filenum=%d, files=%d, plots=%d, rot_angle=%d,\n'],filepath,filenumstart,files,plots,rot_angle);
    % %     fclose(f_log);
end

end