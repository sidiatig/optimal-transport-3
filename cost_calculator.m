function [cost_total,unc_total,bad,output,config_out] = cost_calculator(choice,startfile,fopts)

%Calculates the cost function for M-LOOP optimisation
% Updated
% 27-03-2017 - improved readability and update on cost function st all 3 dim var sum
% 30-03-2017 - TODO: pass options as argument to increase modularity
%            - fixed fft plotting - was plotting index vs amplitude
%% INPUT
%	choice: 'slosh'
%	startfilefile: file ID to analyse - file must be named 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d' + ID
%   fopts: struct with optional fields to configure analysis:
%       filepath: path to data file including data file TOKEN (i.e. '../d')
%       xlim: 1x2 array of x lims
%       ylim: 1x3 array of y lims
%       t_0: first pulse time of arrival
%       dt
%       n_pulse
%       min_count
%       log: log to file
%       graphics: do plots/figures/etc.
%
%% OUTPUT
%	cost_total: objective function output for this data
%	unc_total: cost uncertainty
%	bad: flag raised when encountered with any errorneous behaviour
%   output: all sorts of other things evaluated about the shot "hidden"
%   in the cost function
%   config_out: configs used in the analysis - i.e. fopts used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%slosh
%   output: struct with fields
%       osc: Nx3 array of X,Y,Z oscillations
%       osc_std: standard deviation of trap oscillation; 1x3 array
%       fft: 1x3 cell array of fft (fft format is Nx2 freq vs amplitude)
%       trap_freq: trap frequency; 1x3 array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file
%to improve
%time to wait should be user var

%% analysis
switch choice
    case 'slosh'
        %% Parse inputs
        % A bit long but - could make a cell of varnames<-->defaults and
        % check by loop
        if ~exist('fopts','var')
            %default options
            fopts.filepath='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d';
            fopts.xlim=[-40e-3,35e-3];
            fopts.ylim=[-15e-3,20e-3];
            fopts.t_0=1.9977;
            fopts.dt=5e-3;
            fopts.n_pulse=250;
%             fopts.min_count=3000;
            fopts.log=true;
            fopts.graphics=true;
            fopts.num_win_penalty=30;
            fopts.min_window_pass=10;
        else % missing fields in fopts
            if ~isfield(fopts,'filepath')
                fopts.filepath='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d';
            end
            if ~isfield(fopts,'xlim')
                fopts.xlim=[-40e-3,35e-3];
            end
            if ~isfield(fopts,'ylim')
                fopts.ylim=[-15e-3,20e-3];
            end
            if ~isfield(fopts,'t_0')
                fopts.t_0=1.9977;
            end
            if ~isfield(fopts,'dt')
                fopts.dt=5e-3;
            end
            if ~isfield(fopts,'n_pulse')
                fopts.n_pulse=250;
            end
%             if ~isfield(fopts,'min_count')
%                 fopts.min_count=3000;
%             end
            if ~isfield(fopts,'log')
                fopts.log=true;
            end
            if ~isfield(fopts,'graphics')
                fopts.graphics=true;
            end
            if ~isfield(fopts,'num_win_penalty')
                fopts.num_win_penalty=30;
            end
            if ~isfield(fopts,'min_window_pass')
                fopts.min_window_pass=10;
            end
        end
        
        %% START User Controls
		% NOTE: tweak this code for different profile/sequences
        files=1;
        %DLD OUTPUT
		filepath=fopts.filepath;
        check=0;
        file_check=strcat(filepath,'_txy_forc',num2str(startfile),'.txt');
        n=0;
        outerr=false;
        while (check ~= 2) && n<20 %should make this a user var
            %need to limit loop and send bad out
            check = exist(file_check,'file');
            pause(0.5);
            n=n+1;
        end
        if n==20
            outerr=true;
            if fopts.log
                f_log=fopen(path_log,'a');  % append to log-file
                fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' cost_calculator: Could not find txy_file. ',...
                    'choice=%s, startfile=%d \n'],choice,startfile);
                fclose(f_log);
            end
        end
		
		% window params - atom pulse picking        
        xmin=fopts.xlim(1);
        xmax=fopts.xlim(2);
        ymin=fopts.ylim(1);
        ymax=fopts.ylim(2);
        
        % RF atom laser settings: (n_pulse more than a few thousand takes forever)
        t_0=fopts.t_0;      % DLD resolved time (s) for the first pulse to analyse
        dt=fopts.dt;
        n_pulse=fopts.n_pulse;
        num_win_penalty=fopts.num_win_penalty;
        min_window_pass=fopts.min_window_pass;
        
%         % 3ms/530 pulses
%         dt=0.003;             % pulse-picking spacing in time (s)
%         n_pulse=500;          % number of windows/pulses to pick out

      % 5ms/318 pulse
%         dt=0.0050;
%         n_pulse=250;        
        
%         %10ms/120 pulses
%         dt=0.010;
%         n_pulse=110;        

        dw=dt;          % width of pulse picking window in time
        
		% pass/fail
%         min_count=fopts.min_count;      % minimum counts captured in analysis to pass

		% misc
		v_tof=9.81*0.416;	% z-velocity of atom at detection event
		dld_xy_rot=0.61;    % angle to rotate DLD XY coords
        
        %END User Controls
        
        %% slosh analysis
        %%initalize
        window=zeros(n_pulse,6);
        xyzAvgSD=zeros(files,6);
        
        %create the windows
        for n=1:n_pulse
            tmin=t_0+dt*n-dw/2;		% t edges are centered
            tmax=t_0+dt*n+dw/2;
            window(n,:)=[tmin tmax xmin xmax ymin ymax];
        end
        
        if ~outerr
            [outval,outerr]=PosAnal(filepath,startfile,files,window,false,dld_xy_rot,false);
        end
        
        if outerr
            cost_total = 100;
            unc_total = 100;
            bad='True';
            output=[];
            config_out=fopts;
            
%             %DEBUG - can improve PosAnal
%             size(outval)
%             disp('outerr || size of outval error');
        elseif sum(~isnan(outval(:,2)))<min_window_pass
            % not enough good windows for oscillation analysis - penalise with large cost
            cost_total = 1;
            unc_total = 0.005;
            bad='False';
            output=[];
            config_out=fopts;
            
        else
			% parse output from pulse analysis code
            t_cent=outval(:,1);		% T-centre of window
            x_avg=outval(:,2);		% avg X-coord of counts in window
            y_avg=outval(:,3);
            t_avg=outval(:,4);
            x_std=outval(:,5);		% std of X-coord of counts in window
            y_std=outval(:,6);
            t_std=outval(:,7);
            n_capt=outval(:,8);		% number of counts captured in the defined window
            
            numcounts=sum(n_capt);  % total number
            numpulses=size(t_cent,1);
            
            %number of counts
            output.num_in_win=n_capt;
            output.numcounts=numcounts;
            
            
            %frequency value below which peaks are ignored (in kHz)
            lower_freq_bound = 3;
            
            %% X-oscillations
            meanx=mean(x_avg,'omitnan');
            
            %fourier of average x position
            fftx=FFTxt(t_cent,x_avg(~isnan(x_avg))-meanx);
            
            %eliminate low level frequencies
            fftx_restricted = fftx(fftx(:,1)>lower_freq_bound,:);
            
            % Calculate cost factor from X
            output.osc_std(1)=std(x_avg,'omitnan');   % get oscillation
            costx=output.osc_std(1);
            %costx=max(fftx_restricted(:,2));
            
            penalty_x=penalty_pkpk(x_avg,0.9*(xmax-xmin));
            costx=costx+penalty_x;
            
            % See standard amplitude uncert in fft is sqrt(sum(unc(x_n)^2))/sqrt(num_points) https://stackoverflow.com/questions/27529166/calculate-uncertainty-in-fft-amplitude
            
            %uncx = costx*mean(1./outval(:,2).*outval(:,5)./sqrt(outval(:,8)));
            uncx = sqrt(sum((x_std./sqrt(n_capt)).^2,'omitnan'))/numpulses;
            uncx=uncx+0.05*costx;
            
            %% Y-oscillations
            meany=mean(y_avg,'omitnan');
            %fourier of average y position
            ffty=FFTxt(t_cent,y_avg(~isnan(y_avg))-meany);
            %eliminate low level frequencies
            ffty_restricted = ffty(ffty(:,1)>lower_freq_bound,:);
            
            % Calculate cost factor from Y
            output.osc_std(2)=std(y_avg,'omitnan');   % get oscillation
            %costy=std(y_avg)*sqrt(2); 	% TODO why is Y-axis oscillation divided by sqrt 2?
            costy=output.osc_std(2);
            %costy=max(ffty_restricted);
            
            penalty_y=penalty_pkpk(y_avg,0.9*(ymax-ymin));
            costy=costy+penalty_y;
            
            %uncy = costy*mean(1./outval(:,3).*outval(:,6)./sqrt(outval(:,8)));
            uncy = sqrt(sum((y_std./sqrt(n_capt)).^2,'omitnan'))/numpulses;
            uncy=uncy+0.05*costy;
            
            %% Z-oscillations
            % T - Z conversion
            dz_avg=v_tof*(t_avg-t_cent);		% z-deviation from window centre
            dz_std=v_tof*t_std;
            
            meanz=mean(dz_avg,'omitnan');				% average z-deviation from window centre
            
            % Fourier of avg dZ
            fftz=FFTxt(t_cent,dz_avg(~isnan(dz_avg))-meanz);
            
            %eliminate low level frequencies
            fftz_restricted = fftz(fftz(:,1)>lower_freq_bound,:);
            
            % Calculate cost factor from dZ
            output.osc_std(3)=std(dz_avg,'omitnan');   % get oscillation
            costz=output.osc_std(3);	% cost measure from dZ
            
            penalty_z=penalty_pkpk(dz_avg,0.95*dw*v_tof);
            
            costz=costz+penalty_z;
            
            uncz = sqrt(sum((dz_std./sqrt(n_capt)).^2,'omitnan'))/numpulses;
            uncz=uncz+0.05*costz;
            
            %% Get outputs
            output.osc=[x_avg,y_avg,dz_avg];
            output.fft={fftx_restricted,ffty_restricted,fftz_restricted};
            
            % get trap frequencies from fft
            for ii=1:3
                output.trap_freq(ii)=output.fft{ii}(output.fft{ii}(:,2)==max(output.fft{ii}(:,2)),1);
            end
            
            %% Evaluate cost function
            % X,Y,Z oscillations summed in quadrature
            cost_osc_total = sqrt(sum([costx,costy,costz].^2));
            unc_osc_total = sqrt(sum([uncx,uncy,uncz].^2));
            
            %check for nan oscillation cost - from a very BAD PosAnal results - all bad windows
            if isnan(cost_osc_total)
                cost_osc_total=200e-3;  %penalise as oscillations amplitude at 10cm scale
            end
            if isnan(unc_osc_total)
                unc_osc_total=0.05*cost_osc_total;
            end
            
            %penalty for low counts at detector
            %number of counts
            output.num_in_win=n_capt;
            output.numcounts=numcounts;
            
            %calculate new cost function
            cost_total=cost_osc_total+penalty_num(output.num_in_win,n_pulse*num_win_penalty);
            unc_total=unc_osc_total;
            
            bad = 'False';
            
            %% Plot: Fourier transform
            if fopts.graphics
                figure(1);
                clf;
                set(gcf,'Color',[1 1 1]);
                
                % X
                subplot(3,1,1);
                plot(fftx_restricted(:,1),fftx_restricted(:,2));
                hold on;
                xlim = get( gca, 'Xlim' );
                plot( xlim, [costx costx] )
                plot( xlim, [costx-uncx costx-uncx] )
                plot( xlim, [costx+uncx costx+uncx] )
                annotation('textbox',...
                    [0.6 0.75 0.3 0.15],...
                    'String',{['cost of file ',num2str(startfile)],[num2str(cost_total*10^3,'%5.3f'),'±',num2str(unc_total*10^3,'%5.3f'),'mm'],...
                    [num2str(numpulses),'w ',num2str(size(window,1)),'w ',num2str(numcounts),'c']},...
                    'FontSize',13,...
                    'EdgeColor',[1 1 1]);
                
                hold off
                xlabel('Freq (Hz)');
                ylabel('X_{amp} (m)')
                
                % Y
                subplot(3,1,2)
                plot(ffty_restricted(:,1),ffty_restricted(:,2));
                hold on
                xlim = get( gca, 'Xlim' );
                plot( xlim, [costy costy] )
                plot( xlim, [costy-uncy costy-uncy] )
                plot( xlim, [costy+uncy costy+uncy] )
                hold off
                xlabel('Freq (Hz)');
                ylabel('Y_{amp} (m)')
                
                % Z
                subplot(3,1,3);
                plot(fftz_restricted(:,1),fftz_restricted(:,2));
                hold on;
                xlim=get(gca,'Xlim');
                plot(xlim,costz*[1,1]);
                plot(xlim,(costz-uncz)*[1,1]);
                plot(xlim,(costz+uncz)*[1,1]);
                hold off;
                xlabel('Freq (Hz)');
                ylabel('dZ_{amp} (m)');
                
                %% Plot: Oscillations in X,Y,Z
                figure(2);
                clf;
                set(gcf,'Color',[1 1 1]);
                
                % X
                subplot(3,1,1);
                plot(t_cent,x_avg);		% plot the pulsed atom laser X-oscillations
                hold on;
                xlim = get( gca, 'Xlim' );
                plot( xlim, [costx+meanx costx+meanx] );
                %             plot( xlim, [costx+meanx+uncx costx+meanx+uncx],'r' )
                %             plot( xlim, [costx+meanx-uncx costx+meanx-uncx],'g' )
                plot( xlim, [-costx+meanx -costx+meanx] );
                hold off;
                xlabel('time (s)');
                ylabel('X (m)');
                
                % Y
                subplot(3,1,2);
                plot(t_cent,y_avg);
                hold on;
                xlim = get( gca, 'Xlim' );
                plot( xlim, [costy+meany costy+meany] );
                plot( xlim, [-costy+meany -costy+meany] );
                hold off;
                xlabel('time (s)');
                ylabel('Y (m)');
                
                % Z
                subplot(3,1,3);
                plot(t_cent,dz_avg);
                hold on;
                xlim = get(gca,'Xlim');
                plot(xlim,(costz+meanz)*[1,1]);
                plot(xlim,(-costz+meanz)*[1,1]);
                hold off;
                xlabel('time (s)');
                ylabel('dZ (m)');
                
                %print('BarPlot','-dpng')
            end
        end%if outerr
        if fopts.log
            f_log=fopen(path_log,'a');  % append to log-file
            fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' cost_calculator: finished ok.',...
                'choice=%s, startfile=%d bad=%s cost=%f\n'],choice,startfile,bad,cost_total);
            fclose(f_log);
        end
    % Any other cost-function than trap oscillations    
    otherwise
        cost_total = 100;
        unc_total = 100;
        bad = 'True';
        
end %switch analysis

config_out=fopts;   %return config options to user

