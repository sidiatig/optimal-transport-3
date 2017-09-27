
clear all; clc; close all

%% Configs

dirmain='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\ml_shunt_continued';
dirdata={'insitu_14_150_linear_20170926','insitu_exp_150_best_20170927','insitu_14_150_first_20170926','insitu_14_150_best_20170926'};
f_id=1:100;
win_t=[0.8,1.2];
win_x=[-40,10]*1e-3;
win_y=[-20,30]*1e-3;
window={win_t,win_x,win_y};
mincount=0;
maxcount=Inf;
rot_angle=0.61;     % default
build_txy=0;
verbose=0;
visual=0;

dataname={'linear','exponential','first 14','ML-14'};

%% load data
ndata=numel(dirdata);

txy=cell(1,ndata);
for ii=1:ndata
    txy{ii}=load_txy(fullfile(dirmain,dirdata{ii},'\d'),f_id,window,mincount,maxcount,rot_angle,build_txy,verbose,visual);
end

% collate shots
for ii=1:ndata
    txy{ii}=vertcat(txy{ii}{:});
end

% plot
figure(1);
plot_zxy(txy);

%% density images in X-T
% build 2D grid
n_x=200;
n_t=300;
x_ed=linspace(win_x(1),win_x(2),n_x);
t_ed=linspace(win_t(1),win_t(2),n_t);
x_c=x_ed(1:end-1)+0.5*diff(x_ed);
t_c=t_ed(1:end-1)+0.5*diff(t_ed);
% ed
[xx,tt]=meshgrid(x_c,t_c);

% get x-t array
tx=cell(1,ndata);
for ii=1:ndata
    tx{ii}=txy{ii}(:,[1,2]);
end

% 2D histogram
for ii=1:ndata
    nn{ii}=nhist(tx{ii},{t_ed,x_ed});
    % log transform
    nn{ii}=log(nn{ii});
    % normalise to max=1
    nn{ii}=nn{ii}/max(max(nn{ii}));
end

% plot
hfig=zeros(ndata,1);
for ii=1:ndata
    hfig(ii)=figure('Name',dataname{ii});
%     subplot(1,ndata,ii);
    imagesc(t_ed,x_ed,nn{ii}');
    set(gca,'YDir','normal');
    colormap(inferno);
%     colorbar();
end

%%% Annotations
% shunt times
t_shunt=[0.87,1.02];

for ii=1:ndata
    figure(hfig(ii));
    for jj=1:2
        line(t_shunt(jj)*[1,1],win_x,'Color','w','LineWidth',2);
    end
end

% axes
axcolor='w';
axlinewidth=1;
axfontsize=14;
axplotboxaspectratio=[1,0.5,1];

for ii=1:ndata
    figure(hfig(ii));
    ax=gca;
    ax.XColor=axcolor;
    ax.YColor=axcolor;
    ax.LineWidth=axlinewidth;
    ax.TickLabelInterpreter='latex';
    ax.FontSize=axfontsize;
    xlabel('T [s]');
    ylabel('X [m]');
    
    ax.PlotBoxAspectRatio=axplotboxaspectratio;
end

%% Save figs
% dir
pathdir='C:\Users\HE BEC\Desktop';

for ii=1:ndata
    print(hfig(ii),fullfile(pathdir,dataname{ii}),'-dsvg');
end