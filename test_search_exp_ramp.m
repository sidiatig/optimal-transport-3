%TEST PARAM SCANNER

%%Independent exponential ramp scan for various time constants (2D)
%shunt boundaries (start and end points for exponential ramp design)
Vquad_boundary=[3.4,0.14];
Vshunt_boundary=[0,0.87];

ntau_lim=[2,10];        %param lims - number of exponential constants
n_search=12;            %interp points in lim to search for param
n_shot_avg=10;          %number of shots to take per param set
n_interp=7;             %14[7/7] param shunt
ntau=linspace(ntau_lim(1),ntau_lim(2),n_search);    %single-dim param in linear space

%create all permutations of param_values
ntau_perm=transpose(combvec(ntau,ntau));

%build param_values
n_total_perm=length(ntau_perm);
param_values=zeros(n_total_perm,n_interp*2);    %n_total_perm X 14 (7/7 ramp)
for ii=1:n_total_perm
    param_values(ii,:)=[exp_ramp(Vquad_boundary(1),Vquad_boundary(2),n_interp,ntau_perm(ii,1)),...
        exp_ramp(Vshunt_boundary(1),Vshunt_boundary(2),n_interp,ntau_perm(ii,2))];
end

hfig_exp_ramps=figure();
hold on;
for ii=1:n_total_perm
    plot(0:n_interp+1,[Vquad_boundary(1),param_values(ii,1:n_interp),Vquad_boundary(2)],...
        0:n_interp+1,[Vshunt_boundary(1),param_values(ii,n_interp+1:end),Vshunt_boundary(2)]);
end