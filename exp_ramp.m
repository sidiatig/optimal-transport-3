function ramp=exp_ramp(xi,xf,n_interp,ntau)
% Creates an exponential ramp joining xi and xf with n_interp points (ends exclusive)
% K determines how many "life-times" OR exponent | default is 3

if ~exist('ntau','var')
    ntau=3.5;
end
K=(ntau/n_interp);    % ~3 time constants to reach xf

n_i=1:n_interp;
ramp=xi+(xf-xi)*(1-exp(-K*n_i));

figure();
plot([xi,ramp,xf]);