function y_interp = exp_ramp_cont(y,n_tau,n_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_INTERP = exp_ramp_cont(Y,N_TAU,N_INTERP)
%
% Continuous "exponential" curve interpolation between two points
%   exp form: Y=C*exp(-t/tau)+C0
%
% INPUT
% * y: 1x2 vector of start and final value
% * n_tau: number of exponential time constants taken in interpolation
% * n_interp: number of linearly spaced interpolated points returned (excl)
%
% OUTPUT
% * y_interp: 1Xn_interp points to interpolate the boundaries (excl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1=y(1);
Y2=y(2);

C=(Y2-Y1)/(exp(-n_tau)-1);
C0=(Y2-exp(-n_tau)*Y1)/(1-exp(-n_tau));

ii=0:(n_interp+1);  %including start and end points

y_interp=C*exp(-ii*n_tau/(n_interp+1))+C0;

% %plot
% figure();
% hold on;
% scatter([0,n_interp+1],y,...
%     'bo','filled');
% plot(ii,y_interp,...
%     'r^--');
% box on; grid on;

y_interp=y_interp(2:end-1); %return the interp pts only (ends excluded)
end