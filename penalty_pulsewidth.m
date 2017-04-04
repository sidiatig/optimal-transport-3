function penalty=penalty_pulsewidth(pulsewidth,saturation,penalty_at_saturation)
% penalty=penalty_pulsewidth(pulsewidth,saturation)
%
% INPUT:
%   * pulsewidth: N_pulses X 3 array of pulse width (in std) in X,Y,Z coordinates
%   * saturation: 1 X 3 array of critical pulse widths in X,Y,Z to incur penalty_at_saturation
%   * (optional) penalty_at_saturation: user defined penalty incurred when pulsewidth==saturation
% OUTPUT:
%   * penalty: penalty to apply for observed pulse train's pulse widths
%
% Methods:
% penalty_pulsewidth is a piecewise defined function: narrow pulse width (i.e. momentum spread in-trap)
% incur no penalty while large widths are penalised by a power law
%
% TODO: see test_penalty_pulsewidth.m for behaviour
%

if ~exist('penalty_at_saturation','var')
    penalty_at_saturation=0.1;   % 0.1 (m) penalty at saturation if undefined 
end

ndim=size(pulsewidth,2);    %number of dims supplied for pulse width

pulsewidth_avg=mean(pulsewidth,1);      %get average pulse width in each dim [m]

val=pulsewidth_avg./saturation;         %normalised pulsewidths (dimless)

threshold=2/3;    %point to begin penalty - some appropriate fraction
penalty=0;  %initialise penalty
%evaluate penalty per dim
for i=1:ndim
    % piecewise defined penalty function
    if val(i)<threshold
        penalty_temp=0;
    else
        %evaluate penalty scaling such that penalty_at_saturation is met
        penalty_scaler=penalty_at_saturation/((1-threshold)^4);
        
        %calculate penalty
        penalty_temp=penalty_scaler*(val(i)-threshold)^4;
    end
    
    %update penalty
    penalty=penalty+penalty_temp;
end