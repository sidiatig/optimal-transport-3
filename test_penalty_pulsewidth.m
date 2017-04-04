%Script to test penalty_pulsewidth behaviour
%
%since pulse width penalty is arithmetically summed in each dimension, here we test the penalty
%function for a single dimension
%

nsamp=1000;

pulse_width_saturation=0.008;   %pulse width saturation in standard deviation
penalty_width_saturation=0.1;   %penalty cost at saturation
pulse_width=linspace(0,1.5*pulse_width_saturation,nsamp);

%evaluate penalties at each 1D pulse width (avg for shot)
penalty_pulse_width=zeros(1,nsamp);
for i=1:nsamp
    penalty_pulse_width(i)=penalty_pulsewidth(pulse_width(i),pulse_width_saturation,penalty_width_saturation);
end

%Plot
figure();
plot(pulse_width/pulse_width_saturation,penalty_pulse_width,'LineWidth',2);
titlestr=sprintf('penaltywidthsaturation=%0.3g',...
    penalty_width_saturation);
title(titlestr);
xlabel('width/width_{sat}');
ylabel('penalty cost [m]');

figure();
plot(pulse_width,penalty_pulse_width,'LineWidth',2);
titlestr=sprintf('pulsewidthsaturation=%0.3g; penaltywidthsaturation=%0.3g',...
    pulse_width_saturation,penalty_width_saturation);
title(titlestr);
xlabel('width [m]');
ylabel('penalty cost [m]');