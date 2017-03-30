function[out]=FFTxt(t,x)

T = (t(2)-t(1));             % Sampling period
%disp(num2str(T));
Fs=1/T;
L = size(x,1);             % Length of signal
%disp(num2str(L))
%t = (0:L-1)*T;        % Time vector
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
f=transpose(f);
out=[f P1];

end
