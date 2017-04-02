%%characterise
saturation=30*110;

n_samp=1000;

numcounts=linspace(0,100*110,n_samp);

penalty=zeros(1,n_samp);

hfig=figure();
for i=1:n_samp
    %calculate penalty
    penalty(i)=penalty_num(numcounts(i),saturation);
end

%%Plot
hfig;
plot(numcounts,penalty,...
    'LineWidth',3);

grid on;
title('Penalty cost for total number captured from pulsed atom laser');
xlabel('Total number');
ylabel('penalty cost');