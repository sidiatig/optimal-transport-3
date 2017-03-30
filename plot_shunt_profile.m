%Plot shunt params: 14 params (7 QUAD /7 SHUNT )
function hfig=plot_shunt_profile(params,boundary,hfig)
% boundary: [QUAD_i, QUAD_f, SHUNT_i, SHUNT_f]

if exist('hfig','var')
    figure(hfig);
else
    hfig=figure();
end

quad=[boundary(1),params(1:7),boundary(2)];
shunt=[boundary(3),params(8:end),boundary(4)];

plot(quad,'bo-','LineWidth',2);
hold on;
plot(shunt,'go-','LineWidth',2);
end