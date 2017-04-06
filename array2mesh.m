%3D array --> mesh

function [Xq,Yq,Zq] = array2mesh(x,y,z,n_interp)
if ~exist('n_interp','var')
    n_interp=10;    %default interpolating points
end
x_interp=linspace(min(x),max(x),n_interp);
y_interp=linspace(min(y),max(y),n_interp);

%build meshgrid
[Xq,Yq]=meshgrid(x_interp,y_interp);

%call griddata to interpolate
Zq=griddata(x,y,z,Xq,Yq);