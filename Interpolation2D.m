    
%% Section - 1 - Gaussian Radial Basis Functions
phi =@(X1,X2,Y1,Y2,alpha) exp(-alpha^2*((X1-X2).^2+(Y1-Y2).^2));
phi_x =@(X1,X2,Y1,Y2,alpha) 2*alpha^2*(X1-X2).*exp(-alpha^2*((X1-X2).^2+(Y1-Y2).^2));
phi_y =@(X1,X2,Y1,Y2,alpha) 2*alpha^2*(Y1-Y2).*exp(-alpha^2*((X1-X2).^2+(Y1-Y2).^2));
phi_xx =@(X1,X2,Y1,Y2,alpha) (-2*alpha^2+4*(alpha^4)*(X1-X2)^2).*exp(-alpha^2*((X1-X2).^2+(Y1-Y2).^2));
phi_yy =@(X1,X2,Y1,Y2,alpha) (-2*alpha^2+4*(alpha^4)*(Y1-Y2)^2).*exp(-alpha^2*((X1-X2).^2+(Y1-Y2).^2));
phi_lap =@(X1,X2,Y1,Y2,alpha) (-4*alpha^2+4*alpha^2*((X1-X2).^2+(Y1-Y2).^2)).*exp(-alpha^2*((X1-X2).^2+(Y1-Y2).^2));
%% Section - 2 - 2D node distribution-Cartesian - RBF Centers.
Nx=100;
Ny=100;
xd = linspace(-2,2,Nx);
yd = linspace(-2,2,Ny);
[xd,yd]=meshgrid(xd,yd);
x = xd(:);
y = yd(:);
%% Section - 3 - Shape Paratemer
alpha=2;
%% Section - 4 - 2D Evaluation Points - Uniformly Distributed Random Points.
% xhat=-2+4*rand(1000,1);
% yhat=-2+4*rand(1000,1);
xd1 = linspace(-2,2,Nx/2);
yd1 = linspace(-2,2,Ny/2);
[xd1,yd1]=meshgrid(xd1,yd1);
xhat = xd1(:);
yhat = yd1(:);
%% Section - 5 - Function that we interpolate
f=x.*exp(-x.^2-y.^2);
%plot3(x,y,f,'.')
%% Section - 6 - Interpolation
A = zeros(length(x)); % Collocation Matrix
Ahat = zeros(length(x),length(xhat));
[XT1,XT2]=meshgrid(x);
[YT1,YT2]=meshgrid(y);
A = phi(XT1,XT2,YT1,YT2,alpha);
[XTe1,XTe2]=meshgrid(x,xhat);
[YTe1,YTe2]=meshgrid(y,yhat);
Ahat = phi(XTe1,XTe2,YTe1,YTe2,alpha);
%% Section - 7 - Find the interpolated function
fhat= Ahat*(A\f);
    