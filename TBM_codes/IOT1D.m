function [f,df,u,du]=IOT1D(J0,J1);
% The goal is to find the mass preserving transform of a 1D pdf into another one
%                       df(x)*J1(f(x))=J0(x)
% Input: 
%               J0,J1= PDF functions (Normalized 1D functions of the same
%               size). Note that J0 is the template signal.
% Output:       f    = The mass preserving mapping
%               df   = The Jacobian (in 1D is just the gradient) of f. 
% Written by Soheil Kolouri June 4, 2013.
% Center for BioImage Informatics, Carnegie Mellon University.

if length(J0) ~= length(J1)
    error('Signals must have same length!')
end
if (sum(J0<0) ~= 0 || sum(J1<0) ~= 0)
    error('Signals must be non-negative!')
end
if (size(J0,2) == 1)
    J0=J0';
end
if (size(J1,2) == 1)
    J1=J1';
end
M = 50;               % 1/M is the interpolation resolution. 

xk = 1:length(J0);    % xk is the domain of the original signal.

x = 1:1/M:length(J0); % x is the domain of the interpolated signal, 
                      % higher resolution.
J0 = interp1(xk,J0,x,'spline',0); % Interpolate J0 (Upsample J0).

J1 = interp1(xk,J1,x,'spline',0); % Interpolate J1 (Upsample J1).

J0 = J0/sum(J0); % Normalize interpolated J0

J1 = J1/sum(J1); % Normalize interpolated J1

cJ0 = cumsum(J0);% cJ0 is the CDF of J0 
cJ1 = cumsum(J1);% cJ1 is the CDF of J1 


xtilde = linspace(0,1,length(J0)); % A Grid for CDF with the same size as J0

XJ0 = interp1(cJ0,x,xtilde,'spline');% Find the xs which match the grided
                                     % CDF for cJ0
XJ1 = interp1(cJ1,x,xtilde,'spline');% Find the xs which match the grided 
                                     % CDF for cJ1
%% Plot the cdfs and the gridded points to make sure everything is working 
%  correctly.
% subplot(1,2,1)
% plot(x,cJ0);hold on;plot(x,cJ1,'r')
% axis([1 length(xk) 0 1])
% grid on
% subplot(1,2,2)
% plot(xtilde,XJ0);hold on;plot(xtilde,XJ1,'r')
% axis([0 1 1 length(xk)])
% grid on
%%
u = (XJ0-XJ1); % Calculate the translation function u(XJ0)=XJ0-f(XJ0)
u = interp1(XJ0,u,x,'spline');% Regrid u to get u(x)=x-f(x)
du = gradient(u,1/M);         % Calculate gradient of u
u = interp1(linspace(1,length(xk),length(xtilde)),u,xk,'spline',0); %Downsample u
du = interp1(linspace(1,length(xk),length(xtilde)),du,xk,'spline',0); %Downsample du
f = xk-u; %get the mapping f=x-u
df = 1-du;%get the gradient of f
end
