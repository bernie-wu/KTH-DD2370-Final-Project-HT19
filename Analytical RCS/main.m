close all;
clear;
clc;

c       = 299792458;
eps0    = 8.85418782e-012;
mu0     = 1.25663706e-006;
global lambda eta 
lambda= 1; %Validation can be done by setting Lambda = 1 and eta = 0.99 which should converge to PEC case
eta = 1.00;
eps1 = lambda*eps0;
imp = sqrt(mu0/eps0);
mu1 = mu0;
c1 = 1/sqrt(eps1*mu1);
lmax = 60;
resolution = 200;
maxkb = 4.5;
% ka = linspace(0,20*eta,200);
ka = logspace(0,1,resolution);
ka = maxkb*eta*(ka-1)/(max(ka)-1);
kb = ka/eta;
sigma = zeros(size(ka));
for l=1:lmax
    temp = 2*l*(l+1)*(abs(a(l,ka)).^2 + abs(b(l,ka)).^2);%Change ad(l,ka) to a(l,ka) and bd(l,ka) to b(l,ka) for PEC case. 
%     plot(ka,abs(b(l,ka))./ka);
%     plot(abs(b(l,ka)));
%     plot(temp);
    hold on;
    sigma = sigma+temp;
end
sigma = sigma./(pi*kb.^2);
plot(kb,sigma,'LineWidth',2);