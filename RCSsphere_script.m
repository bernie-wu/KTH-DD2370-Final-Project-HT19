%%
close all;
clear;
clc;

load('analyticalresults')

n = 51;
r = linspace(0.5,20,n) * 1.0e-2;
k = zeros(1,n);
TC = zeros(1,n);
TU = zeros(1,n);

for i = 1:n
    [k(i),TU(i),TC(i)] = RWG_Sphere_v2(1,r(i),1e9,0,0,0,0,0);
    TU(i) = TU(i)/(pi*(1*r(i))^2);
    disp(i);
end
%%
close all;
%%
figure(1);
plot(k.*r , TC , 'LineWidth',2);
xlabel('kr')
ylabel('Current')
figure(2);
plot(k.*r , TU , 'r', 'LineWidth',2);
%title('RCS/(\pi a^2)')
ylabel('\sigma_t_o_t /(\pi a^2)')
xlabel('ka')
%saveas(gcf,'RCSKA5cm4.png')
hold on;
%%
plot(kb(1:5:200),sigma(1:5:200),'--k','LineWidth',2)
%%
legend( 'MoM' , 'Analytical' , 'Location' , 'southeast')





