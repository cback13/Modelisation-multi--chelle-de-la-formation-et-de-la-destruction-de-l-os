clear all
clc


%valeurs des paramètres pour test1
% k1 = 2.5;
% n_ob1 = 4;
% alpha1 =0.1;
% k2 =1.5;
% n_ob2 = 4;
% KB =10;
% KC =10;
% teta1 = 1.5;
% k3 =0.9;
% n_oc1 = 2;
% n_oc2 =2;
% alpha2 = 1;
% k4 =0.4 ; 
% Eymax = 33;
% teta2 = 0.9;
% alpha3 = 0.9;
% k5 =1.5;
% n_y1 = 5;
% n_y2 = 5;
% k6 = 1.5;
% teta4 = 0.9; 
% teta0 = 1.9;
% teta5 = 1.4;

%paramètres test 2 avec Ey en terme soustraction 
% Ob(k+1) = Ob(k) + dt* ((alpha1 *(k1^n_ob1/(teta0^n_ob1+Ey(k)^n_ob1))*Ob(k)*(1-Ob(k)/KB)) - (Ob(k)*k2*(Ey(k)^n_ob2/(teta1^n_ob2+Ey(k)^n_ob2))));
%      Oc(k+1) = Oc(k) + dt* ((k3 *(Ey(k)^n_oc1/(teta2^n_oc1+Ey(k)^n_oc1))*(1-Oc(k)/KC)*Oc(k)) - (Oc(k)*alpha2*(k4^n_oc2/(k4^n_oc2+Ey(k)^n_oc2))));
%      Ey(k+1) = Ey(k) + dt*((Ey(k)*(1-Ey(k)/Eymax)*alpha3*(k5^n_y1/(k5^n_y1+Ob(k)^n_y1))) - (Ey(k)*k6*(Oc(k)^n_y2/(teta3^n_y2+Oc(k)^n_y2))));
% k1 = 1.5;
% n_ob1 = 4;
% alpha1 =1;
% k2 =0.3;
% n_ob2 = 4;
% KB =10;
% KC =10;
% teta1 = 1.5;
% k3 =1;
% n_oc1 = 2;
% n_oc2 =2;
% alpha2 = 1;
% k4 =1 ; 
% Eymax = 33;
% teta2 = 1;
% alpha3 = 1;
% k5 =1;
% n_y1 = 4;
% n_y2 = 4;
% k6 = 1;
% teta3 = 1; 
% teta0 = 2;

k1 = 1.8;
n_ob1 = 0.1;
alpha1 =1.1;
k2 =4;
n_ob2 = 0.1;
KB =10;
KC =10;
teta1 = 1;
k3 =1.5;
n_oc1 = 0.3;
n_oc2 = 0.3;
alpha2 = 1;
k4 =1 ; 
Eymax = 33;
teta2 = 1.5;
alpha3 = 0.5;
k5 =1.5;
n_y1 = 0.1;
n_y2 = 0.1;
k6 = 0.06;
teta3 = 0.05; 
teta0 = 3;
krb = 0.01;
kry = 0.1;
krc = 0.5;

%temps
dt = 0.01 ;
T=3000;t=0:dt:T;tn=length(t);

%initialisation variables
Ob(1) = 0.00002 ;
Oc(1) = 0.00001;
Ey(1) = 0.000015;

for k = 1:tn-1
     Ob(k+1) = Ob(k) + dt* ((alpha1 *(k1^n_ob1/(teta0^n_ob1+Ey(k)^n_ob1+krb*Ey(k)))*Ob(k)*(1-Ob(k)/KB)) - (Ob(k)*k2*(Oc(k)^n_ob2/(teta1^n_ob2+Oc(k)^n_ob2))));
     Oc(k+1) = Oc(k) + dt* ((k3 *(Ey(k)^n_oc1/(teta2^n_oc1+Ey(k)^n_oc1+krc*Ey(k)))*(1-Oc(k)/KC)*Oc(k)) - (Oc(k)*alpha2*(k4^n_oc2/(k4^n_oc2+Ey(k)^n_oc2+krc*Ey(k)))));
     Ey(k+1) = Ey(k) + dt*((Ey(k)*(1-Ey(k)/Eymax)*alpha3*(Ob(k)^n_y1/(k5^n_y1+Ob(k)^n_y1))) - (Ey(k)*k6*(Oc(k)^n_y2/(teta3^n_y2+Oc(k)^n_y2))));
 end

 save('Ob.mat', 'Ob');
%-lbd2*(Ey(k)/(teta1+Ey(k)))*Ob(k));
%-lbd4*(1/(teta5+Ey(k)))*Oc(k)));
%-lbd6*(Oc(k)/(teta4+Oc(k))));
 subplot(1,3,1);
 plot(t, Ob,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('Ob(t)','FontSize', 12);title("Evolution of Osteoblasts",'FontSize', 14);hold on;
 grid on;
 subplot(1,3,2);
 plot(t, Oc,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('Oc(t)','FontSize', 12);title("Evolution of Osteoclasts",'FontSize', 14);hold on;
 grid on;
 subplot(1,3,3);
 plot(t, Ey,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('Ey(t)','FontSize', 12);title("Evolution of the Stiffness",'FontSize', 14);hold on;
 grid on;