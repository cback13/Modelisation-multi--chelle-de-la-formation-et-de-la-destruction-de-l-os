clear all
clc
%activation = palier max ob ey =  10jours 
%resorbtion = action ostéoclastes 21 jours 
% reversal = palier à 0 Ob ey = 5 jours
%formation et minéralisation 90 jours à 1 an , 90 jours collagen, 1 an
%minéraux ? 
% en jour total = 126 jours 
% en heure = 120+2160+240+ 504
%modèle Ey F Oc Ob P
k1 = 1.38;
n_ob1 = 4;
alpha1 =0.7;
k2 =4;
n_ob2 = 4;
KB =10;
KC =10;
teta1 = 1;
k3 =0.9;
n_oc1 = 4;
n_oc2 = 4;
alpha2 = 1.3;
k4 =0.09 ; 
Eymax = 33;
teta2 = 1;
alpha3 = 0.02;
k5 =0.3;
n_y1 = 6;
n_y2 = 6;
k6 = 0.4;
teta3 = 0.5; 
teta0 = 2.3;
phi=1;

%temps
dt = 0.01 ;
T=12500;t=0:dt:T;tn=length(t);

%age
da=0.1;
A=80;a=0:da:A;an=length(a);
%seuil maturité
alpha=41.7;
alphav=0:da:alpha;alphan=length(alphav);

%initialisation variables
Ob(1) = 0.00002 ;
Oc(1) = 0.00001;
Ey(1) = 0.000015;
% condition au bord c
for k=1:tn
    c(k,1)=0;
end
% condition initial c
for j=2:an
  c(1,j)=(1/50)*exp(-(j*da-A/2)^2 / 10);
end
%condition initial C1 C2 Ctot F
s1=0;
s2=0;
for i=2:an
    s1=s1+c(1,i)*da;  
end
Ctot(1)=s1;
%mature collagen 
for i=alphan+1:an
    s2=s2+c(1,i)*da;  
end
C2(1)=s2;
%immature collagen
C1(1)=s1-s2;
F(1)=C2(1)/Ctot(1);


for k = 1:tn-1
     Ob(k+1) = Ob(k) + dt* ((alpha1 *(k1^n_ob1/(teta0^n_ob1+Ey(k)^n_ob1))*Ob(k)*(1-Ob(k)/KB)) - (Ob(k)*k2*(Oc(k)^n_ob2/(teta1^n_ob2+Oc(k)^n_ob2))));
     Oc(k+1) = Oc(k) + dt* ((k3 *(Ey(k)^n_oc1/(teta2^n_oc1+Ey(k)^n_oc1))*(1-Oc(k)/KC)*Oc(k)) - (Oc(k)*alpha2*(k4^n_oc2/(k4^n_oc2+Ey(k)^n_oc2))));
     Ey(k+1) = Ey(k) + dt*((Ey(k)*(1-Ey(k)/Eymax)*alpha3*(F(k)^n_y1/(k5^n_y1+F(k)^n_y1))) - (Ey(k)*k6*(Oc(k)^n_y2/(teta3^n_y2+Oc(k)^n_y2))));
    %calcul de c
     s1=0;
     s2=0;
     if round(Ob(k+1),1)==round(Ob(k),1) && Ob(k+1)<KB/2
        for j=2:an
           c(k+1,j)= (1/50)*exp(-(j*da-A/2)^2 / 10);
        end
     elseif round(Ob(k+1),1)==round(Ob(k),1) && Ob(k+1)>KB/2
         for j=2:an
           c(k+1,j)= c(k,j);  
        end
     elseif Ob(k)<Ob(k+1)
        for j=2:an 
        % Séchama upwind
            c(k+1,j)= c(k,j)-14*(dt /da)*(c(k,j)-c(k,j-1)) ;
        end
     elseif Ob(k)>Ob(k+1) 
        for j=2:an 
            c(k+1,j)= (1/50)*exp(-(j*da-A/2)^2 / 10);
        end
     end
     for j=2:an 
        s1=s1+c(k+1,j)*da;
     end
     for i=alphan+1:an
        s2=s2+c(k+1,i)*da;  
     end
     Ctot(k+1)=s1;
     %mature collagen
     C2(k+1)=s2;
     %immature collagen
     C1(k+1)=s1-s2;
     F(k+1)=C2(k+1)/Ctot(k+1);
end

 subplot(1,4,1);
 plot(t, Ob,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('Ob(t)','FontSize', 12);title("Evolution of Osteoblasts",'FontSize', 10);hold on;
 grid on;
 subplot(1,4,2);
 plot(t, Oc,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('Oc(t)','FontSize', 12);title("Evolution of Osteoclasts",'FontSize', 10);hold on;
 grid on;
 subplot(1,4,3);
 plot(t, Ey,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('Ey(t)','FontSize', 12);title("Evolution of the Stiffness",'FontSize', 10);hold on;
 grid on;
 subplot(1,4,4);
 plot(t, F,'r','LineWidth', 2);xlabel('time','FontSize', 12);ylabel('F(t)','FontSize', 12);title("Evolution of the Collagen maturity",'FontSize', 10);hold on;
 grid on;