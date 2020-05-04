% ******************************************************************************************************
% Morphodynamic equilibrium of a convergent tidal channel: fluvial case (Seminara et al., 2012)
% Authors: Niccolò Ragno
% Modified on: 27-April-2020
% Open Source code, distributed under GNU General Public Licence (GPLv3)
% NB:
%   *Research code, developed for academic use only!
%   *Not fully tested: please report possible issues at niccolo.ragno@edu.unige.it
% ******************************************************************************************************

close all
clear all
clc

addpath('Functions')

%% Input parameters

Q=500;          %Water discharge [m^3/s]
S=0.001;     %Downstream gradient [m/m]
d50=0.0001;        %Representative grain size [m]
Wa=150;          %Main channel width [m]
Delta=1.65;      %Relative submerged density of sediment [-]
a=0.5;            %Tidal amplitude [m]
g=9.81;             % gravity acceleration [m/s2]
i=sqrt(-1);         % Imaginary unit
ks=35;              %Gauckler-Strickler coefficient 
T=43200;            % Tide period [s]
omega=2*pi/T;         % Tide frequency [s^-1]

Da=(Q/(Wa*sqrt(S)*ks))^(3/5);           %Water depth [m]
C=ks*Da^(1/6)/sqrt(g);              %Chezy coefficient [-]

beta=Wa/(2*Da);                             % Aspect ratio [-]
theta=Q^2/(Wa^2*C^2*g*Delta*d50);           % Shields parameter [-]
ds=d50/Da;                                  % Relative grain size[-]
epsilon=a/Da;        % Small parameter tied to the magnitude of forcing tide relative to uniform flow depth

lambda=(omega*Da^2*Wa)/(S*Q);
mu=5*(1+sqrt(1+18*lambda*i/25))/3;           % Defined parameter of the differential problem


N=10000;
x_list=linspace(0,1e2,N);

a_list=[0,0.25,0.5,0.75,1];
Na=length(a_list);
Q_list=[2000,1000,500,250,100];
NQ=length(Q_list);

for j=1:NQ
    Q=Q_list(j);
    Da=(Q/(Wa*sqrt(S)*ks))^(3/5);           %Water depth [m]
    C=ks*Da^(1/6)/sqrt(g);              %Chezy coefficient [-]
    
    beta=Wa/(2*Da);                             % Aspect ratio [-]
    theta=Q^2/(Wa^2*C^2*g*Delta*d50);           % Shields parameter [-]
    ds=d50/Da;                                  % Relative grain size[-]
    epsilon=a/Da        % Small parameter tied to the magnitude of forcing tide relative to uniform flow depth
    
    lambda=(omega*Da^2*Wa)/(S*Q)
    mu=5*(1+sqrt(1+18*lambda*i/25))/3;           % Defined parameter of the differential problem
    
    for k=1:N
        x=x_list(k);
        x_star(k)=x*Da./S;
        %% Zero Order
        
        q_0(j,k)=-1;
        D_0(j,k)=1;   % Modulus of the complex number (q0=-1 -> D0=1 Uniform flow in absence of channel convergence)
        h_0(j,k)=x;
        
        %% First order
        
        q_11(j,k)=(i*lambda)./(2*mu);
        D_11(j,k)=0.5;
        h_11(j,k)=D_11(j,k);
        
        h_10(j,k)=0;
        
        %% Second order
        
        D_20(j,k)=(+10*lambda.^2./(11*abs(mu).^2) + 143/88 + (5/2)*lambda*imag(mu)./(abs(mu).^2))*exp(-(mu+conj(mu))*x);
        h_20(j,k)=((167*lambda.^2 ./ (66*abs(mu).^2))+ 65/36 + 5*lambda*imag(mu)./(abs(mu).^2))*(exp(-(mu+conj(mu))*x)-1)./(mu+conj(mu));
        eta_20(j,k)=h_20(j,k)-D_20(j,k);
        %% Complete solution
        
        h(j,k)=h_0(j,k)+epsilon.*(h_10(j,k))+epsilon.^2.*h_20(j,k);
        D(j,k)=D_0(j,k)+epsilon.^2*D_20(j,k);
        eta(j,k)=(h_0(j,k)-D_0(j,k))+epsilon*(h_10(j,k))+epsilon.^2*(eta_20(j,k));
    end
end



%% Plot

figure('Name','Bed surface profile')
pl1=plot(x_list,h(1,:),'k','Linewidth',1);
hold on
pl2=plot(x_list,eta(1,:),'k','Linewidth',1);
pl3=plot(x_list,h(2,:),'b-','Linewidth',1);
pl4=plot(x_list,eta(2,:),'b-','Linewidth',1);
pl5=plot(x_list,h(3,:),'g-.','Linewidth',1);
pl6=plot(x_list,eta(3,:),'g-.','Linewidth',1);
pl7=plot(x_list,h(4,:),'r--','Linewidth',1);
pl8=plot(x_list,eta(4,:),'r--','Linewidth',1);
pl9=plot(x_list,h(5,:),'m-.','Linewidth',1);
pl10=plot(x_list,eta(5,:),'m-.','Linewidth',1);
grid on
xlim([0 1])
ylim([-3 1])
% legend([pl1(1),pl3(1),pl5(1),pl7(1),pl9(1)],{'$\theta=0.5$','$\theta=1$','$\theta=1.5$','$\theta=2$','$\theta=2.5$'},'FontSize',10,'Location','SE','NumColumns',1,'Interpreter','latex');
% legend([pl1(1),pl3(1),pl5(1),pl7(1),pl9(1)],{'$\epsilon=0$','$\epsilon=0.1$','$\epsilon=0.2$','$\epsilon=0.3$','$\epsilon=0.4$'},'FontSize',10,'Location','SE','NumColumns',1,'Interpreter','latex');
xlabel('Longitudinal coordinate $x$ ','Interpreter','latex')
ylabel('Dimensionless elevation $(\eta,H)$','Interpreter','latex')
ylim([-2 1])
