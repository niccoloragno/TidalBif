% ******************************************************************************************************
% Morphodynamic equilibrium of a convergent tidal channel: fluvial case (Seminara et al., 2012)
% Authors: Niccolò Ragno
% Modified on: 1-February-2020
% Open Source code, distributed under GNU General Public Licence (GPLv3)
% ******************************************************************************************************

close all
clear all
clc

%% For a given set of three dimesionless parameters, namely the Shields strees (theta_u), the aspect ratio (beta_u),
% and relative roughness (dsu), the free water surface and bed profile at
% equilibrium are computed.

epsilon=0.15;       % Small parameter tied to the magnitude of forcing tide relative to uniform flow depth
beta_u=30;          % Width-to-depth ratio
theta_u=5;          % Shields stress
dsu=0.00008;        % Relative grain size
Bu_star=100;        % Width (inlet/channel) [m]
g=9.81;             % gravity acceleration [m/s2]
Delta=1.65;         % Submerged density of sediments
i=sqrt(-1);         % Imaginary unit

Du_star=Bu_star/beta_u          % Water depth [m]
Su=theta_u*dsu*Delta            % Channel slope

a0_star=epsilon*Du_star;        % Tide amplitude [m]
T_star=43200;                   % Tide period [s]
omega_star=2*pi/T_star;         % Tide frequency [s^-1]
lambda=(omega_star*sqrt(Bu_star))/(Su^(3/2)*sqrt(beta_u*g)) % Effect of local inertia relative to convective transport
mu=5*(1+sqrt(1+18*lambda*i/25))/3;           % Defined parameter of the differential problem

N=1e5;                          % Number of x points
x_list=linspace(0,1e2,N);

epsilon_list=[0,0.1,0.2,0.3,0.4];           % Set of relative tidal amplitude to investigate
Nepsilon=length(epsilon_list);

for j=1:Nepsilon
    epsilon=epsilon_list(j);
    for k=1:N
        x=x_list(k);
        x_star(k)=x*Du_star./Su;
        
        %% Leading Order
        
        q_0(j,k)=-1;
        D_0(j,k)=1;
        h_0(j,k)=x;
        
        %% First order
        
        q_11(j,k)=(i*lambda)./(2*mu);
        h_10(j,k)=0;
        
        %% Second order
        
        D_20(j,k)=(+10*lambda.^2./(11*abs(mu).^2) + 143/88 + (5/2)*lambda*imag(mu)./(abs(mu).^2))*exp(-(mu+conj(mu))*x);
        h_20(j,k)=((167*lambda.^2 ./ (66*abs(mu).^2))+ 65/36 + 5*lambda*imag(mu)./(abs(mu).^2))*(exp(-(mu+conj(mu))*x)-1)./(mu+conj(mu));
        eta_20(j,k)=h_20(j,k)-D_20(j,k);
        
        %% Complete solution
        
        h(j,k)=h_0(j,k)+epsilon.*(h_10(j,k))+epsilon.^2.*h_20(j,k);         % Water surface level
        D(j,k)=D_0(j,k)+epsilon.^2*D_20(j,k);                               % Water flow depth
        eta(j,k)=(h_0(j,k)-D_0(j,k))+epsilon*(h_10(j,k))+epsilon.^2*(eta_20(j,k));      % Bed level
    end
end

%% Plot

figure('Name','Bed surface profile')
pl1=plot(x_list,h(1,:),'k-','Linewidth',1.2);
hold on
pl2=plot(x_list,eta(1,:),'k-','Linewidth',1.2);
pl3=plot(x_list,h(2,:),'b-','Linewidth',1);
pl4=plot(x_list,eta(2,:),'b-','Linewidth',1);
pl5=plot(x_list,h(3,:),'g-.','Linewidth',1);
pl6=plot(x_list,eta(3,:),'g-.','Linewidth',1);
pl7=plot(x_list,h(4,:),'--r','Linewidth',1);
pl8=plot(x_list,eta(4,:),'--r','Linewidth',1);
pl9=plot(x_list,h(5,:),'m-.','Linewidth',1);
pl10=plot(x_list,eta(5,:),'m-.','Linewidth',1);
grid on
xlim([0 1])
ylim([-3 1])
legend([pl1(1),pl3(1),pl5(1),pl7(1),pl9(1)],{'$\epsilon=0$','$\epsilon=0.1$','$\epsilon=0.2$','$\epsilon=0.3$','$\epsilon=0.4$'},'FontSize',10,'Location','SE','NumColumns',1,'Interpreter','latex');
xlabel('Longitudinal coordinate $x$ ','Interpreter','latex')
ylabel('Dimensionless elevation $(\eta,H)$','Interpreter','latex')
ylim([-2 1])