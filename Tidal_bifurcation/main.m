% ******************************************************************************************************
%% Morphodynamic equilibrium of a tidal bifurcation with prismatic-like channels
% Authors: Niccolò Ragno
% Modified on: 19-February-2020
% Open Source code, distributed under GNU General Public Licence (GPLv3)
% ******************************************************************************************************

close all
clear all
clc

%% Model parameters

g=9.81;             % gravity acceleration [m/s2]
Delta=1.65;         % Submerged density of sediments
i=sqrt(-1);         % Imaginary unit
r=0.5;               % Ikeda coefficient
alpha=4;            % Dimensionless length of the upstream cells

%% Basic flow parameters & Boundary conditions

beta_a=15       % Half width-to-depth ratio
theta_a=1       % Shields stress
Da=1;           % Dimensionless upstream depth
qa=1;           % Dimensionless uniform upstream discharge
qsa=1;          % Upstream sediment discharge which feed each cell
epsilon=0.3;            % Scaled tidal forcing in comparison to Da_star
epsilon_b=0.2;         % Scaled tidal amplitude for channel b
epsilon_c=0.2;         % Scaled tidal amplitude for channel c
x=0.25;          % Dimensionless longitudinal coordinate  x_star*Sa/(Da_star)
x_b=0.25;       % Dimensionless longitudinal coordinate of channel b
x_c=0.25;       % Dimensionless longitudinal coordinate of channel c
zeta=-x;                % Definition of a dimensionless coordinate zeta=-x ()
zeta_b=-x_b;
zeta_c=-x_c;
T_star=43200;                   % Tide period [s] M2 tide (12h)
omega_star=2*pi/T_star;         % Tide frequency [s^-1]
lambda_a=10;

ra=1;    % Ratio between main channel width and sum of downstream widths (Wa/(Wb+Wc))
rb=0.5;  % Width of channel B relative to the sum of downstream widths (Wb/(Wb+Wc))
rc=1-rb;

%% Solution of the non-linear system

N=100;
beta_min=4;
beta_max=30;
beta_list=linspace(beta_min,beta_max,N);

theta_min=0.5;
theta_max=5;
theta_list=linspace(theta_min,theta_max,N);

zeta_list=[-1,-0.75,-0.5,-0.25,-0.1];
epsilon_list=[0,0.1,0.2,0.3,0.4];
Neps=length(epsilon_list);
Nzeta=length(zeta_list);

% for s=1:Neps
%     epsilon_b=epsilon_list(s);
%     epsilon_c=epsilon_list(s);
%     zeta_b=zeta_list(s);
%     zeta_c=zeta_list(s);

% Set of initial conditions

%     IC(:,1)=[0.5*qa 1.5*qa 0.5*Da 1.5*Da 0.5*Da 1.5*Da];
%     IC(:,2)=[qa qa Da Da Da Da];
%     IC(:,3)=[1.5*qa 0.5*qa 1.5*Da 0.5*Da 1.5*Da 0.5*Da];
%     IC(:,4)=[1.1*qa 0.9*qa 1.1 0.9 1.1 0.9];
%     IC(:,3)=[0.5*qa 1.1*qa 0.9 1.1 0.9 1.1];
%
%     for j=1:size(IC,2)
%         for k=1:N
%             beta_a=beta_list(k);
%
%             param.g=g;
%             param.alpha=alpha;
%             param.r=r;
%             param.i=i;
%             param.zeta=zeta;
%             param.zeta_b=zeta_b;
%             param.zeta_c=zeta_c;
%             param.theta_a=theta_a;
%             param.beta_a=beta_a;
%             param.epsilon=epsilon;
%             param.epsilon_b=epsilon_b;
%             param.epsilon_c=epsilon_c;
%             param.Delta=Delta;
%             param.lambda_a=lambda_a;
%             param.omega_star=omega_star;
%             param.qa=qa;
%             param.qsa=qsa;
%             param.ra=ra;
%             param.rb=rb;

% Solution of the non-linear system

%             % y(1)=qb
%             % y(2)=qc
%             % y(3)=Dbu
%             % y(4)=Dcu
%             % y(5)=Db
%             % y(6)=Dc

%             f=@(y) Func(y,param);
%             options = optimoptions('fsolve', 'Algorithm', 'trust-region', 'FiniteDifferenceType','central');
%             [Sol feval] = fsolve(f, IC(:,j),options);
%
%             Sol=real(Sol);
%             qb(k,s)=Sol(1);
%             qc(k,j,s)=Sol(2);
%             Dbu(k,s)=Sol(3);
%             Dcu(k,s)=Sol(4);
%             Db(k,j,s,t)=Sol(5);
%             Dc(k,j,s,t)=Sol(6);
%             DeltaQ(k,j,s)=1-2*qc(k,j,s)/qa*rc/ra;          % Tidally averaged discharge asymmetry

%         end
%     end

%% Computation of the critical aspect ratio

for t=1:Nzeta
    epsilon=epsilon_list(t);
    for k=1:N
        theta_a=theta_list(k);
        param.g=g;
        param.alpha=alpha;
        param.r=r;
        param.i=i;
        param.zeta=zeta;
        param.theta_a=theta_a;
        param.beta_a=beta_a;
        param.epsilon=epsilon;
        param.Delta=Delta;
        param.lambda_a=lambda_a;
        param.omega_star=omega_star;
        param.qa=qa;
        param.qsa=qsa;
        param.ra=ra;
        param.rb=rb;
        beta_C(k,t)=beta_crit(param);
    end
end

%% Computation of the critical length

% for t=1:Ntheta
%     theta_a=theta_list(t);
%     for k=1:N
%         beta_a=beta_list(k);
%         param.g=g;
%         param.alpha=alpha;
%         param.r=r;
%         param.i=i;
%         param.beta_a=beta_a;
%         param.theta_a=theta_a;
%         param.epsilon=epsilon;
%         param.Delta=Delta;
%         param.lambda_a=lambda_a;
%         param.omega_star=omega_star;
%
%         toll=1e-7;              % Tolerance (max norm of the residuals)
%         Nmax=300;               % Maximum number of iterations
%         y0=[-0.1];               % Same bed elevation
%
%         f=@(y) length_crit(param,y);
%         options = optimoptions('fsolve', 'Algorithm', 'trust-region', 'FiniteDifferenceType','central');
%         [Sol feval] = fsolve(f, y0, options);
%          Length_C(t,k)=-Sol;             % The '-' is to convert coordinate zeta to x
%
%     end
%          diffL=diff(Length_C(t,:));             % Exlude the values of beta<beta_cr (asymptotic solution)
%          [valmax,indmax]=max(diffL);
%          Length_C(t,1:indmax)=NaN;
% end

%% Plot

figure('Name','Critical aspect ratio')
pl1=plot(theta_list,beta_C(:,1),'k-','Linewidth',1.2);
hold on
pl2=plot(theta_list,beta_C(:,2),'b-','Linewidth',1);
pl3=plot(theta_list,beta_C(:,3),'g-.','Linewidth',1);
pl4=plot(theta_list,beta_C(:,4),'--r','Linewidth',1);
pl5=plot(theta_list,beta_C(:,5),'m-.','Linewidth',1);
grid on
legend([pl1(1),pl2(1),pl3(1),pl4(1),pl5(1)],{'$\epsilon=0$','$\epsilon=0.1$','$\epsilon=0.2$','$\epsilon=0.3$','$\epsilon=0.4$'},'FontSize',10,'Location','NE','NumColumns',1,'Interpreter','latex');
xlabel('Shields stress $\theta_a$','Interpreter','latex')
ylabel('Critical aspect ratio $\beta_{cr}$','Interpreter','latex')