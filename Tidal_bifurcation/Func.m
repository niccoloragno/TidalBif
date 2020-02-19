%****************************************************************************
% Authors: Niccolò Ragno
% Modified on: 19-February-2020
% Open Source code, distributed under GNU General Public Licence (GPLv3)
%****************************************************************************

%% System of equations to solve in qb,qc,qsb,qsc,Dbu,Dcu,Db,Dc

function z=Func(y,param)

% Set di parametri

alpha=param.alpha;
r=param.r;
i=param.i;
zeta=param.zeta;
zeta_b=param.zeta_b;
zeta_c=param.zeta_c;
theta_a=param.theta_a;
beta_a=param.beta_a;
epsilon=param.epsilon;
epsilon_b=param.epsilon_b;
epsilon_c=param.epsilon_c;
Delta=param.Delta;
lambda_a=param.lambda_a;
omega_star=param.omega_star;
qa=param.qa;
qsa=param.qsa;
ra=param.ra;
rb=param.rb;
rc=1-rb;

% System coefficients

mu_a=(5/3)*(1+sqrt(1+18*i*lambda_a/25));

lambda_b=lambda_a*y(3)^(16/3)/(y(1)^3);
lambda_c=lambda_a*y(4)^(16/3)/(y(2)^3);
mu_b=(5/3)*(1+sqrt(1+18*i*lambda_b/25));
mu_c=(5/3)*(1+sqrt(1+18*i*lambda_c/25));

A_b=(mu_b + conj(mu_b))*y(1)^2 * y(3)^(-13/3);
A_c=(mu_c + conj(mu_c))*y(2)^2 * y(4)^(-13/3);
Flambda_b=143/88 + 5*lambda_b*imag(mu_b)/(2*abs(mu_b)^2) + 10*lambda_b^2/(11*abs(mu_b^2));
Flambda_c=143/88 + 5*lambda_c*imag(mu_c)/(2*abs(mu_c)^2) + 10*lambda_c^2/(11*abs(mu_c^2));
% Flambda_b=0;
% Flambda_c=0;
Glambda_b=(1/(mu_b + conj(mu_b))) * (65/36 + 5*lambda_b*imag(mu_b)/(abs(mu_b)^2) + 167*lambda_b^2/(66*abs(mu_b)^2));
Glambda_c=(1/(mu_c + conj(mu_c))) * (65/36 + 5*lambda_c*imag(mu_c)/(abs(mu_c)^2) + 167*lambda_c^2/(66*abs(mu_c)^2));
%Glambda_b=0;
%Glambda_c=0;

%% Ist mode
% System of equations  

% y(1)=qb
% y(2)=qc
% y(3)=Dbu
% y(4)=Dcu
% y(5)=Db
% y(6)=Dc

q_y =y(1)*rb/ra - qa*rb;
qs_y=y(1)^5*y(3)^(-11/2)*rb/ra - qsa*rb;

z(1)=y(1)*rb/ra+y(2)*rc/ra - qa;                      % Flow continuity
z(2)=y(1)^5*y(3)^(-11/2)*rb/ra + y(2)^5*y(4)^(-11/2)*rc/ra - qsa;                    % Exner
z(3)= qs_y/qsa - (q_y/qa - r*alpha*(y(6)-y(5))*ra/(beta_a * sqrt(theta_a)));     % Nodal
z(4)=-y(5)+y(3)+(epsilon_b^2)*(1/y(3))*Flambda_b*exp(A_b*zeta_b) ;
z(5)=-y(6)+y(4)+(epsilon_c^2)*(1/y(4))*Flambda_c*exp(A_c*zeta_c) ;
z(6)=(-zeta_b)*(y(1)^2)*(y(3)^(-10/3)) + (1/y(3))*(epsilon_b^2)*Glambda_b*(exp(A_b*zeta_b)-1) + ...
    zeta_c*(y(2)^2)*(y(4)^(-10/3)) - (1/y(4))*(epsilon_c^2)*Glambda_c*(exp(A_c*zeta_c)-1);
end