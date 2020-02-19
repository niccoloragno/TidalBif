% ******************************************************************************************************
% Code for the computation of the critical length (Linear analysis)
% Authors: Niccolò Ragno
% Modified on: 19-February-2020
% Open Source code, distributed under GNU General Public Licence (GPLv3)
% ******************************************************************************************************

function Length_C=length_crit(param,y)

%% Input parameters

%syms y              % Critical length

g=param.g;
alpha=param.alpha;
r=param.r;
i=param.i;
beta_a=param.beta_a;
theta_a=param.theta_a;
epsilon=param.epsilon;
Delta=param.Delta;
lambda_a=param.lambda_a;
omega_star=param.omega_star;

%% Parameters

mu_a=(5/3)*(1+sqrt(1+18*i*lambda_a/25));
mu_1=3*i/(3*mu_a - 5);
A0 = 2*real(mu_a);
A1_q = 4*real(mu_a) - 6*real(mu_1)*lambda_a;
A1_D= real(mu_1)*lambda_a*32/3 - 26*real(mu_a)/3 ;
N0= 143/88 + 5*lambda_a*imag(mu_a)/(2*abs(mu_a)^2) + 10*lambda_a^2/(11*abs(mu_a)^2);
N1= 20*lambda_a/(11*abs(mu_a)^2) - 10*(mu_1*conj(mu_a)+conj(mu_1)*mu_a)*lambda_a^2/(11*abs(mu_a)^4) + ...
    5*imag(mu_a)*lambda_a/(2*abs(mu_a)^2) + 5*lambda_a*imag(mu_1)/(2*abs(mu_a)^2) - ...
    5*lambda_a*(mu_1*conj(mu_a)+conj(mu_1)*mu_a)*imag(mu_a)/(2*abs(mu_a)^4);
G0=(65/36 + 5*lambda_a*imag(mu_a)/(abs(mu_a)^2) + 167*lambda_a^2/(66*abs(mu_a)^2))/A0;
G1=(1/(2*real(mu_a)))*(-G0*(mu_1 + conj(mu_1)) +  167*lambda_a/(33*abs(mu_a)^2) -...
    167*lambda_a^2/(66*abs(mu_a)^4)*(mu_1*conj(mu_a)+conj(mu_1)*(mu_a))) + 5*imag(mu_a)/(abs(mu_a)^2) + ...
    5*lambda_a*imag(mu_1)/(abs(mu_a)^2) - 5*lambda_a*(mu_1*conj(mu_a)+conj(mu_1)*(mu_a))*imag(mu_a)/(abs(mu_a)^4) ;

psi_11=4 -4*(r*alpha)/(beta_a*sqrt(theta_a))*((epsilon^2)*exp(A0*y)*(-3*lambda_a*N1 + N0*A1_q*y));
psi_12=-11/2 - 4*(r*alpha)/(beta_a*sqrt(theta_a))*(1+(epsilon^2)*exp(A0*y)*(-N0 + 16*N1*lambda_a/3 + N0*A1_D*y));
psi_21=-2*y + (epsilon^2)*(-3*G1*lambda_a*(exp(A0*y)-1) + G0*exp(A0*y)*A1_q*y);
psi_22=10*y/3 + (epsilon^2)*((16/3)*G1*lambda_a*(exp(A0*y)-1) + G0*exp(A0*y)*A1_D*y - G0*(exp(A0*y)-1));

A=[psi_11 psi_12
    psi_21 psi_22];

Length_C=det(A);

%% Approximate form of L_cr, neglecting the fourth-order terms of epsilon

% r_hat=4*r*alpha/(beta_a*sqrt(theta_a));
%
% Length_C=40*y/3 + 64*epsilon^2*G1*lambda_a*exp(A0*y)/3 - (64/3)*epsilon^2*G1*lambda_a + 4*epsilon^2*G0*A1_D*y*exp(A0*y) - ...
%     4*epsilon^2*G0*exp(A0*y) + 4*epsilon^2*G0 + (10/3)*y*r_hat*epsilon^2*3*lambda_a*N1*exp(A0*y) -...
%     (10/3)*y^2*r_hat*epsilon^2*exp(A0*y)*N0*A1_q - 11*y - (33/2)*epsilon^2*G1*lambda_a*exp(A0*y) +...
%     (33/2)*epsilon^2*G1*lambda_a + 5.5*epsilon^2*G0*exp(A0*y)*A1_q*y - 2*r_hat - 3*r_hat*epsilon^2*G1*lambda_a*exp(A0*y) +...
%     3*r_hat*epsilon^2*G1*lambda_a + r_hat*epsilon^2*G0*exp(A0*y)*A1_q*y + 2*r_hat*y*epsilon^2*exp(A0*y)*N0 -...
%     2*y*r_hat*epsilon^2*exp(A0*y)*16*lambda_a*N1/3 + 2*r_hat*y^2*epsilon^2*exp(A0*y)*N0*A1_D ;

% Length_C_symb=solve(LC,y);                        % Symbolic solution of the 2x2 matrix
% Length_C_compl=double(Length_C_symb);          % Complete solution (2 roots)
% Length_C=abs(Length_C_symb);
% Length_C=Length_C_compl(1)  ;    % Higher root

return
