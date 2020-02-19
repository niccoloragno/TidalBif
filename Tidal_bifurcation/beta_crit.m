% ******************************************************************************************************
% Function for the computation of the critical aspect ratio beta_cr (linear analysis)
% Authors: Niccolò Ragno
% Modified on: 19-February-2020
% Open Source code, distributed under GNU General Public Licence (GPLv3)
% ******************************************************************************************************

function beta_C=beta_crit(param)

%% Input parameters

syms y              % Critical aspect-ratio

g=param.g;
alpha=param.alpha;
r=param.r;
i=param.i;
zeta=param.zeta;
theta_a=param.theta_a;
epsilon=param.epsilon;
Delta=param.Delta;
lambda_a=param.lambda_a;
omega_star=param.omega_star;

%% Parameters of the linear analysis

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

%% Coefficients of the linear system Delta_{jk}, {j,k}={1,2}

D_11=exp(A0*zeta)*(-3*lambda_a*N1 + N0*A1_q*zeta);
D_12=exp(A0*zeta)*(N0 + 16*N1*lambda_a/3 + N0*A1_D*zeta);
D_21=-3*G1*lambda_a*(exp(A0*zeta)-1) + G0*exp(A0*zeta)*A1_q*zeta;
D_22=(16/3)*G1*lambda_a*(exp(A0*zeta)-1) + G0*exp(A0*zeta)*A1_D*zeta - G0*(exp(A0*zeta)-1);

%% Beta cr (closed form)

beta_C0 = alpha*r*24/(7*sqrt(theta_a));         % Critical aspect ratio in absence of tide -> Solution of Bolla Pittaluga et al. (2015)

beta_C= (alpha*r/(sqrt(theta_a)))*(8*zeta + epsilon^2*(40/3*zeta*D_11 + 8*zeta*D_12 - 4*D_21) + epsilon^4*(4*D_11*D_22 - 4*D_12*D_21))/...
    (7/3*zeta + 4*epsilon^2*D_22 + 11*epsilon^2*D_21/2);

return
