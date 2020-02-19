% Func derivative->Jacobian
% Central finite difference method

%% Numerical derivative of a vector function

function dz = DerFunc(x,Func,eps)

N = length(x);
dz = zeros(N,N);

for j=[1:N]
    xp = x;
    xp(j) = xp(j) + eps;
    xm = x;
    xm(j) = xm(j) - eps;
    dz(:,j) = (Func(xp)-Func(xm))/(2*eps);
end

end