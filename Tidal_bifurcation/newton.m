%% Globally convergent Newton method for general N x N systems 

% La funzione di Newton prende in input un vettore x0, una tolleranza, un
% numero massimo di iterazioni, la funzione (vettoriale) con la sua
% derivata, più epsilon, il "dx" delle differenze finite, il metodo
% utilizzato per il calcolo della derivata

function [x,exitflag]=newton(x0,toll,Nmax,Func,DerFunc,eps)

x=x0;
exitflag=false;                 % Exitflag mi dice se il ciclo si interrompe perchè arriva a convergenza o perchè raggiunge Nmax

for i=1:Nmax
    f=Func(x);                  % Evaluate the nonlinear function 
    res=norm(f);                % Compute the residual of the nonlinear function
    if(res<toll)                % If tolerance has been reached, stop iterations 
        exitflag=true;
        break
    end
    df=DerFunc(x,Func,eps);     % Compute the derivative of the function 
    dx=-f/df;                   % Solve the linear equation system 
    delta = 1;                  % Always try an entire Newton step
    for j=1:6
        if(norm(Func(x+delta*dx))<norm(Func(x)))
            break;
        else
            delta = delta*0.1;  %If the residual does not decrease
        end
    end
    x=x+dx;

  
end

return
