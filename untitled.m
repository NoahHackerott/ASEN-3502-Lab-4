%% 3.1 Adam Function implementation
% These are given from teh first week parameters
params = struct();
params.rho = 1.225; %kg/m^2
params.CL0 = 0.247;
params.CM0 = 0.1;
params.DCLp0 = 0.6;
params.V = 30; %m/s
params.S = 16.3; %m^2
params.W = 11346; %N
params.CLal = 6.0; %rad^-1
params.CMal = -1.3; %rad^-1
params.alp = 0.05; %rad
params.CLdel = 0.43; %rad^-1
params.CMdel = -1.28; %rad^-1


x = [0.2;0.1];
l = 0.02;
beta1 = 0.3;
beta2 = 0.99;
epsilon = 10e-8;
tol = 10e-5;
maxit=1000;
f = @(alpha,deltaE) ...
    ((0.5*params.rho*params.V^2*params.S .* ...
    (params.CL0 + params.CLal.*alpha + params.CLdel.*deltaE ...
    + params.DCLp0.*exp(-(alpha./params.alp).^2)) ...
    - params.W)./params.W).^2 ...
    + (params.CM0 - params.CMal.*alpha ...
    - params.CMdel.*deltaE).^2;
% Function used for trim gradient
g = @(x) trim_gradient(x(1), x(2), params);

[position path] = Adam(g,x,tol,maxit,l,beta2,beta1,epsilon);
fprintf('The final position is %d\n',position)


plot_contours(1,f);

hold on;
plot(path(1,:)*180/pi,path(2,:)*180/pi,'r-o','LineWidth',2);
plot(position(1)*180/pi,position(2)*180/pi,'bo','MarkerFaceColor','b');
legend('Objective Function','ADAM Path','Final Position');
hold off;



%% Functions
function G = trim_gradient(alpha, deltaE,params)
rho = params.rho; %kg/m^2
CL0 = params.CL0;
CM0 = params.CM0;
DCLp0 = params.DCLp0;
V = params.V; %m/s
S = params.S; %m^2
W = params.W; %N
CLal = params.CLal; %rad^-1
CMal = params.CMal; %rad^-1
alp = params.alp; %rad
CLdel = params.CLdel; %rad^-1
CMdel = params.CMdel; %rad^-1

syms a de real 
f = ((0.5*rho*V^2*S * (CL0 + CLal*a + CLdel*de + CLal*exp(-(a/alp)^2)) - W)/W)^2 + (CM0 - CMal*a - CMdel*de)^2;

Gsym = [diff(f,a); diff(f,de)];


    G = double(subs(Gsym, [a, de], [alpha, deltaE]));
end


function [position, path] = Adam(g,x,tol,maxit,l,gamv,gams,epsilon)

    x = x(:);

    m = zeros(2,1);
    v = zeros(2,1);

    k = 0;

    path = x;

    while k < maxit

        k = k + 1;

        grad = g(x);

        m = gams*m + (1-gams)*grad;

        v = gamv*v + (1-gamv)*grad.^2;

        mhat = m/(1-gams^k);

        vhat = v/(1-gamv^k);

        xnew = x - l*mhat./(sqrt(vhat)+epsilon);

        path(:,end+1) = xnew;

         if norm(g(xnew)) < tol
           x = xnew;
            fprintf('%d iterations occurred\n',k)
        break;
        end
        
        x = xnew;

    end

    position = x;

end