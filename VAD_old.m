close all

% User defined parameters:
r = 0.05;           % Risk free interest rate
sigma = 0.2;        % Volatility
Startvalue = 8;     % Starting value for the underlying asset at time 0.
T = 10;             % Time horizon
N = 100;            % # simulation points on [0,T];
K = 10;             % Strike price
nbr_MC = 1000;      % # of Monte Carlo simulations
nbrMC_z = 100;      % # of samples over the barrier

% Declaring nedded variables:
h = T/(N-1);
t = linspace(0,T,N+1);
X = zeros(N,nbr_MC);
X(1,:) = Startvalue*ones(1,nbr_MC);

Y_delta = zeros(N,nbr_MC);
Y_delta(1,:) = ones(1,nbr_MC);

Y_vega = zeros(N,nbr_MC);

X_end = zeros(1,nbr_MC);
Delta_vec = zeros(1,nbr_MC);
Vega_vec = zeros(1,nbr_MC);

Delta_pathwise = zeros(N,nbr_MC);
Vega_pathwise = zeros(N,nbr_MC);

for i=1:nbr_MC
    %W = [0; cumsum(randn(N,1))];
    Z = randn(N,1);
    for n=2:N
        % Asset
        X(n,i) = X(n-1,i) + r*h*X(n-1,i) + sigma*sqrt(h)*X(n-1,i)*Z(n);

        % Tangent Processes:
        % Delta: theta = X0:
        Y_delta(n,i) = Y_delta(n-1,i) + r*h*Y_delta(n-1,i) + 0*X(n-1,i) ...
            + (Y_delta(n-1,i)*sigma*sqrt(h)+ 0*X(n-1,i))*Z(n); 
        
        % Vega: heta = sigma 
        Y_vega(n,i) = Y_vega(n-1,i) + r*h*Y_vega(n-1,i) + 0*X(n-1,i) ...
            + (Y_vega(n-1,i)*sigma*sqrt(h)+ sqrt(h)*X(n-1,i))*Z(n); 
    end
    % Generate random nubers for the step over the barrier.
    Z = randn(1,nbrMC_z);
    
    % Obtain valued for the process and derivatives on the other side of
    % the barrier:
    Delta = 0;
    Vega = 0;
    for k=1:nbrMC_z
    
        X_end(i) = X_end(i) + X(end,i)*(1+r*h + sigma*sqrt(h)*Z(k));
      
        % First derivatives:
        X_Tplus = X(end,i)+r*h*X(end,i) + sigma*X(end,i)*sqrt(h)*Z(k);
        X_Tminus = X(end,i)+r*h*X(end,i) - sigma*X(end,i)*sqrt(h)*Z(k);
        X_Tdot = X(end,i)+r*h*X(end,i);

        V_Tplus = payoff(X_Tplus,K);
        V_Tminus = payoff(X_Tminus,K);
        V_Tdot = payoff(X_Tdot,K);

        % Delta
        dmu_dtheta = Y_delta(end,i)*(1+r*h) + X(end,i)*0;
        dsig_dtheta = Y_delta(end,i)*sigma*sqrt(h) + X(end,i)*0;

        Delta = Delta + dmu_dtheta*(1/2)*(V_Tplus-V_Tminus)*(Z(k)/(X(end,i)*sigma*sqrt(h)))...
            + dsig_dtheta*(1/2)*( V_Tplus-2*V_Tdot+V_Tminus )*((Z(k)^2-1)/(X(end,i)*sigma*sqrt(h)));

        % Vega
        dmu_dtheta = Y_vega(end,i)*(1+r*h) + X(end,i)*0;
        dsig_dtheta = Y_vega(end,i)*sigma*sqrt(h) + X(end,i)*sqrt(h);

        Vega = Vega + dmu_dtheta*(1/2)*(V_Tplus-V_Tminus)*(Z(k)/(X(end,i)*sigma*sqrt(h)))...
            + dsig_dtheta*(1/2)*( V_Tplus-2*V_Tdot+V_Tminus )*((Z(k)^2-1)/(X(end,i)*sigma*sqrt(h)));
    end
    
    % Here: apply an AD method on Delta, Vega etc. to obtain 2nd order
    % derivatives.
    
    % Averaging over nbrMC_z:
    X_end(i) = (1/nbrMC_z)*X_end(i);
    Delta_vec(i) = (1/nbrMC_z)*Delta;
    Vega_vec(i) = (1/nbrMC_z)*Vega;
    
    % Place current value of Delta and vega in the corresponding matricies
    Delta_pathwise(:,i) = Y_delta(:,i).*(payoff_delta(X(:,i),K))';
    Vega_pathwise(:,i) = Y_vega(:,i).*(payoff_vega(X(:,i),K))';
    
end


% Calculate and Average Pathwise Delta/Vega
Delta_pathwise = (1/nbr_MC)*sum(Delta_pathwise');
Vega_pathwise = (1/nbr_MC)*sum(Vega_pathwise');

% Append X-vector with the last time-step:
X = [X; X_end];

% Average the values over all Monte Carlo paths:
X_avg = (1/nbr_MC)*sum(X');

Delta = (1/nbr_MC)*sum(Delta_vec);
Vega = (1/nbr_MC)*sum(Vega_vec);

% Option value at t=T:
V_end = X_avg(end);

%Discounted option value:
t_ = fliplr(t);
V = V_end*exp(-r*t_);

% Plots:
figure
plot(t,X_avg,'r')
hold on
plot(t,V,'c')
plot(t,[Delta_pathwise Delta],'+g')
plot(t,[Vega_pathwise Vega],'+b')
title('Vibrato MC')
xlabel('time')
ylabel('Value')
legend('average of simulated GBM','discounted option value','Delta','Vega')
grid on

figure
plot(t,payoff_eval(X_avg',K),'r')
hold on
plot(t,[Delta_pathwise Delta],'+--g')
plot(t,[Vega_pathwise Vega],'+--b')
plot(t,V-K,'c')

%% Verification:

X_exact = zeros(1,101);
for i=1:10000
X_exact = X_exact + geometric_browninan_exact(N,r,sigma,T)';
end
X_exact = X_exact/10000;

tt = 0;
d1 = (log(K/Startvalue) + (r+0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
d2 = d1 - sigma*sqrt(T-t);

delta_verification = cdf('Normal',d1,0,1);
vega_verification = X_exact.*pdf('Normal',d1,0,1).*sqrt(T-t);
price_verification = X_exact.*cdf('Normal',d1,0,1) - exp(-r*(T-t)).*K.*cdf('Normal',d2,0,1);

plot(t,price_verification,'m')
plot(t,delta_verification,'y')
plot(t,vega_verification,'k')
title('Vibrato MC v.s. analytical expressions')
xlabel('time')
ylabel('Value')
legend('simulated payoff','simulated Delta','simulated Vega','Discounted sumulated option price',...
    'analytical option price','analytical Delta','analytical Vega')
grid on


%-------------------------------------------------------------------------%
function P = payoff(x,K)
    if (x-K)>0
        P = x-K;
    else 
        P = 0;
    end
end

function P = payoff_eval(x,K)
    P = zeros(1,size(x,1));
    for i=1:length(P)
        if (x(i)-K)>0
            P(i) = x(i)-K;
        else 
            P(i) = 0;
        end
    end
end

function P = payoff_delta(x,K)
    P = zeros(1,size(x,1));
    for i=1:length(P)
        if (x(i)-K)>0
            P(i) = 1;
        else 
            P(i) = 0;
        end
    end
end

function P = payoff_vega(x,K)
    P = zeros(1,size(x,1));
    for i=1:length(P)
        if (x(i)-K)>0
            P(i) = 1;
        else 
            P(i) = 0;
        end
    end
end
