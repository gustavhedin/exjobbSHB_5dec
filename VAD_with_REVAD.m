STRIKEPRICE = [];
DELTA = [];
VEGA = [];
RHO = [];
x_mean = [];

for p = 9:9

    % User defined parameters:
    r = 0.05;           % risk free interest rate
    sigma = 0.2;               % volatility
    Startvalue = ADRev(p);            % Starting value for the underlying asset at time 0.
    T = 1;                     % Time horizon
    N = 3;                     % # simulation points on [0,T];
    K = 10;                    % Strike price
    nbr_MC = 1;                % # of Monte Carlo simulations
    nbrMC_z = 1;               % # of samples over the barrier


    % Declaring nedded variables:
    h = T/N;                            % Stepsize
    %X = zeros(N,nbr_MC);                % Path matrix
    X = ADRev(Startvalue);      % Initiate startvalue
    %Y_delta = zeros(N,nbr_MC);         % Tangent process, Delta
    Y_delta = ADRev(0); %ones(1,nbr_MC);       % Initiate Delta to be 1.
    Y_vega = ADRev(0); %zeros(N,nbr_MC);       % Tangent process, Vega
    Y_rho = ADRev(0); %zeros(N,nbr_MC);        % Tangent process, Rho

     for n=2:N
             Z = randn(1,nbr_MC);
            % Asset
            X(n,:) = X(n-1,:) + r*X(n-1,:)*h + sigma*sqrt(h)*X(n-1,:)*Z;

            % Tangent Processes: 
            % Delta: theta = X0:
            Y_delta(n,:) = Y_delta(n-1,:) + r*h*Y_delta(n-1,:) + 0*X(n-1,:) ...
                + (Y_delta(n-1,:)*sigma*sqrt(h)+ 0*X(n-1,:))*Z; 

            % Vega: theta = sigma 
            Y_vega(n,:) = Y_vega(n-1,:) + r*h*Y_vega(n-1,:) + 0*X(n-1,:) ...
                + (Y_vega(n-1,:)*sigma*sqrt(h)+ sqrt(h)*X(n-1,:))*Z; 

            % Rho: theta = sigma 
            Y_rho(n,:) = Y_rho(n-1,:) + r*h*Y_rho(n-1,:) + h*X(n-1,:) ...
                + (Y_rho(n-1,:)*sigma*sqrt(h)+ 0*X(n-1,:))*Z; 
     end


     % Find price 
     Z = randn(1,nbr_MC);
     Xend = X(end,:)+r*h*X(end,:) + sigma*sqrt(h)*X(end,:)*Z;

     % Check asset at barrier checkpoints ?     
     % payoff = V(end)*prodsum(I1, I2 ...) , Ii = 1 if X > K else 0, at a specific barrier checkingpoint.  
     
     
     Price = (payoff(Xend,K))*exp(-r*T);   

     % Generate random nubers for the step over the barrier.
     Z = randn(nbrMC_z,nbr_MC);

     % Obtain valued for the process and derivatives on the other side of
     % the barrier:

     firstpart = repmat(X(end,:)+r*h*X(end,:),nbrMC_z,1);
     lastpart = repmat(sigma*X(end,:)*sqrt(h),nbrMC_z,1);
     % First derivatives:
        X_Tplus = firstpart + Z*lastpart;
        X_Tminus = firstpart - Z*lastpart;
        X_Tdot = firstpart;

        V_Tplus = payoff(X_Tplus,K);
        V_Tminus = payoff(X_Tminus,K);
        V_Tdot = payoff(X_Tdot,K);

        % Delta
         dmu_dtheta = repmat(Y_delta(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
         dsig_dtheta = repmat(Y_delta(end,:)*sigma*sqrt(h) + X(end,:)*0,nbrMC_z,1);
         divfactor = repmat(1/(X(end,:)*sigma*sqrt(h)),nbrMC_z,1);
         Delta = ((dmu_dtheta*(1/2)*(V_Tplus-V_Tminus)*(Z*divfactor)...
                + dsig_dtheta*(V_Tplus-2*V_Tdot+V_Tminus)*((Z.^2-1)*divfactor)))*exp((-1)*r*T);

        % Vega
         dmu_dtheta =  repmat(Y_vega(end,:)*(1+r*h) + X(end,:)*0,nbrMC_z,1);
         dsig_dtheta =  repmat(Y_vega(end,:)*sigma*sqrt(h) + X(end,:)*sqrt(h),nbrMC_z,1);
         Vega = ((dmu_dtheta*(1/2)*(V_Tplus-V_Tminus)*(Z*divfactor)...
                + dsig_dtheta*(V_Tplus-2*V_Tdot+V_Tminus)*((Z.^2-1)*divfactor)))*exp((-1)*r*T);

        % Rho
         dmu_dtheta = repmat(Y_rho(end,:)*(1+r*h) + X(end,:)*h,nbrMC_z,1);
         dsig_dtheta =  repmat(Y_rho(end,:)*sigma*sqrt(h) + X(end,:)*0,nbrMC_z,1);
         Rho = ((dmu_dtheta*(1/2)*(V_Tplus-V_Tminus)*(Z*divfactor)...
                + dsig_dtheta*(V_Tplus-2*V_Tdot+V_Tminus)*((Z.^2-1)*divfactor)-T*V_Tdot))*exp((-1)*r*T); 

        % Here, apply an AD-tequnique to get 2nd order derivatives. 
        
        X = [X;Xend];
        
        %Delta_ = chainRule(Delta);


    % Verification:------------------------------------------------------------

%    tt = 0;
%    d1 = (ln(Startvalue/K) + (r+0.5*sigma^2)*(T-tt))./(sigma*sqrt(T-tt));
%    d2 = d1 - sigma*sqrt(T-tt);

%    delta_verification = cdf('Normal',d1,0,1);
%    vega_verification = Startvalue.*pdf('Normal',d1,0,1).*sqrt(T-tt);
%    rho_verification = K*T*exp(-r*T)*cdf('Normal',d2,0,1);
%    price_verification = Startvalue.*cdf('Normal',d1,0,1) - exp(-r*(T-tt)).*K.*cdf('Normal',d2,0,1);

    
   STRIKEPRICE = [STRIKEPRICE Startvalue];
   DELTA = [DELTA Delta];
   VEGA = [VEGA Vega];
   RHO = [RHO Rho];
    
%    MEAN = mean(X.value(end));
%    x_mean = [x_mean MEAN(end)];

end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure
    %plot(STRIKEPRICE,payoff(X.value,K),'y')
    %hold on
    %plot(STRIKEPRICE,DELTA,'o--r')
    %plot(STRIKEPRICE,VEGA,'o--m')
    %plot(STRIKEPRICE,RHO,'o--b')
    %legend('PAYOFF','DELTA','VEGA','RHO')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
function P = payoff(X,K)
   %P=max(x-K,0);
   %P = (x < 12).* max(x-K,0); % Strike K, Barrier at 12
   %P = adr_payoff(X,K)*adr_neg_smoothed(X,K+2,0.01);
   P = adr_sigmoidbarrier_payoff(X,K,14,0.05);
end


