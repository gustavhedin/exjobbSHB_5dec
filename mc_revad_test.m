clear all

%smoothing_parameter = [1 0.5 0.1 0.05 0.01 0.005];
smoothing_parameter = [5 7.5 10 12.5 15 17.5]; % Use this for sigmoid
%smoothing_parameter = [1 0.5 0.25 0.1 0.05 0.01]; % Use this for derivative smoothing.
DELTA_mat = [];
RHO_mat = [];
VEGA_mat = [];
THETA_mat = [];
OPTIONPRICE_mat = [];

for q=1:6

UNDERLYING = [];
OPTIONPRICE = [];
DELTA = [];
VEGA = [];
RHO = [];
THETA = [];

    for s=5:20    
        % inputs:
        % volatility sigma
        % interest rate r
        % maturity time tau
        start_der = 0;
        sigma_der = 0;
        r_der = 0;
        tau_der = 0;

        K = 10; % Strike price
        H = 14; % Barrrier

        % constants:
        N = 12; % # steps
        nbrMC = 10000; % # MC simulations
        X = zeros(N+1,1);

        for i = 1:nbrMC
            startvalue = ADRev(s); 
            sigma      = ADRev(0.2);
            r       =   ADRev(0.05);
            tau      =     ADRev(1);

            eRev_forward = geometric_brownian_adr(N,r,sigma,tau,startvalue,K,H,smoothing_parameter(q));
            reverse_sweep = chainRule(eRev_forward); 
            X = X + reverse_sweep.value(end);
            start_der = start_der + startvalue.derivative(end);
            sigma_der = sigma_der + sigma.derivative(end);
            r_der = r_der + r.derivative(end);
            tau_der = tau_der + tau.derivative(end);
        end

        X = X/nbrMC;
        start_der = start_der/nbrMC;
        sigma_der = sigma_der/nbrMC;
        r_der = r_der/nbrMC;
        tau_der = tau_der/nbrMC;

        optionprice = X(end)*exp(-r.value*tau.value);

        UNDERLYING = [UNDERLYING startvalue.value];
        OPTIONPRICE = [OPTIONPRICE optionprice(1)];
        DELTA = [DELTA start_der];
        VEGA = [VEGA sigma_der];
        RHO = [RHO r_der];
        THETA = [THETA tau_der];

    end

    RHO = (RHO-tau.value*OPTIONPRICE)*exp(-r.value*tau.value);
    DELTA = DELTA*exp(-r.value*tau.value);
    VEGA = VEGA*exp(-r.value*tau.value);
    THETA = (THETA-r.value*OPTIONPRICE)*exp(-r.value*tau.value);
    
    DELTA_mat = [DELTA_mat; DELTA];
    RHO_mat = [RHO_mat; RHO];
    VEGA_mat = [VEGA_mat; VEGA];
    THETA_mat = [THETA_mat; THETA];
    OPTIONPRICE_mat = [OPTIONPRICE_mat; OPTIONPRICE];

end

%%

figure
plot(UNDERLYING,OPTIONPRICE_mat(1,:),'o--r')
hold on
plot(UNDERLYING,OPTIONPRICE_mat(2,:),'o--g')
plot(UNDERLYING,OPTIONPRICE_mat(3,:),'o--b')
plot(UNDERLYING,OPTIONPRICE_mat(4,:),'o--m')
plot(UNDERLYING,OPTIONPRICE_mat(5,:),'o--c')
plot(UNDERLYING,OPTIONPRICE_mat(6,:),'o--k')
title('OPTIONPRICE')
xlabel('Underlying')
ylabel('Value')
legend('smoothing=5','smoothing=7.5', 'smoothing=10', 'smoothing=12.5', 'smoothing=15', 'smoothing=17.5','Location', 'Best')
grid on

figure
plot(UNDERLYING,DELTA_mat(1,:),'o--r')
hold on
plot(UNDERLYING,DELTA_mat(2,:),'o--g')
plot(UNDERLYING,DELTA_mat(3,:),'o--b')
plot(UNDERLYING,DELTA_mat(4,:),'o--m')
plot(UNDERLYING,DELTA_mat(5,:),'o--c')
plot(UNDERLYING,DELTA_mat(6,:),'o--k')
title('DELTA')
xlabel('Underlying')
ylabel('Value')
legend('smoothing=5','smoothing=7.5', 'smoothing=10', 'smoothing=12.5', 'smoothing=15', 'smoothing=17.5','Location', 'Best')
grid on

figure
plot(UNDERLYING,VEGA_mat(1,:),'o--r')
hold on
plot(UNDERLYING,VEGA_mat(2,:),'o--g')
plot(UNDERLYING,VEGA_mat(3,:),'o--b')
plot(UNDERLYING,VEGA_mat(4,:),'o--m')
plot(UNDERLYING,VEGA_mat(5,:),'o--c')
plot(UNDERLYING,VEGA_mat(6,:),'o--k')
title('VEGA')
xlabel('Underlying')
ylabel('Value')
legend('smoothing=5','smoothing=7.5', 'smoothing=10', 'smoothing=12.5', 'smoothing=15', 'smoothing=17.5')
grid on

figure
plot(UNDERLYING,RHO_mat(1,:),'o--r')
hold on
plot(UNDERLYING,RHO_mat(2,:),'o--g')
plot(UNDERLYING,RHO_mat(3,:),'o--b')
plot(UNDERLYING,RHO_mat(4,:),'o--m')
plot(UNDERLYING,RHO_mat(5,:),'o--c')
plot(UNDERLYING,RHO_mat(6,:),'o--k')
title('RHO')
xlabel('Underlying')
ylabel('Value')
legend('smoothing=5','smoothing=7.5', 'smoothing=10', 'smoothing=12.5', 'smoothing=15', 'smoothing=17.5')
grid on

figure
plot(UNDERLYING,THETA_mat(1,:),'o--r')
hold on
plot(UNDERLYING,THETA_mat(2,:),'o--g')
plot(UNDERLYING,THETA_mat(3,:),'o--b')
plot(UNDERLYING,THETA_mat(4,:),'o--m')
plot(UNDERLYING,THETA_mat(5,:),'o--c')
plot(UNDERLYING,THETA_mat(6,:),'o--k')
title('THETA')
xlabel('Underlying')
ylabel('Value')
legend('smoothing=5','smoothing=7.5', 'smoothing=10', 'smoothing=12.5', 'smoothing=15', 'smoothing=17.5')
grid on

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure
%     plot(UNDERLYING,OPTIONPRICE,'o--c')
%     hold on
%     plot(UNDERLYING,DELTA,'o--r')
%     plot(UNDERLYING,VEGA,'o--m')
%     plot(UNDERLYING,RHO,'o--b')
%     plot(UNDERLYING,THETA,'o--g')
%     legend('OPTIONPRICE','DELTA','VEGA','RHO','THETA')
 

%    a = (log(UNDERLYING/K) - (r.value+0.5*sigma.value^2)*(tau.value))./(sigma.value*sqrt(tau.value));
%    b = (log(UNDERLYING/H) - (r.value+0.5*sigma.value^2)*(tau.value))./(sigma.value*sqrt(tau.value));
    
%    price_ver = -exp(-r.value*tau.value)*(-K*(cdf('Normal',b,0,1)-cdf('Normal',a,0,1))...
%         + UNDERLYING.*exp((r.value-0.5*sigma.value^2)*tau.value+sigma.value*sqrt(tau.value).*UNDERLYING).*...
%        (cdf('Normal',b-sigma.value*sqrt(tau.value),0,1)-cdf('Normal',a-sigma.value*sqrt(tau.value),0,1)));
    
%    delta_ver = cdf('Normal',b-sigma.value*sqrt(tau.value),0,1)-cdf('Normal',a-sigma.value*sqrt(tau.value),0,1) - ...
%        pdf('Normal',b-sigma.value*sqrt(tau.value),0,1)*(1/(sigma.value*sqrt(tau.value))) + ...
%        K*exp(-r.value*tau.value)*(1./(UNDERLYING.*sigma.value*sqrt(tau.value))).*pdf('Normal',b,0,1);
    
%    plot(UNDERLYING,price_ver,'k')
%    plot(UNDERLYING,delta_ver,'g')

%    d1 = (log(UNDERLYING./K) + (r.value+0.5*sigma.value^2)*(tau.value))./(sigma.value*sqrt(tau.value));
%    d2 = d1 - sigma.value*sqrt(tau.value);

%    delta_verification = cdf('Normal',d1,0,1);
%    vega_verification = UNDERLYING.*pdf('Normal',d1,0,1).*sqrt(tau.value);
%    rho_verification = K*tau.value*exp(-r.value*tau.value).*cdf('Normal',d2,0,1);
    
%    plot(UNDERLYING,delta_verification,'k')
%    plot(UNDERLYING,vega_verification,'y')
%    plot(UNDERLYING,rho_verification,'g')
%    legend('OPTIONPRICE','DELTA','VEGA','RHO','delta_ver','vega_ver','rho_ver')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

