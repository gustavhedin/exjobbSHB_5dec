
K = 10;
H = 14;
delta = 0.25;

x_values = [];
f_values = [];
derivatives = [];

for x = 1:0.01:20
    X = ADRev(x);
    f = adr_linbarrier_payoff(X,K,H,delta);
    %f = adr_sigmoidbarrier_payoff(X,K,H,delta);
    rev = chainRule(f);
    x_values = [x_values x];
    f_values = [f_values f.value];
    derivatives = [derivatives X.derivative];
end


plot(x_values,f_values,'b')
hold on
plot(x_values,derivatives,'r')