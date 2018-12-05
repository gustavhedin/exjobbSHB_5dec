UNDERLYING = [];
VALUE = [];
DELTA = [];
VEGA = [];
THETA = [];
RHO = [];

for s = 5:0.1:15
    [val,valD,valV,valR,valT,var,varD,varV,varR,varT] = VMC_barr(0.05,0.2,1,s,12,10,100000,100,10,0.7);
    UNDERLYING = [UNDERLYING s];
    VALUE = [VALUE val];
    DELTA = [DELTA valD];
    VEGA = [VEGA valV];
    THETA = [THETA valT];
    RHO = [RHO valR];
end
%% 
figure
plot(UNDERLYING,VALUE,'o--g')
hold on
plot(UNDERLYING,DELTA,'o--r')
plot(UNDERLYING,VEGA,'o--m')
plot(UNDERLYING,THETA,'o--c')
plot(UNDERLYING,RHO,'o--b')
legend('PRICE','DELTA','VEGA','THETA','RHO')
xlabel('Underlying')
ylabel('value')
grid on
