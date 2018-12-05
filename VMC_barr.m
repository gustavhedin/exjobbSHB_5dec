function [val,valD,valV,valR,valT,var,varD,varV,varR,varT]...
= VMC_barr(r,sigma,T,S0,B,K,M,N,d,T1)
% vibrato computes VMC values of 'value','delta','vega' in digital call at
% strike K and expiry T case with discrete up and out barrier at time T1
% Uses antithetic vars
%
% Black-Scholes model, digital call case
% Euler-Maruyama Approximation
%
% r - interest rate
% sigma - volatility
% T - time to expiry
% S0 - asset initial value(s)
% K - strike price(s)
% M - number of samples
% N - number of timesteps
% d - number of final samples for Z
% T1 - intermediate time where barrier constraint is checked
if nargin ~= 10
    error('wrong number of args')
end

h=T/N; %timestep
N1=floor(T1/h);%index of barrier
S=S0*ones(1,M); % array of underlying value
dS=ones(1,M); %array of deltas
vS=zeros(1,M); %array of vegas
rS=zeros(1,M); %array of rhos
tS=zeros(1,M); %array of thetas

for n=1:1:N1-1 % till step N1-1, is normal euler simulation
    dW=sqrt(h)*randn(1,M);
    dS=dS.*(1+r*h+sigma*dW); % delta evol
    vS=vS.*(1+r*h+sigma*dW)+S.*dW; % vega evol
    rS=rS.*(1+r*h+sigma*dW)+h*S; % rho evol
    tS=tS.*(1+r*h+sigma*dW)+S.*(r/N+sigma/(2*N*h)*dW); %theta evol
    S=S.*(1+r*h+sigma*dW); % value evol
end

%storing values @ N1-1 for final evaluation
S1=S;dS1=dS;vS1=vS;rS1=rS;tS1=tS;

%skipping barrier time T1 - size 2*h
dW=sqrt(2*h)*randn(1,M);
dS=dS.*(1+r*2*h+sigma*dW); % delta evol
vS=vS.*(1+r*2*h+sigma*dW)+S.*dW; % vega evol
rS=rS.*(1+r*2*h+sigma*dW)+2*h*S; % rho evol
tS=tS.*(1+r*2*h+sigma*dW)+S.*(r/N+sigma/(2*N*2*h)*dW); %theta evol
S=S.*(1+r*2*h+sigma*dW); % value evol

% Brownian Bridge - drifts, vols @ T1
muw1=(S1+S)/2; %array of S-drifts @T1 (dim M)
sigmaw1=S1*(sigma*sqrt(h/2)); %array of vols - includes sqrt(h/2) term
dmuw1=(dS1+dS)/2; % delta-drifts
dsigmaw1=dS1*(sigma*sqrt(h/2)); % delta-vols
vmuw1=(vS1+vS)/2; % vega-drifts
vsigmaw1=(vS1*(sigma*sqrt(h/2))+S1*sqrt(h/2)); % vega-vols
rmuw1=(rS1+rS)/2; % vega-drifts
rsigmaw1=rS1*(sigma*sqrt(h/2)); % rho-vols
tmuw1=(tS1+tS)/2; % theta-drifts
tsigmaw1=tS1*sigma*sqrt(h/2)+sigma/(2*N*sqrt(h*2))*S1; % theta-vols

for n=N1+1:1:N-1 % again normal euler simulation
    dW=sqrt(h)*randn(1,M);
    dS=dS.*(1+r*h+sigma*dW); % delta evol
    vS=vS.*(1+r*h+sigma*dW)+S.*dW; % vega evol
    rS=rS.*(1+r*h+sigma*dW)+h*S; % rho evol
    tS=tS.*(1+r*h+sigma*dW)+S.*(r/N+sigma/(2*N*h)*dW); % theta evol
    S=S.*(1+r*h+sigma*dW); % value evol
end

%S=S;dS=dS;vS=vS;rS=rS;tS=tS; % stored values at final euler-time N-1

% final vibrato step, N-1 to N
Z1=randn(d,M);%d=d sims, M= M paths
Z=randn(d,M);
% drifts, vol @T
muw=S*(1+r*h); % array of S-drifts @T (dim M)
sigmaw=S*(sigma*sqrt(h)); % array of vols - includes sqrt(h) term
dmuw=dS*(1+r*h); % delta-drifts
dsigmaw=dS*(sigma*sqrt(h)); % delta-vols
vmuw=vS*(1+r*h); % vega-drifts
vsigmaw=(vS*(sigma*sqrt(h))+S*sqrt(h)); % vega-vols
rmuw=rS*(1+r*h)+h*S; % vega-drifts
rsigmaw=rS*(sigma*sqrt(h)); % rho-vols
tmuw=tS*(1+r*h)+r*S/N; % theta-drifts
tsigmaw=tS*sigma*sqrt(h)+sigma/(2*N*sqrt(h))*S; % theta-vols

%value
SNmatP1=ones(d,1)*muw1+(ones(d,1)*sigmaw1).*Z1; %size d*M
SNmatM1=ones(d,1)*muw1-(ones(d,1)*sigmaw1).*Z1;
SNmat1=ones(d,1)*muw1;

SNmatP=ones(d,1)*muw+(ones(d,1)*sigmaw).*Z; %size d*M
SNmatM=ones(d,1)*muw-(ones(d,1)*sigmaw).*Z;
SNmat=ones(d,1)*muw;

VNmatP=payoff_barrier(SNmatP1,SNmatP,B,K,r,T1,T); % size d*M, f(mu+Cz)
VNmatM=payoff_barrier(SNmatM1,SNmatM,B,K,r,T1,T);
VNmatC=payoff_barrier(SNmat1,SNmat,B,K,r,T1,T);
VNmat=1/2*(VNmatP+VNmatM); %used in rho calculation when differentiating discount term
dfdr_mat=-T*VNmat;
%same for theta calculation
dfdT_mat=-r*VNmat;
vals=sum(VNmat)/d; %vector of values averaged wrt Z
val=sum(vals)/M; % average of (average over Z) over all trajectories
varz=1/d*sum(VNmat.^2)-(1/d*sum(VNmat)).^2; % vector of variances wrt Z
% vector of variance wrt W of avgd value (over Z)
varw=1/M*sum(vals.^2)-(val)^2;
var=1/M*varw+1/(M*d)*1/M*sum(varz); %variance of value

%delta
% matrix d*M, Z/sig*1/2(fp-fm)
matmu1=Z1.*(ones(d,1)*(1./sigmaw1))*1/2.*(VNmatP-VNmatM);
% matrix d*M, Z/sig*1/2(fp-fm)
matmu=Z.*(ones(d,1)*(1./sigmaw))*1/2.*(VNmatP-VNmatM);
% matrix d*M, (Z1^2-1)/sig1*1/4*(...
matsi1=(Z1.^2-1).*(ones(d,1)*(1./sigmaw1))*1/4.*(VNmatP-2*VNmatC+VNmatM);
% matrix d*M, (Z^2-1)/sig*1/4*(...
matsi=(Z.^2-1).*(ones(d,1)*(1./sigmaw))*1/4.*(VNmatP-2*VNmatC+VNmatM);
%valsD=1/d*(dmuw1.*sum(matmu1)+dmuw.*sum(matmu))+1/d*...
%(2*sigmaw1.*dsigmaw1.*sum(matsi1)+2*sigmaw.*dsigmaw.*sum(matsi))%array Ez()
%is not possible to compute variance later !
%matrix d*M to be avgd over Z then W to get estimator
mat=(ones(d,1)*dmuw1).*matmu1+(ones(d,1)*dmuw).*matmu+...
(ones(d,1)*(2*sigmaw1.*dsigmaw1)).*matsi1+(ones(d,1)*...
(2*sigmaw.*dsigmaw)).*matsi;
valsD=1/d*sum(mat); % array Ez()
valD=sum(valsD)/M; %delta value =Ew(Ez())
varD=(1/M)*(1/M*sum(valsD.^2)-(valD)^2)...
+1/(M*d)*1/M*sum(sum(mat.^2)/d-valsD.^2); %variance of delta

%vega
%we keep same matmu and matsi - not dependent on sensitivity type
mat=(ones(d,1)*vmuw1).*matmu1+(ones(d,1)*vmuw).*matmu...
+(ones(d,1)*(2*sigmaw1.*vsigmaw1)).*matsi1...
+(ones(d,1)*(2*sigmaw.*vsigmaw)).*matsi; %matrix d*M
valsV=1/d*sum(mat); % array Ez()
valV=1/M*sum(valsV); % vega value =Ew(Ez())
varV=(1/M)*(1/M*sum(valsV.^2)-(valV)^2)...
+1/(M*d)*1/M*sum(sum(mat.^2)/d-valsV.^2); %variance of vega

%rho
mat=(ones(d,1)*rmuw1).*matmu1+(ones(d,1)*rmuw).*matmu...
+(ones(d,1)*(2*sigmaw1.*rsigmaw1)).*matsi1...
+(ones(d,1)*(2*sigmaw.*rsigmaw)).*matsi; %matrix d*M
valsR=1/d*sum(mat+dfdr_mat); % array Ez()
valR=1/M*sum(valsR); % vega value =Ew(Ez())
varR=(1/M)*(1/M*sum(valsR.^2)-(valR)^2)...
+1/(M*d)*1/M*sum(sum((mat+dfdr_mat).^2)/d-valsR.^2); %variance of rho

%theta
mat=(ones(d,1)*tmuw1).*matmu1+(ones(d,1)*tmuw).*matmu...
+(ones(d,1)*(2*sigmaw1.*tsigmaw1)).*matsi1...
+(ones(d,1)*(2*sigmaw.*tsigmaw)).*matsi; %matrix d*M
valsT=1/d*sum(mat+dfdT_mat); % array Ez()
valT=sum(valsT)/M; %delta value =Ew(Ez())
varT=(1/M)*(1/M*sum(valsT.^2)-(valT)^2)...
+1/(M*d)*1/M*sum(sum((mat+dfdT_mat).^2)/d-valsT.^2); %variance of theta
end

function V = payoff_barrier(S1,S,B,K,r,T1,T)
% returns discounted payoff of digital call @(K,T)
% with discrete barrier @(B,T1)
% takes as two first argument values @ T1 and T
SM=max(S1,S);
V=exp(-r*T)*(S>K).*(SM<B);
end