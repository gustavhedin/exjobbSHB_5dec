function [X,Delta,Gamma,Vega,t,x] = barrierpde_solve(xint,tint,grid_size,par,Tcheck,K,H)
%function [X,Delta,Gamma,Vega,t,x] = barrierpde_solve(xint,tint,grid_size,par,Tcheck,K,H)
%
% Calculates Value and Greeks for discretely monitored up out call barrier
% using Crank-Nichols finite difference PDE scheme with boundary
% corrections
%
% ---------------input-------------------------
%
% xint: [lower and upper] end point of space grid (use 0 for lower 1.5*H forupper)
% tint: [start and end point] for time grid (maturity is t_int(2))
% grid_size: [nof space points, nof time points]
% par : Parameter for dynamics [r sigma]
% Tcheck: Monitoring times for barrier excluding maturity time
%         if equal to [] only maturity will be monitored
% 
% K: Strike price
% H: Barrier level
% -------------output--------------------------
%
%  The first four varibles all are given as  
%    Var(:,1)= equal all values at start time  
%    Var(:,end)= equal all values at maturity
% X: Value at all times in time grid for all x in space grid
% Delta: sensitivty derivative for S0
% Gamma: second sensitivty derivative for S0
% Vega:  sensitivty derivative for sigma
%  t: time grid
%  x: space grid
% 
% To vizualize e.g. Value for all time use:
% mesh(x,t,X'), xlabel('Stock') , ylabel('Time'), view(9,14)
%
% Note that the value will be discontinuous at monitoring times
% the other greeks behave badly also at time points just prior to
% monitoring times
%
% (C) Magnus Wiktorsson (2018)


M=grid_size(1); N=grid_size(2);
x_min=xint(1); x_max=xint(2);
t_min=tint(1); t_max=tint(2);
% form space and time grid
x=linspace(x_min,x_max,M)';
t=linspace(t_min,t_max,N);
if ~isempty(Tcheck) % check if any monitoring prior to maturity
 [~,Im]=min(abs(Tcheck-t')); % find index of grid times closest to monitor times
else
  Im=[];  % no monitoring prior to maturity     
end 
dx=(x(2)-x(1)); % space step
dt=(t(2)-t(1)); % time step
% intialize value
X=zeros(M,N);
X(:,N)=max(x-K,0).*(x<H); % final pay off
% build matrices for finite differences 
[Mu,V]=f_ms2_mod(t(1),x,par);   
Ax=diag([ones(1,M-1)].*Mu(1:end-1)',1)+diag([-ones(1,M-1)].*Mu(2:end)',-1);
Ax(1,1:3)=[-3 4 -1]*Mu(1); % boundary correction
Ax(end,end-2:end)=[1 -4 3]*Mu(end); % boundary correction
Axx=diag([1 -2*ones(1,M-2) 1].*V')+diag([-2 ones(1,M-2)].*V(1:end-1)',1)+diag([ones(1,M-2) -2].*V(2:end)',-1);
Axx(1,3)=V(1);   
Axx(end,end-2)=V(end);
A1=(-par(1)*diag(ones(size(x)))+Axx/dx/dx+Ax/dx/2)*dt;
A2=A1;
A2((x>=H),:)=[0]; % fix barrier boundary condition
% Build exponential matrices
G1=expm(A1); % standard cond
G2=expm(A2); % special barrier cond
% loop over time points backward from maturity
for i=(N-1):-1:1
    if any(i==Im) % check if barrier time 
        X(:,i)=G2*X(:,i+1); % update solution with barrier matrix
        X(x>=H,i)=0;   % set solution to zero above barrier
    else
      X(:,i)=G1*X(:,i+1); % update according to non-barrier dynamics
    end  
end
d=diag(ones(M-1,1),1)-eye(M); % matrix for finite diff approx of Delta
% calculate Delta from X with finite diff and spline interpolation
Delta=interp1((x(1:end-1)+x(2:end))/2,d(1:end-1,:)*X/dx,x,'PCHIP');
% calculate Gamma from Delta with finite diff and spline interpolation
Gamma=interp1((x(1:end-1)+x(2:end))/2,d(1:end-1,:)*Delta/dx,x,'PCHIP');
%--------- Code for Vega-----------------
h=0.001; % signa step size
% rebuild matrices with new sigma (= sigma+h)
A1=(-par(1)*diag(ones(size(x)))+Axx/dx/dx/par(2)^2*(par(2)+h)^2+Ax/dx/2)*dt;
A2=A1;
G1=expm(A1);
G2=expm(A2);
A2((x>=H),:)=[0];
% intialize value
Xvp(:,N)=X(:,N);
% loop over time points backward from maturity
for i=(N-1):-1:1
if any(i==Im) % check if barrier time 
 Xvp(:,i)=G2*Xvp(:,i+1); % update solution with barrier matrix
  X(x>=H,i)=0;% set solution to zero above barrier
else
  Xvp(:,i)=G1*Xvp(:,i+1); % update according to non-barrier dynamics
end  
end
Vega=(Xvp-X)/h; % form finite diff approx of vega



function [m,s2]=f_ms2_mod(t,x,par)
% drift and squared diffusion for Black Scholes
% copy and update for other dynamics
 m=par(1)*x;
 s2=x.^2*par(2)^2/2;
