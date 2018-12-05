function X_path = geometric_brownian_adr(N,r,sigma,T,startvalue,K,H,smoothing_parameter)

% Conceptual arlgorithm:
t = (0:1:N)'/N;                   % t is the column vector [0 1/N 2/N ... 1]
W = [0; cumsum(randn(N,1))]/sqrt(N); % S is running sum of N(0,1/N) variables
t = t*T;
W = W*sqrt(T);
Y = (r-(sigma^2)/2)*t + sigma * W;
X = startvalue*exp(Y);

%X_path = adr_payoff(X,K)*adr_neg_indicator(X,K+4,2);
%X_path = adr_linbarrier_payoff(X,K,H,smoothing_parameter);%*adr_neg_payoff(X,K+3);
X_path = adr_sigmoidbarrier_payoff(X,K,14,smoothing_parameter);
%X_path = adr_barrier_dersmoothed_payoff(X,K,H,smoothing_parameter);

%%% Med ADRev-objects:
% t = (0:1:N)'/N;                   % t is the column vector [0 1/N 2/N ... 1]
% W = [0; cumsum(randn(N,1))];%/sqrt(N); % S is running sum of N(0,1/N) variables
% t = adr_mul(t,T);
% h = adr_div(T,N);
% W = adr_mul(W,adr_sqrt(h));
% Y = adr_add(adr_mul(adr_sub(r,adr_div(adr_pow(sigma,2),2)),t) , adr_mul(sigma,W));
% X = adr_mul(startvalue,adr_exp(Y));
% X_path = payoff(X,K);

end

% function P = payoff(X,K)
%     P = X;
%     for i = 1:size(P.value,1)
%         if X.value(i)-K > 0
%             P.value(i) = X.value(i)-K;
%         else
%             P.value(i) = 0;
%         end
%     end
% %    if  X.value(size(P.value,1))-K <= 0
% %        P.derivative = 'zeroFlag';
% %    end
% end