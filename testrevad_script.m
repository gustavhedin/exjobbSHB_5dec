% Initiate a as an ADRev object with value 2
aRev = ADRev(1);
% Initiate b as an ADRev object with value 1
bRev = ADRev(4);

% Evaluate function f(a,b) = (a+b)*(b+1), and in the same time construct 
% the calculation tree:
eRev_forward = function_testrevad(aRev,bRev);

% Use the chainRule to accumulate all derivatives from the top node down to
% the input parameters:
eRev_reverse = chainRule(eRev_forward);

% If you now type eRev_reverse.value you get f(2,1)
%     - // -      aRev.derivative you get df/da evaluated in (a,b) = (2,1)
%     - // -      aRev.derivative you get df/db evaluated in (a,b) = (2,1)

