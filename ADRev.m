classdef ADRev <  handle
    % This is a class defineing ADRev objects, which are used to calculate
    % function value and derivatives using the reverso mode of Automatic
    % differentiation. 
    % 
    % Implemented by Gustav Hedin, 2018. With inspiration from Laksh Gupta.
    % 
    % Methods in this class are: 
    %   addition          adr_sub()
    %   subtraction       adr_sub()
    %   multiplication    adr_mul()
    %   division          adr_div()
    %   natural exponent  adr_exp() 
    %   natural logarithm adr_ln()  
    %   sqrare root       adr_sqrt()
    %   step function     adr_step()

    properties
        value
        derivative
        derivativeOp
        parents
    end
    
    methods
        % Constructor: 
        function obj = ADRev(val,der)
            if nargin == 1
                obj.value = val;
                obj.derivative = 0;
                obj.derivativeOp = @ adr_constD;
                obj.parents = [];
            else
                obj.value = val;
                obj.derivative = der;
                obj.derivativeOp = @ adr_constD;
                obj.parents = [];
            end
        end
        
        function modified_parents = adr_constD(~, ~)
            %adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            %adNodes(2).derivative = adNodes(2).derivative + prevDerivative;
            modified_parents = [];
        end
        
        % Fuctions for defineing the basic calculation operations for
        % ADRev-objects:
        
        % Addition
        function result = plus(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value + y.value);
            result.derivativeOp = @ adr_addD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_addD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            adNodes(2).derivative = adNodes(2).derivative + prevDerivative;
            modified_parents = adNodes;
        end
        
        % Subtraction
        function result = minus(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value - y.value);
            result.derivativeOp = @ adr_subD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_subD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            adNodes(2).derivative = adNodes(2).derivative - prevDerivative;
            modified_parents = adNodes;
        end
        
        % Multiplication:
        function result = mtimes(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value .* y.value);
            result.derivativeOp = @ adr_mulD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_mulD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*adNodes(2).value;
            adNodes(2).derivative = adNodes(2).derivative + prevDerivative.*adNodes(1).value;
            modified_parents = adNodes;
        end
        
        % Division
        function result = mrdivide(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value ./ y.value);
            result.derivativeOp = @ adr_divD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_divD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative*(1/adNodes(2).value);
            adNodes(2).derivative = adNodes(2).derivative + ...
                prevDerivative.*(-adNodes(1).value./adNodes(2).value.^2);
            modified_parents = adNodes;
        end
        
        % Natural exponent
        function result = exp(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(exp(x.value));
            result.derivativeOp = @ adr_expD;
            result.parents = x;
        end
        
        function modified_parents = adr_expD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*exp(adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % Natural logarithm
        function result = ln(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(log(x.value));
            result.derivativeOp = @ adr_lnD;
            result.parents = x;
        end
        
        function modified_parents = adr_lnD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*(1/adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % Power
        function result = mpower(x,p)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(x.value.^p);
            switch p
                case 2
                    result.derivativeOp = @ adr_powD_2;
                case 3
                    result.derivativeOp = @ adr_powD_3;
                case 4
                    result.derivativeOp = @ adr_powD_4;
                otherwise
                    error('x^p for p<4 not implemented yet. use a lower p or extend implementation in ADRev')
            end
            result.parents = x;
        end
        
        function modified_parents = adr_powD_2(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + 2*prevDerivative.*adNodes(1).value;
            modified_parents = adNodes;
        end
        
        function modified_parents = adr_powD_3(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + 3*prevDerivative.*adNodes(1).value^2;
            modified_parents = adNodes;
        end
        
        function modified_parents = adr_powD_4(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + 4*prevDerivative.*adNodes(1).value^3;
            modified_parents = adNodes;
        end
        
        % sine
        function result = sin(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(sin(x.value));
            result.derivativeOp = @ adr_sinD;
            result.parents = x;
        end
        
        function modified_parents = adr_sinD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*cos(adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % cosine
        function result = cos(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(cos(x.value));
            result.derivativeOp = @ adr_cosD;
            result.parents = x;
        end
        
        function modified_parents = adr_cosD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*sin(adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % square root
        function result = sqrt(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(sqrt(x.value));
            result.derivativeOp = @ adr_sqrtD;
            result.parents = x;
        end
        
        function modified_parents = adr_sqrtD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*1/(2*sqrt(adNodes(1).value));
            modified_parents = adNodes;
        end
        
        
        % step-function
        function result = step(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            v = x.value;
            unistep = v>=0;
            result = ADRev(unistep);
            result.derivativeOp = @ adr_stepD;
            result.parents = x;
        end
        
        function modified_parents = adr_stepD(prevDerivative, adNodes)
            a = 0.001; % S?tta denna parameter n?gon annanstans? i properties?
            new_derivative = myDirac(adNodes(1).value);
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*new_derivative;
            modified_parents = adNodes;
        end
        
        % Payoff for European vanilla call-option
        function result = adr_payoff(X,K)
            P = zeros(size(X.value,1),size(X.value,2));
            
            if size(P,1) > size(P,2)
                d = size(P,1);
            else
                d = size(P,2);
            end
            
            for i=d:d
                P(i) = max(0,X.value(i)-K);
            end
            result = ADRev(P);
            result.derivativeOp = @ adr_payoffD;
            result.parents = [X K];
        end
        
        function modified_parents = adr_payoffD(prevDerivative, adNodes)
            P = zeros(size(adNodes(1).value,1),size(adNodes(1).value,2));
            % fin out if rowvector or colonvector
                if size(P,1) > size(P,2)
                    d = size(P,1);
                else
                    d = size(P,2);
                end
            %--------------------------------------%
            for i=d:d
                if adNodes(1).value(i)-adNodes(2).value > 0                
                    P(i) = 1;
                else
                    P(i) = 0;
                end
            end
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        %-----------------------------------------------------------------%
        % Payoff for European barrier call-option
        function result = adr_linbarrier_payoff(X,K,H,delta)
            %P = zeros(size(X.value,1),size(X.value,2));
            if X.value(end)-K < 0
                P = 0;
            elseif X.value(end) < ((K+2)*delta+H)/(delta+1)
                P = X.value(end)-K;
            elseif X.value(end)-K < H-K+2*delta
                %P = 12*(1/delta) - (3/delta)*(X.value(end)-K);
                P = 2+(H/delta) - (1/delta)*(X.value(end));
            else
                P = 0;
            end

            result = ADRev(P);
            result.derivativeOp = @ adr_linbarrier_payoffD;
            result.parents = [X K H delta];
        end
        
        function modified_parents = adr_linbarrier_payoffD(prevDerivative, adNodes)
           %P = zeros(size(adNodes(1).value,1),size(adNodes(1).value,2));
           if adNodes(1).value(end)-adNodes(2).value < 0
               P = 0;
           elseif adNodes(1).value(end) < ((adNodes(2).value+2)*adNodes(4).value+adNodes(3).value)/(adNodes(4).value+1) 
               P = 1;
           elseif adNodes(1).value(end)-adNodes(2).value < adNodes(3).value-adNodes(2).value+2*adNodes(4).value
               P = -1/adNodes(4).value;
           else
               P = 0;
           end
           adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
           modified_parents = adNodes;
        end
        
        %-----------------------------------------------------------------%
        % Payoff for European barrier call-option
        function result = adr_curvedbarrier_payoff(X,K,H,d)
            %P = zeros(size(X.value,1),size(X.value,2));
            %for i=1:max(size(X.value,1),size(X.value,2))
                P = max(0,X.value(end)-K)*0.5*erfc((X.value(end)-H)/(d*sqrt(2)));
            %end
            result = ADRev(P);
            result.derivativeOp = @ adr_curvedbarrier_payoffD;
            result.parents = [X K H d];
        end
        
        function modified_parents = adr_curvedbarrier_payoffD(prevDerivative, adNodes)
           %P = zeros(size(adNodes(1).value,1),size(adNodes(1).value,2));
            if adNodes(1).value(end)-adNodes(2).value < 0
                P = 0;
            else
                P = 0.5*erfc((adNodes(1).value-adNodes(3).value)/(adNodes(4).value*sqrt(2)))-(adNodes(1).value(end)-adNodes(2).value)*2/(sqrt(pi)*adNodes(4).value)*exp(-1*(adNodes(1).value-adNodes(3).value).^2/(2*adNodes(4).value));
            end
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        %-----------------------------------------------------------------%
        % Payoff for European barrier call-option (sigmoid smoothing)
        function result = adr_sigmoidbarrier_payoff(X,K,H,d)
            if X.value(end)-K < 0
                P = 0;
            else
                sigmoid = 1/(1+exp(d*(X.value(end)-H)));
                P = (X.value(end)-K)*sigmoid;
            end
            result = ADRev(P);
            result.derivativeOp = @ adr_sigmoidbarrier_payoffD;
            result.parents = [X K H d];
        end
        
        function modified_parents = adr_sigmoidbarrier_payoffD(prevDerivative, adNodes)
           if adNodes(1).value(end)-adNodes(2).value < 0
               P = 0;
           else
               sigmoid = 1/(1+exp(adNodes(4).value*(adNodes(1).value(end)-adNodes(3).value)));
               Dsigmoid = (-adNodes(4).value*exp(adNodes(4).value*(adNodes(1).value(end)...
                   -adNodes(3).value)))/(1+exp(adNodes(4).value*(adNodes(1).value(end)-adNodes(3).value)))^2;
               P = sigmoid + (adNodes(1).value(end)-adNodes(2).value)*Dsigmoid;
           end
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        %-----------------------------------------------------------------%
        % Payoff for European barrier call-option (only derivative is smoothed)
        function result = adr_barrier_dersmoothed_payoff(X,K,H,d)
            if X.value(end)-K < 0
                P = 0;
            elseif X.value(end) < H 
                P = (X.value(end)-K);
            else 
                P = 0;
            end
            result = ADRev(P);
            result.derivativeOp = @ adr_barrier_dersmoothed_payoffD;
            result.parents = [X K H d];
        end
        
        function modified_parents = adr_barrier_dersmoothed_payoffD(prevDerivative, adNodes)
           if adNodes(1).value(end)-adNodes(2).value < 0
               P = 0;
           elseif adNodes(1).value(end)-adNodes(2).value < adNodes(3).value - adNodes(4).value
               P = 1;
           elseif adNodes(1).value(end)-adNodes(2).value < adNodes(3).value + adNodes(4).value
               P = -4;
           else
               P = 0;
           end
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        
    end
   
end

% --- Help functions -----------------------------------------------------%
function output = myDirac(x)
output = x;
for i=1:length(x)
    if x(i)==0
        output(i) = 100; % 100 or how big? 
    else
        output(i) = 0;
    end
end
end


