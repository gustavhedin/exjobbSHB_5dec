%-----------------------------------------------------------------%        
        % normal-smoothed indicator function with variance sigma:
        function result = adr_smoothed(X,K,sigma)
            P = 0.5*(1+erf((X.value-K)/(sigma*sqrt(2))));
            result = ADRev(P);
            result.derivativeOp = @ adr_smoothedD;
            result.parents = [X K sigma];
        end
        
        function modified_parents = adr_smoothedD(prevDerivative, adNodes)
            P = 1/(2*pi*adNodes(3).value)*exp(-1*(adNodes(1).value-adNodes(2).value).^2/(2*adNodes(3).value));
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        % normal-smoothed reverse indicator function with variance sigma:
        function result = adr_neg_smoothed(X,K,sigma)
            P = 0.5*erfc((X.value-K)/(sigma*sqrt(2)));
            result = ADRev(P);
            result.derivativeOp = @ adr_neg_smoothedD;
            result.parents = [X K sigma];
        end
        
        function modified_parents = adr_neg_smoothedD(prevDerivative, adNodes)
            P = -1/(2*pi*adNodes(3).value)*exp(-1*(adNodes(1).value-adNodes(2).value).^2/(2*adNodes(3).value));
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        % Smoothed indicator function, w/ smoothing parameter delta
        function result = adr_pos_indicator(X,K,delta)
            P = zeros(size(X.value,1),size(X.value,2));
            if size(P,1) > size(P,2)
                d = size(P,1);
            else
                d = size(P,2);
            end
            for i=1:d
                if ( X.value(i) < K - delta) 
                    P(i) = 0;
                elseif (X.value(i) > K + delta)
                    P(i) = 1;
                else
                    temp = min((X.value(i)-(K-delta))/(2*delta),1);
                    P(i) = max(0,temp);
                end
            end
            result = ADRev(P);
            result.derivativeOp = @ adr_pos_indicatorD;
            result.parents = [X K delta];
        end
        
        function modified_parents = adr_pos_indicatorD(prevDerivative, adNodes)
            P = zeros(size(adNodes(1).value,1),size(adNodes(1).value,2));
            % fin out if rowvector or colonvector
                if size(P,1) > size(P,2)
                    d = size(P,1);
                else
                    d = size(P,2);
                end
            %--------------------------------------%
            for i=1:d
                if (adNodes(2).value-adNodes(3).value < adNodes(1).value(i)) && (adNodes(2).value+adNodes(3).value > adNodes(1).value(i))               
                    P(i) = 1/(2*adNodes(3).value);
                else
                    P(i) = 0;
                end
            end
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end
        
        % Smoothed reverse indicator function, w/ smoothing parameter delta
        function result = adr_neg_indicator(X,K,delta)
            P = ones(size(X.value,1),size(X.value,2));
            if size(P,1) > size(P,2)
                d = size(P,1);
            else
                d = size(P,2);
            end
            for i=d:d
                if ( X.value(i) < K - delta) 
                    P(i) = 1;
                elseif (X.value(i) > K + delta)
                    P(i) = 0;
                else
                    temp = 1-(X.value(i)-(K-delta))/(2*delta);
                    P(i) = max(0,temp);
                end
            end
            result = ADRev(P);
            result.derivativeOp = @ adr_neg_indicatorD;
            result.parents = [X K delta];
        end
        
        function modified_parents = adr_neg_indicatorD(prevDerivative, adNodes)
            P = zeros(size(adNodes(1).value,1),size(adNodes(1).value,2));
            % find out if rowvector or colonvector
                if size(P,1) > size(P,2)
                    d = size(P,1);
                else
                    d = size(P,2);
                end
            %--------------------------------------%
            for i=1:d
                if (adNodes(2).value-adNodes(3).value < adNodes(1).value(i)) && (adNodes(2).value+adNodes(3).value > adNodes(1).value(i))               
                    P(i) = -1/(2*adNodes(3).value);
                else
                    P(i) = 0;
                end
            end
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*P;
            modified_parents = adNodes;
        end