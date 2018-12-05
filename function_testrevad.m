function output = function_testrevad(x,y)
    %output = adr_mul(adr_add(x,y),adr_add(y,ADRev(1))); % times and plus
    %output = adr_div(x,y); % division
    %output = adr_sub(x,y); % sub
    %output = adr_mul(x,y); % mult
    %output = adr_mul(adr_exp(x),y); % mult. and exp.
    %output = adr_ln(adr_mul(x,y));
    %output = adr_pow(x,2);
    %output = adr_sin(x);
    %output = adr_step(x);
    
    %output = adr_mul(adr_add(adr_add(adr_pow(x,2),adr_mul(x,2)),1),y);
    %output = adr_mul(adr_add(x,3),adr_add(y,2));
    output = x*y;
end
