function mu = damageFactor(x)
    global S0 S1
    if x < S0(2) && x > -1
        mu = 1;
    elseif x > S0(2) && x < S1(2)
        mu = (S1(2) - x)/(S1(2) - S0(2));
    else
        mu = 0;
    end
end