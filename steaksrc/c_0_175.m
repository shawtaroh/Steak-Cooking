%
% c_0_175(P,t) returns c_0 as a function of time (non-dimensional)
%

function c_0 = c_0_175(P,t)
    b1 = 2.928;
    b2 = -0.08363;
    b3 = 26.55;
    c_0 = b3/(1+exp(b1+b2*(t*P.t_0/60)));
    c_0 = c_0/100;
end