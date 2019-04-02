%
% c_0_225(P,t) returns c_0 as a function of time (non-dimensional)
%

function c_0 = c_0_225(P,t)
    b1 = 3.939-1;
    %b2 = -0.2262; %median of 95% CI
    b2 = -0.3915; % lower bound 95% CI
    %b3 = 31.94;
    b3 = 38.22;
    
    c_0 = b3*ones(size(t))./(1+exp(b1+b2*(t*P.t_0/60)));
    c_0 = c_0/100;
end