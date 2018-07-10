%
%   mu.m - Returns mu(T)
%
%
function y = mu(T)

y = 2.414e-5 * 10.^(247.8*ones(size(T))./(T-140*ones(size(T))));
%y = 2.414e-5 * 10^(247.8/(T-140));

end