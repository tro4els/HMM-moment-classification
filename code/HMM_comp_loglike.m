function [d] = HMM_comp_loglike(seq,P123_model)
%
% Composite likelihood of third order moment
%
% inputs:
%   seq: vector of observed symbols
%   P123_model: [K x K x K x maxNN] tensor with minus the logarithm of third order moments at maxNN time steps (the markov chain is assumed stationary after maxNN steps)
%
% ouputs:
%   d: composite likelihood of model (P123_model) given the sequence (seq)
%
% Rasmus Troelsgaard, rast@dtu.dk
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
%

l=length(seq);

e = zeros(l-2,1);

maxNN = size(P123_model,4); % number of burn-in steps

% using sub2ind might be faster than the for loops below

if maxNN==1
    for n=1:l-2
        e(n) = P123_model(seq(n),seq(n+1),seq(n+2));
    end
else
    for n=1:min(l-2,maxNN)
        e(n) = P123_model(seq(n),seq(n+1),seq(n+2),n);
    end
    for n=maxNN+1:l-2
        e(n) = P123_model(seq(n),seq(n+1),seq(n+2),maxNN);
    end
end

d = mean(e);


end