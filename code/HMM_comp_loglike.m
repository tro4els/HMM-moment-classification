function [d] = HMM_comp_loglike(seq,P123_model)
%
% composite likelihood of third order moment
%
% inputs:
%   seq: vector of observed symbols
%   P123_model: [K x K x K x maxNN] tensor (markov chain is assumed stationary after maxNN steps)
%
% ouputs:
%   d: composite likelihood of model (P123) given the sequence (seq)
%
%

l=length(seq);

e = zeros(l-2,1);

maxNN = size(P123_model,4); % number of burn-in steps

% using sub2ind might be faster than the for loops below

if maxNN==1
    for n=1:l-2
        e(n) = -log(P123_model(seq(n),seq(n+1),seq(n+2)));
    end
else
    for n=1:min(l-2,maxNN)
        e(n) = -log(P123_model(seq(n),seq(n+1),seq(n+2),n));
    end
    for n=maxNN+1:l-2
        e(n) = -log(P123_model(seq(n),seq(n+1),seq(n+2),maxNN));
    end
end

d = mean(e);


end