function [tbound,pi0] = markovConvergenceTimeBound(T,PI,epsilon)
%
% Calculates upper bound on convergence time given an upper limit on the
% total variation distance from the stationary distribution.
%
% Inputs:
%   T: Transition probability matrix. Size: S x S
%   PI: Initial state probability vector. Size: 1 x S
%
% Outputs:
%   tbound: time step at which convergence criterion is met
%   pi0: stationary distribution of the Markov chain
%
%
% Rasmus Troelsgaard, rast@dtu.dk
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
%

S=size(T,1);
if S==1 % special case: no dynamics, just mixture model
    tbound=0;
    pi0 = 1;
else
    
    % find stationary distribution
    pi0 = (ones(1,S)/(eye(S)-T+ones(S))')';
    
    Ttilde = (T.*(ones(S,1)*pi0')./(pi0*ones(1,S)))';
    
    MP=T*Ttilde;
    bet=sort(eigs(MP,2),1,'descend');
    
    chi0 = sqrt(sum(((PI-pi0).^2)./pi0));
    epsilon = min(min(1,epsilon),chi0/2);
    tbound = 2*(log(2*epsilon)-log(chi0))/log(bet(2));
end
end