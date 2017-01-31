function [tbound,pi0] = markovConvergenceTimeBound(T,PI,epsilon)

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