function [P123_model] = calc_model_moments(O, T, PI, maxNN)
% 
% O: [K x S] observation probability matrix
% T: [S x S] transition probability matrix
% PI: [S x 1] initial state probability vector
% maxNN >= 1 (integer) position (inclusive) at which the Markov chain is assumed stationary

[K,S] = size(O);

% find stationary distribution
pi0 = (ones(1,S)/(eye(S)-T+ones(S))')';

tt = [PI,zeros(S,maxNN-1)];
for i=2:maxNN
    tt(:,i)=T*tt(:,i-1);
end

TOdiags=T'*repmat(eye(S),[1,K]).*(ones(S,1)*reshape(O',[1,S*K])); 

P123_model = zeros(K,K,K,maxNN);
for n=1:maxNN-1
    P123_model(:,:,:,n) = reshape(reshape(permute(reshape(O*diag(tt(:,n))*TOdiags,[K,S,K]),[1,3,2]),[K*K,S])*T'*O',[K,K,K]);
end

% place stationary third order moment at the end
P123_model(:,:,:,maxNN) = reshape(reshape(permute(reshape(O*diag(pi0)*TOdiags,[K,S,K]),[1,3,2]),[K*K,S])*T'*O',[K,K,K]);

end
