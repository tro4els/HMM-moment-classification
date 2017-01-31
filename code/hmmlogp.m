function loglike = hmmlogp(y,A,B,pi1)
%  hmmlogp --> Log-likelihood for given observation sequence and HMM.
%
%  <Synopsis>
%    loglike = hmmlogp(y,A,B,pi1)
%
%  <Description>
%    The function calculates the log-likelihood for a given
%    observation sequence y, and Hidden Markov Model A, B, and  pi1.
%    Here, A is the S-by-S state transition matrix, B is the K-by-S
%    observation probability matrix, and pi1 is the initial state
%    probability vector.
%
%  <See Also>
%    hmmrecog --> HMM based word classifier.

%  <References>
%  [1] J.R Deller, J.G. Proakis and F.H.L. Hansen, "Discrete-Time
%      Processing of Speech Signals", IEEE Press, chapter 12, (2000).
%
%  <Revision>
%    Peter S.K. Hansen, IMM, Technical University of Denmark
%
%    Last revised: September 30, 2000
%-----------------------------------------------------------------------

T = length(y);               % Length of observation sequence.
S = length(A);               % Number of hidden states.

alpha = zeros(S,T);                  % Forward recursion.
scale = zeros(T,1);                  % Scale factors.
alpha(:,1) = pi1.*B(y(1),:)' + eps;  % Alpha for t=1.
scale(1)   = sum(alpha(:,1));        % Scale factor for t=1.
alpha(:,1) = alpha(:,1)/scale(1);    % Scaled alpha for t=1.
for (t = 2:T)
  alpha(:,t) = A*alpha(:,t-1).*B(y(t),:)';
  scale(t)   = sum(alpha(:,t));      % Scale factor for t.
  alpha(:,t) = alpha(:,t)/(scale(t) + eps); % Scaled alpha for t.
end

loglike = sum(log10(scale + eps));

%-----------------------------------------------------------------------
% End of function hmmlogp
%-----------------------------------------------------------------------
