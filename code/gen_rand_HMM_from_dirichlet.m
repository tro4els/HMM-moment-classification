function [T,O,PI] =  gen_rand_HMM_from_dirichlet(K,S,Tdiag,vinter,Osparsity,PI_base)
%
% Generate parameters of an HMM with discrete observations
%
% inputs:
%   K: size of observation space (number of symbols)
%   S: size of state space
%   Tdiag: mean value of diagonal elements
%   vinter: ]0;1] parameter controlling the variance of the dirichlet distributions to
%           sample from. It specifies the interpolation between no variance and
%           maximum variance given Tdiag
%
% Outputs:
%   T: Transition probability matrix. Size: S x S
%   O: Observation probability matrix. Size: K x S
%   PI: Initial state probability vector. Size: 1 x S
%
%
% Rasmus Troelsgaard, rast@dtu.dk
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
%
if nargin<6
   PI_base = ones(1,S)/S; 
end

if nargin<5 || isempty(Osparsity);
    Osparsity=1;
end
if nargin<4
    vinter=1;
end


xmin=S;% this corresponds to minimal scaling parameter equal to S
lam=(xmin+1-vinter)/(xmin*vinter);
v = lam*xmin*[Tdiag,((1-Tdiag)/(S-1))*ones(1,S-1)];

T = drchrnd(v,S)';
for i=2:S
    tmp=T(i,i);
    T(i,i)=T(1,i);
    T(1,i)=tmp;
end


PI = drchrnd(lam*xmin*PI_base,1)';


xmin=1;% this corresponds to minimal scaling parameter equal to 1
lam=(xmin+1-vinter)/(xmin*vinter);

O=drchrnd(lam*xmin*Osparsity*ones(1,K),S)';

end