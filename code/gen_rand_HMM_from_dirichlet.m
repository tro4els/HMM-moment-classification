function [T,O,PI] =  gen_rand_HMM_from_dirichlet(K,S,Tdiag,vinter,Osparsity,PI_base)
%
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