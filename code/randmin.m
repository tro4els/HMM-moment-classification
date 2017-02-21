function [M,i] = randmin(x,y,dim)
% Finds and returns the minimum of a matrix and its index.
% Differs from standard min() in that it chooses randomly among multiple
% identical values. Thus the indices i can differ from run to run.
%
% Rasmus Troelsgaard, rast@dtu.dk
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
%

if nargin < 2
    dim=1;
elseif ~isempty(y)
    error('data in second argument is not supported. Use []')
end

d = size(x);

restdim = setdiff(1:length(d),dim);
if dim ~= 1
    x = permute(x,[dim,restdim]);
end

% make 2D
if length(d)>2
    x=reshape(x,[d(dim) prod(d(restdim))]);
end

[M,I] = min(x);

mask = ones(d(dim),1)*M == x;

idxMod = find(sum(mask)>1);

for ii = idxMod
    idx = find(mask(:,ii));
    I(ii) = idx(randi(length(idx)));
end

%reshape to fit original size
if length(d)>2
    I=reshape(I,[1,d(restdim)]);
    M=reshape(M,[1,d(restdim)]);
end

if dim ~= 1
    i=ipermute(I,[dim,restdim]);
    M=ipermute(M,[dim,restdim]);
else
    i=I;
end
