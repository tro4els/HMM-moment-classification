function y = my_hmmgenerate(seqLengths,T,O,PI)

if any(seqLengths<1)
    error('encountered sequence length shorter than 1 symbol. Exiting');
end
    

N = length(seqLengths);
if size(seqLengths,1)<N
    seqLengths = seqLengths';
end
[K,S] = size(O);

if N/mean(seqLengths) > 10
    
    X=zeros(N,max(seqLengths)); % cut off the rest when done
    Y=X;
    
    % draw initial state for all chains:
    idx = randperm(N);
    counts = mnrnd(N,PI);
    tmpidx = 1;
    for s=1:S
        if counts(s)>0
            X(idx(tmpidx:tmpidx+counts(s)-1),1) = s;
            tmpidx = counts(s)+tmpidx;
        end
    end
    
    % draw rest of states
    for i = 2:max(seqLengths)
        for s=1:S
            idxs = find(X(:,i-1)==s);
            if ~isempty(idxs)
                nexts = mnrnd(length(idxs),T(:,s));
                X(idxs(randperm(length(idxs))),i)=repvec(S,nexts);
            end
        end
    end
    
    % generate observations
    for s=1:S
        idxs = find(X(:)==s);
        if ~isempty(idxs)
            obs = mnrnd(length(idxs),O(:,s));
            Y(idxs(randperm(length(idxs)))) = repvec(K,obs);
        end
    end
    
    %cut off unnessecary observations
    
    y = cellfun(@(x,z) x(1:z),mat2cell(Y,ones(1,N)),mat2cell(seqLengths,ones(1,N)),'uniformOutput',false);
    
else % use Matalb's own implementation
    y = cell(N,1);
    for j = 1:N
        y{j} = hmmgenerate(seqLengths(j),[0,PI';zeros(S,1),T'],[zeros(1,K);O']);
    end
end

end

function reps = repvec(S,r)
x=1:S;
t = r > 0;
a = cumsum(r(t));
b = zeros(1,a(end));
b(a - r(t) + 1) = 1;
x1 = x(t);
reps = x1(cumsum(b));
end