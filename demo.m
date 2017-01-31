% This is a short demo script for the composite likelihood based
% classification of discrete time series based on third order moments
addpath('./code/');


S=6; % state space size
K=15; % number of unique symbols
Nclass = 4; % number of classes

% upper bound on Total Variation Distance to stationary distribution
epsilon = 1e-3;

% number of sequences from each class
N=50;

% relatedness of classes: ]0;1]
c=0.3;

% approximate size of diagonal elements of transition probability matrices: [0;1[
d=0.5;


T=cell(Nclass,1);
O=cell(Nclass,1);
PI=cell(Nclass,1);
P123=cell(Nclass,1);

sequences = cell(N*Nclass,1);
labels = zeros(N*Nclass,1);

% generate classes and sequences
for j=1:Nclass

    [T_,O_,PI_] = gen_rand_HMM_from_dirichlet(K,S,d,c);

    tbound = markovConvergenceTimeBound(T_,PI_,epsilon);
    
    P123_ = calc_model_moments(O_, T_, PI_, ceil(tbound));
    P123{j} = P123_;
    
    
    T{j} = T_;
    O{j} = O_;
    PI{j} = PI_;
    
    NN = poissrnd(ones(N,1)*300);
    
    sequences((j-1)*N+1:j*N) = my_hmmgenerate(NN,T_,O_,PI_);
    labels((j-1)*N+1:j*N) = ones(N,1)*j;
    
end


% classify according to composite log-likelihood
comp_loglike = zeros(Nclass*N,Nclass);
loglike = zeros(Nclass*N,Nclass);
for j=1:Nclass
    for i=1:Nclass*N;
        comp_loglike(i,j) = HMM_comp_loglike(sequences{i},P123{j});
        loglike(i,j) = -hmmlogp(sequences{i},T{j},O{j},PI{j});
    end
end

[~,guess] = randmin(comp_loglike,[],2);
acc_comp = sum(labels==guess)/length(labels);

[~,guess] = randmin(loglike,[],2);
acc_loglike= sum(labels==guess)/length(labels);


fprintf(1,'Classification accuracy:\n');
fprintf(1,'\tComposite log-likelihood method: %.4f\n',acc_comp);
fprintf(1,'\tTrue log-likelihood method: %.4f\n',acc_loglike);

