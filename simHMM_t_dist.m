function [z,y]=simHMM_t_dist(pi,A,mu,sigmasq,nu,T)
%  Simulate an univariate gaussian HMM observation series y and its corresponding state
%  series z of length T.
%  Input:
%    Initial probability---pi, a K by 1 probability vector
%    Transition matrux ---A, a K by K matrix with column sum up to 1
%    Emission mean---mu, the emission mean, i.e., E(y_t|x_t=k)=mu_k
%    Emission variance, a K by 1 vector---sigma, the emission variance, i.e., Var(y_t|x_t=k)=sigmasq_k
%    Length-T, the length of the series
%  Output:
%    z--the state series
%    y--the observation series
z=zeros(T,1);
y=zeros(T,1);
K=length(pi);
%Initialize z[1]
z(1)=randsample(K,1,true,pi);

for i=2:T
    %Simulate the state series z
    z(i)=randsample(K,1,true,A(:,z(i-1)));
end

for i=1:T
    %Simulate the observation series y from z
    y(i) = mu(z(i)) + sqrt(sigmasq(z(i))) * trnd(nu);
end