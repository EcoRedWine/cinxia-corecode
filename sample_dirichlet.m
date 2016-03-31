function theta = sample_dirichlet(alpha, N)
% SAMPLE_DIRICHLET Sample N vectors from Dir(alpha(1), ..., alpha(k))
% theta = sample_dirichlet(alpha, N)
% theta(i,j) = i'th sample of theta_j, where theta ~ Dir

% % We use the method from p. 482 of "Bayesian Data Analysis", Gelman et al.
% 
% k = length(alpha);
% theta = zeros(N, k);
% scale = 1; % arbitrary?
% for i=1:k
%   theta(:,i) = gamrnd(alpha(i), scale, N, 1);
%   %theta(:,i) = sample_gamma(alpha(i), scale, N, 1);
% end
% S = sum(theta,2);
% theta = theta ./ repmat(S, 1, k); 


%The Dirichlet is a vector of unit-scale gamma random variables, normalized by
%their sum. So, with no error checking, this will get you that:

% a = [2.1 3.2 4.3];
% n = 10000;
% r = drchrnd(a,n)
% 
% function r = drchrnd(a,n)
p = length(alpha);
r = gamrnd(repmat(alpha',N,1),1,N,p); % if the shape parameter (A, the first input) is large, it spits out the same number as the shape param (= independent dispersal)
% if A is > 10^10, it spits out virtually the same number as A. =
% completely random dispersal
theta = r ./ repmat(sum(r,2),1,p);