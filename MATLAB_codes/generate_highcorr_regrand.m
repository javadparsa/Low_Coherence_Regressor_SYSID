function [phi,s] = generate_highcorr_regrand(m,n, ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates random matrix with high mutual coherence and random sparse vector. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax:
%               [phi,s] = generate_highcorr_regrand(m,n, ns);             full version
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       m:                              number of observation
%                       n:                              number of parameters
%                       ns:                             number of non-zero parameters
% %%%%%%%%%%%%%%%%%%%% OUTPUTs %%%%%%%%%%%%%%%%%%%%%
%
% The algorithm returns "phi", a regressor with high mutual coherence, together with "s", a sparse random vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Javad Parsa                                                                                                   %%%
%%%     Electrical Engineering School, 
%%%     Division of Decision and Control Systems                                                               %%%
%%%     KTH Royal Institute of Technology,                                                                       %%%
%%%     Stockholm, Sweden                                                                                                         %%%
%%%     E-mail: javadp@kth.se                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dif = fix(n/(ns-1));
as = 1:fix(dif*2):n;  
aa = 0.02;  % parameter to set the mutual coherence
%% regressor with high mutual coherence
[C,~] = randcorr(n,aa,as);
L = chol(C);
X = null([randn(m-n-1,m);ones(1,m)])*L;
phi = normc(X);
%% Sparse random vector
s = zeros(n,1);
for jj=1:dif:n
    s(jj) = 2*randn(1);
end