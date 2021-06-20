function [phi,s] = generate_highcorr_reg(m,n, ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates random matrix with high mutual coherence and random sparse vector. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax:
%               [phi,s] = generate_highcorr_reg(m,n, ns);             full version
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
dis = fix(n/(ns-1));
dif = fix(dis*2);
phi=randn(m,n);
newmat=rand(m,n);
newmat=normc(newmat);
phi=normc(phi);
aw=newmat*newmat';
awn=normc(aw);
s = zeros(n,1);
for j=1:dif:n
    phi(:,j)=awn(:,((j-1)/dif)+1);
end
for jj=1:dis:n
    s(jj) = 2*randn(1);
end