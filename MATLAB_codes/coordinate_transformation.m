function [F,Hp,mu, muHp]=coordinate_transformation(IF, M, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the design of new regressor H, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax:
%               [F,Hp,mu, muHp]=coordinate_transformation(IF, M, muF);      Short version
%               [F,Hp,mu, muHp]=coordinate_transformation(IF, M, muF, No, opts);              full version
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       M:                                  main regressor
%                       No:                                 number of iterations
% Optional:
%
%                        opts.c:                            threshold decaying factor (default: c=0.9)
%                        opts.alpha0:                 factor of the initial value for alpha (default: alpha0=500)
%                        opts.wF :                       weighting constant of F-update (default: wF=0.85)
%                        opts.wQ :                      weighting constant of Q-update (default: wQ=0.85)
%                        opts.Ni :                        number of inner-loop iterations (default: Ni=15)
%                        opts.T :                          number of U-update iterations (default: T=3)
%                                                                                                                          
% %%%%%%%%%%%%%%%%%%%% OUTPUTs %%%%%%%%%%%%%%%%%%%%%
%
% The algorithm returns "H", a regressor with low mutual coherence, together with "muHp and mu", a vector
% containing the mutual coherence values for regressor H and a vector including the undesired effect on the noise along the iterations of the algorithm, respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Javad Parsa                                                                                                   %%%
%%%     Electrical Engineering School, 
%%%     Division of Decision and Control Systems                                                               %%%
%%%     KTH Royal Institute of Technology,                                                                       %%%
%%%     Stockholm, Sweden                                                                                                         %%%
%%%     E-mail: javadp@kth.se                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Not enough input parameters!');
end

% parse input parameters
if (isfield(opts,'No'))
    No=opts.No;
else
    No=200;
end

if (isfield(opts,'c_min'))
    c_min=opts.c_min;
else
    c_min=50;
end

if (isfield(opts,'c_max'))
    c_max=opts.c_max;
else
    c_max=500;
end

if (isfield(opts,'b_min'))
    b_min=opts.b_min;
else
    b_min=50;
end

if (isfield(opts,'b_max'))
    b_max=opts.b_max;
else
    b_max=500;
end

if (isfield(opts,'wF'))
    wF=opts.wF;
else
    wF=0.85;
end

if (isfield(opts,'wQ'))
    wQ=opts.wQ;
else
    wQ=0.85;
end

if (isfield(opts,'T'))
    T=opts.T;
else
    T=4;
end
n = size(M,2);
m = size(M,1);
Hp = normc(rand(n,n));    %initial random regressor H
F = M*Hp;  
mu=zeros(1,(No));
Z0=F'*F-IF;
B0 = Hp*Hp'-eye(n);
Hpi = inv(Hp);
for iter=1:No
        mu(iter)= norm((inv(F'*F)-IF),'fro')/(sqrt((size(F,1))*(size(F,2))));
        muHp(iter)=  max(max(abs(Hpi'*Hpi-eye(n))));
        % Q-update
        alph = learnrate(iter, No, c_min, c_max);
        Z0o=Z0;
        Z2=(inv(F'*F)-IF);
        Z0=Z2+wQ*(Z2-Z0o);
        z0=reshape(Z0,[n^2,1]);
        z1=ProjectOntoL1Ball(z0,alph);
        Z1=reshape(z1,[n,n]);
        Q=Z0-Z1+IF;
        beta = learnrate(iter, No, b_min, b_max);
        lam = learnrate(iter, 150, 0.01, 1);
        %lam = 0.5;
        % P-update
        B0o=B0;
        B2=(Hpi'*Hpi-eye(n));
        B0=B2+wQ*(B2-B0o);
        b0=reshape(B0,[n^2,1]);
        b1=ProjectOntoL1Ball(b0,beta);
        B1=reshape(b1,[n,n]);
        P=B0-B1+eye(n);
        M11 = (((P+eye(n))^(-0.5)))*Hp;
        [s11,~,d11]= svd(M11,'econ');
        V = (d11*s11');
        M1 = (((Q+IF)^(-0.5)))*F';
        [s,~,d]= svd(M1,'econ');
        U = d*s';
        for it=1:T
            Y = U*((((Q+IF))^(-0.5)));
            Z = ((((P+eye(n)))^(-0.5)))*V';
            Hp = (inv(M'*M+lam*eye(n)))*(M'*Y+lam*Z);
            F = M*Hp;
            X = F;
            M1 = (((Q+IF)^(-0.5)))*F';
            [s,~,d]= svd(M1,'econ');
            U = d*s';
            M11 = (((P+eye(n))^(-0.5)))*Hp;
            [sp,~,dp]= svd(M11,'econ');
            V = (dp*sp');
            Xo = X;
        end
           F = real(F);
           Hp = real(Hp);
           Hpi = inv(Hp);
           Hpi = normc(Hpi);
           Hpo(:,:,iter) = Hpi;
           Hp = inv(Hpi);
           F = M*Hp;
           Fo(:,:,iter) = F;
end
iii = min(find(muHp==min(muHp)));
F = Fo(:,:,iii);
Hpi = Hpo(:,:,iii);
Hp = inv(Hpi);
end
