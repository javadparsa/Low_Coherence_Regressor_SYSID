function [D , sp_pomp] = RECOSOMP(y,phi,optsr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements RECOSOMP, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax:
%               [D , sp_pomp] = RECOSOMP(y,phi);      Short version
%               [D , sp_pomp] = RECOSOMP(y,phi,optsr);              full version
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       phi:                                  main regressor
%                       y:                                 observation
% Optional:
%
%                        optsr.lambda:                            threshold decaying factor (default: c=0.9)
%                        optsr.lambdad:                 factor of the initial value for alpha (default: alpha0=500)
%                        optsr.lambdaprim :                       weighting constant of F-update (default: wF=0.85)
%                        optsr.num :                      weighting constant of Q-update (default: wQ=0.85)
%                        optsr.inite :                        number of inner-loop iterations (default: Ni=15)                                                                                                                          
% %%%%%%%%%%%%%%%%%%%% OUTPUTs %%%%%%%%%%%%%%%%%%%%%
%
% The algorithm returns "D", a regressor with low mutual coherence, together with "sp_omp", an estimated sparse parameter vector
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
if (isfield(optsr,'lambda'))
    lambda=optsr.lambda;
else
    lambda=0.001;
end

if (isfield(optsr,'lambdad'))
    lambdad=optsr.lambdad;
else
    lambdad=0.006;
end

if (isfield(optsr,'lambdaprim'))
    lambdaprim=optsr.lambdaprim;
else
    lambdaprim=0.0001;
end

if (isfield(optsr,'num'))
    num=optsr.num;
else
    num=40;
end

if (isfield(optsr,'inite'))
    inite=optsr.inite;
else
    inite=40;
end
alph=lambdaprim;
m = size(phi,1);
D=eye(m);
A=phi;
n1 = size(A,2);
for int=1:inite
      for inl=1:num
         A=A-lambda*((A*(A'*A-eye(n1)))+lambdad*(A-D*phi));
         A=normc(A);
      end
       D=(A*phi'*(pinv(phi*phi')));
       for iterd=1:num
        zx=diag(diag(D'*D));
        D=(D-lambda*((D*phi-A)*(phi')+alph*(D*(D'*D-zx))));
        D=normc(D);
       end
end
 Da=D;
 Dnew=normc(D*phi);
 A=Dnew;
 DtD=D*D';
 n = size(D,2);
 D1=eye(n);
for int1=1:inite
  for inl=1:num
         A=A-lambda*((A*(A'*A-eye(n1)))+lambdad*(A-D1*Dnew));
         A=normc(A);
  end
       D1=(A*Dnew'*(pinv(Dnew*Dnew')));
       for iterd1=1:num
        zx=diag(diag((D1*D)'*(D1*D)));
        D1=(D1-lambda*((D1*Dnew-A)*(Dnew')+alph*((D1*DtD*D1'-zx)*D1*DtD)));
        D1=normc(D1);
       end
end
D=D1*D;
Dnew=normc(D*phi);
A=Dnew;
Da=D;
D1=eye(size(D));
DtD=D*D';
for int1=1:inite
  for inl=1:num
         A=A-lambda*((A*(A'*A-eye(n1)))+lambdad*(A-D1*Dnew));
         A=normc(A);
  end
       D1=(A*Dnew'*(pinv(Dnew*Dnew')));
       for iterd1=1:num
        zx=diag(diag((D1*D)'*(D1*D)));
        D1=(D1-lambda*((D1*Dnew-A)*(Dnew')+alph*((D1*DtD*D1'-zx)*D1*DtD)));
        D1=normc(D1);
       end
end
D=D1*D;
ned=D*phi;
ned=normc(ned);
y_new=D*y;
w1=omp(ned,(y_new),[],optsr.ns);
w=w1;
G1=D*phi;
for qw=1:size(G1,2)
    b(qw,1)=norm(G1(:,qw));
end
for rq=1:size(y,2)
 sp_pomp(:,rq)=(w(:,rq))./b;
end
end