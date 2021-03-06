%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a demo for implemntation of CORED-OMP, CORED-ADMM, OMP, ADMM,
% RECOSOMP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       phi:                                  main regressor
%                       y:                                    observation                                                                                          
% %%%%%%%%%%%%%%%%%%%% OUTPUTs %%%%%%%%%%%%%%%%%%%%%
%
% The demo returns POD and NRMSE of each algorithm.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Javad Parsa                                                                                                   %%%
%%%     Electrical Engineering School, 
%%%     Division of Decision and Control Systems                                                               %%%
%%%     KTH Royal Institute of Technology,                                                                       %%%
%%%     Stockholm, Sweden                                                                                                         %%%
%%%     E-mail: javadp@kth.se                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
numou = 6;   % number of experiments
snrm = [15 20 25 30 40 50 60 80 100];  % different levels of SNR (noise)

%% Initial zeros vectors to report the POD, NRMSE and computational time.
error1=zeros(size(snrm));
error2=zeros(size(snrm));
error3=zeros(size(snrm));
error4=zeros(size(snrm));
error5=zeros(size(snrm));
error6=zeros(size(snrm));
error7=zeros(size(snrm));
error8=zeros(size(snrm));
errormLS=zeros(size(snrm));
errorHF=zeros(size(snrm));
errorCOREDOMP=zeros(size(snrm));
err1=zeros(size(snrm));
err2=zeros(size(snrm));
err3=zeros(size(snrm));
errCOREDADMM=zeros(size(snrm));
err5=zeros(size(snrm));
err6=zeros(size(snrm));
err7=zeros(size(snrm));
err8=zeros(size(snrm));
errmLS=zeros(size(snrm));
err10=zeros(size(snrm));
err11=zeros(size(snrm));
errHF=zeros(size(snrm));
errLS=zeros(size(snrm));
ersp1=zeros(size(snrm));
ersp2=zeros(size(snrm));
ersp3=zeros(size(snrm));
erspCOREDADMM=zeros(size(snrm));
ersp5=zeros(size(snrm));
ersp6=zeros(size(snrm));
ersp7=zeros(size(snrm));
ersp8=zeros(size(snrm));
erspmLS=zeros(size(snrm));
ersp10=zeros(size(snrm));
ersp11=zeros(size(snrm));
erspHF=zeros(size(snrm));
erspCOREDOMP=zeros(size(snrm));
erss1=zeros(size(snrm));
erss2=zeros(size(snrm));
erss3=zeros(size(snrm));
erss4=zeros(size(snrm));
erss5=zeros(size(snrm));
erss6=zeros(size(snrm));
erss7=zeros(size(snrm));
erss8=zeros(size(snrm));
erssmLS=zeros(size(snrm));
erss10=zeros(size(snrm));
erss11=zeros(size(snrm));
erssHF=zeros(size(snrm));
erssLS=zeros(size(snrm));
sreco1=zeros(size(snrm));
sreco2=zeros(size(snrm));
sreco3=zeros(size(snrm));
srecoCOREDADMM=zeros(size(snrm));
sreco5=zeros(size(snrm));
sreco6=zeros(size(snrm));
sreco7=zeros(size(snrm));
sreco8=zeros(size(snrm));
srecomLS=zeros(size(snrm));
sreco10=zeros(size(snrm));
sreco11=zeros(size(snrm));
srecoHF=zeros(size(snrm));
srecoCOREDOMP=zeros(size(snrm));
srec1=zeros(size(snrm));
srec2=zeros(size(snrm));
srec3=zeros(size(snrm));
srecCOREDADMM=zeros(size(snrm));
srec5=zeros(size(snrm));
srec6=zeros(size(snrm));
srec7=zeros(size(snrm));
srec8=zeros(size(snrm));
srecmLS=zeros(size(snrm));
srecHF=zeros(size(snrm));
srecLS=zeros(size(snrm));
tim1=zeros(size(snrm));
timCOREDADMM=zeros(size(snrm));
timtcdc=zeros(size(snrm));
timadmm=zeros(size(snrm));
timHF=zeros(size(snrm));
timLS=zeros(size(snrm));
timf1=zeros(size(snrm));
timfCOREDADMM=zeros(size(snrm));
timftcdc=zeros(size(snrm));
timfadmm=zeros(size(snrm));
timfLS=zeros(size(snrm));
%% Hyperparameters for coordinate transformation
m = 300;     % number of observations
n = 50;      % number of parameters
ns = 7;      % number of non-zero elements in parameter vector
IF = eye(n,n); % identity matrix                     
opts.c=0.95;               % threshold decaying factor (default: c=0.9)
opts.wF=0.85;            % weighting constant of F-update (default: wF=0.85)
opts.wQ=0.6;           % weighting constant of Q-update (default: wQ=0.85)
opts.T=4;                     % number of H-update iterations (default: T=3)
opts.c_min = 30;
opts.c_max = 300;
opts.b_min = 30;
opts.b_max =300;
opts.No = 300;    % number of iterations for coordinate_transformation algorithm
%% hyperparmeters for RECOSOMP
optsr.lambda = 0.0008;
optsr.lambdad = 0.004;
optsr.lambdaprim = 0.0001;
optsr.inite = 25;
optsr.num = 25;
optsr.ns = ns;
for out = 1:numou
[phi,s_true] = generate_highcorr_reg(m,n, ns); % generating sparse random parameter vector and random regressor with high mutual coherence
%[phi,s_true] = generate_highcorr_regrand(m,n, ns);
aw1=phi'*phi;
M = phi;
for ii =1:size(snrm,2)
SNR =snrm(ii);
[Y,r] = generate_observation(phi,s_true,m,SNR); % generating the observation
d1_Y = size(Y,1); 
d2_Y = size(Y,2);
tic;s_OMP = omp(M,Y,[],ns);t1=toc; % OMP sparse estimation method
diff_OMP = norm(s_true-s_OMP)/norm(s_OMP); % NRMSE for OMP
tic;[F,Hp,mu, muHp] = coordinate_transformation(IF, M, opts);tcoortran = toc; % coordinate transformation to design the new regressor H
Hpi = inv(Hp); % H^{-1}
F = M*Hp; 
tic;[D , s_RECOSOMP] = RECOSOMP(Y,phi,optsr);tcdc = toc; % RECOSOMP algorithm
x = (inv(F'*F))*F'*Y; % Lease square method to estimate parameter vector x
tic;options.verbose  = 1;options.relative_tolerance = 1e-2;[wbp,~]= lassoadmm(Hpi, x, 0.14, 1.5, 1.1); sp_COREDADMM=(wbp);t_COREDADMM=toc;t_COREDADMM=tcoortran+t_COREDADMM; % CORED_ADMM algorithm
tic;[weadmm1,~]=lassoadmm(phi,Y, 0.14, 1, 1);weadmm=weadmm1;tadmm=toc; % Lasso via ADMM method
tic;sp_COREDOMP = omp(Hpi,x,[],ns);t_COREDOMP= toc;t_COREDOMP=t_COREDOMP+tcoortran; % CORED_OMP algorithm
tic;aaa = find(sp_COREDOMP~=0);
PP = M(:,aaa);
xx = (inv(PP'*PP))*PP'*Y;
s_COREDOMP = zeros(n,1);
s_LS = phi\Y;
s_COREDOMP(aaa) = xx;t_LS = toc;t_LS = t_LS+tcoortran;
errsHF = (norm(s_true-sp_COREDOMP))/(norm(s_true));
diff_COREDOMP = norm(s_true-s_COREDOMP)/norm(s_true);
diff_RECOSOMP = norm(s_true-s_RECOSOMP)/norm(s_true);
y=Y;
%% error of each algorithm which is normalized by the dimensions of observation
error1(ii)=(norm((Y-phi*(s_OMP)),'fro')/(sqrt(d1_Y*d2_Y)))+error1(ii);
error3(ii)=(norm((Y-phi*(s_RECOSOMP)),'fro')/(sqrt(d1_Y*d2_Y)))+error3(ii);
error2(ii)=(norm((Y-phi*(weadmm))/(sqrt(d1_Y*d2_Y)),'fro'))+error2(ii);
error4(ii)=(norm((y-phi*(sp_COREDADMM)),'fro')/(sqrt(d1_Y*d2_Y)))+error4(ii);
errormLS(ii)=(norm((y-phi*(s_LS)),'fro')/(sqrt(d1_Y*d2_Y)))+errormLS(ii);
errorCOREDOMP(ii)=(norm((y-phi*(s_COREDOMP)),'fro')/(sqrt(d1_Y*d2_Y)))+errorCOREDOMP(ii);
%% NRMSE of each algorithm
erspCOREDADMM(ii)=norm((s_true-(sp_COREDADMM)),'fro')/norm(s_true,'fro')+erspCOREDADMM(ii);
ersp1(ii)=norm((s_true-s_OMP),'fro')/norm(s_true,'fro')+ersp1(ii);
ersp3(ii)=norm((s_true-(s_RECOSOMP)),'fro')/norm(s_true,'fro')+ersp3(ii);
ersp2(ii)=norm((s_true-(weadmm)),'fro')/norm(s_true,'fro')+ersp2(ii);
erspmLS(ii)=norm((s_true-(s_LS)),'fro')/norm(s_true,'fro')+erspmLS(ii);
erspCOREDOMP(ii)=norm((s_true-(s_COREDOMP)),'fro')/norm(s_true,'fro')+erspCOREDOMP(ii);
%% POD of each algorithm
sreco1(ii)=(nnz((s_true.*s_OMP)))/ns+sreco1(ii);
sreco3(ii)=(nnz((s_true.*(s_RECOSOMP))))/ns+sreco3(ii);
sreco2(ii)=(nnz((s_true.*(weadmm))))/ns+sreco2(ii);
srecoCOREDADMM(ii)=(nnz((s_true.*(sp_COREDADMM))))/ns+srecoCOREDADMM(ii);
srecomLS(ii)=(nnz((s_true.*(s_LS))))/ns+srecomLS(ii);
srecoHF(ii)=(nnz((s_true.*(sp_COREDOMP))))/ns+srecoHF(ii);
srecoCOREDOMP(ii)=(nnz((s_true.*(s_COREDOMP))))/ns+srecoCOREDOMP(ii);
%% Computational time
tim1(ii)=t1+tim1(ii);
timCOREDADMM(ii)=t_COREDADMM+timCOREDADMM(ii);
timtcdc(ii)=tcdc+timtcdc(ii);
timadmm(ii)=tadmm+timadmm(ii);
timLS(ii)=t_LS+timLS(ii);
end
out
end
for isu=1:numou
    err1=error1/numou;
    err2=error2/numou;
    err3=error3/numou;
    errCOREDADMM=error4/numou;
    err5=error5/numou;
    err6=error6/numou;
    err7=error7/numou;
    err8=error8/numou;
    errmLS=errormLS/numou;
    errHF=errorHF/numou;
    errLS=errorCOREDOMP/numou;
    erss1=ersp1/numou;
    erss2=ersp2/numou;
    erss3=ersp3/numou;
    erssCOREDADMM=erspCOREDADMM/numou;
    erss5=ersp5/numou;
    erss6=ersp6/numou;
    erss7=ersp7/numou;
    erss8=ersp8/numou;
    erssmLS=erspmLS/numou;
    erssHF=erspHF/numou;
    erssLS=erspCOREDOMP/numou;
    srec1=sreco1/numou;
    srec2=sreco2/numou;
    srec3=sreco3/numou;
    srecCOREDADMM=srecoCOREDADMM/numou;
    srec5=sreco5/numou;
    srec6=sreco6/numou;
    srec7=sreco7/numou;
    srec8=sreco8/numou;
    srecmLS=srecomLS/numou;
    srecHF=srecoHF/numou;
    srecLS=srecoCOREDOMP/numou;
    timf1=tim1/numou;
    timfCOREDADMM=timCOREDADMM/numou;
    timftcdc=timtcdc/numou;
    timfadmm=timadmm/numou;
    timfHF=timHF/numou;
    timfLS=timLS/numou;
end
%% plot figures
figure(1);
plot(snrm,errLS,'-o',snrm,err1,'-s',snrm,err3,'-p',snrm,err2,'-d',snrm,errCOREDADMM,'-^','linewidth',2.5,'MarkerSize',7);
leg=legend('CORED-OMP','OMP','RECOSOMP','ADMM','CORED-ADMM');set(leg,'fontsize',15);
set(gca,'fontsize',15);
xlabel('SNR','FontName','Times New Roman','fontsize',15);
ylabel('$\frac{\|{y-\Phi\theta}\|_2}{\sqrt{m}}$','interpreter','latex','fontsize',15);
xlim([min(snrm) max(snrm)])
figure(2);
plot(snrm,erssLS,'-o',snrm,erss1,'-s',snrm,erss3,'-p',snrm,erss2,'-d',snrm,erssCOREDADMM,'-^','linewidth',2.5,'MarkerSize',7);
leg=legend('CORED-OMP','OMP','RECOSOMP','ADMM','CORED-ADMM');set(leg,'fontsize',15);
set(gca,'fontsize',15);
xlabel('SNR','FontName','Times New Roman','fontsize',15);
ylabel('$\frac{\|{\hat{\theta}-\theta_0}\|_2}{\|{\theta_0}\|_2}$','interpreter','latex','fontsize',15);
xlim([min(snrm) max(snrm)])
figure(3);
plot(snrm,srecLS,'-o',snrm,srec1,'-s',snrm,srec3,'-p',snrm,srec2,'-d',snrm,srecCOREDADMM,'-^','linewidth',2.5,'MarkerSize',7);
leg=legend('CORED-OMP','OMP','RECOSOMP','ADMM','CORED-ADMM');set(leg,'fontsize',15);
set(gca,'fontsize',15);
xlabel('SNR','FontName','Times New Roman','fontsize',15);
ylabel('POD','interpreter','latex','fontsize',15);
xlim([min(snrm) max(snrm)])