function [Y,r] = generate_observation(phi,s,m,SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates observations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax:
%               [Y,r] = generate_observation(phi,s,m,SNR);             full version
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       m:                              number of observation
%                       phi:                            main regressor
%                       s:                              sparse parameter vector
%                       SNR:                            level of noise
% %%%%%%%%%%%%%%%%%%%% OUTPUTs %%%%%%%%%%%%%%%%%%%%%
%
% The algorithm returns "Y" observations, together with "r", level of noise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Javad Parsa                                                                                                   %%%
%%%     Electrical Engineering School, 
%%%     Division of Decision and Control Systems                                                               %%%
%%%     KTH Royal Institute of Technology,                                                                       %%%
%%%     Stockholm, Sweden                                                                                                         %%%
%%%     E-mail: javadp@kth.se                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y1 = phi*s;
YN = randn(m,1);
le = (10^((20*log10(norm(Y1))-SNR)/20))/norm(YN);
Y = Y1+le*YN;
r = snr(Y1,Y-Y1);
end