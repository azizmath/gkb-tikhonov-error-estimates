% clc;clf; 
% clear;
% delta = 1e-4;
% lambda = 1e-2;
% noise_operator = 1;
% %[ke,xexact,h,y,y_delta] = Baart_GKB(delta,lambda,noise_operator);
% %[ke,xexact,h,y,y_delta] = Foxgood_GKB(delta,lambda,noise_operator);
% %[ke,xexact,h,y,y_delta] = Gravity_GKB(delta,lambda,noise_operator);
% [ke,xexact,h,y,y_delta] = Shaw_GKB(delta,lambda,noise_operator);
% 
% ell = 30; iv = 2;
% tau = 1e-12;
% [C,U,V,gamma_ell] = continuous_GKB(ke,y_delta,xexact,iv,tau,ell);
%  T = C(1:ell+1,1:ell)*C(1:ell,1:ell)';
%  dT = diag(T,0);
%  sT = diag(T,-1);
%  s = svd(ke);
%  ss = s.^2;
%  rhs = zeros(ell,1);
%  lhs = zeros(ell,1);
%  for i= 1: ell
%      rhs(i,1) = prod(ss(1:i));
%      lhs(i,1) = prod(sT(1:i));
%  end
%  figure (1);semilogy(lhs,'r-');hold on;semilogy(rhs,'b-');
%  hold off;legend('lhs','rhs');
%  figure(2);semilogy(dT,'r-');hold on;semilogy(sT,'b-');
%  hold off;legend('diag','off-diag');
 clc
for i = 1:10 
    fprintf(['(%e,%e)%%\n'],i,rhs(i,1));
 end
% %  
% %  
% %   