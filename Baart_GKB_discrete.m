function [A,xexact,h,y,y_delta]= Baart_GKB_discrete(n,delta,noise_operator)
disp('Baart')
%BAART Test problem: Fredholm integral equation of the first kind. 
% 
% [A,b,x] = baart_alt(n) 
% 
% Discretization by a Nystrom method based on the trapezoidal rule
% of a first-kind Fredholm integral equation with kernel K and 
% right-hand side g given by 
%
%    K(s,t) = exp(s*cos(t)),  g(s) = 2*sinh(s)/s, 
%
% and with integration intervals  s in [0,pi/2],  t in [0,pi]. 
% The solution is given by 
%
%    f(t) = sin(t). 
%
% The n by n matrix A is a discretization of the integral operator,
% the n vector b a discretization of g, and x a discretization of
% f. Note that Ax differs from b. 
  
% Reference: M. L. Baart, "The use of auto-correlation for pseudo- 
% rank determination in noisy ill-conditioned linear least-squares 
% problems", IMA J. Numer. Anal. 2 (1982), 241-247. 

% Check input. 
if (n < 1), error('The order n must be positive'), end 
 randn('seed',11);
% Generate the matrix. 
hs = pi/(2*(n-1)); ht = pi/(n-1); A = zeros(n,n);

for i = 1:n
    hsi = (i-1)*hs;
    for j = 1:n
        htj = (j-1)*ht;
        A(i,j) = ht*exp(hsi*cos(htj));
    end
end 
for i = 1:n
    A(i,1)=A(i,1)/2;
    A(i,n)=A(i,n)/2;
end
% Generate the right-hand side. 
y = zeros(n,1);
 y(1)=2;
for i = 2:n
    hsi = (i-1)*hs;
    y(i) = 2*sinh(hsi)/hsi;
end
% Generate the solution. 
xexact = zeros(n,1);
for j = 1:n
    htj = (j-1)*ht;
    xexact(j) = sin(htj);
end
if noise_operator==1
    nA=norm(A);
    noise = randn(n,n);       % generate noise
    noise = delta*noise*nA/norm(noise);   % adjust norm of the noise
    A_delta = A + noise;
    h=norm(A-A_delta);
    A=A_delta;
else
    h=0;
end
ny=norm(y);
noise = randn(n,1);       % generate noise
noise = delta*noise*ny/norm(noise);   % adjust norm of the noise
y_delta = y + noise; 
end
