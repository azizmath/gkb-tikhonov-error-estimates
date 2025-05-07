function [A,xexact,h,y,y_delta]= Gravity_GKB_discrete(n,delta,noise_operator)
disp('Gravity')
%Gravity Test problem: Fredholm integral equation of the first kind. 
% 
% 
% Discretization by a Nystrom method based on the trapezoidal rule
% of a first-kind Fredholm integral equation with kernel K and 
% right-hand side g given by 
 
%randn('seed',11);
% Check input. 
if (n < 1), error('The order n must be positive'), end 
 
% Generate the matrix. 
hs = 1/(n-1); ht = 1/(n-1); A = zeros(n,n);

for i = 1:n
    hsi = (i-1)*hs;
    for j = 1:n
        htj = (j-1)*ht;
        d=0.25;
        A(i,j) = ht*d*(d^2+(hsi-htj)^2)^(-3/2);
    end
end 
for i = 1:n
    A(i,1)=A(i,1)/2;
    A(i,n)=A(i,n)/2;
end

% Generate the solution. 
xexact = zeros(n,1);
for j = 1:n
    htj = (j-1)*ht;
    xexact(j,1) = sin(pi*htj)+(1/2)*sin(2*pi*htj);
end
% Generate the right-hand side. 
y=A*xexact;
%%%
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
