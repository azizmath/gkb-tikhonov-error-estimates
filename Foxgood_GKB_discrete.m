function [A,xexact,h,y,y_delta]= Foxgood_GKB_discrete(n,delta,noise_operator)
disp('Foxgood')
%Foxgood Test problem: Fredholm integral equation of the first kind. 
% 
% 
% Discretization by a Nystrom method based on the trapezoidal rule
% of a first-kind Fredholm integral equation with kernel K and 
% right-hand side g given by 
 
randn('seed',11);
% Check input. 
if (n < 1), error('The order n must be positive'), end 
 
% Generate the matrix. 
h = 1/(n-1);  A = zeros(n,n);

for i = 1:n
    hsi = (i-1)*h;
    for j = 1:n
        htj = (j-1)*h;
        A(i,j) = h*sqrt(hsi.^2+htj.^2);
    end
end 
for i = 1:n
    A(i,1)=A(i,1)/2;
    A(i,n)=A(i,n)/2;
end
% Generate the right-hand side. 
y = zeros(n,1);
for i = 1:n
    hsi = (i-1)*h;
    y(i,1) = (1/3)*((1+hsi^2)^(1.5)-hsi^3);
end
% Generate the solution. 
xexact = zeros(n,1);
for j = 1:n
    htj = (j-1)*h;
    xexact(j,1) = htj;
end
%y=A*xexact;
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
