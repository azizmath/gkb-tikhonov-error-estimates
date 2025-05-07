function [A,xexact,h,y,y_delta]= Shaw_GKB_discrete(n,delta,noise_operator)
disp('Shaw')
% Check input.
if (n < 1), error('The order n must be positive'), end
randn('seed',12);
% Initialization.
h = pi/(n-1); A = zeros(n,n);

% Compute the matrix A.
cos_ = cos(-pi/2 + [0:(n-1)]*h);
u_ = pi*sin(-pi/2 + [0:(n-1)]*h);
for i=1:n
    for j=i:n
        if (j ~= (n-i+1))
            cos_sum = cos_(i) + cos_(j);
            u = u_(i) + u_(j);
            A(i,j) = (cos_sum*sin(u)/u)^2;
            A(n-j+1,n-i+1) = A(i,j);
        end
    end
    A(i,n-i+1) = (2*cos_(i))^2;
end
A = A + triu(A,1)'; A = A*h;
for i = 1:n
    A(i,1)=A(i,1)/2;
    A(i,n)=A(i,n)/2;
end

% Compute the vectors x and b.
a1 = 2; c1 = 6; t1 =  .8;
a2 = 1; c2 = 2; t2 = -.5;
if (nargout>1)
  xexact =   a1*exp(-c1*(-pi/2 + [0:(n-1)]'*h - t1).^2) ...
      + a2*exp(-c2*(-pi/2 + [0:(n-1)]'*h - t2).^2);
  y = A*xexact;
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