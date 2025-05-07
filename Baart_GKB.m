function [ke,xexact,h,y,y_delta]= Baart_GKB(delta,lambda,noise_operator)
    disp('Baart')
    
    ke = chebfun2(@(s,t) exp(s*cos(t)) ,[0 pi/2 0 pi],'eps',1e-16,'vectorize');% Kernel
     no_points_per_col = 51;
     l = linspace(ke.domain(1),ke.domain(2),no_points_per_col);
     q = linspace(ke.domain(3),ke.domain(4),no_points_per_col);
%          for ii = l
%              for jj = q
%                  fprintf('(%e,%e,%e)%%\n',ii,jj,ke(ii,jj));
%              end
%           end
%          zmax=max(max(ke))
%          zmin=min(min(ke))
        
    xexact=chebfun(@(t) sin(t),[0 pi],'eps',1e-16,'vectorize'); % Exact solution
    y = chebfun(@(s) 2*sinh(s)/s,[0 pi/2],'eps',1e-16,'vectorize'); %Error free rhs function
    s1 = ke.domain(1); s2 = ke.domain(2); 
    t1 = ke.domain(3); t2 = ke.domain(4);
    %figure(1);plot(ke);
    if noise_operator==1
        nke=norm(ke);
        noise = randnfun2(lambda,[s1 s2 t1 t2]); % 375s
        noise=delta*noise*nke/norm(noise);
        ke_delta = ke+noise;
        h=norm(sum(transpose(ke)*xexact,2)-sum(transpose(ke_delta)*xexact,2));
        ke=ke_delta;
    else
        h=0;
    end
   %figure(2);plot(ke);
ny=norm(y);
noise = randnfun(lambda,[s1 s2]);       % generate noise
noise = delta*noise*ny/norm(noise);   % adjust norm of the noise
y_delta = y + noise; 
   
end