function [ke,xexact,h,y,y_delta]= Gravity_GKB(delta,lambda,noise_operator)
disp('Gravity')
        d=0.25;
        ke = chebfun2(@(s,t) d*(d^2 + (s-t)^2)^(-3/2),[0 1 0 1],'eps',1e-16,'vectorize','splitting','on');
         %figure(1);plot(ke);
         %no_points_per_col = 51;
         %l = linspace(ke.domain(1),ke.domain(2),no_points_per_col);
         %q = linspace(ke.domain(3),ke.domain(4),no_points_per_col);
         %for ii = l
          %  for jj = q
           %     fprintf('    (%e,%e,%e)%%\n',ii,jj,ke(ii,jj));
           % end
         %end
        xexact = chebfun(@(t) sin(pi*t) + 0.5*sin(2*pi*t),[0 1],'eps',1e-16,'vectorize','splitting','on');
        y = sum(transpose(ke)*xexact,2);
        s1 = ke.domain(1); s2 = ke.domain(2); 
        t1 = ke.domain(3); t2 = ke.domain(4);
        if noise_operator==1
            nke=norm(ke);
            noise = randnfun2(lambda,[s1 s2 t1 t2]); % 375s
            noise=delta*noise*nke/norm(noise);
            ke_delta = ke+noise;
            
            h=norm(sum(transpose(ke)*xexact,2)-sum(transpose(ke_delta)*xexact,2));
            ke=ke_delta;
            %figure(2);plot(ke);
        else
            h=0;
        end
        ny=norm(y);
        noise = randnfun(lambda,[s1 s2]);       % generate noise
        noise = delta*noise*ny/norm(noise);   % adjust norm of the noise
        y_delta = y + noise; 
    end