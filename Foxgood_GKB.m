function [ke,xexact,h,y,y_delta]= Foxgood_GKB(delta,lambda,noise_operator)
disp('Foxgood')
        ke = chebfun2(@(s,t) (s^2+t^2)^(1/2) ,[0 1 0 1],'eps',1e-16,'vectorize');
%         no_points_per_col = 51;
%         l = linspace(ke.domain(1),ke.domain(2),no_points_per_col);
%         q = linspace(ke.domain(3),ke.domain(4),no_points_per_col);
%         for ii = l
%             	for jj = q
%                 		fprintf('(%e,%e,%e)%%\n',ii,jj,ke(ii,jj));
%             	end
%         end
%         zmax=max(max(ke))
%         zmin=min(min(ke)) 
        xexact= chebfun(@(t) t,[0 1],'eps',1e-16,'vectorize');
        y =chebfun(@(s) 1/3*((1+s^2)^(3/2)-s^3) ,[0 1],'eps',1e-16,'vectorize');
        s1 = ke.domain(1); s2 = ke.domain(2); 
        t1 = ke.domain(3); t2 = ke.domain(4);
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
        ny=norm(y);
        noise = randnfun(lambda,[s1 s2]);       % generate noise
        noise = delta*noise*ny/norm(noise);   % adjust norm of the noise
        y_delta = y + noise; 
    end