function [ke,xexact,h,y,y_delta]= Shaw_GKB(delta,lambda,noise_operator)
   disp('Shaw')
        u = chebfun2(@(s,t) pi()*(sin(s)+sin(t)),[-pi()/2 pi()/2 -pi()/2 ...
        pi()/2],'eps',1e-16,'splitting','on');
        q = chebfun(@(x) (sin(x)/x).^2, [-2,2],'vectorize');
        qu = chebfun2(@(s,t) q(u(s,t)), ...
        [-pi()/2 pi()/2 -pi()/2 pi()/2],'eps',1e-16,'splitting','on');
        cs = chebfun2(@(s,t) (cos(s)+cos(t)),[-pi()/2 pi()/2 -pi()/2 pi()/2],'eps',1e-16);
        ke = cs.*qu;
%          no_points_per_col = 51;
%          l = linspace(ke.domain(1),ke.domain(2),no_points_per_col);
%          q = linspace(ke.domain(3),ke.domain(4),no_points_per_col);
%         for ii = l
%             	for jj = q
%                 		fprintf('(%e,%e,%e)%%\n',ii,jj,ke(ii,jj));
%             	end
%         end
%         zmax=max(max(ke))
%         zmin=min(min(ke)) 
        a1 = 2;a2 = 1;c1 = 6; c2 = 2; t1 =  .8;t2 = -.5;
        xexact = chebfun(@(t) a1*exp(-c1*(t-t1)^2)+a2*exp(-c2*(t-t2)^2),[-pi()/2 pi()/2],'eps',1e-16,'vectorize');
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
        else
            h=0;
        end
        ny=norm(y);
        noise = randnfun(lambda,[s1 s2]);       % generate noise
        noise = delta*noise*ny/norm(noise);   % adjust norm of the noise
        y_delta = y + noise; 
    end