function [C,U,V,gamma_ell] = continuous_GKB(ke,y_delta,xexact,iv,tau,ell)
    %%% GKB process
    Beta=zeros(ell+1,1);
    Alpha=zeros(ell+1,1);
    U=[];
    V=[];
    C = zeros(ell+1,ell+1);
    %%%% Define the operators A_h and A_h'
     A_h = @(u) sum(transpose(ke)*u,2) ;
     TA_h = @(u) sum(ke*u,2) ;
    %%%%%%%%%%%%%%
    if iv == 1
        Beta(1,1)=norm(y_delta);
        U=[y_delta./Beta(1,1)];
    else
       yy = TA_h(y_delta);
       yy=A_h(yy);
       Beta(1,1)=norm(yy);
       U=[yy./Beta(1,1)];
    end
    v=TA_h(U(:,1));%  v=sum(ke*U(:,1),2);
    Alpha(1,1)=norm(v);
    V=[v./Alpha(1,1)];
    j=1;
    while (abs(Beta(j,1))>tau)&& (Alpha(j,1)>tau) && (j<ell+1)
          j=j+1;
          %u=sum(transpose(ke)*V(:,j-1),2)-Alpha(j-1,1)*U(:,j-1);
          u=A_h(V(:,j-1))-Alpha(j-1,1).*U(:,j-1);
          % Reorthogonalization step:
          sqrteps=sqrt(eps);
            for i=1:j-1
                C1=U(:,1:i)'*u;
                u=u-U(:,1:i)*C1;
                while norm(C1)> sqrteps
                    C1=U(:,1:i)'*u;
                    u=u-U(:,1:i)*C1;
                end
            end
        Beta(j,1)=norm(u);
        U=[U u./Beta(j,1)];
      %v=sum(ke*U(:,j),2)-Beta(j,1)*V(:,j-1);
       v=TA_h(U(:,j))-Beta(j,1).*V(:,j-1); 
        % Reorthogonalization step:
        for i=1:j-1
            C1=V(:,1:i)'*v;
            v=v-V(:,1:i)*C1;
            while norm(C1)> sqrteps
                C1=V(:,1:i)'*v;
                v=v-V(:,1:i)*C1;
            end
         end
         Alpha(j,1)=norm(v);
         V=[V v./Alpha(j,1)];
    end
   if j < ell
       warning(sprintf('GKB breakdown happened at it. %d',j));
       ell=j-1;
       C = zeros(ell+1,ell+1);
   end
   for i=1:ell+1
        C(i,i)=Alpha(i,1);
   end
   for i = 1:ell
        C(i+1,i)=Beta(i+1,1);
   end
   A_h_ell=U*C*transpose(V);
   gamma_ell=norm(A_h(xexact)-A_h_ell*xexact);
   %gamma_ell=norm(sum(transpose(ke)*xexact,2)-A_h_ell*xexact);
   %% 
%    T = C(1:ell+1,1:ell)*C(1:ell,1:ell)';
%    tt  = diag(T);
%    ttt = diag(T,-1);
%    for i= 1: ell
%        fprintf(['(%e,%e)%%\n'],i,tt(i,1));
%    end
%    for i= 1: ell
%        fprintf(['(%e,%e)%%\n'],i,ttt(i,1));
%    end
end


