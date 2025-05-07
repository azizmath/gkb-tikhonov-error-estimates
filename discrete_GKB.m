function [C,U,V,A_hell,gamma_ell] = discete_GKB(A_h,y_delta,iv,tau,ell)
    %GKB
    [m,n] = size(A_h);
    Beta = zeros(ell+1,1);
    Alpha = zeros(ell+1,1);
    U = zeros(m,ell+1);
    V = zeros(n,ell+1);
    if iv == 1
        Beta(1,1) = norm(y_delta);
        U(:,1) = y_delta/Beta(1,1);
    else
        yy=A_h*(A_h'*y_delta);
        Beta(1,1) = norm(yy);
        U(:,1) = yy/Beta(1,1);
    end
    v = A_h'*U(:,1);
    Alpha(1,1) = norm(v);
    V(:,1) = v/Alpha(1,1);
    j = 1;
    while (abs(Beta(j,1))>tau)&& (abs(Alpha(j,1))>tau) && (j<ell+1)
        j = j+1;
        u = A_h*V(:,j-1) - Alpha(j-1,1)*U(:,j-1);
        % Reorthogonalization step:
        sqrteps = sqrt(eps);
        for i = 1:j-1
            C1 = U(:,1:i)'*u;
            u = u - U(:,1:i)*C1;
            while norm(C1)> sqrteps
                  C1 = U(:,1:i)'*u;
                  u = u- U(:,1:i)*C1;
            end
        end
        Beta(j,1) = norm(u);
        U(:,j) = u/Beta(j,1);
        v = A_h'*U(:,j) - Beta(j,1)*V(:,j-1);
        % Reorthogonalization step:
        for i = 1:j-1
            C1 = V(:,1:i)'*v;
            v = v-V(:,1:i)*C1;
            while norm(C1)> sqrteps
                  C1 = V(:,1:i)'*v;
                  v = v - V(:,1:i)*C1;
            end
        end
        Alpha(j,1) = norm(v);
        V(:,j) = v/Alpha(j,1);
    end
    if j < ell
        warning(sprintf('GKB breakdown happened at it. %d',j));
        ell = j-1;
        U = U(:,1:ell+1);
    end
    C = zeros(ell+1,ell+1);
    for i = 1:ell+1
        C(i,i) = Alpha(i,1);
    end
    for i = 1:ell
        C(i+1,i)=Beta(i+1,1);
    end
    A_hell = U*C*V';
    gamma_ell = norm(A_h-A_hell);
end
    