% This function helps us to compute f(alpha) and f'(alpha) for applying
% Newton method to find a good value of the regulrization parameter in 
% tikhonov regularization  
function [f_alpha,fprime_alpha]= fofalpha(mu,y_delta,U,U_tilde,Sigma,ell,rhs)
    ss=diag(Sigma);
    s=ss.^2;
    D=zeros(ell+1);
    for i=1:ell+1
        D(i,i)=(mu^3)./(s(i,1)+mu).^3;
    end
    %f_alpha=transpose(y_delta)*U*Z*Iq*D*Iq*Z'*transpose(U)*y_delta-rhs;
    f_alpha=transpose(y_delta)*U*U_tilde*D*U_tilde'*transpose(U)*y_delta-rhs;
    if nargout > 1
        D1=zeros(ell+1);
        for i=1:ell+1
           % D1(i,i)=((3*mu^2)./((s(i,1)+mu).^3))-((3*mu.^3)./((s(i,1)+mu).^4));
%             ff = 3*(mu^(2)*s(i,1)+ mu^(3)-mu^(2));
%             D1(i,i)= ff./((s(i,1)+alpha).^(4));
            D1(i,i) = (3*mu^2)./((s(i,1)+mu).^4);

        end
        %fprime_alpha=transpose(y_delta)*U*Z*Iq*D1*Iq*Z'*transpose(U)*y_delta;
        fprime_alpha=transpose(y_delta)*U*U_tilde*D1*U_tilde'*transpose(U)*y_delta;

    end
end