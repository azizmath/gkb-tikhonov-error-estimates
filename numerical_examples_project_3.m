clear;clf;close all;clc;
format shorte;
disp('%------------------------------------------------------------%')
disp('% Error estimates for Golub-Kahan bidiagonalization–Tikhonov %')
disp('%------------------------------------------------------------%')
addpath('C:\Users\azhay\OneDrive\Desktop\Research\2nd_project_Chebfub research\chebfun');
format shortE
whichexamples = input('Select example: 1)baart,2)foxgood 3)gravity, 4)shaw, 5)wing: ');
deltas = input('Enter the level of noise: ');
noise_operator = input('Do you want to have a noise in the operator A? 0) no, 1) yes Ans: ');
disc_point = input('How many number of discretization points would you like to take? : '); %
lambda = 1e-2; %  specifies the minimal wavelength.
tau = 1e-12; % the break down paramter
ellsteps = input('How many steps of GKB process would you like to take? Ans:'); 
iv = input('1/initial vector y_delta, 2/initial vector Ay_delta: ');
noruns = input('How many runs would you like to take? Ans:');%results are averaged
whichmethods = input('Which method would you like to use? 1) Chebfun, 2) Discerte, [1,2] both? Ans:');
%mu_0 = input('Enter an initial guess for zero finder iteration:');
for method = whichmethods
    if method == 1
        disp('%---------------------------------------------------------%')
        disp('% continuous GKB-tikhonov with identity as regularization %')
        disp('%---------------------------------------------------------%')
        for exp = whichexamples
            disp('%****************************************************%')
            for delta = deltas
                disp('%---------------------------------------------------%')
                if exp == 1
                    [ke,xexact,h,y,y_delta] = Baart_GKB(delta,lambda,noise_operator);
                elseif exp==2
                    [ke,xexact,h,y,y_delta] = Foxgood_GKB(delta,lambda,noise_operator);
                elseif exp==3
                    [ke,xexact,h,y,y_delta] = Gravity_GKB(delta,lambda,noise_operator);
                elseif exp==4
                    [ke,xexact,h,y,y_delta]= Shaw_GKB(delta,lambda,noise_operator);
                else 
                    [ke,xexact,h,y,y_delta] = Wing_GKB(delta,lambda,noise_operator);
                end
                j =0;
                for ell = ellsteps
                    j=j+1;
                    gamma_ell = zeros(noruns,1);
                    mu = zeros(noruns,1);
                    RE = zeros(noruns,1);
                    for runs = 1:noruns
                        [C,U,V,gamma_ell(runs,1)] = continuous_GKB(ke,y_delta,xexact,iv,tau,ell);
                        A_h_ell = U*C*V';
                        [WW,Sigma,ZZ] = svd(C); % Compute the SVD of the bidiagonal matrix C
                        %%% Computing the regulriaztion paramter alpha using Eq 3.22
                        Qell = U*transpose(U);
                        % RHS of the eq
                        E = norm(xexact);
                        CC = 1;
                        delta_new = norm(y-y_delta); % use this delta for determing value of alpha
                        rhs = (E*gamma_ell(runs,1)+CC*delta_new)^2;
                        mu_0 = delta;
                        mu(runs,1) = newton_project3(mu_0,y_delta,U,WW,Sigma,ell,rhs); % use Newton method
                        SS = zeros(ell+1,ell+1);
                        for i=1:ell+1
                            SS(i,i) = Sigma(i,i)./(Sigma(i,i).^2+mu(runs,1));
                        end
                        z = ZZ*SS*transpose(WW)*transpose(U)*y_delta;
                        x_GKB_Tik = transpose(V)\z;
                        %figure(j);plot(xexact,'-b');hold on;plot(x_GKB_Tik,'-r');
                       % x_GKB_Tik1= pinv(transpose(V))*z;
                        RE(runs,1) = norm(xexact-x_GKB_Tik)/norm(xexact);
%                       x_GKB_Tik = VV*V_tilde*SS*U_tilde'*transpose(U)*y_delta;
%                       RE(runs,1) = norm(xexact-x_GKB_Tik)/norm(xexact);
                    end
                    gamma_ell = sum(gamma_ell)/noruns;
                    mu = sum(mu)/noruns;
                    RE = sum(RE)/noruns;
                    temp1 = floor(log10(gamma_ell));
                    temp2 = floor(log10(mu));
                    temp3 = floor(log10(RE));
fprintf('&&&$%d$&$%6.4f\\cdot 10^{%+d}$&$%6.4f\\cdot 10^{%+d}$&$%6.4f\\cdot 10^{%+d}$\\\\\n',...
 ell, gamma_ell*10^(-temp1), temp1,mu*10^(-temp2),temp2,RE*10^(-temp3),temp3);
                    %[delta,h,ell,gamma_ell, mu,RE]
                       % how to set yy= a: (b-a)/2000:b;
%                          yy=transpose(0:pi/1999:pi);
%                          for ii = 1:length(yy)
%                              fprintf(['(%e,%e)%%\n'],ii,x_GKB_Tik(yy(ii)));
%                          end
                end
            end
        end
    else
          disp('%---------------------------------------------------------%')
        disp('% Discrete GKB-tikhonov with identity as regularization %')
        disp('%---------------------------------------------------------%')
        for exp=whichexamples
        disp('%****************************************************%')
            for delta = deltas
                disp('---------------------------------------------------')
                for n = disc_point
                    disp('---------------------------------------------------')
                    % generate linear discrete ill-posed problem 
                    if exp == 1
                        [A,xexact,h,y,y_delta] = Baart_GKB_discrete(n,delta,noise_operator); 
                    elseif exp==2
                        [A,xexact,h,y,y_delta] = Foxgood_GKB_discrete(n,delta,noise_operator);
                    elseif exp==3
                        [A,xexact,h,y,y_delta] = Gravity_GKB_discrete(n,delta,noise_operator); 
                    elseif exp==4
                        [A,xexact,h,y,y_delta] = Shaw_GKB_discrete(n,delta,noise_operator);
                    else 
                        [ke,xexact,h,y,y_delta] = Wing_GKB(delta,lambda,noise_operator);
                    end
                    delta
                    temp4 = floor(log10(h));
                    fprintf('$%d$&$%6.4f\\cdot 10^{%+d}$\\\\\n',n,h*10^(-temp4),temp4);
                    for ell = ellsteps
                        gamma_ell = zeros(noruns,1);
                        mu = zeros(noruns,1);
                        RE = zeros(noruns,1);
                        for runs = 1:noruns
                            [C,U,V,Appx_operaotr,gamma_ell(runs,1)] = discete_GKB(A,y_delta,iv,tau,ell);
                            A_h_ell = U*C*V';
                            [WW,Sigma,ZZ] = svd(C); % Compute the SVD of the bidiagonal matrix C
                            %%% Computing the regulriaztion paramter alpha using Eq 3.22
                            %%% find minimal reg_parameter using Newton
                            delta_new = norm(y-y_delta); % use this delta for determing value of alpha
                            D = norm(xexact);
                            CC = 1;
                            rhs = (D*gamma_ell(runs,1)+CC*delta_new)^2;
                            mu(runs,1) = newton_project3(mu_0,y_delta,U,WW,Sigma,ell,rhs);
                            SS = zeros(ell+1,ell+1);
                            for i=1:ell+1
                                SS(i,i) = Sigma(i,i)./(Sigma(i,i).^2+mu(runs,1));
                            end
                            z = ZZ*SS*transpose(WW)*transpose(U)*y_delta;
                            %x_GKB_Tik = transpose(V)\z;
                            x_GKB_Tik= pinv(V')*z;
                            RE(runs,1) = norm(xexact-x_GKB_Tik)/norm(xexact);
                        end
                        gamma_ell = sum(gamma_ell)/noruns;
                        mu = sum(mu)/noruns;
                        RE = sum(RE)/noruns;
                        temp1 = floor(log10(gamma_ell));
                        temp2 = floor(log10(mu));
                        temp3 = floor(log10(RE));
fprintf('&&&$%d$&$%6.4f\\cdot 10^{%+d}$&$%6.4f\\cdot 10^{%+d}$&$%6.4f\\cdot 10^{%+d}$\\\\\n',...
ell, gamma_ell*10^(-temp1), temp1,mu*10^(-temp2),temp2,RE*10^(-temp3),temp3);
                        %[delta,n,h,ell,gamma_ell,mu,RE]
%                         yy=transpose(0:(1)/1999:1);
%                         for ii = 1:length(yy)
%                             fprintf(['(%e,%e)%%\n'],ii,x_GKB_Tik((ii),1));
%                         end
%                         clf;plot(xexact);hold on;plot(x_GKB_Tik);
                   
                    end
                end
            end
        end
    end
end
