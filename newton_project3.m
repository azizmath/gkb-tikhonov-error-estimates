% We are going to implement Newton's method to choose a good value for
% regularization parameter alpha in Tikhonov regularization
function alpha = newton_project3(alpha_0,y_delta,U,U_tilde,Sigma,ell,rhs)
% Set defaults.
thr = sqrt(eps);  % Relative stopping criterion.
it_max = 5000;      % Max number of iterations.
% Initialization.
if (alpha_0 < 0)
  error('Initial guess alpha_0 must be nonnegative')
end
alpha = alpha_0; 
step = 1; 
it = 0;
while (abs(step) > thr*alpha && abs(step) > thr && it < it_max)
    it = it+1;
    [f_alpha,fprime_alpah]= fofalpha(alpha,y_delta,U,U_tilde,Sigma,ell,rhs);
    step=(f_alpha/fprime_alpah);
    alpha = alpha-step; 
    % If alpha < 0 then restart with smaller initial guess.
  if (alpha < 0) 
      alpha = 0.5*alpha_0; 
      alpha_0 = 0.5*alpha_0;
  end   
end
% Terminate with an error if too many iterations.
if (abs(step) > thr*alpha && abs(step) > thr)
  %error(['Max. number of iterations (',num2str(it_max),') reached'])
end
it
end