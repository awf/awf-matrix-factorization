function [x, f, log_data] = awf_levmarq(x, func, opts)

% AWF_LEVMARQ   Home-grown LM with line search
%               [x, J, OPTS] = awf_levmarq(x, f)
%               x:  Initial estimate
%               func:  Function to be called.

if nargin == 1 && strcmp(x, 'opts')
  opts.MaxIter = 100;
  opts.MaxFunEvals = 1000;
  opts.Display = 'iter';
  opts.DECOMP_LU = 0;  % Use backslash or PCG?
  opts.LAMBDA_MIN = 1e-10;
  opts.LAMBDA_DECREASE = 2;
  opts.LAMBDA_MAX = 1e8;
  opts.LAMBDA_INCREASE_BASE = 10;
  opts.USE_LINMIN = 1;
  opts.CHECK_DERIVATIVES = 1;
  opts.USE_JTJ = 0;
  x = opts;
  return
end

if nargin < 3
  opts = awf_levmarq('opts');
end

switch opts.Display
  case {'none','off'}
    VERBOSE = 0;
  case 'iter'
    VERBOSE = 3;
  case 'final'
    VERBOSE = 1;
  otherwise
    error(['Display [' opts.Display '] not recognized']);
end

sumsq = @(x) sum(x.^2);

% Call function
x = x(:);
[e,J] = func(x);
f = sum(e.^2);

nparams = length(x);
nresiduals = numel(e);

%% Optionally check the Jacobian
if opts.CHECK_DERIVATIVES
  check_derivatives(func,x,J);
end

if VERBOSE >= 1
  fprintf('awf_levmarq: Beginning LM params %d, residuals %d, error = %g\n', ...
    nparams, nresiduals, f)
end

%% Begin the optimization

if VERBOSE >= 3
  fprintf('awf_levmarq: Initial f %g\n', f);
end

% set up options

lm_lambda = 1e-6;
log_data = [lm_lambda f 0 1];
funevals = 0;
iter = 0;
Id = speye(nparams);
while true
  % This outer loop is called for each Jacobian computation
  % We assume that computing J is an expensive operation so 
  % try to squeeze as much out of each J as possible.
  
  % Estimate Jacobian norm and preconditioner for PCG
  diag_JtJ = sum(J.^2,1)';
  norm_estimate = full(mean(diag_JtJ));
  if ~opts.DECOMP_LU
    Preconditioner = spdiags(diag_JtJ,0,nparams,nparams);
  end

  % Do a variety of line-searches using this Jacobian,
  % varying lambda each time, as well as LAMBDA_INCREASE 
  % This attempts to get to large lambda as fast as possible
  % if we are rejecting, so that the only convergence 
  % criterion is large lambda (i.e. we stop when gradient
  % descent with a tiny step cannot reduce the function)

  % Reset LAMBDA_INCREASE
  opts.LAMBDA_INCREASE = opts.LAMBDA_INCREASE_BASE;
  
  found_reduction = 0;
  while ~found_reduction
    if VERBOSE >= 3
      fprintf('awf_levmarq: iter %2d, lambda=%6.2e, norm=%4.2e, ', iter, lm_lambda, norm_estimate);
      fprintf('solve ');
    end
    
    % The computations below assume J is quite portrait, so
    % tend to work with J'*J -- may need to be changed for other
    % shapes/sparsity patterns.  For example, replacing
    %   (J' J + mu I) \ J' e 
    % by
    %   [J ; sqrt(mu) I] \ [e; O]
    if opts.DECOMP_LU
      N = size(J,2);
      if iter == 1 && size(J,1) < size(J,2)
        fprintf(2, 'awf_levmarq: warning: J not portrait\n');
      end
      
      if opts.USE_JTJ
        JtJ = J'*J;
        %Sometimes this may be faster:
        %JtJ = J'*J;
        %indices = awf_eyeind(nparams);
        %JtJ(indices) = JtJ(indices) + (lm_lambda*norm_estimate);
        dx = -(JtJ + lm_lambda*Id) \ (J'*e);
      else
        dx = -[J; lm_lambda*Id]\[e; zeros(nparams,1)];
      end
    else
      % pcg
      AugmentedJtJ = J'*J + (lm_lambda)*speye(nparams);
      PCG_ITERS = 20;
      [dx, ~] = pcg(AugmentedJtJ, -(J'*e), 1e-7 * norm_estimate, PCG_ITERS); %;, Preconditioner);
      if ~all(isfinite(dx))
        error('infinite dx...');
      end
    end

    if opts.USE_LINMIN
        % linmin: Line search along dx
        % Define 1D function for linmin
        f1d = @(t) sumsq(func(x + t * dx));
        % Set options for linmin
        fminbnd_options = optimset('fminbnd');
        fminbnd_options.TolX = 1e-5;  % We don't need this search to be very precise.
        fminbnd_options.MaxFunEvals = min(120, opts.MaxFunEvals - funevals);
        fminbnd_options.Display = 'off';
        
        % Execute linmin
        [t,f_test,~,fminbnd_output] = fminbnd(f1d, .2, 10, fminbnd_options);
        x_test = x + t * dx;
        funevals = funevals + fminbnd_output.funcCount;
        if VERBOSE >= 3, fprintf('linmin [t=%4.2f], f=%d, ', t, funevals); end
        % record log data for both rejections and acceptances,
        log_data = [log_data; lm_lambda, f_test, t, funevals];

    else
        x_test = x + dx;
        f_test = sumsq(func(x_test));
        funevals = funevals + 1;
        % record log data for both rejections and acceptances,
        log_data = [log_data; lm_lambda, f_test, 1, funevals];

    end

    if f_test < f
      if VERBOSE >= 3, drawnow; fprintf('Accept %g\n', f_test); end
      if lm_lambda > opts.LAMBDA_MIN
        lm_lambda = lm_lambda / opts.LAMBDA_DECREASE;
      end
      f = f_test;
      break
    else
      if VERBOSE >= 3, drawnow; fprintf('**rej* %g\n', f_test); end
      % This means linmin failed -- the whole direction is wrong, so LAMBDA_MAX should be large
      if lm_lambda > opts.LAMBDA_MAX
        endmsg = 'exceeded lambda_max';
        found_reduction = 1;
        break
      else
        lm_lambda = lm_lambda * opts.LAMBDA_INCREASE;

        % successive rejects should be exponential, so now LAMBDA_INCREASE
        % becomes 1e1, 1e2, 1e4, 1e8 etc
        opts.LAMBDA_INCREASE = opts.LAMBDA_INCREASE^2;
        % continue
      end
    end
  end

  % terminate?
  if found_reduction
    break
  end

  if funevals > opts.MaxFunEvals
    endmsg = '>MaxFunEvals';
    break
  end
  
  % use log_data to check cgce
  if size(log_data,1) > 10
    last_f = log_data(end, 2);
    mid_f = log_data(ceil(end/2), 2);
    if abs(last_f - mid_f) < 1e-8*last_f
      endmsg = 'flatlined';
      break
    end
  end

  %% Need a new iter.

  % Record new x
  x = x_test;

  % Recompute jacobian
  [e,J] = func(x);
  funevals = funevals + 1;
  iter = iter + 1;

  % Record new f
  f = f_test;

  % Just check that it's the same as above
  if abs(sumsq(e) - f) > 1e-8
    keyboard
  end

  if (VERBOSE > 1) && (VERBOSE < 3)
    fprintf(' %g', f);
  end
end
if (VERBOSE >= 1) && (VERBOSE < 3)
  if VERBOSE == 1, fprintf('awf_levmarq:'); end
  fprintf(' done [%s], rms=%g, %d iterations, %d evals\n', endmsg, rms(e), iter, funevals);
end

%%
function r = rms(e)
r = sqrt(mean(e.^2));

%%
function check_derivatives(f,x,J)
delta=1e-4;
scale = 1/2/delta;
[n,p] = size(J);
assert(p == numel(x))
% Check derivatives OK
fdJ=0*J;
for k=1:p
  e = zeros(size(x));
  e(k) = delta; 
  fdJ(:,k) = (f(x+e) - f(x-e))*scale;
end
if max((fdJ(:) - J(:)).^2) > 1e-7
  error('derivatives wrong')
end
