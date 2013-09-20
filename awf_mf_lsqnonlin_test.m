function allres = awf_mf_lsqnonlin_test(tests, save_file)
% AWF_MF_LSQNONLIN_TEST  Run some factorization tests
%      Input parameter TESTS is a struct array with fields
%             M       Measurement matrix
%             W       Valid elements mask
%             A0      Initial estimate for A
%             B0      Initial estimate for B
%             regularizer_lambda  Times norm([A B], 'fro')
%             gauge_fix_weight    Times norm(A(1:r,1:r) - eye(r), 'fro')
%             alg                 Algorithm, e.g. 'lm', 'lsq', 'altdn01a'
%             lsopts              Options for inner solver
%      Input parameter SAVE_FILE is a filename to which to save results
%
%      Output is TESTS augmented with e.g.
%             rms       Measurement error
%             time      Wall-clock time in seconds
%             lsout     Inner solver diagnostics
%

if nargin == 0
  evalin('base', 'run_tests');
  return;
end

opts = awf_mf_lsqnonlin('opts');

for k=1:length(tests)
  alg=tests(k).alg;
  opts.regularizer_lambda = tests(k).regularizer_lambda;
  opts.gauge_fix_weight = tests(k).gauge_fix_weight;
  
  M = tests(k).M;
  W = tests(k).W;
  A0 = tests(k).A0;
  B0 = tests(k).B0;
  
  fprintf('%s, reg_lambda=%g, gauge_wt=%g... ', alg, ...
    opts.regularizer_lambda, opts.gauge_fix_weight);
  tic
  switch alg
    case 'lsq'
      opts.alg = 'lsqnonlin';
      opts.lsopts.Algorithm = 'trust-region-reflective';
      [A,B,lsout] = awf_mf_lsqnonlin(W,M,A0,B0, opts);
    case 'lm'
      opts.alg = 'lsqnonlin';
      opts.lsopts.Algorithm = 'levenberg-marquardt';
      [A,B,lsout] = awf_mf_lsqnonlin(W,M,A0,B0, opts);
    case 'awflm'
      opts.alg = 'awf';
      opts.awopts.MaxFunEvals = 10000;
      opts.awopts.USE_LINMIN = 0;
      opts.awopts.Display = 'none';
      [A,B,lsout] = awf_mf_lsqnonlin(W,M,A0,B0, opts);
    case 'scg'
      opts.alg = 'scg';
      opts.lsopts.Algorithm = 'scg';
      [A,B,lsout] = awf_mf_lsqnonlin(W,M,A0,B0, opts);
    case 'dw'
      % [U, V] = damped_wiberg(Y, H, r, Vini) uses Vini as initial values of V
      % instead of random initial values.
      maxiter = tests(k).wiberg_iters;
      [A,B,~,iter] = damped_wiberg(M,W,size(A0,2), B0, 1e-6, maxiter, 0);
      clear lsout
      lsout.iterations = iter;
    case 'wl'
      % 'wiberg-lm'
      [Aw,Bw] = damped_wiberg(M,W,size(A0,2), B0, 1e-6, tests(k).wiberg_iters, 0);
      opts.lsopts.Algorithm = 'levenberg-marquardt';
      [A,B,lsout] = awf_mf_lsqnonlin(W,M,Aw,Bw, opts);
    otherwise
      [A,B] = amb_mfwmd(A0,B0,M,double(W), ...
        'algorithm', alg,...
        'L', opts.regularizer_lambda, ...  % total regularization
        'graph', 0, ...
        'iterations', 5000, ...
        'tolerance', 1e-8);
  end
  t=toc;
  rms = sqrt(norm((M - A*B').*W, 'fro')^2/nnz(W));
  fprintf('rms=%g, time=%g sec\n', rms, t);
  
  out = tests(k);
  out.rms = rms;
  out.time = t;
  out.A = A;
  out.B = B;
  out.lsout = lsout;
  
  struct_push_back allres out
  
  save(save_file,'allres');
  
  % Plot output
  Mrecon = A*B';
  hold off
  plot(Mrecon(1:2:end,:), Mrecon(2:2:end,:),'k-')
  hold on
  plot(M(1:2:end,:), M(2:2:end,:),'r.')
  title(sprintf('%s, \\lambda=%g, GW=%g, rms=%g, time=%g', ...
    alg, out.regularizer_lambda, out.gauge_fix_weight, out.rms, out.time));
  drawnow
end
