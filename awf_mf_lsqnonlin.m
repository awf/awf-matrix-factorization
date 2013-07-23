function [Aout,Bout,output] = awf_mf_lsqnonlin(W,M,A0,B0, opts)

if 0
  %%
  vec = @(x) x(:);
  rms = @(Ap,Bp) sqrt(sum(vec(W.*(M-Ap*Bp')).^2)/nnz(W));

  load_dino
  %M = M(1:36,1:150);
  %W = W(1:36,1:150);
  rng(1);
  A0 = randn(size(M,1), 4);
  B0 = randn(size(M,2), 4);
  
  profile on
  opts = awf_mf_lsqnonlin;
  opts.regularizer_lambda = 0;
  opts.gauge_fix_weight = 0;
  
  opts.lsopts.Algorithm = 'awf';
  opts.awopts.MaxFunEvals = 4000;
  opts.awopts.Display = 'final';
  opts.awopts.USE_LINMIN = 0;

  opts.awopts.SCHUR_SPLIT = 0;%4*72;
  opts.awopts.USE_JTJ = 1;
  opts.awopts.DECOMP_LU = 0;
  tic
  [A1,B1,out2] = awf_mf_lsqnonlin(W,M,A0,B0, opts);
  fprintf('awf rms = %g\n', rms(A1,B1));
  toc

  
  clf
  loglog(9+(1:size(out1.log_data,1)), out1.log_data(:,2))
  hold on
  loglog(9+(1:size(out2.log_data,1)), out2.log_data(:,2))
  
  profile viewer
  
  tic
  opts = awf_mf_lsqnonlin;
  % opts.lsopts.Display = 'iter';
  opts.lsopts.MaxIter = 800;
  [A,B,out] = awf_mf_lsqnonlin(W,M,A0,B0, opts);
  fprintf('lsq rms = %g\n', rms(A,B));
  toc

  out
  
  %%
  clf
  plot_tracks(M,W,'k.');
  hold on
  %set(plot_tracks(A*B','r-'), 'color', 'r');
  set(plot_tracks(A1*B1','b-'), 'color', 'b');
end

if nargin == 0
  opts.regularizer_lambda = 1e-5;
  opts.gauge_fix_weight = 1e-5;
  opts.lsopts = optimset('lsqnonlin');
  opts.lsopts.Algorithm = 'levenberg-marquardt';
  opts.lsopts.LargeScale = 'on';
  opts.lsopts.Display = 'off';
  opts.lsopts.MaxIter = 800;
  opts.lsopts.TolFun = 1e-8;
  opts.awopts = awf_levmarq('opts');
  Aout = opts;
  return
end

[i,j] = find(W);
Mvalues = M(W~=0)';

% Condition
scale = std(Mvalues);
Mvalues = Mvalues / scale;

nresiduals = length(Mvalues);

r = size(A0,2);
[m,n] = size(M);
mlp_assert_equal m size(A0,1)
mlp_assert_equal n size(B0,1)
mlp_assert_equal r size(B0,2)

% M - A'*B  (transposed formulation as matlab is column-major)
A0 = A0';
B0 = B0';

% Pack into vector
v0 = pack(A0,B0);

% Check Jacobian
if 0
  A = A0;
  B = B0;
  disp('testing jacobian');
  tic;
  Jslow = zeros(nresiduals, (m+n)*r);
  for k=1:nresiduals
    Jslow(k,(i(k)-1)*r + (1:r)) = -B(:,j(k));
    Jslow(k,(j(k)-1)*r + (1:r) + m*r) = -A(:,i(k));
  end
  toc
  tic
  [~,Jfast] = f_data(pack(A,B));
  toc
  if rand <.1
    disp('J\x');
    tic
    Jfast \ randn(size(Jfast,1),1);
    toc
  end
  mlp_assert_equal Jfast sparse(Jslow) 1e-8
  
  [A0test,B0test] = unpack(v0);
  mlp_assert_equal A0test A0
  mlp_assert_equal B0test B0
end

gauge_fix_target = v0(1:r*r);

if strcmp(opts.lsopts.Algorithm, 'awf')
  % Call awf_levmarq
  awopts = opts.awopts;
  awopts.CHECK_DERIVATIVES = 0;
  [v, ~, log_data] = awf_levmarq(v0, @f, awopts);
  output.exitflag = 1;
  output.log_data = log_data;
else
  % Call lsqnonlin
  lsopts = opts.lsopts;
  lsopts.Jacobian = 'on';
  lsopts.PrecondBandwidth = r;
  
  [v,~,~, exitflag, lmout] = lsqnonlin(@f, v0, [], [], lsopts);
  output = lmout;
  output.exitflag = exitflag;
end

[Aout,Bout] = unpack(v);
Aout = Aout' * scale;
Bout = Bout';

%%
  function v = pack(A,B)
    v = [A(:); B(:)];
  end

  function [A,B] = unpack(v)
    % Unpack
    A = reshape(v(1:m*r), r,m);
    B = reshape(v((1:n*r)+m*r), r,n);
  end

  function [residuals, Jacobian] = f(v)
    if nargout == 1
      rd = f_data(v);
      if opts.regularizer_lambda > 0
        rr = f_reg(v);
      else
        rr=[];
      end
      if opts.gauge_fix_weight > 0
        ru = f_gauge(v);
      else
        ru=[];
      end
      residuals = [
        rd
        opts.regularizer_lambda/scale * rr
        opts.gauge_fix_weight/scale * ru
        ];
    else
      % Want f and Jacobian
      [rd, Jd] = f_data(v);
      if opts.regularizer_lambda > 0
        [rr, Jr] = f_reg(v);
      else
        rr=[];Jr=[];
      end
      if opts.gauge_fix_weight > 0
        [ru, Ju] = f_gauge(v);
      else
        ru=[];Ju=[];
      end
      
      residuals = [
        rd
        opts.regularizer_lambda/scale * rr
        opts.gauge_fix_weight/scale * ru
        ];
      Jacobian = [
        Jd
        opts.regularizer_lambda/scale  * Jr
        opts.gauge_fix_weight/scale  * Ju
        ];
    end
  end

  function [residuals, Jacobian] = f_data(v)
    [A,B] = unpack(v);
    
    % Compute residual
    residuals = Mvalues - sum(A(:,i).*B(:,j));
    
    residuals = residuals';
    
    if nargout == 1
      return
    end
    
    % Gradient of kth residual:
    % rk = M(i(k), j(k)) - A(:,i(k))'*B(:,j(k))
    % drk/dAp = [p=i(k)]*(-B(:,j(k)))
    % drk/dBp = [p=j(k)]*(-A(:,i(k)))
    
    k_reps = (1:nresiduals)' * ones(1,r);
    i_reps = bsxfun(@plus, (i-1)*r, 1:r);
    j_reps = bsxfun(@plus, (j-1)*r, 1:r);
    JA = sparse(k_reps', i_reps', -B(:,j), nresiduals, m*r);
    JB = sparse(k_reps', j_reps', -A(:,i), nresiduals, n*r);
    Jacobian = [JA JB];
  end

  function [residuals, Jacobian] = f_reg(v)
    residuals = v.^2;
    if nargout > 1
      Jacobian = speye((m+n)*r);
    end
  end

  function [residuals, Jacobian] = f_gauge(v)
    residuals = v(1:r*r) - gauge_fix_target;
    if nargout > 1
      Jacobian = speye(r*r, (m+n)*r);
    end
  end

end
