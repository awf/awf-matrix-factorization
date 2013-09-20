
clear tests

for sd = 1:100
  rng(sd)
  for dataset = {'dino', 'giraffe'}
    switch dataset{1}
      case 'giraffe'
        load_giraffe
        r = 6;
      case 'dino'
        [M,W, Mx, My] = load_dino;
        r = 4;
    end
    
    if 0
      crop = 30;
      M = M(1:crop,1:2*crop);
      W = W(1:crop,1:2*crop);
    end
    
    if 0
      ok_rows = sum(W,2) > 2;
      ok_pairs = ok_rows(2:2:end) & ok_rows(1:2:end-1);
      ok_rows(1:2:end)= ok_pairs;
      ok_rows(2:2:end)= ok_pairs;
      M=M(ok_rows,:);
      W=W(ok_rows,:);
    end
    
    [m,n] = size(M);
    
    A0 = randn(m, r);
    A0(1:r,1:r) = eye(r);
    B0 = randn(n, r);
    clear test
    test.M = M;
    test.W = W;
    test.A0 = A0;
    test.B0 = B0;
    test.dataset = dataset;
    test.regularizer_lambda = 0;
    test.gauge_fix_weight =  0;
    test.wiberg_iters = 100;
    
    nu = 1e-5;
    test_init = {
       100    0 'dw'
    %   0 0 'lm'
       0 0 'awflm'
    %   0 0 'lsq'
    %   0 nu 'lm'
    %   nu 0 'lm'
    %   nu nu 'lm'
   %    30 nu 'wl'
   %    9 nu 'wl'
   %    100    0 'dw'
     %  0  0 'damped newton'
     % nu  0 'damped newton'
     %  0  0 'lm'
     % nu nu 'lm'
      };
    
    for k=1:size(test_init,1)
      test.alg = test_init{k,3};
      switch test.alg
        case {'lm', 'lsq', 'awflm', 'scg'}
          test.regularizer_lambda = test_init{k,1};
          test.gauge_fix_weight =  test_init{k,2};
        case {'dw', 'wl'}
          test.wiberg_iters = test_init{k,1};
      end
      
      struct_push_back tests test
    end
  end
end

tests

%%
save_file = 'c:\tmp\mf_tests_aug20.mat';

if exist(save_file, 'file')
  error(['!del ' save_file]);
end
res = awf_mf_lsqnonlin_test(tests, save_file);
