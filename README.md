Matrix Factorization
====================

A simple implementation of low-rank matrix factorization using 
MATLAB's built-in Levenberg Marquard (or the simple LM from awful.codeplex.com).

*Usage*: 

%1. Obtain M, W, e.g.

        [M,W] = load_giraffe;

%2. Choose initial estimates A0, B0

        A0 = rand(size(M,1), 6);
        B0 = rand(size(M,2), 6);

%3. Call 

        opts = awf_mf_lsqnonlin('opts');
        opts.lsopts.Display = 'iter';
        [A,B] = awf_mf_lsqnonlin(W,M,A0,B0,opts);
        tic
        opts.lsopts.Display = 'none';
        [A,B] = awf_mf_lsqnonlin(W,M,A0,B0,opts);
        toc % will take about 1 minute 
		
%4. or, with prior and/or Gauge fixing

		opts.regularizer_lambda = 1e-5;
    	opts.gauge_fix_weight = 1e-5;
		[A,B] = awf_mf_lsqnonlin(W,M,A0,B0,opts);

5. or, with au_levmarq

        opts = awf_mf_lsqnonlin('opts');
		opts.Algorithm = 'awf';
        tic
		[A,B] = awf_mf_lsqnonlin(W,M,A0,B0,opts);
        toc


Sources
=======
damped_wiberg.zip from [http://www.vision.is.tohoku.ac.jp/us/download]
md5 hash 6E72535A32F06045C9C7EC626A38D543

Data from [http://www.robots.ox.ac.uk/~amb]

