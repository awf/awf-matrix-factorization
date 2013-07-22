Matrix Factorization
====================

A simple implementation of low-rank matrix factorization using MATLAB's built-in Levenberg Marquardt.

*Usage*: 

1. Obtain M, W

2. Choose initial estimates A0, B0

3. Call 

		[A,B] = awf_mfwmd_lsqnonlin(W,M,A0,B0);

4. or, with prior and/or Gauge fixing

		opts.regularizer_lambda = 1e-5;
    	opts.gauge_fix_weight = 1e-5;
		[A,B] = awf_mfwmd_lsqnonlin(W,M,A0,B0);


Sources
=======
damped_wiberg.zip from [[http://www.vision.is.tohoku.ac.jp/us/download]]
md5 hash 6E72535A32F06045C9C7EC626A38D543

Data from [[http://www.robots.ox.ac.uk/~amb]]

