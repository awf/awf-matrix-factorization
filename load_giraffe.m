%script load_giraffe

M = load('data/giraffe.M.txt');
M = M(:,2:end);
W = isfinite(M);
M(~W) = 0;

