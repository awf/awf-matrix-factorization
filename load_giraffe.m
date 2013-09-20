function [M,W] = load_giraffe(small)
% LOAD_GIRAFFE Load giraffe data
%              [M,W] = load_giraffe

M = load('data/giraffe.M.txt');
M = M(:,2:end);
W = isfinite(M);
M(~W) = 0;

if nargin >= 1 && any(small)
    M=M(1:40,1:40);
    W=W(1:40,1:40);
end
