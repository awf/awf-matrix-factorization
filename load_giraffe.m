function [M,W] = load_giraffe
% LOAD_GIRAFFE Load giraffe data
%              [M,W] = load_giraffe


M = load('data/giraffe.M.txt');
M = M(:,2:end);
W = isfinite(M);
M(~W) = 0;

