function [M,W, Mx, My] = load_dino
% LOAD_DINO  Load data for "rotating dinosaur" dataset
%  [M,W, Mx, My] = load_dino


M = load('data/rotdino.M_inliers.txt');

W = M~=0;

if nargout > 2
  Mx = M(1:2:end,:);
  My = M(2:2:end,:);
  Mx(W(1:2:end,:)==0) = nan;
end
