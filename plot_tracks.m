function h=plot_tracks(M, W, varargin)
% PLOT_TRACKS  Plot point trajectories from measurement and weight matrix
%         PLOT_TRACKS(M, W, ...)

x = M(1:2:end-1,:);
y = M(2:2:end,:);
x(~W(1:2:end,:)) = nan;
hh=plot(x,y,varargin{:});
axis([0 640 0 480])
axis ij

if nargout > 0
  h=hh;
end
