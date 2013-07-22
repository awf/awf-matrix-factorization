function logout = awf_levmarq_test

% Constant in Rosenbrock function
K = 10;

% Starting point, around the corner from global min
x=[-1.9, 2];

%% Set go = 0 to do without the graphical display
GX.go = 1;
if GX.go
  %% plot function
  range = linspace(-4,4,60);
  [xx,yy] = meshgrid(range);
  zz = xx;
  for i=1:numel(xx)
    e = rosenbrock([xx(i) yy(i)], K, struct('go',0));
    zz(i) = sum(e.^2);
  end
  hold off
  contour(xx,yy,zz,2.^[1 4 6:18])
  hold on
  axis xy
  drawnow
  % GX.ax = get(h,'parent');
end

% define optimization function
f = @(x) rosenbrock(x, K, GX);

% check it at some values
%[e,J] = f(x);
%f([1,1]);

% call LM
opts = awf_levmarq('opts');
opts.Display = 'final';
opts.USE_LINMIN = 1;
[~,~,log] = awf_levmarq(x,f, opts);

%% Define the Rosenbrock function
function [e,J] = rosenbrock(x, K, GX)
if nargout == 1
  e = awf_rosenbrock(x, K);
else
  [e,J] = awf_rosenbrock(x, K);
end

% Graphics
if GX.go
  h = plot(x(1),x(2), 'b.');
  if nargout > 1
    set(h,'marker', 'o', 'markersize', 10, 'color', 'r');
  end
  drawnow
end
