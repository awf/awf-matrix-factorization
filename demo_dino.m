function demo_dino

p = load('demo_dino_data');
p = p.p;

opts = awf_mf_lsqnonlin('opts');
opts.awopts.Display = 'iter';
[m,n] = size(p.M);
r = size(p.A0,2);
scale = std(p.M(p.M~=0));
rms = @(Ap,Bp) sqrt(sum(vec(p.W.*(p.M-Ap*Bp')).^2)/nnz(p.W));

opts.awopts.IterStartFcn = @plotfn;

[A,B] = awf_mf_lsqnonlin(p.W, p.M, p.A0, p.B0, opts);
%plot_tracks(A*B', true);

  function x = plotfn(x)
    A = reshape(x(1:m*r), r,m)';
    B = reshape(x((1:n*r)+m*r), r,n)';
    

    plot_tracks(A*B'*scale, true);
    axis([-100 700 -100 600])
    axis ij
    rat = rms(A*scale,B)/1.084673;
    title(sprintf('rms=%f ratio=%f', rms(A*scale,B), rat));
    drawnow
    if rat <1.01
      title(sprintf('rms=%f ratio=%f [WAITING]', rms(A*scale,B), rat));
      
      waitforbuttonpress
    end
    
  end
end
