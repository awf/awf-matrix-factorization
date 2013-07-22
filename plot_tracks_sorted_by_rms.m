
[~,i] = sort([allres.rms]);
for k=1:length(i)
  r = allres(i(k)); 
  hold off; 
  plot_tracks(r.A*r.B', true); 
  hold on; 
  plot_tracks(r.M,r.W, 'k.'); 
  title(sprintf('k=%d, rms=%.6f, alg=%s',k,r.rms,r.alg));
  waitforbuttonpress; 
end
