
% results from awf_mf_lsqnonlin_test
%load c:\tmp\mf_tests2
load c:\tmp\mf_tests_jul16

% restrict to 'dino'
dataset = cellfun(@(x) x{1}, {allres.dataset}, 'uni', 0);
allres = allres(strcmp('dino', dataset));

%%
alg_is_lm = strcmp('lm', {allres.alg});
alg_is_dw = strcmp('dw', {allres.alg});
alg_is_wl = strcmp('wl', {allres.alg});

L = [allres.regularizer_lambda];
G = [allres.gauge_fix_weight];
WC = [allres.wiberg_iters];

reg00 = (L == 0) & (G == 0);

selectors = {
  'lm', alg_is_lm & reg00, 'kx'
  'lm_reg', alg_is_lm & (L ~= 0), 'ko'
  'lm_gauge', alg_is_lm & (G ~= 0) & (L == 0), 'k.'
  'dw', alg_is_dw & (WC > 9), 'rx'
  'dw9', alg_is_dw & (WC <= 9), 'r+'
  'wl', alg_is_wl & (WC ~= 9), 'bo'
  'wl9', alg_is_wl & (WC == 9), 'b.'
  };

minrms = min([allres.rms]);
clf
hold on

%%
for subset_index = 1:size(selectors, 1)
  tag = selectors{subset_index,1};
  mask = selectors{subset_index,2};
  style = selectors{subset_index,3};
  %tag = selectors{subset_index,1};
  dtag = strrep(tag, '_', '\_');
 
  results = allres(mask);

  n = length(results);
  if n > 0
    n_min = sum([results.rms] - minrms < 1e-6);
    
    h = plot([results.rms], [results.time], style);
    set(h, 'displayname', sprintf('%2d/%-2d %s', n_min, n, dtag));
  end
end

ylabel('time(sec)')
xlabel('rms')
title(sprintf('minrms = %.6f', minrms))

%axis([1.0 5 0 200]);

legend('Location', 'east')

vline(min([allres.rms]))
