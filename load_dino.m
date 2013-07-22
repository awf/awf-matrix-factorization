M = load('data/rotdino.M_inliers.txt');

W = M~=0;

Mx = M(1:2:end,:); 
My = M(2:2:end,:); 
Mx(W(1:2:end,:)==0) = nan; 


