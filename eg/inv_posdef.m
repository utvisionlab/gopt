function Minv = inv_posdef(M)
R = chol(M);
Rinv = R \ eye(size(R,1)); % Faster version of inv_triu(R); 
Minv = (Rinv * Rinv.');
