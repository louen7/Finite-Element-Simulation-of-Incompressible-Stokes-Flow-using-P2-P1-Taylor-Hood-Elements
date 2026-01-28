function Mel = matM_elemP1_P1(S1, S2, S3)
% matM_elem : Code fait au TD1
% calcul la matrices de masse elementaire en P1 lagrange
%
% SYNOPSIS Mel = matM_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exact (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par
%          les coordonnees barycentriques
%


% On fixe les sommets du traingle
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l''aire d''un triangle est nulle!!!');
end

% On calcule de la matrice de masse
Mel = zeros(3,3);
coef = abs(D)/24;   % = |T|/12
for i = 1:3
  for j = 1:3
    % diag: 2*coef (=|T|/6), hors diag: 1*coef (=|T|/12)
    Mel(i,j) = coef * (1 + (i==j));
  end
end
