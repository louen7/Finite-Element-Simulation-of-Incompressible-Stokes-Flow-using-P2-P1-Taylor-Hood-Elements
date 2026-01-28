function [Kel] = matk_elemP1_P1(S1, S2, S3)

% mat_elem P1-P1
% calcul la matrice en P1 lagrange fait au TD1
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2]; % associee a S1 (arete S2-S3)
norm(2, :) = [y3-y1, x1-x3]; % associee a S2 (arete S3-S1)
norm(3, :) = [y1-y2, x2-x1]; % associee a S3 (arete S1-S2)

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l''aire d''un triangle est nulle !!!');
end

% calcul de la matrice de rigidité
Kel = zeros(3,3);
for i = 1:3
  for j = 1:3
    % K_ij = (n_i · n_j) / (2*|D|)
    Kel(i,j) = dot(norm(i,:), norm(j,:)) / (2*abs(D));
  end
end
