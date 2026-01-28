function [Mel] = matM_elem_p2(S1, S2, S3)

% Matrice de masse élémentaire P2 (triangle 2D) par quadrature d'ordre 4
% Entrée : S1,S2,S3 = [x y] des 3 sommets
% Sortie : Mel (6x6)


% On récupère les coordonnées des sommets du triangle réel
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% On applique la transformation : F(M̂) = S1 + B * M̂
B = [x2-x1, x3-x1; ...
     y2-y1, y3-y1];


detB = abs(det(B));

% Points de quadrature d'ordre 4 (6 points) sur le triangle de référence T̂

S_hat = [0.0915762135098, 0.0915762135098;  ...
         0.8168475729805, 0.0915762135098;  ...
         0.0915762135098, 0.8168475729805;  ...
         0.1081030181681, 0.4459484909160;  ...
         0.4459484909160, 0.1081030181681;  ...
         0.4459484909160, 0.4459484909160];

% Poids associés à chacun des 6 points de quadrature
poids = [0.05497587183, 0.05497587183, 0.05497587183, ...
         0.1116907948,  0.1116907948,  0.1116907948];

% On initialise la matrice de masse (6 noeuds P2)
Mel = zeros(6,6);

% Boucle sur les points de quadrature
for q = 1:6
    % Coordonnées du point de quadrature q dans le triangle de référence
    xi  = S_hat(q,1); 
    eta = S_hat(q,2);

    % Coordonnées barycentriques (L1, L2, L3) associées à (xi, eta)
    % On a L1 + L2 + L3 = 1
    L1 = 1 - xi - eta;  
    L2 = xi;  
    L3 = eta;

    % Valeurs des 6 fonctions de forme P2 au point de quadrature
    N = zeros(6,1);
    N(1) = L1*(2*L1 - 1);
    N(2) = L2*(2*L2 - 1);
    N(3) = L3*(2*L3 - 1);
    N(4) = 4*L1*L2;
    N(5) = 4*L2*L3;
    N(6) = 4*L3*L1;

    % On applique la formule de quadrature totale
    Mel = Mel + poids(q) * (N * N.') * detB;
end

end

