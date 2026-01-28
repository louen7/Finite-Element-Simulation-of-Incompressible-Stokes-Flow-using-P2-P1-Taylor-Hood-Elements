function [Kel] = matK_elem_p2(S1, S2, S3)

% Matrice de rigidité élémentaire P2 (triangle 2D) par quadrature d'ordre 2
% Entrée : S1,S2,S3 = [x y] des 3 sommets
% Sortie : Kel (6x6)


% On fixe les sommets du triangle réel
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Matrice de transformation : F(M̂) = S1 + B * M̂
B    = [x2-x1, x3-x1; ...
        y2-y1, y3-y1];

detB = abs(det(B));

% Points de quadrature d'ordre 2 sur le triangle de référence T̂
% (3 points de Gauss, chacun avec le même poids)
S_hat = [1/6, 1/6;  ...  
         2/3, 1/6;  ...  
         1/6, 2/3];      

% Poids associés aux 3 points
poids = [1/6, 1/6, 1/6];

% Initialisation de la matrice de rigidité 
Kel = zeros(6,6);


for q = 1:3
    xi  = S_hat(q,1); 
    eta = S_hat(q,2);

    % Coordonnées barycentriques (L1, L2, L3)
    % L1 + L2 + L3 = 1
    L1 = 1 - xi - eta;  
    L2 = xi;  
    L3 = eta;

    % Dérivées des 6 fonctions de forme P2 par rapport à (xi, eta)
    % dN_dxi(i)  = dNi/dxi
    % dN_deta(i) = dNi/deta
    dN_dxi  = zeros(6,1);
    dN_deta = zeros(6,1);

    dN_dxi(1)  = -4*L1 + 1;      
    dN_deta(1) = -4*L1 + 1;

   
    dN_dxi(2)  =  4*L2 - 1;      
    dN_deta(2) =  0;

  
    dN_dxi(3)  =  0;             
    dN_deta(3) =  4*L3 - 1;
    
    dN_dxi(4)  =  4*(L1 - L2);   
    dN_deta(4) = -4*L2;

    dN_dxi(5)  =  4*L3;          
    dN_deta(5) =  4*L2;

    dN_dxi(6)  = -4*L3;          
    dN_deta(6) =  4*(L1 - L3);

    % On regroupe les dérivées en un tableau 2x6 :
    
    Grad_h = [dN_dxi.'; dN_deta.'];

    % Passage au triangle réel :

    Grad = (B.') \ Grad_h;  

    %On utilise la formule de quadrature pour obtenir la matrice
    Kel = Kel + poids(q) * (Grad.' * Grad) * detB;
end

end

