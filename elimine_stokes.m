function [A_mod, L_mod] = elimine_stokes(A, L, Refneu, Coorneu, numsommets)
% Pseudo-elimination pour le probleme de Stokes :
% - on impose u = g sur Gamma_D (Dirichlet) pour la vitesse
% - on fixe p = 0 en un sommet de pression sur x = 3
%
% Entrées :
%   A         : matrice globale du système Stokes
%   L         : second membre global
%   Refneu    : références des nœuds (non utilisé ici mais en général utile)
%   Coorneu   : coordonnées des nœuds (Nbpt x 2)
%   numsommets: indices globaux des sommets (pour la pression P1)
%
% Sorties :
%   A_mod : matrice modifiée qui tient compte des CL
%   L_mod : second membre modifié

A_mod = A;   % on part de A et L, puis on les modifie
L_mod = L;

Nbpt = size(Coorneu,1); % on calcule le nombre total de nœuds P2 (pour u1,u2)
Ntot = size(A,1);       % taille totale du système (u1,u2,p)

% On récupère les coordonnées x,y de tous les nœuds
x = Coorneu(:,1);
y = Coorneu(:,2);
tol = 1e-8;  % on fixe une tolérance pour étudier les comparaisons sur les bords


% On repère les nœuds Dirichlet pour la vitesse avec les conditions aux
% limites
neoud_derriere = abs(y - 0) < tol;   % bord du bas (y=0)
neoud_dessus   = abs(y - 2) < tol;   % bord du haut (y=2)
neoud_gauche   = abs(x - 0) < tol;   % bord gauche (x=0)

% On regroupe tous ces nœuds dans un seul vecteur
noeuds_Dirichlet = find(neoud_derriere | neoud_dessus | neoud_gauche);


% On calcule la valeur de la condition g = (g1,g2) sur les nœuds de Dirichlet

xD = Coorneu(noeuds_Dirichlet,1); % coordonnées x des nœuds Dirichlet
yD = Coorneu(noeuds_Dirichlet,2); % coordonnées y des nœuds Dirichlet

% g1(x,y) et g2(x,y) 
g1_vals = g1(xD, yD); 
g2_vals = g2(xD, yD);  

% On construit les indices associés aux nœuds de bord

indice_Dirichlet = zeros(2*length(noeuds_Dirichlet),1); % tous les indices de vitesse sont imposés
g_vec  = zeros(Ntot,1);               % vecteur qui contient g (u1,u2, et 0 pour p)

for k = 1:length(noeuds_Dirichlet)
    n = noeuds_Dirichlet(k);

    id1 = n;         % on fixe l'indice global pour u1 au nœud n
    id2 = Nbpt + n;  % on fixe l'indice global pour u2 au nœud n

    % On stocke ces indices dans indice_Dirichlet
    indice_Dirichlet(2*k-1) = id1;
    indice_Dirichlet(2*k)   = id2;

    % On met les valeurs de g1,g2 dans le vecteur g_vec
    g_vec(id1) = g1_vals(k);
    g_vec(id2) = g2_vals(k);
end


% On corrige le second membre pour les points intérieurs (pour les points
% de bord)



I = setdiff(1:Ntot, indice_Dirichlet);  

L_mod(I) = L_mod(I) - A_mod(I, indice_Dirichlet) * g_vec(indice_Dirichlet);

% On annule les lignes et les colonnes des indices de Dirichlet
A_mod(indice_Dirichlet,:) = 0;
A_mod(:,indice_Dirichlet) = 0;

% On met diag = 1 pour ces les indices de dirichlet
A_mod(sub2ind([Ntot Ntot], indice_Dirichlet, indice_Dirichlet)) = 1;

% On impose la valeur U = g sur ces indices de dirichlet
L_mod(indice_Dirichlet) = g_vec(indice_Dirichlet);


% On Fixe la constante de pression
%    Pour avoir une solution unique, on fixe p = 0 en un point. on choisit
%    x = 3

x_sommets = Coorneu(numsommets,1);               % x des sommets de pression
indice_droit = find(abs(x_sommets - 3) < tol, 1);   % on prend le premier sommet avec x ≈ 3
ip = 2*Nbpt + indice_droit;  

% On remplace la ligne et la colonne par 0, sauf un 1 sur la diagonale
A_mod(ip,:) = 0;
A_mod(:,ip) = 0;
A_mod(ip,ip) = 1;

% On impose p = 0
L_mod(ip) = 0; 

end
