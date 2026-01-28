

% Stokes Taylor-Hood (P2-P1) dans un canal avec cavité ("marche")
%
% Domaine défini dans geomRectangle_partie3-5.geo :
%   - Canal principal : x ∈ [0,3], y ∈ [0,1]
%   - Cavité          : x ∈ [0,1], y ∈ [-0.5, 0]
%
% Éléments finis utilisés :
%   - Vitesse  : P2 (6 nœuds par triangle : 3 sommets + 3 milieux)
%   - Pression : P1 (pression définie aux sommets uniquement)
%
% Conditions aux limites :
%   - Entrée (bord gauche, y>=0) :
%         u1 = (1 - y)*y,  u2 = 0  (profil de Poiseuille)
%   - Cavité + parois haut/bas :
%         u = 0 (vitesse nulle)
%   - Sortie (bord droit) :
%         condition naturelle (Neumann)
%   - Pression :
%         on fixe p = 0 en un point pour lever l'indétermination
%
% Résultats affichés :
%   - u1 (composante horizontale de la vitesse)
%   - u2 (composante verticale de la vitesse)
%   - p  (pression)
%   - on affiche le champ de vitesse (quiver) sur les sommets
%   - on affiche les différentes erreurs u1-u1h, u2-u2h, p-ph
% 

clear; close all; clc;


% Paramètre physique

nu = 1;  % viscosité 

% Lecture du maillage 
nom_maillage = 'geomRectangle_partie3-5.msh'; % le nouveau maillage avec une géometrie différenre
%on lit le maillage
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, ~, ~, ~, ~] = ...
    lecture_msh_ordre2(nom_maillage);

fprintf('Maillage : %s\n', nom_maillage);
fprintf('  Nbpt  = %d\n', Nbpt);
fprintf('  Nbtri = %d\n', Nbtri);


MM     = sparse(Nbpt, Nbpt);
KK     = sparse(Nbpt, Nbpt);
EE     = zeros(Nbpt, Nbpt);
FF     = zeros(Nbpt, Nbpt);
LL = zeros(Nbpt, 1);

% On parcourt tous les triangles
for l = 1:Nbtri

    % Coordonnées des trois sommets (P1) du triangle l
    S1 = Coorneu(Numtri(l,1),:);
    S2 = Coorneu(Numtri(l,2),:);
    S3 = Coorneu(Numtri(l,3),:);

    % On marque les sommets du triangle comme "sommets de pression P1"
    LL(Numtri(l,1)) = 1;
    LL(Numtri(l,2)) = 1;
    LL(Numtri(l,3)) = 1;

 
    Mel = matM_elem_p2(S1, S2, S3);
    Kel = matK_elem_p2(S1, S2, S3);
    Eel = matE_elem(S1, S2, S3);
    Fel = matF_elem(S1, S2, S3);

   
    indP2 = Numtri(l,:);       % 6 nœuds pour la vitesse
    indP1 = Numtri(l,1:3);     % 3 sommets pour la pression

    % Assemblage dans les matrices globales de vitesse
    MM(indP2,indP2) = MM(indP2,indP2) + Mel;
    KK(indP2,indP2) = KK(indP2,indP2) + Kel;

    % Assemblage dans les blocs pression-vitesse que l'on va ensuite
    % modifier
    EE(indP2,indP1) = EE(indP2,indP1) + Eel;
    FF(indP2,indP1) = FF(indP2,indP1) + Fel;
end


% Nous gerons les sommets P1 (pour la pression)

% On récupère les indices globaux des nœuds qui sont des sommets
numsommets = find(LL == 1);
Ns = length(numsommets);   % nombre de sommets (pression P1)

% On extrait seulement les colonnes associées aux sommets
% pour construire les blocs E et F de taille (Nbpt x Ns)
EE = EE(:, numsommets);  % bloc (p, u1)
FF = FF(:, numsommets);  % bloc (p, u2)


% Construction de la matrice par blocs A

AA = sparse(2*Nbpt + Ns, 2*Nbpt + Ns);

% Matrices vitesse-vitesse (u1-u1 et u2-u2)
AA(1:Nbpt, 1:Nbpt) = nu * KK;                    % bloc pour u1
AA(Nbpt+1:2*Nbpt, Nbpt+1:2*Nbpt) = nu * KK;      % bloc pour u2

% Parties de la matrice vitesse-pression (u1-p et u2-p)
AA(1:Nbpt,           2*Nbpt+1:2*Nbpt+Ns) = EE;   % bloc (u1, p)
AA(Nbpt+1:2*Nbpt,    2*Nbpt+1:2*Nbpt+Ns) = FF;   % bloc (u2, p)

% Parties de la matrice pression-vitesse (p-u1 et p-u2), ce sont les transposées (question 3.2) 
AA(2*Nbpt+1:2*Nbpt+Ns, 1:Nbpt)        = EE';
AA(2*Nbpt+1:2*Nbpt+Ns, Nbpt+1:2*Nbpt) = FF';

% Second membre (on prend f = 0)
LL = zeros(2*Nbpt + Ns, 1);

[tilde_AA, tilde_LL] = elimine_stokes1(AA, LL, Refneu, Coorneu, numsommets);


% Résolution du système linéaire
UU = tilde_AA \ tilde_LL;

U1 = UU(1:Nbpt);
U2 = UU(Nbpt+1:2*Nbpt);
P  = UU(2*Nbpt+1 : 2*Nbpt + Ns);


% Affichages des résultats


% Affichage de u1 (P2)
affiche_ordre2(U1, Numtri, Coorneu, ...
    sprintf('u_1 - %s', nom_maillage));

% Affichage de u2 (P2)
affiche_ordre2(U2, Numtri, Coorneu, ...
    sprintf('u_2 - %s', nom_maillage));

% Pression (sur sommets P1 uniquement)
% On construit un tableau SOMMET qui associe à chaque
% nœud P2 un numéro de sommet P1 ou 0 si ce n'est pas un sommet.
SOMMET = zeros(Nbpt, 1);
SOMMET(numsommets) = 1:Ns;

% Numtri(:,1:3) contient les 3 sommets P1 de chaque triangle
% On remplace les numéros de nœuds P2 par les numéros de sommets P1
Numtrichangement = SOMMET(Numtri(:,1:3));

affiche_ordre1(P, Numtrichangement, Coorneu(numsommets,:), ...
    sprintf('pression - %s', nom_maillage));

% Champ de vitesse (quiver) sur les sommets pour plus de lisibilité
Ux_som = U1(numsommets);            % u1 aux sommets
Uy_som = U2(numsommets);            % u2 aux sommets
X_som  = Coorneu(numsommets,1);     % x des sommets
Y_som  = Coorneu(numsommets,2);     % y des sommets

figure;
quiver(X_som, Y_som, Ux_som, Uy_som, ...
       'AutoScale','on','AutoScaleFactor',0.8);
axis equal;
xlabel('x'); ylabel('y');
title(sprintf('Champ de vitesse (sommets) - %s', nom_maillage));
grid on;


 