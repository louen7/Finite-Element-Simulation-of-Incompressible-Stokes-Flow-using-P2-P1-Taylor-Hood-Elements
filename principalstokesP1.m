% Stokes 2D - Écoulement de Poiseuille
% Schéma P1-P1, test sur UN seul maillage (h = 0.1)

clear; close all; clc;

nu = 1;   % viscosité

% Fichier maillage ordre 1 
nom_maillage = 'geomCarreP1_P1.msh';  % maillage du TP1 d'ordre 1 
h = 0.1;

fprintf('\n=== Maillage %s (h = %.4f) ===\n', nom_maillage, h);




% Lecture du maillage d'ordre 1

[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = ...
    lecture_msh(nom_maillage);


MM = sparse(Nbpt, Nbpt);   % matrice de masse (P1) qui est néanmoins inutile dans ce code
KK = sparse(Nbpt, Nbpt);   % matrice de rigidité (P1)
EE = sparse(Nbpt, Nbpt);   % bloc (u1,p) P1-P1
FF = sparse(Nbpt, Nbpt);   % bloc (u2,p) P1-P1

for l = 1:Nbtri

    % coordonnées des 3 sommets P1 du triangle
    S1 = Coorneu(Numtri(l,1),:);
    S2 = Coorneu(Numtri(l,2),:);
    S3 = Coorneu(Numtri(l,3),:);

    % matrices élémentaires P1/P1 programées à la fois dans le TP1 et TP2
    Mel = matM_elemP1_P1(S1, S2, S3);        % masse P1
    Kel = matk_elemP1_P1(S1, S2, S3);        % rigidité P1
    Eel = matE_elemP1_P1(S1, S2, S3);     % (p, v1)
    Fel = matF_elemP1_P1(S1, S2, S3);     % (p, v2)

    % indices locaux P1 (3 noeuds)
    ind = Numtri(l,:);   % 1x3

    % assemblage des MM et KK 
    MM(ind, ind) = MM(ind, ind) + Mel;
    KK(ind, ind) = KK(ind, ind) + Kel;

    % assemblage vitesse-pression
    EE(ind, ind) = EE(ind, ind) + Eel;   % (u1, p)
    FF(ind, ind) = FF(ind, ind) + Fel;   % (u2, p)
end


numsommets = (1:Nbpt).';
Ns         = Nbpt;


AA = sparse(2*Nbpt + Ns, 2*Nbpt + Ns);

AA(1:Nbpt, 1:Nbpt) = nu * KK;                    
AA(Nbpt+1:2*Nbpt, Nbpt+1:2*Nbpt) = nu * KK;      


AA(1:Nbpt,           2*Nbpt+1:2*Nbpt+Ns) = EE;   
AA(Nbpt+1:2*Nbpt,    2*Nbpt+1:2*Nbpt+Ns) = FF;   

AA(2*Nbpt+1:2*Nbpt+Ns, 1:Nbpt)         = EE';
AA(2*Nbpt+1:2*Nbpt+Ns, Nbpt+1:2*Nbpt)  = FF';

% second membre global (f = 0)
LLg = zeros(2*Nbpt + Ns, 1);

[tilde_AA, tilde_LL] = elimine_stokes(AA, LLg, Refneu, Coorneu, numsommets);


% Résolution du système linéaire

UU = tilde_AA \ tilde_LL;

% on sépare les composantes
U1 = UU(1:Nbpt);                % vitesse u1 (P1)
U2 = UU(Nbpt+1 : 2*Nbpt);       % vitesse u2 (P1)
P  = UU(2*Nbpt+1 : 2*Nbpt+Ns);  % pression (P1)

% Calcul des erreurs (vitesse)
x = Coorneu(:,1);
y = Coorneu(:,2);

U1_exact = (2 - y) .* y;
U2_exact = zeros(Nbpt,1);

e1 = U1_exact - U1;
e2 = U2_exact - U2;

% norme L2(u) 
eL2_sq = e1' * (MM * e1) + e2' * (MM * e2);
eL2    = sqrt(eL2_sq);

% norme H1(u) 
eH1_sq = e1' * (KK * e1) + e2' * (KK * e2);
eH1    = sqrt(eH1_sq);

% normes exactes de u pour erreurs relatives
nL2_sq = U1_exact' * (MM * U1_exact);   % U2_exact = 0
nH1_sq = U1_exact' * (KK * U1_exact);
nL2    = sqrt(nL2_sq);
nH1    = sqrt(nH1_sq);

errL2u = eL2 / nL2;
errH1u = eH1 / nH1;

fprintf('Erreur relative L2(u) : %.3e\n', errL2u);
fprintf('Erreur relative H1(u) : %.3e\n', errH1u);

% Visualisations

affiche_ordre1(U1, Numtri, Coorneu, ...
    sprintf('U1 (P1-P1) - %s', nom_maillage));
affiche_ordre1(U2, Numtri, Coorneu, ...
    sprintf('U2 (P1-P1) - %s', nom_maillage));

P_exact = -2 * (x - 3);
affiche_ordre1(P, Numtri, Coorneu, ...
    sprintf('Pression approx (P1-P1) - %s', nom_maillage));

% Erreurs
errU1 = U1_exact - U1;
errU2 = U2_exact - U2;
errP  = P_exact - P;

affiche_ordre1(errU1, Numtri, Coorneu, ...
    sprintf('Erreur u1 (P1-P1) - %s', nom_maillage));
affiche_ordre1(errU2, Numtri, Coorneu, ...
    sprintf('Erreur u2 (P1-P1) - %s', nom_maillage));
affiche_ordre1(errP, Numtri, Coorneu, ...
    sprintf('Erreur pression (P1-P1) - %s', nom_maillage));

% champ de vitesse (quiver) sur tous les noeuds
figure;
quiver(x, y, U1, U2, ...
       'AutoScale','on','AutoScaleFactor',0.8);
axis equal;
xlabel('x'); ylabel('y');
title(sprintf('Champ de vitesse (P1-P1) - %s', nom_maillage));
grid on;
