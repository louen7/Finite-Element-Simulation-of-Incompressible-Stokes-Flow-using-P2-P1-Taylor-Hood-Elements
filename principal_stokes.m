% =====================================================
% Projet EF Stokes 2D - Écoulement de Poiseuille
%
% Problème de Stokes stationnaire dans Omega = ]0,3[ x ]0,2[
%
%   - nu * Delta u + Grad p = f    dans Omega
%                    div u  = 0    dans Omega
%                          u = g   sur Gamma_D
%            nu * du/dn - p n = 0  sur Gamma_N
%
% Solution exacte de Poiseuille :
%   u1(x,y) = (2 - y) * y
%   u2(x,y) = 0
%   p(x,y)  = -2 * (x - 3)
%
% Conditions aux limites :
%   - y = 0 et y = 2  : u = (0,0) (Dirichlet homogène, parois)
%   - x = 0          : u = ((2 - y)*y, 0) (profil imposé à l'entrée)
%   - x = 3          : Neumann homogène (sortie libre)
%
clear; close all; clc;

% Paramètres et liste des maillages

nu = 1;   % viscosité

mesh_files = { ...
    'geomRectangle_partie3.msh', ...
    'geomRectangle_partie3-2.msh', ...
    'geomRectangle_partie3-3.msh', ...
    'geomRectangle_partie3-4.msh' ... 
};

H = [0.2; 0.1; 0.05; 0.0375]; %valeurs de h

nm = numel(mesh_files); %nombre de fichiers mesh_files

% tableaux pour les erreurs 
errL2u  = zeros(nm,1);
errH1u  = zeros(nm,1);




% Débit exact de l'écoulement de Poiseuille :
% Q_exact 4/3
Q_exact  = 4/3;
Qh       = zeros(nm,1);  % débits numériques
errQ_rel = zeros(nm,1);  % erreurs  sur Q

x0   = 0.0;      % section choisie : x = 0 (entrée)
tolerance = 1e-10;    % on fixe une tolérance pour x proche de x0 pour la quadrature des trapezes pour le debit

% on boucle sur les quatre maillages
for k = 1:nm

    nom_maillage = mesh_files{k};
    fprintf('\n=== Maillage %s (h = %.4f) ===\n', nom_maillage, H(k));

    % Lecture du maillage d'ordre 2 (Gmsh)
    
    [Nbpt, Nbtri, Coorneu, Refneu, Numtri, ~, ~, ~, ~] = ...
        lecture_msh_ordre2(nom_maillage);

    
    % Assemblage des matrices EF pour Stokes
    % (Taylor-Hood P2-P1)
    MM        = sparse(Nbpt, Nbpt);  % masse 
    KK        = sparse(Nbpt, Nbpt);  % rigidité 
    EEtmp     = zeros(Nbpt, Nbpt);   % bloc (p, v1) 
    FFtmp     = zeros(Nbpt, Nbpt);   % bloc (p, v2)
    LL = zeros(Nbpt, 1);      % repère les sommets (P1)

    for l = 1:Nbtri

        % coordonnées des 3 sommets P1 du triangle
        S1 = Coorneu(Numtri(l,1),:);
        S2 = Coorneu(Numtri(l,2),:);
        S3 = Coorneu(Numtri(l,3),:);

        % on marque les sommets P1
        LL(Numtri(l,1)) = 1;
        LL(Numtri(l,2)) = 1;
        LL(Numtri(l,3)) = 1;

        % matrices élémentaires P2/P1
        Mel = matM_elem_p2(S1, S2, S3);   % masse 
        Kel = matK_elem_p2(S1, S2, S3);   % rigidité 
        Eel = matE_elem(S1, S2, S3);      % (p, v1) 
        Fel = matF_elem(S1, S2, S3);      % (p, v2) 

        % indices P2 (6 nœuds) et P1 (3 sommets)
        indP2 = Numtri(l,:);       % 6 nœuds P2 pour la vitesse
        indP1 = Numtri(l,1:3);     % 3 sommets P1 pour la pression

        % assemblage globaux vitesse
        MM(indP2, indP2) = MM(indP2, indP2) + Mel;
        KK(indP2, indP2) = KK(indP2, indP2) + Kel;

        % assemblage blocs pression-vitesse (temporaires)
        EEtmp(indP2, indP1) = EEtmp(indP2, indP1) + Eel;
        FFtmp(indP2, indP1) = FFtmp(indP2, indP1) + Fel;
    end
    
    % in extrait tous sommets P1

    numsommets = find(LL == 1);   % indices globaux des sommets
    Ns         = length(numsommets);     % nombre de sommets

    EE = EEtmp(:, numsommets);           % bloc (p, v1) 
    FF = FFtmp(:, numsommets);           % bloc (p, v2) 

    
 
    AA = sparse(2*Nbpt + Ns, 2*Nbpt + Ns);

    % partie de la matrice vitesse-vitesse
    AA(1:Nbpt, 1:Nbpt) = nu * KK;                    % (u1,u1)
    AA(Nbpt+1:2*Nbpt, Nbpt+1:2*Nbpt) = nu * KK;      % (u2,u2)

    % partie de la matrice vitesse-pression
    AA(1:Nbpt,           2*Nbpt+1:2*Nbpt+Ns) = EE;   % (u1,p)
    AA(Nbpt+1:2*Nbpt,    2*Nbpt+1:2*Nbpt+Ns) = FF;   % (u2,p)

    % partie de la matrice pression-vitesse
    AA(2*Nbpt+1:2*Nbpt+Ns, 1:Nbpt)         = EE';
    AA(2*Nbpt+1:2*Nbpt+Ns, Nbpt+1:2*Nbpt)  = FF';

    % second membre global (ici f = 0)
    LL = zeros(2*Nbpt + Ns, 1);

    % Pseudo-élimination :
    [tilde_AA, tilde_LL] = elimine_stokes(AA, LL, Refneu, Coorneu, numsommets);

    % Résolution du système linéaire
    UU = tilde_AA \ tilde_LL;

    % séparation des composantes
    U1 = UU(1:Nbpt);
    U2 = UU(Nbpt+1 : 2*Nbpt);
    P  = UU(2*Nbpt+1 : 2*Nbpt+Ns);

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

    errL2u(k) = eL2 / nL2;
    errH1u(k) = eH1 / nH1;

    fprintf('Erreur relative L2(u) : %.3e\n', errL2u(k));
    fprintf('Erreur relative H1(u) : %.3e\n', errH1u(k));

   
    % Débit approximatif Q_h sur la section x = x0
    %   (nous avons utilisé de la quadrature 1D par trapèzes sur les nœuds x ≈ x0 avec la fonction trapz)
    
    x_entier = Coorneu(:,1);
    y_entier = Coorneu(:,2);

    ind_x0= find(abs(x_entier - x0) < tolerance);   % nœuds  x = x0

    y_x0  = y_entier(ind_x0);
    u1_x0 = U1(ind_x0);

    % on trie les points de la section par y croissant
    [y_x0_sorted, perm] = sort(y_x0);
    u1_x0_sorted        = u1_x0(perm);

    % fonction trapz pour la quadrature pour l'approximation du débit
    Qh(k) = trapz(y_x0_sorted, u1_x0_sorted);
    errQ_rel(k) = abs(Qh(k) - Q_exact) / abs(Q_exact);

    fprintf('Débit : Q_exact = %.10f, Q_h = %.10f, erreur relative Q = %.3e\n', ...
            Q_exact, Qh(k), errQ_rel(k));


    % Visualisation (uniquement pour le maillage h=0.1)
    if k == 2
        % u1 et u2 en champs scalaires
        affiche_ordre2(U1, Numtri, Coorneu, sprintf('U1 - %s', nom_maillage));
        affiche_ordre2(U2, Numtri, Coorneu, sprintf('U2 - %s', nom_maillage));

        % pression P (sur les sommets seulement)
        SOMMETS = zeros(Nbpt,1);
        SOMMETS(numsommets) = 1:Ns;
        Numtrichangement = SOMMETS(Numtri(:,1:3));

        x_sommets = Coorneu(numsommets,1);
        P_exact   = -2 * (x_sommets - 3);

        affiche_ordre1(P, Numtrichangement, Coorneu(numsommets,:), ...
            sprintf('Pression approximée - %s', nom_maillage));
        % Erreur sur u1 : u1_exact - u1_h
        errU1 = U1_exact - U1;   % vecteur de taille Nbpt
        affiche_ordre2(errU1, Numtri, Coorneu, ...
            sprintf('Erreur u1 (u1_{exact} - u1_h) - %s', nom_maillage));

        % Erreur sur u2 : u2_exact - u2_h
        errU2 = U2_exact - U2;   % ici U2_exact = 0
        affiche_ordre2(errU2, Numtri, Coorneu, ...
            sprintf('Erreur u2 (u2_{exact} - u2_h) - %s', nom_maillage));

        % Erreur sur p : p_exact - p_h (sur les sommets seulement)
        errP = P_exact - P;      % taille Ns
        affiche_ordre1(errP, Numtrichangement, Coorneu(numsommets,:), ...
            sprintf('Erreur pression (p_{exact} - p_h) - %s', nom_maillage));
       
        % champ de vitesse (quiver) sur les sommets
        Ux_som = U1(numsommets);
        Uy_som = U2(numsommets);
        X_som  = Coorneu(numsommets,1);
        Y_som  = Coorneu(numsommets,2);

        figure;
        quiver(X_som, Y_som, Ux_som, Uy_som, ...
               'AutoScale','on','AutoScaleFactor',0.8);
        axis equal;
        xlabel('x'); ylabel('y');
        title(sprintf('Champ de vitesse (approx, sommets) - %s', nom_maillage));
        grid on;
    end

end 


% Récapitulatif des résultats numériques
fprintf('\n================= RÉCAPITULATIF =================\n');
fprintf('Maillage\t h\t  errL2(u)\t  errH1(u)\t  errQ_rel\n');
for k = 1:nm
    fprintf('%d\t\t %.4f\t %.3e\t %.3e\t %.3e\n', ...
        k, H(k), errL2u(k), errH1u(k), errQ_rel(k));
end
fprintf('=================================================\n');

