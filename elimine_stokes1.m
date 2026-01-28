function [A_mod, L_mod] = elimine_stokes1(A, L, Refneu, Coorneu, numsommets)
% Pseudo-élimination pour le problème de Stokes (canal 0<x<3, 0<y<2)
%On s'inspire de elimine mais on modifie légerement les C.L
%
% CL utilisées ici (d’après la géométrie) :
%   -  x = xmin  -> u1 = (2 - y)*y, u2 = 0      (profil de Poiseuille)
%   -  x = xmax  -> Neumann naturel (on n’impose rien)
%   - Parois (haut/bas) : u = 0   (vitesse nulle)
%
% Entrées :
%   A         : matrice globale du système (taille Nt x N)
%   L         : second membre global (taille Ntot x 1)
%   Refneu    : références des nœuds (N_bpt x 1), >0 si noeud de bord
%   Coorneu   : coordonnées des nœuds (N_bpt x 2)
%   numsommets: indices globaux des sommets (pour la pression P1, pas
%               directement utilisé ici, mais passé pour cohérence)
%
% Sorties :
%   A_mod : matrice A modifiée avec les CL
%   L_mod : second membre modifié avec les CL

    A_mod = A;
    L_mod = L;

    % Nombre de nœuds P2 (pour la vitesse)
    nb_noeuds = size(Coorneu,1);
    % Taille totale du système (u1, u2, p)
    nb_dirichlet = size(A,1);

    % On récupère les coordonnées (x,y) de tous les nœuds
    x_noeuds = Coorneu(:,1);
    y_noeuds = Coorneu(:,2);

    % On calcule xmin et xmax du domaine (sur les nœuds)
    x_min = min(x_noeuds);
    x_max = max(x_noeuds);
    tolerance = 1e-10; %on met une tolerance pour les conditions aux bords

    % Nœuds de bord (Refneu > 0)
    est_bord = (Refneu > 0);

    % On distingue les différents types de bord avec la géométrie du pb
    est_entree = est_bord & (abs(x_noeuds - x_min) < tolerance);   % 
    est_sortie = est_bord & (abs(x_noeuds - x_max) < tolerance);   % 
    est_paroi  = est_bord & ~(est_entree | est_sortie);            % 

    % On impose Dirichlet sur :
    %   - l’entrée (est_entree)
    %   - les parois (est_paroi)
    % On met une condition de Neumann sur la paroie de sortie.
    indices_noeuds_D = find(est_entree | est_paroi);

    % Valeurs de g1,g2 (vitesse imposée) sur ces nœuds
    g1_valeurs = zeros(length(indices_noeuds_D),1);  % composante u1
    g2_valeurs = zeros(length(indices_noeuds_D),1);  % composante u2 (ici 0 partout)

    for k = 1:length(indices_noeuds_D)
        i_noeud = indices_noeuds_D(k);

        if est_entree(i_noeud)
            % Profil de Poiseuille sur l’entrée : u1 = (2 - y)*y
            g1_valeurs(k) = (2 - y_noeuds(i_noeud)) * y_noeuds(i_noeud);
        else
            % Condition aux limites sur les parois u1 = 0
            g1_valeurs(k) = 0;
        end

        %u2 = 0 sur les bords Dirichlet
        g2_valeurs(k) = 0;
    end


    
    indice_dirichlet = zeros(2*length(indices_noeuds_D),1);
    g_global       = zeros(nb_dirichlet,1);  % vecteur qui stocke g (u1,u2,p)

    for k = 1:length(indices_noeuds_D)
        n = indices_noeuds_D(k);  % numéro de nœud P2

       
        id_u1 = n;
        id_u2 = nb_noeuds + n;

        % On stocke ces indices dans la liste des DOF Dirichlet
        indice_dirichlet(2*k-1) = id_u1;
        indice_dirichlet(2*k)   = id_u2;

        % On stocke aussi les valeurs de la condition g dans g_global
        g_global(id_u1) = g1_valeurs(k);
        g_global(id_u2) = g2_valeurs(k);
    end

   
    indices_interieurs = setdiff(1:nb_dirichlet, indice_dirichlet);

   
    L_mod(indices_interieurs) = L_mod(indices_interieurs) ...
                              - A_mod(indices_interieurs, indice_dirichlet) * g_global(indice_dirichlet);

 
    A_mod(indice_dirichlet,:) = 0;
    A_mod(:,indice_dirichlet) = 0;

    A_mod(sub2ind([nb_dirichlet nb_dirichlet], indice_dirichlet, indice_dirichlet)) = 1;
    L_mod(indice_dirichlet) = g_global(indice_dirichlet);

    
    indice_p_fixe = 2*nb_noeuds + 1;

    A_mod(indice_p_fixe,:)  = 0;
    A_mod(:,indice_p_fixe)  = 0;
    A_mod(indice_p_fixe,indice_p_fixe) = 1;
    L_mod(indice_p_fixe)    = 0;   % on impose p = 0 en ce point

end

