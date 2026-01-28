function [A_mod, L_mod] = elimine_stokes1(A, L, Refneu, Coorneu, numsommets)
% Pseudo-élmination pour le problème de Stokes (canal 0<x<8, 0<y<2)
% Notre but est d'adapter le "elimine" trouvé pour stokes pour ce nouveau
% problème avec une autre géométrie

    A_mod = A;
    L_mod = L;

    Nbpt = size(Coorneu,1);
    Ntot = size(A,1);

   
    % 1) Nœuds de Dirichlet (vitesse)
   
    isDir    = (Refneu == 1) | (Refneu == 2);
    noeuds_D = find(isDir);

    % Coordonnées de ces nœuds
    xD = Coorneu(noeuds_D,1);
    yD = Coorneu(noeuds_D,2);

    % g1,g2 sur ces nœuds :
    %   - si Refneu = 1 -> profil -4(y-2)(y-1)
    %   - si Refneu = 2 -> 0
    g1_vals = zeros(length(noeuds_D),1);
    g2_vals = zeros(length(noeuds_D),1);  % u2 = 0 partout

    for k = 1:length(noeuds_D)
        n   = noeuds_D(k);
        ref = Refneu(n);
        yk  = Coorneu(n,2);

        if ref == 1
            % Γ_left : profil d'écoulement entrant
            g1_vals(k) = -4*(yk - 2)*(yk - 1);
        else
            % Γ_down et le reste des parois : u1 = 0
            g1_vals(k) = 0;
        end

        % composante verticale : toujours 0
        g2_vals(k) = 0;
    end

    % On construit la liste des indices de Dirichlet
    indices_Dirichlet = zeros(2*length(noeuds_D),1);
    g_vec  = zeros(Ntot,1);

    for k = 1:length(noeuds_D)
        n = noeuds_D(k);

        id1 = n;           % indice de dirichlet correspondant à u1 en ce nœud
        id2 = Nbpt + n;    % indice de dirichlet correspondant à u2 en ce nœud

        indices_Dirichlet(2*k-1) = id1;
        indices_Dirichlet(2*k)   = id2;

        g_vec(id1) = g1_vals(k);
        g_vec(id2) = g2_vals(k);
    end

    % Indices des neouds intérieurs (tout sauf les indices de Dirichlet)
    I = setdiff(1:Ntot, indices_Dirichlet);

    % Nous corrigeons le second membre pour les équations intérieures
    L_mod(I) = L_mod(I) - A_mod(I, indices_Dirichlet) * g_vec(indices_Dirichlet);

    % Annuler lignes et colonnes des indices de Dirichlet
    A_mod(indices_Dirichlet,:) = 0;
    A_mod(:,indices_Dirichlet) = 0;

    % Nous imposons diag = 1 et RHS = g sur ces indices de bord
    A_mod(sub2ind([Ntot Ntot], indices_Dirichlet, indices_Dirichlet)) = 1;
    L_mod(indices_Dirichlet) = g_vec(indices_Dirichlet);

   
    % 2) On fixe la constante de pression à un point pour avoir l'unicité
  
    % On fixe p = 0 sur un sommet de pression (par ex. le premier)
    %
    % numsommets = indices des sommets P1 dans Coorneu

    ip = 2*Nbpt + 1;  

    A_mod(ip,:)  = 0;
    A_mod(ip,ip) = 1;
    L_mod(ip)    = 0;

end

