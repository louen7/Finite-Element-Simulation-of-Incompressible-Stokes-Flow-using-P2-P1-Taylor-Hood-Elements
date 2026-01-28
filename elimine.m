function [Ae, Le] = elimine(AA, LL, Refneu, Coorneu)
% ELIMINE : pseudo-élimination des nœuds de Dirichlet
% (éventuellement non homogène)
%
% On veut résoudre un problème du type :
%       AA * U = LL
% avec la condition :
%       U = g sur le bord Dirichlet.
%
% Cette fonction modifie :
%   - la matrice AA  → Ae
%   - le second membre LL → Le
% pour que la condition U = g soit imposée.
%
% Entrées :
%   AA      : matrice globale N×N 
%   LL      : second membre N×1
%   Refneu  : vecteur N×1
%             Refneu(i) = 1 si le nœud i est sur le bord Dirichlet
%   Coorneu : coordonnées des nœuds (N×2),
%
% pour g homogène elimine est plus simple et elle est codée dans elimine2

    % Taille du système (N = nombre total de nœuds)
    N  = size(AA,1);

    % On copie AA et LL dans Ae et Le pour les modifier
    Ae = AA;
    Le = LL;

    
    %On repère les nœuds de bord Dirichlet
  
    b = find(Refneu == 1);  % indices des nœuds sur le bord Dirichlet

    % S'il n'y a pas de nœud de Dirichlet :
    if isempty(b)
        return; 
    end

    % I = indices des nœuds intérieurs (tout sauf les nœuds de bord)
    I = setdiff(1:N, b); 
    % (setdiff renvoie tous les indices de 1 à N qui ne sont pas dans b)

    % On calcule la valeur de g aux nœuds de Dirichlet
    xB = Coorneu(b,1); 
    yB = Coorneu(b,2);  

    % On évalue la condition de Dirichlet avec la fonction g.
    gB = g(xB, yB);  

   
    % On applique le programme de pseudo élimination
   
    Le(I) = Le(I) - Ae(I,b) * gB;

   
    % 4) On met à zéro les lignes et les colonnes correspondant aux nœuds
    %    de bord, puis on impose Ae(b,b) = 1 et Le(b) = gB.
    
    Ae(b, :) = 0;
    Ae(:, b) = 0;

    % On met des 1 sur la diagonale Ae(b,b)
    Ae(sub2ind([N N], b, b)) = 1;

    % On impose la valeur de g
    Le(b) = gB;
end

