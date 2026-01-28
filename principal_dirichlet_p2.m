clear; clc; close all;

% Convergence EF P2 Lagrange pour :
%   -Delta u + u = f  dans Ω
%   u = g             sur ∂Ω  (Dirichlet non homogène)
%
% Solution exacte choisie :
%   u_ex(x,y) = 1 + 3*sin(pi*x)*sin(2*pi*y)
%
% f1(x,y) vérifie :  f = -Delta u_ex + u_ex
% g(x,y)  = u_ex(x,y) sur le bord
%
% Objectif :
%   - Nous avons étudié la convergence en normes L2 et H1 en fonction de h
%   - Nous avons tracé log(1/h) -> log(erreur relative)
%   - Nous avons affiché sur les graphes des erreurs le coefficient directeur
%     (ordre de convergence) de la droite de régression

mesh_files = { ...
    'geomRectangle.msh', ...
    'geomRectangle2.msh', ...
    'geomRectangle3.msh', ...
    'geomRectangle4.msh' ...
};

H = [0.2; 0.1; 0.05;0.025] ;   % valeurs de h

nm     = numel(mesh_files); %nombre de fichiers 
errL2r = zeros(nm,1);
errH1r = zeros(nm,1);

% Quadrature de degré 4 (6 points) sur triangle de référence
S_hat = [ ...
    0.0915762135098, 0.0915762135098;
    0.8168475729805, 0.0915762135098;
    0.0915762135098, 0.8168475729805;
    0.1081030181681, 0.4459484909160;
    0.4459484909160, 0.1081030181681;
    0.4459484909160, 0.4459484909160];

poids = [ ...
    0.05497587183;
    0.05497587183;
    0.05497587183;
    0.1116907948;
    0.1116907948;
    0.1116907948];


% BOUCLE SUR LES MAILLAGES

for k = 1:nm

    nom_maillage = mesh_files{k};
    fprintf('\n=== Maillage %s (h = %.3f) ===\n', nom_maillage, H(k));

    %  Lecture du maillage
    [Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = ...
        lecture_msh_ordre2(nom_maillage); 

    % Matrices globales 
    KK = sparse(Nbpt, Nbpt);   % matrice de rigidité
    MM = sparse(Nbpt, Nbpt);   % matrice de masse
    LL = zeros(Nbpt, 1);       % second membre

    % Assemblage élément par élément pour tous les triangles
    for l = 1:Nbtri

        % Nœuds P2 locaux 
        loc = Numtri(l,:);
        S1  = Coorneu(loc(1),:);
        S2  = Coorneu(loc(2),:);
        S3  = Coorneu(loc(3),:);

        % Matrices élémentaires P2
        Kel = matK_elem_p2(S1, S2, S3);
        Mel = matM_elem_p2(S1, S2, S3);

        for i=1:6
            I = Numtri(l,i);
            for j=1:6
                J = Numtri(l,j);
                MM(I,J) = MM(I,J) + Mel(i,j);
                KK(I,J) = KK(I,J) + Kel(i,j);
            end
        end

       
    end 

      x = Coorneu(:,1);
      y = Coorneu(:,2);
      FF = f1(x,y);     % Appel à ta fonction f1.m
      LL = MM * FF;     

    AA = MM + KK;

    %  Pseudo-ÉLIMINATION (Dirichlet non homogène)
    % elimine.m applique u = g(x,y) sur le bord
    [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu); %on peut mettre elimine2 pour tester le cas homogène

    % Résolution
    UU = tilde_AA \ tilde_LL;

    % Calcul des erreurs 
    % Solution exacte prise :
    %U_exact = sin(pix)*sin(piy)
    U_exact = 3 * cos(pi*Coorneu(:,1)) .* cos(2*pi*Coorneu(:,2));
    % Erreurs approchées via M et K
    e   = U_exact - UU;
    eL2 = sqrt( e' * (MM * e) );
    eH1 = sqrt( e' * (KK * e) );

    nL2 = sqrt( U_exact' * (MM * U_exact) );
    nH1 = sqrt( U_exact' * (KK * U_exact) );

    errL2r(k) = eL2 / nL2;
    errH1r(k) = eH1 / nH1;

    fprintf('Erreur relative L2 : %.3e\n', errL2r(k));
    fprintf('Erreur relative H1 : %.3e\n', errH1r(k));
    affiche_ordre2(UU,Numtri, Coorneu,sprintf('Dirichlet - %s',nom_maillage)) ;
end % boucle sur les maillages



% Tracés log-log + régression linéaire

X  = log(1./H);      % abscisse : log(1/h) 
Y1 = log(errL2r);    % ordonnée : log(erreur L2 relative)
Y2 = log(errH1r);    % ordonnée : log(erreur H1 relative)

% Régressions linéaires 
pL2   = polyfit(X, Y1, 1);
pH1   = polyfit(X, Y2, 1);
Y1fit = polyval(pL2, X);
Y2fit = polyval(pH1, X);

% pentes "observées" (ordre de convergence ≈ -pente)
penteL2 = -pL2(1);
penteH1 = -pH1(1);

fprintf('\n===== Pentes observées (ordre de convergence) =====\n');
fprintf('L2 : pente = %.4f  => ordre ≈ %.4f\n', pL2(1), penteL2);
fprintf('H1 : pente = %.4f  => ordre ≈ %.4f\n', pH1(1), penteH1);

% Graphe L2
figure; hold on; grid on; box on;
plot(X, Y1, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
     'DisplayName','log(erreur L^2 relative)');
plot(X, Y1fit, 'k--', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Régression (p = %.3f)', penteL2));
xlabel('log(1/h)');
ylabel('log(||u - u_h||_{L^2} / ||u||_{L^2})');
title('Convergence en norme L^2 (Dirichlet non homogène, P2)');
legend('Location','best');

% Graphe H1
figure; hold on; grid on; box on;
plot(X, Y2, 's-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
     'DisplayName','log(erreur H^1 relative)');
plot(X, Y2fit, 'k--', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Régression (p = %.3f)', penteH1));
xlabel('log(1/h)');
ylabel('log(|u - u_h|_{H^1} / |u|_{H^1})');
title('Convergence en semi-norme H^1 (Dirichlet non homogène, P2)');
legend('Location','best');
