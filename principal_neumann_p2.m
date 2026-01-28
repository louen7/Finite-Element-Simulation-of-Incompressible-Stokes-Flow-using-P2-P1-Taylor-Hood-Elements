clear; clc; close all;

% =====================================================
% Convergence EF P2 pour le problème :
%   -Delta u + u = f  dans Ω
%   ∂u/∂n = g         sur ∂Ω   (Neumann)
%
% Solution exacte :
%   u(x,y) = 3 cos(pi*x) cos(2*pi*y)
%
% Objectifs :
%   - Nous avons étudié l’erreur en fonction du pas de maille h
%   - Nous avons établi les ordres de convergence (L2 et H1)
%   - Nous avons tracé des coupes 1D de la solution
%   - Nous avons affiché sur les graphes des erreurs la pente de la droite de régression




% Définition des fichiers de maillage
mesh_files = { ...
    'geomrectangle.msh', ...
    'geomrectangle2.msh', ...
    'geomrectangle3.msh', ...
    'geomrectangle4.msh', ...
};

H = [0.2; 0.1;0.05;0.025];   % valeurs de h

nm     = numel(mesh_files); % nombre de fichiers de maillage
errL2r = zeros(nm,1);
errH1r = zeros(nm,1);

% Pour les coupes
Xligne        = cell(nm,1);
Uligne        = cell(nm,1);
y_mid_global = [];

% Quadrature du second membre (ordre 4)
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

    % Lecture du maillage P2
    [Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = ...
        lecture_msh_ordre2(nom_maillage); 

    % Initialisation des matrices globales
    KK = sparse(Nbpt, Nbpt); % matrice de rigidité
    MM = sparse(Nbpt, Nbpt); % matrice de masse
    LL = zeros(Nbpt, 1);     % second membre

    % Assemblage triangle par triangle
    for l = 1:Nbtri

        loc = Numtri(l,:);    % les 6 nœuds locaux
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
    end % fin triangles
    % Second membre 
     x = Coorneu(:,1);
     y = Coorneu(:,2);
     FF = f(x,y);
     LL = MM * FF;

    % Résolution du système linéaire
    UU = (MM + KK) \ LL;

    % Solution exacte aux nœuds
    U_exact = 3*cos(pi*Coorneu(:,1)) .* cos(2*pi*Coorneu(:,2));

    % erreurs
    e   = U_exact - UU;
    eL2 = sqrt( e' * (MM * e) );
    eH1 = sqrt( e' * (KK * e) );

    nL2 = sqrt( U_exact' * (MM * U_exact) );
    nH1 = sqrt( U_exact' * (KK * U_exact) );

    errL2r(k) = eL2 / nL2;
    errH1r(k) = eH1 / nH1;

    fprintf('Erreur relative L2 : %.3e\n', errL2r(k));
    fprintf('Erreur relative H1 : %.3e\n', errH1r(k));
    affiche_ordre2(UU,Numtri, Coorneu,sprintf('Neumann - %s',nom_maillage)) ;

  
    % Coupe 1D à y = y_mid
    
    x_all = Coorneu(:,1);
    y_all = Coorneu(:,2);

    y_min = min(y_all);
    y_max = max(y_all);
    y_mid = 0.5*(y_min + y_max); %on prend la médiane pour que cela soit le plus représentatif possible

    if k == 1
        y_mid_global = y_mid;
    end

    ind_line = find(abs(y_all - y_mid) < 1e-6);

    x_line = x_all(ind_line);
    u_line = UU(ind_line);

    [x_line_sorted, idx] = sort(x_line);
    u_line_sorted        = u_line(idx);

    Xligne{k} = x_line_sorted;
    Uligne{k} = u_line_sorted;

end 




% Graphes convergence log-log

X  = log(1./H);
Y1 = log(errL2r);
Y2 = log(errH1r);

pL2   = polyfit(X, Y1, 1);
pH1   = polyfit(X, Y2, 1);
Y1fit = polyval(pL2, X);
Y2fit = polyval(pH1, X);

penteL2 = -pL2(1);
penteH1 = -pH1(1);

fprintf('\n===== Pentes observées =====\n');
fprintf('L2 : %.4f\n', penteL2);
fprintf('H1 : %.4f\n', penteH1);

% Convergence L2
figure; hold on; grid on; box on;
plot(X, Y1, 'o-', 'LineWidth',1.5);
plot(X, Y1fit, 'k--', 'LineWidth',1.5);
xlabel('log(1/h)'); ylabel('log(err L2)');
title('Convergence L2 (P2)');
legend('Erreur', sprintf('Régression (p = %.3f)', penteL2), ...
       'Location','SouthWest');

% Convergence H1
figure; hold on; grid on; box on;
plot(X, Y2, 's-', 'LineWidth',1.5);
plot(X, Y2fit, 'k--', 'LineWidth',1.5);
xlabel('log(1/h)'); ylabel('log(err H1)');
title('Convergence H1 (P2)');
legend('Erreur', sprintf('Régression (p = %.3f)', penteH1), ...
       'Location','SouthWest');




% Coupes 1D de la solution

x_min = inf; x_max = -inf;
for k = 1:nm
    x_min = min(x_min, min(Xligne{k}));
    x_max = max(x_max, max(Xligne{k}));
end

x_exact       = linspace(x_min, x_max, 400);
u_exact_line  = 3*cos(pi*x_exact) .* cos(2*pi*y_mid_global);

figure; hold on; grid on; box on;
plot(x_exact, u_exact_line, 'k-', 'LineWidth',2);

for k = 1:nm
    plot(Xligne{k}, Uligne{k}, 'o-');
end

xlabel('x');
ylabel(sprintf('u(x, y_{mid}=%.3f)', y_mid_global));
title('Coupes 1D de la solution');
legend('Exacte','u_h h=0.2','u_h h=0.1','u_h h=0.05','u_h h=0.025');
