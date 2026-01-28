function val = f2(x, y)
% F2  : second membre f(x,y) = (1 + 2*pi^2) * sin(pi*x)*sin(pi*y)
%
% Entrées :
%     x, y : coordonnées (peuvent être scalaires, vecteurs ou matrices)
%
% Sortie :
%     val : valeur(s) de f(x,y)

    val = (1 + 2*pi^2) .* sin(pi*x) .* sin(pi*y);

end
