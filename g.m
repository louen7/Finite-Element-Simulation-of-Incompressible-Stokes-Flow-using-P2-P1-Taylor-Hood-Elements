function val = g(x, y)
% g : condition de Dirichlet non homogène sur le bord
% ici on choisit : u_ex(x,y) = 1 + 3*sin(pi*x)*sin(2*pi*y)
%val = 1 + 3*sin(pi*x).*sin(2*pi*y);
val = 3 * cos(pi*x) .* cos(2*pi*y); %uexact = 3 * cos(pi*x) .* cos(2*pi*y)
end
