function val = f(x,y)
% Donnée pour u = 3*cos(pi*x)*cos(2*pi*y) avec -Δu + u = f
% Ici, Δu = -5*pi^2*u  =>  f = (1 + 5*pi^2) * u
val = 3*(1 + 5*pi^2) .* cos(pi*x) .* cos(2*pi*y);
end
