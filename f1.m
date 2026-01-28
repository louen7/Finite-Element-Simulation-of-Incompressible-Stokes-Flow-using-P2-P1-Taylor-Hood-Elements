function val = f1(x, y)
% f1 : second membre pour u_ex = 1 + 3*sin(pi*x)*sin(2*pi*y)
% -Delta u_ex + u_ex = 1 + (15*pi^2 + 3) sin(pi x) sin(2 pi y)

%s = sin(pi*x).*sin(2*pi*y);
%val = 1 + (15*pi^2 + 3) * s;
val = 3 * (5*pi^2 + 1) * cos(pi*x) .* cos(2*pi*y);
end
