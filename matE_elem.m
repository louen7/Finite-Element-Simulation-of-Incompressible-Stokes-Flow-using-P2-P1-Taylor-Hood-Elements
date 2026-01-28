function [Eel] = matE_elem(S1, S2, S3)

% calcul de la matrice elementaire du bloc rectangulaire (p, dv1/dx)
%
% SYNOPSIS [Eel] = matE_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Eel matrice elementaire rectangulaire (matrice 6x3)
%

x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% matrice J de la transformation 
J  = [x2-x1, x3-x1;
      y2-y1, y3-y1];
detJ  = abs(det(J));
invJT = inv(J)';   % (J^{-T})

% quadrature d'ordre 2 sur le triangle de reference
% nous notons les points (xi,eta) et poids correspondants
gp = [ 1/6, 1/6;
       2/3, 1/6;
       1/6, 2/3 ];
w  = [ 1/6, 1/6, 1/6 ];  

Eel = zeros(6,3);

for k = 1:3
    xi  = gp(k,1);
    eta = gp(k,2);

    % coordonnées barycentriques
    l1 = 1 - xi - eta;
    l2 = xi;
    l3 = eta;

   
    %  P1 : fonctions de base pression
    Np = [l1; l2; l3];   
    % dérivées par rapport à (xi,eta)
    dN_dxi  = zeros(6,1);
    dN_deta = zeros(6,1);

    dN_dxi(1)  = -(4*l1 - 1);
    dN_deta(1) = -(4*l1 - 1);

    dN_dxi(2)  =  (4*l2 - 1);
    dN_deta(2) =  0;

    dN_dxi(3)  =  0;
    dN_deta(3) =  (4*l3 - 1);

    dN_dxi(4)  = 4*(l1 - l2);
    dN_deta(4) = -4*l2;

    dN_dxi(5)  = 4*l3;
    dN_deta(5) = 4*l2;

    % N6 = 4 l1 l3
    dN_dxi(6)  = -4*l3;
    dN_deta(6) =  4*(l1 - l3);

    % passage à (x,y) par la transformation
    dN_dx = zeros(6,1);
    for i = 1:6
        grad_ref   = [dN_dxi(i); dN_deta(i)];
        grad_h  = invJT * grad_ref;
        dN_dx(i)   = grad_h(1);   % on ne garde que dN_i/dx
    end

   
    for i = 1:6
        for j = 1:3
            Eel(i,j) = Eel(i,j) - detJ * w(k) * Np(j) * dN_dx(i);
        end
    end
end
end
