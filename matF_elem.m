function [Fel] = matF_elem(S1, S2, S3)
% calcul de la matrice elementaire du bloc rectangulaire (p, dv2/dy)
%
% SYNOPSIS [Fel] = matF_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Fel matrice elementaire rectangulaire (matrice 6x3)

%On fixe les sommets du traingle
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Transformation
J    = [x2-x1, x3-x1;
        y2-y1, y3-y1];
detJ  = abs(det(J));
invJT = inv(J)'; 

% quadrature d'ordre 2 sur le triangle de reference idendité
gp = [ 1/6, 1/6;
       2/3, 1/6;
       1/6, 2/3 ];
w  = [ 1/6, 1/6, 1/6 ];

Fel = zeros(6,3);

for k = 1:3
    xi  = gp(k,1);
    eta = gp(k,2);

    % coordonées barycentriques
    l1 = 1 - xi - eta;
    l2 = xi;
    l3 = eta;

    % fonctions de base pression
    Np = [l1; l2; l3];

    %  derivees des fonctions comme d'habitude
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

    dN_dxi(6)  = -4*l3;
    dN_deta(6) =  4*(l1 - l3);


    dN_dy = zeros(6,1);
    for i = 1:6
        grad_ref   = [dN_dxi(i); dN_deta(i)];
        grad_h  = invJT * grad_ref;
        dN_dy(i)   = grad_h(2);  
    end
 
    % on utilise la formule de quadrature
    for i = 1:6
        for j = 1:3
            Fel(i,j) = Fel(i,j) - detJ * w(k) * Np(j) * dN_dy(i);
        end
    end
end

end
