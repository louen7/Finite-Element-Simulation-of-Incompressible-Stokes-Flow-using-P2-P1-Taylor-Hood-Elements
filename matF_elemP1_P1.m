function [Fel] = matF_elemP1_P1(S1, S2, S3)

% calcul de la matrice elementaire F^l pour le couplage (p, v2) en P1-P1
% On utilise les formules déterminées sur papier
% (Fl)I1 = -(2/3) [x1 - x3 ; x2 - x1] · grad_hat(lambda_I)
% (Fl)I2 = -(1/6) [x1 - x3 ; x2 - x1] · grad_hat(lambda_I)
% (Fl)I3 = -(1/6) [x1 - x3 ; x2 - x1] · grad_hat(lambda_I)
% INPUT  * S1, S2, S3 : coordonnees des 3 sommets du triangle [x y]
% OUTPUT - Fel : matrice elementaire 3x3

% coordonnees des sommets
x1 = S1(1);  y1 = S1(2);
x2 = S2(1);  y2 = S2(2);
x3 = S3(1);  y3 = S3(2);

% vecteur venant de la 2eme ligne de (B_l^T)^(-1) * |det B_l|
vF = [x1 - x3;
      x2 - x1];

% gradients barycentriques sur le triangle de reference "identité"
grads_h = [-1, -1;   
              1,  0;   
              0,  1];  

Fel = zeros(3,3);

for I = 1:3
    gI   = grads_h(I, :).';       
    prod = vF.' * gI;               


    Fel(I,1) = -(2/3) * prod;       
    Fel(I,2) = -(1/6) * prod;       
    Fel(I,3) = -(1/6) * prod;      
end

end
