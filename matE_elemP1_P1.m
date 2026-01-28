function [Eel] = matE_elemP1_P1(S1, S2, S3)

% calcul de la matrice elementaire E^l pour le couplage (p, v1) en P1-P1
% selon les formules trouvées sur papier à l'aide la question 2 
%
% (El)I1 = -(2/3) [y3 - y1 ; y1 - y2] · grad_hat(lambda_I)
% (El)I2 = -(1/6) [y3 - y1 ; y1 - y2] · grad_hat(lambda_I)
% (El)I3 = -(1/6) [y3 - y1 ; y1 - y2] · grad_hat(lambda_I)
%
% SYNOPSIS [Eel] = matE_elem_P1_P1(S1, S2, S3)
%
% INPUT  * S1, S2, S3 : coordonnees des 3 sommets du triangle [x y]
% OUTPUT - Eel : matrice elementaire 3x3


% On fixe coordonnees des sommets
x1 = S1(1);  y1 = S1(2);
x2 = S2(1);  y2 = S2(2);
x3 = S3(1);  y3 = S3(2);

% vecteur de la transformation
vE = [y3 - y1;
      y1 - y2];

% gradients barycentriques sur le triangle de reference
grads_hat = [-1, -1;   
              1,  0;  
              0,  1];  

Eel = zeros(3,3);

for I = 1:3
    gI   = grads_hat(I, :).';       
    prod = vE.' * gI;               

    Eel(I,1) = -(2/3) * prod;      
    Eel(I,2) = -(1/6) * prod;       
    Eel(I,3) = -(1/6) * prod;       
end

end


