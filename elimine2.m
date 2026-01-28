function [AA_tilde, LL_tilde] = elimine2(AA, LL, Refneu)
% =====================================================
% Réalise la pseudo-élimination pour conditions de Dirichlet homogènes :
% Pour chaque nœud M sur le bord (Refneu(M) = 1),
% on impose :  AA(M,:) = 0, AA(:,M) = 0, AA(M,M) = 1, LL(M) = 0.
% on a pris une constante ρ = 1.
% =====================================================

AA_tilde = AA;
LL_tilde = LL;

for I = 1:length(Refneu)
    if Refneu(I) == 1  % noeud sur le bord
        AA_tilde(I, :) = 0;      % ligne nulle
        AA_tilde(:, I) = 0;      % colonne nulle
        AA_tilde(I, I) = 1;      % imposer 1 sur la diagonale
        LL_tilde(I) = 0;         % valeur imposée (u = 0)
    end
end

end