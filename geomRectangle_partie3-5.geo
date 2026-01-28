// ====================================================
// geomRectangle_partie3-5.geo
// Domaine type "marche" / cavité en bas à gauche
//
// Géométrie :
//   - Canal principal : x ∈ [0,3], y ∈ [0,1]
//   - Cavité : x ∈ [0,1], y ∈ [-0.5, 0]
//   -> domaine en forme de marche, avec une cuvette à gauche en bas.
//
// Bords physiques :
//   Gamma_down : fond de la cavité
//   Bottom     : bas du canal principal
//   Right      : paroi droite (sortie potentielle)
//   Top        : haut du canal
//   Left       : paroi gauche (verticale)
//
// Maillage :
//   - Élément d'ordre 2 (P2)
//   - Format .msh 2.2 ASCII (=> ligne MeshFormat "2.2 0 8")
// ====================================================

h = 0.1;    // taille de maille cible

// ---------------------------
// Points (numérotation anti-horaire)
// ---------------------------
// Cavité (en bas à gauche)
Point(1) = {0,  -0.5, 0, h};  // A : bas gauche cavité
Point(2) = {1,  -0.5, 0, h};  // B : bas droite cavité
Point(3) = {1,   0.0, 0, h};  // C : coin marche

// Canal principal
Point(4) = {3,   0.0, 0, h};  // D : bas droit canal
Point(5) = {3,   1.0, 0, h};  // E : haut droit canal
Point(6) = {0,   1.0, 0, h};  // F : haut gauche canal
Point(7) = {0,   0.0, 0, h};  // G : haut cavité / bas gauche canal

// ---------------------------
// Lignes de bord
// ---------------------------
// Fond de cavité (Gamma_down)
Line(1) = {1, 2};   // y = -0.5

// Paroi verticale entre cavité et canal (marche)
Line(2) = {2, 3};   // x = 1

// Bas du canal principal
Line(3) = {3, 4};   // y = 0, x de 1 à 3

// Paroi droite
Line(4) = {4, 5};   // x = 3

// Haut du canal
Line(5) = {5, 6};   // y = 1

// Paroi gauche (partie haute)
Line(6) = {6, 7};   // x = 0, y de 1 à 0

// Paroi gauche (partie basse, côté cavité)
Line(7) = {7, 1};   // x = 0, y de 0 à -0.5

// ---------------------------
// Surface fluide
// ---------------------------
Line Loop(10) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(20) = {10};

// ---------------------------
// Groupes physiques
// ---------------------------
// Surface fluide
Physical Surface("Fluide") = {20};

// Fond de cavité
Physical Line("Gamma_down") = {1};

// Bas du canal principal
Physical Line("Bottom")     = {3};

// Paroi droite
Physical Line("Right")      = {4};

// Haut du canal
Physical Line("Top")        = {5};

// Paroi gauche (toute la verticale x=0)
Physical Line("Left")       = {6, 7};

// (Optionnel : la marche verticale interne comme paroi aussi)
// Physical Line("Step")   = {2};

// ---------------------------
// Options de maillage pour Matlab (P2 + .msh 2.2)
// ---------------------------
Mesh.ElementOrder           = 2;    // éléments d'ordre 2 (P2)
Mesh.SecondOrderIncomplete  = 0;    // éléments quadratiques complets
Mesh.MshFileVersion         = 2.2;  // format .msh 2.2 (ligne : "2.2 0 8")
