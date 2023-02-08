function stress = Stress(noelmn_t, conecs, coords, displacement, failStat, E)
% This function calculates the stress of each element.

coords_x = coords(:, 1);
coords_y = coords(:, 2);

stress=zeros(noelmn_t, 1);
for i = 1:noelmn_t
    if failStat(i) == 0
        indice = conecs(i, :);
        elmndof = [2*indice(1)-1 2*indice(1) 2*indice(2)-1 2*indice(2)];   % calculation of global degrees of freedom from the nodes' numbers
        lngth_x = coords_x(indice(2)) - coords_x(indice(1));
        lngth_y = coords_y(indice(2)) - coords_y(indice(1));
        elmnlngth = sqrt(lngth_x*lngth_x + lngth_y*lngth_y);
        C = lngth_x / elmnlngth;                                           % cosine
        S = lngth_y / elmnlngth;                                           % sinus
        stress(i) = (E/elmnlngth)*[-C -S C S]*(displacement(elmndof))';    % calculation of the stress of the element which is obtained by multiplying the modulus of elasticity (E) by the strain
    end
end

end