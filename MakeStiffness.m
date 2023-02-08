function stiffness = MakeStiffness(u, failStat, nodof_t, noelmn_t, conecs, coords, E, A, alpha)
% This function calculates the stiffness matrix of the whole system.

coords_x = coords(:, 1);
coords_y = coords(:, 2);

stiffness = zeros(nodof_t);
for i = 1:noelmn_t
    if failStat(i) == 0
        indice = conecs(i, :);
        elmndof = [2*indice(1)-1 2*indice(1) 2*indice(2)-1 2*indice(2)];   % calculation of global degrees of freedom from the nodes' numbers
        lngth_x = coords_x(indice(2)) - coords_x(indice(1));
        lngth_y = coords_y(indice(2)) - coords_y(indice(1));
        elmlngth = sqrt(lngth_x*lngth_x + lngth_y*lngth_y);
        C = lngth_x / elmlngth;                                            % cosine
        S = lngth_y / elmlngth;                                            % sinus
        uu = sqrt( (u(2*indice(2)-1) - u(2*indice(1)-1))^2 + (u(2*indice(2)) - u(2*indice(1)))^2 );
        Kelmn = (E*A(i)/elmlngth) * 2e-8 * (5000*uu + alpha*uu^3) *...
            [C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;-C*C -C*S C*C C*S;-C*S -S*S C*S S*S];      % calculation stiffness matrix of element
        stiffness(elmndof, elmndof) = stiffness(elmndof, elmndof) + Kelmn;       % assembling elemental stiffness matrix in global stiffness matrix
    end
end

end