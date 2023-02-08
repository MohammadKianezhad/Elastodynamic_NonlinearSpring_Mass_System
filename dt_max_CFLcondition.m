function dt_max = dt_max_CFLcondition(coords, displacement, failStat, conecs, E, Rho, SF)
% This function calculates the largest time step with which the simulation will be stable.
% To enforce stable results, the time step size is limited so that a stress wave cannot
% travel farther than the smallest element characteristic length in a single time step.
% This is called the Courant Friedrichs Lewy (CFL) condition.
% SF is a safety factor, usually equal to or smaller than 1.

nounfailedElmn = length(find(failStat == 0));

coords_ult = coords + displacement;

C = sqrt(E/Rho);                                                           % the material’s sound speed wave

elmnlngth_t = zeros(nounfailedElmn, 1);
for i = 1:nounfailedElmn
    indice = conecs(i, :);
    lngth_x = coords_ult(indice(2), 1) - coords_ult(indice(1), 1);
    lngth_y = coords_ult(indice(2), 2) - coords_ult(indice(1), 2);
    elmnlngth_t(i) = sqrt(lngth_x^2 + lngth_y^2);
end
h = min(elmnlngth_t);
dt_max = SF*(h/C);

end