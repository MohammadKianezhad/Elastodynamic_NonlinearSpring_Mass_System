function Drawdef_Stress(stress, stress_max, stress_min, coords, conecs, displacement, failStat, A_t)
% This function draws the structure in its deformed configuration along with
% the stress contour on elements.



stress1 = stress - stress_min;
stress_norm = stress1 / stress_max;

noelmn_t = size(conecs, 1);

coords_ult = coords + displacement;

colorMap = jet(256);
colormap(colorMap)
C = colorbar;
C.Ticks = linspace(0, 1, 8);
C.TickLabels = num2cell(round(linspace(stress_min/1e9, stress_max/1e9, 8), 1));
C.Title.String = 'Stress (GPa)';
C.FontName = 'Times';

for i = 1:noelmn_t
    if failStat(i) ==0
        x1elmn = coords_ult(conecs(i, 1), 1);
        y1elmn = coords_ult(conecs(i, 1), 2);
        x2elmn = coords_ult(conecs(i, 2), 1);
        y2elmn = coords_ult(conecs(i, 2), 2);
        line([x1elmn, x2elmn], [y1elmn, y2elmn], 'color',...
            colorMap(ceil(stress_norm(i)*255 + 1e-8), :), 'linestyle', '-', 'linewidth', A_t(i)*150);
    end
end


end
