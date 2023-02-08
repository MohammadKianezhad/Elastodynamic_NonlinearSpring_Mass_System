function Drawdef(coords, conecs, displacement, failStat, A_t)
% This function draws the structure in its deformed configuration.


noelmn_t = size(conecs, 1);

coords_ult = coords + displacement;

axis equal
for i = 1:noelmn_t
    if failStat(i) == 0
        x1elmn = coords_ult(conecs(i, 1), 1);
        y1elmn = coords_ult(conecs(i, 1), 2);
        x2elmn = coords_ult(conecs(i, 2), 1);
        y2elmn = coords_ult(conecs(i, 2), 2);
        line([x1elmn, x2elmn], [y1elmn, y2elmn], 'color', 'r', 'linestyle', '-', 'linewidth', A_t(i)*150);
    end
end

hold on

scatter(coords_ult(:, 1), coords_ult(:, 2), 40 * ones(1, size(coords_ult, 1)), 'o', 'MarkerEdgeColor','k',...
              'MarkerFaceColor', [0 0.4470 0.7410], 'LineWidth', 0.2)

end
