function VideoMaker_Stress(coords, conecs, u, failStat, stress, Tsteps, OutputsName, FrameRate, A_t)
% This function makes a video of the structure's deformation along with its
% stress contour throughout the analysis.

nframes = size(u, 1);

coords_ult = zeros(size(u));
for i = 1:nframes
    coords_ult(i, :) = u(i, :) + reshape(coords', 1, '');
end

indice = conecs(1, :);
lngth_x = coords(indice(2), 1) - coords(indice(1), 1);
lngth_y = coords(indice(2), 2) - coords(indice(1), 2);
elmlngth = sqrt(lngth_x^2 + lngth_y^2);
edge = 0.5 * elmlngth;

figure
set(gcf, 'units', 'points', 'position', [500, 150, 600, 600])
v = VideoWriter(['Outputs\\', OutputsName,'.avi']);
v.FrameRate = FrameRate;
open(v)
for i=1:nframes
    clf
    axis([min(min(coords_ult(1:2:end)))-edge max(max(coords_ult(1:2:end)))+edge...
        min(min(coords_ult(2:2:end)))-edge max(max(coords_ult(2:2:end)))+edge])
%     coords_istep = reshape(coords_ult(i, :), 2, '')';
%     axis([min(coords_istep(:, 1))-edge max(coords_istep(:, 1))+edge...
%         min(coords_istep(:, 2))-edge max(coords_istep(:, 2))+edge])
    Drawdef_Stress(stress(i, :), max(max(stress - min(min(stress)))),...
        min(min(stress)), coords, conecs, reshape(u(i, :), 2, '')', failStat(i, :), A_t)
    axis off
    title(sprintf('%s,  T = %.4f (s)', OutputsName, Tsteps(i)),...
        'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'bold')
    xlabel('X (m)', 'FontName', 'Times', 'FontWeight', 'normal')
    ylabel('Y (m)', 'FontName', 'Times', 'FontWeight', 'normal')
    set(gca, 'Box', 'on', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Layer', 'top')
    grid on
    drawnow();
    writeVideo(v, getframe(gcf));
end

end
