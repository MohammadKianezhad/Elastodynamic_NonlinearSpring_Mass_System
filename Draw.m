function Draw(coords, conecs)
% This function draws the structure in its initial configuration.


noelmn_t = size(conecs, 1);

for i = 1:noelmn_t
    x1elmn = coords(conecs(i, 1), 1);
    y1elmn = coords(conecs(i, 1), 2);
    x2elmn = coords(conecs(i, 2), 1);
    y2elmn = coords(conecs(i, 2), 2);
    line([x1elmn, x2elmn], [y1elmn, y2elmn], 'color', [0.5 0.5 0.9], 'linestyle', '-.', 'linewidth', 1.4);
end

end
