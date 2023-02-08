function ddeltau = Solver(Kt, si, impdof)
% function for 0-1 method Solver

Kt(:, impdof) = 0;
Kt(impdof, :) = 0;
for i = 1:length(impdof)
    Kt(impdof(i), impdof(i)) = 1;
end

si(impdof) = 0;

ddeltau = Kt \ si;                                                            % d(deltaU)=invers(Ktotal)*(-Psi)

end