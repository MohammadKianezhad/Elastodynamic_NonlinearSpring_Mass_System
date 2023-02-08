function [coords, conecs, nonod_t, noelem_t, ndtx, eltx] = Mesh_Square(rank_x, rank_y, lngth, mode)
% This function produces a mesh with "rank" iterations of square unitcellin X- and Y-direction.
%
% The shape of the cell unit and the overall mesh for mode "NoDiag" is as follows:
%      _____               _____ _____ _____ _____
%     |     |             |     |     |     |     |
%     |     |             |     |     |     |     |
%     |_____|             |_____|_____|_____|_____|
%     unitcell            |     |     |     |     |
%                         |     |     |     |     |
%                         |_____|_____|_____|_____|
%                         |     |     |     |     |
%                         |     |     |     |     |
%                         |_____|_____|_____|_____|
%                   A "NoDiag" mode mesh with rank 4 in the
%                  X-direction and rank 3 in the Y-direction
%
% The shape of the cell unit and the overall mesh for mode "1Diag" is as follows:
%      _____               _____ _____ _____ _____
%     |    /|             |    /|    /|    /|    /|
%     |  /  |             |  /  |  /  |  /  |  /  |
%     |/____|             |/____|/____|/____|/____|
%     unitcell            |    /|    /|    /|    /|
%                         |  /  |  /  |  /  |  /  |
%                         |/____|/____|/____|/____|
%                         |    /|    /|    /|    /|
%                         |  /  |  /  |  /  |  /  |
%                         |/____|/____|/____|/____|
%                   A "1Diag" mode mesh with rank 4 in the
%                  X-direction and rank 3 in the Y-direction
%
% The shape of the cell unit and the overall mesh for mode "2Diag" is as follows:
%      ____                ____ ____ ____ ____
%     |\  /|              |\  /|\  /|\  /|\  /|
%     | \/ |              | \/ | \/ | \/ | \/ |
%     | /\ |              | /\ | /\ | /\ | /\ |
%     |/__\|              |/__\|/__\|/__\|/__\|
%     unitcell            |\  /|\  /|\  /|\  /|
%                         | \/ | \/ | \/ | \/ |
%                         | /\ | /\ | /\ | /\ |
%                         |/__\|/__\|/__\|/__\|
%                         |\  /|\  /|\  /|\  /|
%                         | \/ | \/ | \/ | \/ |
%                         | /\ | /\ | /\ | /\ |
%                         |/__\|/__\|/__\|/__\|
%                   A "2Diag" mode mesh with rank 4 in the
%                  X-direction and rank 3 in the Y-direction


nonod = 2;                                                                 % number of nodes in each element
nonod_t = (rank_x+1) * (rank_y+1);
coords = zeros(nonod_t, 2);
ndtx = zeros(nonod_t, 3);

inode = 0;
for irank = 1:rank_x+1                                                     % numbering the nodes and calculating their coordinates
    xpoin = (irank-1)*lngth;
    for jrank = 1:rank_y+1
        ypoin = (jrank-1)*lngth;
        inode = inode + 1;
        coords(inode, :) = [xpoin, ypoin];
        ndtx(inode, :) = [inode, xpoin, ypoin];
    end
end

if strcmpi(mode, 'NoDiag')
    noelem_t = 2*rank_x*rank_y + rank_x + rank_y;
    conecs = zeros(noelem_t, nonod);
    eltx = zeros(noelem_t, 3);
    
    ielmn = 0;
    for irank = 1:rank_x                                                       % element numbering and element connectivity
        for jrank = 1:rank_y
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+1];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+1];
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+1];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+1];
        end
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y+1];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y+1];
    end
    for jrank = 1:rank_y
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [rank_x*(rank_y+1) + jrank, rank_x*(rank_y+1) + jrank+1];
        eltx(ielmn, :) = [ielmn, rank_x*(rank_y+1) + jrank, rank_x*(rank_y+1) + jrank+1];
    end
elseif strcmpi(mode, '1Diag')
    noelem_t = 3*rank_x*rank_y + rank_x + rank_y;
    conecs = zeros(noelem_t, nonod);
    eltx = zeros(noelem_t, 3);
    
    ielmn = 0;
    for irank = 1:rank_x                                                       % element numbering and element connectivity
        for jrank = 1:rank_y
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+1];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+1];
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+1];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+1];
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+2];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+2];
        end
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y+1];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y+1];
    end
    for jrank = 1:rank_y
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [rank_x*(rank_y+1) + jrank, rank_x*(rank_y+1) + jrank+1];
        eltx(ielmn, :) = [ielmn, rank_x*(rank_y+1) + jrank, rank_x*(rank_y+1) + jrank+1];
    end
elseif strcmpi(mode, '2diag')
    noelem_t = 4*rank_x*rank_y + rank_x + rank_y;
    conecs = zeros(noelem_t, nonod);
    eltx = zeros(noelem_t, 3);
    
    ielmn = 0;
    for irank = 1:rank_x                                                       % element numbering and element connectivity
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + 1, (irank-1)*(rank_y+1) + 2];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + 1, (irank-1)*(rank_y+1) + 2];
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + 1, (irank-1)*(rank_y+1) + 1+rank_y+1];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + 1, (irank-1)*(rank_y+1) + 1+rank_y+1];
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + 1, (irank-1)*(rank_y+1) + 1+rank_y+2];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + 1, (irank-1)*(rank_y+1) + 1+rank_y+2];
        for jrank = 2:rank_y
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+1];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+1];
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+1];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+1];
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+2];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y+2];
            ielmn = ielmn + 1;
            conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y];
            eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank, (irank-1)*(rank_y+1) + jrank+rank_y];
        end
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y+1];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y+1];
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [(irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y];
        eltx(ielmn, :) = [ielmn, (irank-1)*(rank_y+1) + jrank+1, (irank-1)*(rank_y+1) + jrank+1+rank_y];
    end
    for jrank = 1:rank_y
        ielmn = ielmn + 1;
        conecs(ielmn, :) = [rank_x*(rank_y+1) + jrank, rank_x*(rank_y+1) + jrank+1];
        eltx(ielmn, :) = [ielmn, rank_x*(rank_y+1) + jrank, rank_x*(rank_y+1) + jrank+1];
    end
else
    error('The selected mode is not among the available modes. The mode can only be "NoDiag", "1Diag", or "2Diag".')
end


fid = fopen('Outputs\\mesh.dat', 'w');
%fprintf(fid, '%s %6s %8s\n', 'node', 'x', 'y');
%fprintf(fid, '%2d %10.3f %8.3f\n', ndtx');
fprintf(fid, '%f %f %f\n', ndtx');                                         % writing the number of nodes along with their coordinates in the mesh output file
fprintf(fid, '   \n');
%fprintf(fid, '%s %2s %2s\n', 'element', 'node1', 'node2');
%fprintf(fid, '%4d %6d %5d\n', eltx');
fprintf(fid, '%d %d %d\n', eltx');                                         % Writing the number of elements along with their connectivity in the mesh output file
fclose(fid);


end