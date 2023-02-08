function [coords, conecs, nonod_t, noelem_t, ndtx, eltx] = Mesh_Triangle(rank, lngth)
% This function produces a mesh with "rank" iterations of triangle unitcell.
% The shape of the cell unit and the overall mesh is as follows:
%
%     |\                  |\
%     | \                 | \
%     |__\                |__\
%   unitcell              |\  |\
%                         | \ | \
%                         |__\|__\
%                         |\  |\  |\
%                         | \ | \ | \
%                         |__\|__\|__\
%                       a third rank mesh

nonod = 2;                                                                 % number of nodes in each element
nonod_t = (rank+1)*(rank+2)/2;
noelem_t = (rank)*(rank+1)/2*3;
coords = zeros(nonod_t, 2);
conecs = zeros(noelem_t, nonod);
ndtx = zeros(nonod_t, 3);
eltx = zeros(noelem_t, 3);

inode = 0;
for irank = 1:rank+1                                                       % numbering the nodes and calculating their coordinates
    ypoin = (rank-irank+1)*lngth;
    for jrank = 1:irank
        xpoin = (jrank-1)*lngth;
        inode = inode+1;
        coords(inode, :) = [xpoin, ypoin];
        ndtx(inode, :) = [inode, xpoin, ypoin];
    end
end

kelmn = 3;
ielem = 0;
for irank = 1:rank                                                         % element numbering and element connectivity
    inods = zeros(3, 1);
    
    for jrank = 1:irank
        inods(1) = ((irank-1)^2+(irank-1)+2)/2+jrank-1;
        inods(2) = ((irank)^2+(irank)+2)/2+jrank-1;
        inods(3) = inods(2)+1;
        
        for kelem = 1:kelmn
            ielem = ielem+1;
            inod1 = kelem;
            inod2 = mod(kelem, kelmn)+1;
            conecs(ielem, :) = [inods(inod1), inods(inod2)];
            eltx(ielem, :) = [ielem, inods(inod1), inods(inod2)];
        end
        
    end
    
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