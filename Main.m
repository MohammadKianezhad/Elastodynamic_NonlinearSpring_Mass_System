%% Elasto-dynamic Simulation of Nonlinear Spring-Mass System with Explicit Time Integration

% Notes:
%   - The formulation used in this code is taken from book "DYNAMICS OF
%     STRUCTURES" written by Anil K. Chopra.
%   - in this program the units are in meter (m), Second (s), Kilogram(Kg), and
%     Newtons (N).
%   - in this program "no*" stands for "number of *".
%   - in this program "elmn" stands for "element".
%   - in this program "dof" stands for "degree of freedom".
%   - in this program "u" stands for "displacement".
%   - in this program "ud" stands for "velocity" as the first derivative of "u".
%   - in this program "udd" stands for "accelaration" as the Second derivative of "u".


%% Freshen Up

fclose all;
close all;
clear;
clc;
tic;                                                                       % start of run time calculation


%% Input Information

force = 4e2;                                                               % applied force in N
T = 4;                                                                     % tutal duration of the analysis
rank_x = 30;                                                               % number of lattice in X-direction
rank_y = 30;                                                               % number of lattice in Y-direction
lngth_x = 1;                                                               % the length of the side of the lattice
lngth_y = 1;                                                               % the length of the side of the lattice
E = 2e11;                                                                  % Elasticity module in N/m^2 or Pa
Rho = 7850;                                                                % density
M = 10;                                                                    % mass of the nodes in Kg
C = 0;                                                                     % damping coefficient in N*s/m
sigmaY_t = 9e9;                                                            % yield stress in tension
sigmaY_c = -9e9;                                                           % yield stress in compression
nodof = 2;                                                                 % number of degree of freedom for each node (number of dimension)
alpha = 50;                                                                % a coefficient in the stiffness formulation (5000*u + alpha*u^3). A positive alpha results in Hardening spring and a negative alpha results in Softening spring
noframes = 200;                                                            % the number of video frames of structural deformation and stress contour
A = 0.01;                                                                  % the area of the main cross-sectional area of the elements(the cross-sectional area of the elements will randomly change a little in the future in order to include randomness into the network)
NewMatProp = 'yes';                                                        % if this option is equal to "yes", the random numbers related to the change of the cross-sectional area of the elements will be reproduced and stored in the "Outputs" folder. Otherwise, the previously produced numbers in the "Outputs" folder are used.
DeformVideo = 'yes';                                                       % produces an animation of the structure deformation during the simulation
StressVideo = 'yes';                                                       % produces an animation of the structure deformation along with the stress contour of the elements during the simulation
Display = 'yes';                                                           % if this option is equal to "yes", a report of the analysis process will be displayed in the Command Windows throughout the run.
Defect = 'no';                                                             % if this option is equal to "yes", then the defect will be created in the structure.
    Center = [15, 15];                                                     % The defect will be created in such a way that all the nodes inside the area to the center "Center"
    Radius = 2;                                                            % and to the radius "Radius" along with all the elements connected to them will be deleted.


%% Preprocessing

if ~isfolder('Outputs')
    mkdir('Outputs')
end

fid_log = fopen('Outputs\\log.txt', 'w');                                  % a report of the analysis process will be written in the log.txt file in the Outputs folder.
fprintf(fid_log, 'Elasto-Dynamic Simulation of Nonlinear Spring-Mass System with Explicit Time Integration\n');
fprintf(fid_log, '\n%s\n\n\n', datetime(now,'ConvertFrom','datenum'));

% [coords, conecs, nonod_t, noelmn_t] = Mesh_Triangle(rank_x, lngth);
[coords, conecs, nonod_t, noelmn_t] = Mesh_Square(rank_x, rank_y, lngth_x, lngth_y, '2Diag');  % mesh production

if strcmpi(NewMatProp, 'yes')
    Arand = randi([-400 400], noelmn_t, 1) / 1e5;
    fid = fopen('Outputs\\RandomnessOfArea.txt', 'w');
    fprintf(fid, '%f\n', Arand);
    fclose(fid);
else
    if ~isfile('Outputs\\RandomnessOfArea.txt')
        error(['There is no file related to the random numbers of cross-sectional ',...
            'area in the "Outputs" folder. Please change the option for NewMatProp to "Yes".'])
    end
end
fid2 = fopen('Outputs\\RandomnessOfArea.txt', 'r');
Arand = fscanf(fid2, '%f\n', [1 inf]);
fclose(fid2);
A_t = A + Arand';

defectElmnIDs = [];
if strcmpi(Defect, 'yes')
    defectNodeIDs = find(((coords(:, 1)-Center(1)).^2 + (coords(:, 2)-Center(2)).^2) < Radius^2);
    for i = 1:length(defectNodeIDs)
        defectElmnIDs = [defectElmnIDs; find((conecs(:, 1) == defectNodeIDs(i)) | (conecs(:, 2) == defectNodeIDs(i)))];
    end
    coords(defectNodeIDs, :) = NaN * ones(length(defectNodeIDs), 2);
    conecs(defectElmnIDs, :) = [];
    nonod_t = size(coords, 1);
    noelmn_t = size(conecs, 1);
    A_t(defectElmnIDs) = [];
end

nodof_t = nodof * nonod_t;                                                 % number of total degrees of freedom

dt = dt_max_CFLcondition(coords, zeros(nonod_t, nodof), zeros(noelmn_t, 1), conecs, E, Rho, 0.7);   % calculation of the largest time step for stable simulation, calculated from the CFL condition method
notstep = round(T/dt);

m = M * eye(nodof_t);
c = C * eye(nodof_t);
u0 = zeros(nodof_t, 1);
ud0 = zeros(nodof_t, 1);
u = zeros(notstep+1, nodof_t);
ud = zeros(notstep, nodof_t);
udd = zeros(notstep, nodof_t);
stress = zeros(notstep+1, noelmn_t);
Pn = zeros(notstep+1, noelmn_t);
failStat = zeros(notstep, noelmn_t);                                        % it works as a flag for broken elements. The index of failed elements would be equal to 1.
Tsteps = zeros(notstep, 1);

impdof = [1 2 2*(rank_x*(rank_y+1)+1)];


%% Main Analysis

% Initial calculations
Kt0 = MakeStiffness(u0, zeros(1, noelmn_t), nodof_t, noelmn_t, conecs, coords, E, A_t, alpha);
fs0 = Kt0 * u0;
p0 = zeros(nodof_t, 1);
p0(2*(rank_y+1)-1, 1) = force;
% p0(2*(rank_y+1), 1) = -force;
udd0 = m \ (p0 - c*ud0 - fs0);

u_1 = u0 - dt*ud0 + dt^2*udd0/2;
kcirc = m/(dt^2) + c/(2*dt);
a = m/(dt^2) - c/(2*dt);
b = -2*m / (dt^2);

p = zeros(notstep+1, nodof_t);
p(:, 2*(rank_y+1)-1) = force;                                               % excitation force

% Calculations for the first time step, i = 0
pcirc = p0 - a*u_1 - b*u0 - fs0;
u(1, :) = Solver(kcirc, pcirc, impdof);
Tsteps(1) = dt;
stress(1, :) = Stress(noelmn_t, conecs, coords, u(1, :), zeros(noelmn_t, 1), E);
failElmnID_t = find(stress(1, :) > sigmaY_t);
failElmnID_c = find(stress(1, :) < sigmaY_c);
failStat(:, [failElmnID_t failElmnID_c]) = 1;
Pn(1, :) = stress(1, :) * A_t;                                             % calculation of the force of the elements
fprintf(fid_log, ['In the time step of 1 and time of %.4f (s), the analysis',...
    ' of the structure has been done. Also, in this time step, %d elements were broken.\n'],...
    [dt, length(failElmnID_t) + length(failElmnID_c)]);
if strcmpi(Display, 'yes')
    fprintf(['In the time step of 1 and time of %.4f (s), the analysis of ',...
    'the structure has been done. Also, in this time step, %d elements were broken.\n'],...
    [dt, length(failElmnID_t) + length(failElmnID_c)]);
end

% Calculations for the second time step, i = 1
Kti = MakeStiffness(u(1, :), failStat(1, :), nodof_t, noelmn_t, conecs, coords, E, A_t, alpha);
fsi = Kti * u(1, :)';
pcirc = p(1, :)' - a*u0 - b*u(1, :)' - fsi;
u(2, :) = Solver(kcirc, pcirc, impdof);
dt_max = dt_max_CFLcondition(coords, reshape(u(2, :), 2, '')', failStat(1, :), conecs, E, Rho, 1);
if dt_max < dt
    u(2, :) = zeros(1, nodof_t);
    dt = dt_max;
    kcirc = m/(dt^2) + c/(2*dt);
    a = m/(dt^2) - c/(2*dt);
    b = -2*m / (dt^2);
    pcirc = p(1, :)' - a*u0 - b*u(1, :)' - fsi;
    u(2, :) = Solver(kcirc, pcirc, impdof);
    fprintf(fid_log, ['\n  ** In time step 2, according to method (CFL) condition, the time step was greater',...
        ' than the critical value. Therefore, the time step\n     was changed to the value of %.9f',...
        ' and the analysis in this step will be performed again according to the new time step. **\n'], dt);
    if strcmpi(Display, 'yes')
        fprintf(['\n  ** In time step 2, according to method (CFL) condition, the time step was greater',...
            ' than the critical value. Therefore, the time step\n     was changed to the value of %.9f',...
            ' and the analysis in this step will be performed again according to the new time step. **\n'], dt);
    end
end
Tsteps(2) = Tsteps(1) + dt;
stress(2, :) = Stress(noelmn_t, conecs, coords, u(2, :), failStat(1, :), E);
failElmnID_t = find(stress(2, :) > sigmaY_t);
failElmnID_c = find(stress(2, :) < sigmaY_c);
failStat(2:end, [failElmnID_t failElmnID_c]) = 1;
Pn(2, :) = stress(2, :) * A_t;
ud(1, :) = (u(2, :) - u0') / (2*dt);
udd(1, :) = (u(2, :) - 2*u(1, :) + u0') / (dt^2);
fprintf(fid_log, ['In the time step of 2 and time of %.4f (s), the analysis',...
    ' of the structure has been done. Also, in this time step, %d elements were broken.\n'],...
    [Tsteps(2), length(failElmnID_t) + length(failElmnID_c)]);
if strcmpi(Display, 'yes')
    fprintf(['In the time step of 2 and time of %.4f (s), the analysis of ',...
    'the structure has been done. Also, in this time step, %d elements were broken.\n'],...
    [Tsteps(2), length(failElmnID_t) + length(failElmnID_c)]);
end

% Calculations for each time step, i = 2, 3, 4, ...
% "for "Loop is not used here because considering that the critical time step obtained
% from the CFL condition method in each step may change due to the change in the length
% of the elements, the number of time steps required to reach the desired analysis length (T)
% may change. Therefore, the analysis continues until the analysis time reaches the desired time (T).
istep = 2;
while true
    Kti = MakeStiffness(u(istep, :), failStat(istep, :), nodof_t, noelmn_t, conecs, coords, E, A_t, alpha);
    fsi = Kti * u(istep, :)';
    pcirc = p(istep, :)' - a*u(istep-1, :)' - b*u(istep, :)' - fsi;
    u(istep+1, :) = Solver(kcirc, pcirc, impdof);
    dt_max = dt_max_CFLcondition(coords, reshape(u(istep+1, :), 2, '')', failStat(istep, :), conecs, E, Rho, 1);
    if dt_max < dt
        u(istep+1, :) = zeros(1, nodof_t);
        dt = dt_max;
        kcirc = m/(dt^2) + c/(2*dt);
        a = m/(dt^2) - c/(2*dt);
        b = -2*m / (dt^2);
        pcirc = p(istep, :)' - a*u(istep-1, :)' - b*u(istep, :)' - fsi;
        u(istep+1, :) = Solver(kcirc, pcirc, impdof);
        fprintf(fid_log, ['\n  ** In time step %d, according to method (CFL) condition, the time step was greater',...
            ' than the critical value. Therefore, the time step\n     was changed to the value of %.9f',...
            ' and the analysis in this step will be performed again according to the new time step. **\n'], [istep+1, dt]);
        if strcmpi(Display, 'yes')
            fprintf(['  ** In time step %d, according to method (CFL) condition, the time step was greater',...
                '\n than the critical value. Therefore, the time step\n     was changed to the value of %.9f',...
                ' and the analysis in this step will be performed again according to the new time step. **\n'], [istep+1, dt]);
        end
    end
    Tsteps(istep+1) = Tsteps(istep) + dt;
    stress(istep+1, :) = Stress(noelmn_t, conecs, coords, u(istep+1, :), failStat(istep, :), E);
    failElmnID_t = find(stress(istep+1, :) > sigmaY_t);
    failElmnID_c = find(stress(istep+1, :) < sigmaY_c);
    if istep >= notstep
        failStat = [failStat; failStat(end, :)];
        p = [p; p(end, :)];
    end
    failStat(istep+1:end, [failElmnID_t failElmnID_c]) = 1;
    Pn(istep+1, :) = stress(istep+1, :) * A_t;
    ud(istep, :) = (u(istep+1, :) - u(istep-1, :)) / (2*dt);
    udd(istep, :) = (u(istep+1, :) - 2*u(istep, :) + u(istep-1, :)) / (dt^2);
    fprintf(fid_log, ['In the time step of %d and time of %.4f (s), the analysis',...
        ' of the structure has been done. Also, in this time step, %d elements were broken.\n'],...
        [istep+1, Tsteps(istep+1), length(failElmnID_t) + length(failElmnID_c)]);
    if strcmpi(Display, 'yes')
        fprintf(['In the time step of %d and time of %.4f (s), the analysis of ',...
            'the structure has been done. Also, in this time step, %d elements were broken.\n'],...
            [istep+1, Tsteps(istep+1), length(failElmnID_t) + length(failElmnID_c)]);
    end
    if Tsteps(istep+1) >= T                                                % it checks that the simulation ends when the analysis time reaches the desired time (T)
        break
    end
    istep = istep + 1;
end


%% Postprocessing

notstep = istep;
if noframes > notstep
    frameID = 1:notstep;
else
    frameID = floor(linspace(1, notstep, noframes));
end
if strcmpi(DeformVideo, 'yes')
    VideoMaker(coords, conecs, u(frameID, :), failStat(frameID, :), Tsteps(frameID),...
        'Deformation of the Structure', 10, A_t)
end
if strcmpi(StressVideo, 'yes')
    VideoMaker_Stress(coords, conecs, u(frameID, :), failStat(frameID, :),...
        stress(frameID, :), Tsteps(frameID), 'Stress Contour of the Structure', 10, A_t)
end

elmnDeform = zeros(notstep, noelmn_t);                                     % changing the length of the elements during the analysis
for i = 1:noelmn_t
    indice = conecs(i, :);
    elmnDeform(:, i) = sqrt( (u(:, 2*indice(2)-1) - u(:, 2*indice(1)-1)).^2 +...
        (u(:, 2*indice(2)) - u(:, 2*indice(1))).^2 );
end

C_t = zeros(noelmn_t, 1);
S_t = zeros(noelmn_t, 1);
for i = 1:noelmn_t
    indice = conecs(i, :);
    elementdof = [2*indice(1)-1 2*indice(1) 2*indice(2)-1 2*indice(2)];
    lngth_x = coords(indice(2), 1) - coords(indice(1), 1);
    lngth_y = coords(indice(2), 2) - coords(indice(1), 2);
    elmlngth = sqrt(lngth_x*lngth_x + lngth_y*lngth_y);
    C_t(i) = lngth_x/elmlngth;                                             % the cosine of the angle between the element and the X-axis
    S_t(i) = lngth_y/elmlngth;                                             % The sinus of the angle between the element and the X-axis
end
Fnodal = zeros(notstep, 2*nonod_t);                                        % the resultant force applied to the nodes during the analysis
for i = 1:nonod_t
    elmncnct = find((conecs(:, 1) == i) | (conecs(:, 2) == i));
    Fnodal(:, 2*i-1) = Pn(:, elmncnct) * C_t(elmncnct);
    Fnodal(:, 2*i) = Pn(:, elmncnct) * S_t(elmncnct);
end


fileNameForSave = sprintf("Outputs\\Results_ExplicitDynamicsAnalysis.mat");
save (fileNameForSave, '-v7.3')

runtime = toc;                                                             % end of run time calculation
fprintf(fid_log, '\n    Run Time => %2d:%2d\n    END\n',[floor(runtime/60),round(mod(runtime,60))]);               % run-time report
if strcmpi(Display, 'yes')
    fprintf('\n    Run Time => %2d:%2d\n    END\n',[floor(runtime/60),round(mod(runtime,60))])               % run-time report
end
fclose(fid_log);
