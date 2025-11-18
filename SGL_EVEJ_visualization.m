%% SGL_EVEJ_visualization_2D.m
% 2D schematic E–V–E–J + Solar Oberth visualization
% with animated spacecraft and SGL alignment in the x–y plane only.
% Exports an animated GIF for PowerPoint.

clear; close all; clc;

%% Basic distances (AU)
r_venus   = 0.723;
r_earth   = 1.000;
r_jupiter = 5.203;
r_p       = 0.046;

%% RA/Dec of Proxima Centauri b -> 3D unit vector
RA_h   = 14 + 29/60 + 42.94613/3600;          % hours
RA_deg = 15 * RA_h;                           % deg
Dec_deg = -(62 + 40/60 + 46.16468/3600);      % deg

RA  = RA_deg  * pi/180;
Dec = Dec_deg * pi/180;

u_prox_3D = [cos(Dec)*cos(RA), ...
             cos(Dec)*sin(RA), ...
             sin(Dec)];

% Project onto x–y plane (ignore z) and renormalize
u_prox_2D = u_prox_3D(1:2);
u_prox_2D = u_prox_2D / norm(u_prox_2D);

% Outbound SGL direction: exactly opposite in 2D
u_sgl_2D  = -u_prox_2D;

angle_deg = acosd(dot(u_prox_2D, u_sgl_2D) / ...
                  (norm(u_prox_2D)*norm(u_sgl_2D)));
fprintf('Angle between outbound SGL direction and Proxima vector (2D): %.2f deg\n', angle_deg);

%% Circular orbits in x–y plane
theta = linspace(0, 2*pi, 500);

xV_orb = r_venus   * cos(theta); 
yV_orb = r_venus   * sin(theta);

xE_orb = r_earth   * cos(theta); 
yE_orb = r_earth   * sin(theta);

xJ_orb = r_jupiter * cos(theta); 
yJ_orb = r_jupiter * sin(theta);

xP     = r_p       * cos(theta); 
yP     = r_p       * sin(theta);

%% Schematic planet positions along their orbits (2D)
theta_E1 = 0   * pi/180;   % Earth departure
theta_V1 = -40 * pi/180;   % Venus GA
theta_E2 =  30 * pi/180;   % Earth GA
theta_J  = 110 * pi/180;   % Jupiter GA

rE1 = [r_earth*cos(theta_E1), r_earth*sin(theta_E1)]; % Earth dep
rV1 = [r_venus*cos(theta_V1), r_venus*sin(theta_V1)]; % Venus GA
rE2 = [r_earth*cos(theta_E2), r_earth*sin(theta_E2)]; % Earth GA
rJ  = [r_jupiter*cos(theta_J),r_jupiter*sin(theta_J)];% Jupiter GA

% Perihelion point near the Sun, along line from Jupiter to Sun (2D)
dir_J_to_Sun_2D = -rJ / norm(rJ);
rP_vec_2D       = r_p * dir_J_to_Sun_2D;

%% Trajectory segments (straight lines, schematic, all 2D)
Nseg = 80; 
s    = linspace(0,1,Nseg)';

EV_traj = (1-s).*rE1 + s.*rV1;       % Earth dep -> Venus
VE_traj = (1-s).*rV1 + s.*rE2;       % Venus -> Earth GA
EJ_traj = (1-s).*rE2 + s.*rJ;        % Earth GA -> Jupiter
JS_traj = (1-s).*rJ  + s.*rP_vec_2D; % Jupiter -> perihelion

%% Outbound SGL leg (starting at perihelion, along 2D SGL direction)
L_vec   = 4;                               % 4 AU for compact view
t_line  = linspace(0, L_vec, Nseg)';       % same #points for simplicity
SGL_traj_2D = rP_vec_2D + t_line .* u_sgl_2D;  % perihelion outward in 2D

%% Proxima vector (from Sun, 2D)
t_prox     = linspace(0, L_vec, Nseg)';
Prox_traj_2D  = t_prox .* u_prox_2D;

%% Combine path for animation (avoid duplicating endpoints)
traj_all = [EV_traj;
            VE_traj(2:end,:);
            EJ_traj(2:end,:);
            JS_traj(2:end,:);
            SGL_traj_2D(2:end,:)];

%% Plot static elements (2D only)
fig = figure('Color','w', 'WindowState','maximized');
hold on; grid on; axis equal;

h_orbV = plot(xV_orb, yV_orb, '--', 'Color', [0.7 0.3 0.7], 'LineWidth', 1);
h_orbE = plot(xE_orb, yE_orb, '--', 'Color', [0.3 0.3 0.9], 'LineWidth', 1);
h_orbJ = plot(xJ_orb, yJ_orb, '--', 'Color', [0.3 0.7 0.3], 'LineWidth', 1);
h_rp   = plot(xP,     yP,     ':',  'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% Sun
plot(0,0,'yo','MarkerFaceColor','y','MarkerSize',8);

% Planet markers
plot(rE1(1), rE1(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
text(rE1(1), rE1(2), '  Earth dep', 'Color', 'b');

plot(rV1(1), rV1(2), 'md', 'MarkerFaceColor', 'm', 'MarkerSize', 7);
text(rV1(1), rV1(2), '  Venus GA', 'Color', 'm');

plot(rE2(1), rE2(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
text(rE2(1), rE2(2), '  Earth GA', 'Color', 'b');

plot(rJ(1), rJ(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
text(rJ(1), rJ(2), '  Jupiter GA', 'Color', 'g');

% Perihelion point
plot(rP_vec_2D(1), rP_vec_2D(2), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Static Proxima and SGL direction vectors (2D)
h_proxVec = plot(Prox_traj_2D(:,1), Prox_traj_2D(:,2), 'r', 'LineWidth', 2);
text(Prox_traj_2D(end,1), Prox_traj_2D(end,2), ...
    '  Towards Proxima Centauri b', 'Color', 'r');

h_sglVec = plot(SGL_traj_2D(:,1), SGL_traj_2D(:,2), 'g', 'LineWidth', 2);
text(SGL_traj_2D(end,1), SGL_traj_2D(end,2), ...
    '  Final outbound trajectory (SGL)', 'Color', 'g');

xlabel('x [AU]');
ylabel('y [AU]');
title('2D E–V–E–J + Solar Oberth Trajectory and SGL Alignment');

xlim([-6 6]); ylim([-6 6]);

% Animated path + spacecraft marker (2D)
h_path = plot(NaN,NaN,'c','LineWidth',2);
h_sc   = plot(traj_all(1,1), traj_all(1,2), 'ko', ...
              'MarkerFaceColor','c','MarkerSize',6);


%% Animate spacecraft along full trajectory AND save as 2D GIF
nSteps   = size(traj_all,1);
filename = 'SGL_EVEJ_trajectory_2D.gif';

for k = 1:nSteps
    % update path and spacecraft marker
    set(h_path, 'XData', traj_all(1:k,1), ...
                'YData', traj_all(1:k,2));
    set(h_sc,   'XData', traj_all(k,1), ...
                'YData', traj_all(k,2));
    drawnow;

    % capture current figure as an image
    frame = getframe(gcf);
    im    = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % write or append to GIF
    if k == 1
        imwrite(imind, cm, filename, ...
                'gif', 'Loopcount', inf, 'DelayTime', 0.03);
    else
        imwrite(imind, cm, filename, ...
                'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
    end
end

disp(['Saved animated GIF: ', filename]);
