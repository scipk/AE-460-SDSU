%% Solar Gravity Lens Mission - Gravity Assist Delta-V Analysis
% Calculates delta-V requirements for E-V-E-J sequence + Solar Oberth
% Starting from GEO
% by Parham Khodadi
clear; close all; clc;

%% Conversions
yr_to_sec = 31556952;       % Years to seconds
AU_to_m = 149597870700;     % Astronomical units to meters

%% Gravitational Parameters [m^3/s^2]
muSun     = 1.32712440042e20;
muEarth   = 3.986004418e14;
muVenus   = 3.2485859e14;
muJupiter = 1.26686534e17;

%% Planetary Orbital Data
r_earth_orbit   = 1.000 * AU_to_m;      % Earth orbital radius [m]
r_venus_orbit   = 0.723 * AU_to_m;      % Venus orbital radius [m]
r_jupiter_orbit = 5.203 * AU_to_m;      % Jupiter orbital radius [m]

v_earth_orbit   = sqrt(muSun/r_earth_orbit);    % Earth orbital velocity [m/s]
v_venus_orbit   = sqrt(muSun/r_venus_orbit);    % Venus orbital velocity [m/s]
v_jupiter_orbit = sqrt(muSun/r_jupiter_orbit);  % Jupiter orbital velocity [m/s]

%% Planetary Physical Data
R_earth   = 6378e3;         % Earth radius [m]
R_venus   = 6052e3;         % Venus radius [m]
R_jupiter = 71492e3;        % Jupiter radius [m]
r_GEO     = 42164e3;        % GEO radius [m]

%% Mission Inputs
rp_AU   = 0.14;                        % Solar perihelion [AU]
rf_AU   = 650;                          % Target distance [AU]
T_yr    = 25;                           % Cruise time [years]

%% Unit Conversions
rp      = rp_AU * AU_to_m;              % Perihelion [m]
rf      = rf_AU * AU_to_m;              % Target distance [m]
T       = T_yr * yr_to_sec;             % Cruise time [s]

%% Calculate Required v_infinity at Target
a_from_vinf   = @(vinf) -muSun./(vinf.^2);              % (Brown, Eq. 3.55)
e_from_vinf   = @(vinf) 1 + rp./(-a_from_vinf(vinf));   % (Hale, Eq. 2-8-9)
F_from_vinf   = @(vinf) acosh((rf./(-a_from_vinf(vinf)) + 1)./e_from_vinf(vinf));   % (Hale, Eq. 4-4-6)
tof_from_vinf = @(vinf) sqrt((-a_from_vinf(vinf)).^3/muSun) .* ...
                        (e_from_vinf(vinf).*sinh(F_from_vinf(vinf)) - F_from_vinf(vinf));   % (Hale, Eq. 4-4-3 and 4-4-4)

f = @(vinf) tof_from_vinf(vinf) - T;
vinf_required = fzero(f, rf/T);         % Required heliocentric v_inf [m/s]

v_escape_at_rp = sqrt(2*muSun/rp);      % (Brown, Eq. 3.15)
v_perihelion_required = sqrt(vinf_required^2 + v_escape_at_rp^2);   % (Brown, Eq. 3.57)

fprintf('========================================\n');
fprintf('MISSION REQUIREMENTS\n');
fprintf('========================================\n');
fprintf('Target: %.0f AU in %.0f years\n', rf_AU, T_yr);
fprintf('Solar perihelion: %.3f AU\n', rp_AU);
fprintf('Required v_inf: %.2f km/s\n', vinf_required/1e3);
fprintf('Required v at perihelion: %.2f km/s\n\n', v_perihelion_required/1e3);

%% Patched Conic 1: GEO to Earth Escape
fprintf('========================================\n');
fprintf('PATCHED CONIC 1: GEO TO EARTH ESCAPE\n');
fprintf('========================================\n');

% Initial assumption: we want ~10 km/s geocentric v_HE after Earth escape
% This will be refined based on trajectory optimization
v_HE_earth_departure = 10e3;           % geocentric v_HE leaving Earth [m/s]

% Circular velocity at GEO
v_circular_GEO = sqrt(muEarth/r_GEO);   % (Brown, Eq. 3.6)

% Velocity needed at GEO for hyperbolic escape with v_inf
v_escape_GEO = sqrt(v_HE_earth_departure^2 + 2*muEarth/r_GEO);    % (Brown, Eq. 3.15)

% Delta-V from GEO
dv_GEO = v_escape_GEO - v_circular_GEO;

fprintf('GEO altitude: %.0f km\n', (r_GEO - R_earth)/1e3);
fprintf('v_circular at GEO: %.2f km/s\n', v_circular_GEO/1e3);
fprintf('Desired geocentric v_inf, aka v_HE: %.2f km/s\n', v_HE_earth_departure/1e3);
fprintf('v needed at GEO: %.2f km/s\n', v_escape_GEO/1e3);
fprintf('Delta-V (from GEO): %.2f km/s\n\n', dv_GEO/1e3);

%% Patched Conic 2: Venus Flyby
fprintf('========================================\n');
fprintf('PATCHED CONIC 2: VENUS FLYBY\n');
fprintf('========================================\n');

% Approximate v_inf at Venus arrival
% This depends on the transfer orbit from Earth to Venus
v_inf_venus_arrival = 5.5e3;            % Approximate [m/s]

h_venus = 300e3;                        % Flyby altitude [m]
r_venus_periapsis = R_venus + h_venus;

% Unpowered gravity assist
% Calculate maximum deflection angle
delta_max_venus = 2 * asin(1/(1 + r_venus_periapsis*v_inf_venus_arrival^2/muVenus));

% After optimal deflection, heliocentric v_inf increases
v_inf_venus_departure = 9e3;            % Heliocentric v_inf after Venus [m/s]

fprintf('v_inf going in (Venus-relative): %.2f km/s\n', v_inf_venus_arrival/1e3);
fprintf('Periapsis altitude: %.0f km\n', h_venus/1e3);
fprintf('Max deflection angle: %.1f deg\n', rad2deg(delta_max_venus));
fprintf('Delta-V (unpowered GA): 0.00 km/s\n');
fprintf('v_inf out (heliocentric): %.2f km/s\n\n', v_inf_venus_departure/1e3);

%% Patched Conic 3: Earth Flyby
fprintf('========================================\n');
fprintf('PATCHED CONIC 3: EARTH FLYBY\n');
fprintf('========================================\n');

v_inf_earth_arrival = 6e3;              % Arrival v_inf at Earth [m/s]
h_earth = 300e3;                        % Flyby altitude [m]
r_earth_periapsis = R_earth + h_earth;

% Calculate maximum deflection angle
delta_max_earth = 2 * asin(1/(1 + r_earth_periapsis*v_inf_earth_arrival^2/muEarth));

% Unpowered gravity assist
v_inf_earth_departure_2 = 14e3;         % Heliocentric v_inf after Earth [m/s]

fprintf('v_inf going in (Earth-relative): %.2f km/s\n', v_inf_earth_arrival/1e3);
fprintf('Periapsis altitude: %.0f km\n', h_earth/1e3);
fprintf('Max deflection angle: %.1f deg\n', rad2deg(delta_max_earth));
fprintf('Delta-V (unpowered GA): 0.00 km/s\n');
fprintf('v_inf out (heliocentric): %.2f km/s\n\n', v_inf_earth_departure_2/1e3);

%% Patched Conic 4: Jupiter Powered Flyby
fprintf('========================================\n');
fprintf('PATCHED CONIC 4: JUPITER POWERED FLYBY\n');
fprintf('========================================\n');

% v_inf at Jupiter arrival from Earth-Jupiter transfer
v_inf_jupiter_arrival = 10e3;           % Jupiter-relative v_inf [m/s]

r_jupiter_periapsis = 3 * R_jupiter;    % Periapsis at 3 R_J [m]
dv_jupiter = 2e3;                       % Powered burn [m/s]

% Velocity at periapsis before burn
v_p_jupiter_before = sqrt(v_inf_jupiter_arrival^2 + 2*muJupiter/r_jupiter_periapsis);    % (Hale, Eq. 5-4-1)

% Powered gravity assist formula (Oberth effect)
v_inf_jupiter_departure = sqrt(v_inf_jupiter_arrival^2 + ...
    2*v_p_jupiter_before*dv_jupiter + dv_jupiter^2);

fprintf('v_inf in (Jupiter-relative): %.2f km/s\n', v_inf_jupiter_arrival/1e3);
fprintf('Periapsis: %.1f R_J (%.0f km)\n', r_jupiter_periapsis/R_jupiter, r_jupiter_periapsis/1e3);
fprintf('v at periapsis (before): %.2f km/s\n', v_p_jupiter_before/1e3);
fprintf('Delta-V (at periapsis): %.2f km/s\n', dv_jupiter/1e3);
fprintf('v_inf out (heliocentric): %.2f km/s\n\n', v_inf_jupiter_departure/1e3);

%% Patched Conic 5: Solar Oberth Maneuver
fprintf('========================================\n');
fprintf('PATCHED CONIC 5: SOLAR OBERTH\n');
fprintf('========================================\n');

% Velocity at perihelion before burn
v_perihelion_before = sqrt(v_inf_jupiter_departure^2 + v_escape_at_rp^2); % (Brown, Eq. 3.57)

% Required delta-v at solar perihelion
dv_solar_oberth = v_perihelion_required - v_perihelion_before;

fprintf('v_inf from Jupiter: %.2f km/s\n', v_inf_jupiter_departure/1e3);
fprintf('v_escape at %.3f AU: %.2f km/s\n', rp_AU, v_escape_at_rp/1e3);
fprintf('v at perihelion (before): %.2f km/s\n', v_perihelion_before/1e3);
fprintf('v at perihelion (required): %.2f km/s\n', v_perihelion_required/1e3);
fprintf('Delta-V (at perihelion): %.2f km/s\n\n', dv_solar_oberth/1e3);

%% Total Mission Delta-V Summary
fprintf('========================================\n');
fprintf('TOTAL MISSION DELTA-V\n');
fprintf('========================================\n');

total_dv = dv_GEO + dv_jupiter + dv_solar_oberth;

fprintf('GEO escape:         %.2f km/s\n', dv_GEO/1e3);
fprintf('Venus flyby:        %.2f km/s\n', 0.0);
fprintf('Earth flyby:        %.2f km/s\n', 0.0);
fprintf('Jupiter flyby:      %.2f km/s\n', dv_jupiter/1e3);
fprintf('Solar Oberth:       %.2f km/s\n', dv_solar_oberth/1e3);
fprintf('                    ----------------\n');
fprintf('TOTAL:              %.2f km/s\n\n', total_dv/1e3);

%% ========================================
%  FIGURE Delta-V Budget Bar Chart
%  Uses dv_GEO, dv_jupiter, dv_solar_oberth from main script
%  (all in [m/s])
%  =========================================

dv_labels = categorical({'GEO escape','Jupiter powered GA','Solar Oberth'});
dv_labels = reordercats(dv_labels, ...
    {'GEO escape','Jupiter powered GA','Solar Oberth'});

% Convert from m/s to km/s
dv_vals_kms = [dv_GEO, dv_jupiter, dv_solar_oberth] / 1e3;

figure;
bar(dv_labels, dv_vals_kms);
ylabel('\DeltaV [km/s]');
title(sprintf('Solar Gravity Lens Mission \\DeltaV Budget (r_p = %.3f AU)', rp_AU));
grid on;

%% ========================================
%  FIGURE Delta-V vs Perihelion Radius
%  r_p from 0.04 AU to 0.50 AU in steps of 0.01 AU
%  =========================================

% --- Jupiter powered flyby assumptions ---
v_inf_jup_arrival = 10e3;          % [m/s]
r_jup_periapsis   = 3 * R_jupiter; % [m]

% Velocity at Jupiter periapsis before burn
v_p_jup_before = sqrt(v_inf_jup_arrival^2 + 2*muJupiter/r_jup_periapsis);

% v_infinity leaving Jupiter (heliocentric), assumed independent of r_p
v_inf_jup_departure = sqrt(v_inf_jup_arrival^2 + ...
                           2*v_p_jup_before*dv_jupiter + dv_jupiter^2);

% Earth-departure Delta-V (constant term in total Delta-V)
dv_GEO_kms = dv_GEO/1e3;   % [km/s] from patch 1
dv_J_kms   = dv_jupiter/1e3;

% --- Sweep over perihelion radii ---
rp_vals_AU   = 0.04:0.01:0.50;   % [AU]
n            = numel(rp_vals_AU);
dv_sun_kms   = zeros(1,n);
dv_total_kms = zeros(1,n);

for i = 1:n
    rp_AU_i = rp_vals_AU(i);
    rp_i    = rp_AU_i * AU_to_m;   % [m]
    
    % Hyperbolic orbit parameters as functions of v_inf (for this rp_i)
    a_from_vinf = @(vinf) -muSun./(vinf.^2);                     % hyperbolic a
    e_from_vinf = @(vinf) 1 + rp_i./(-a_from_vinf(vinf));        % eccentricity
    F_from_vinf = @(vinf) acosh( (rf./(-a_from_vinf(vinf)) + 1) ./ ...
                                 e_from_vinf(vinf) );            % hyperbolic anomaly
    
    tof_from_vinf = @(vinf) sqrt((-a_from_vinf(vinf)).^3/muSun) .* ...
        ( e_from_vinf(vinf).*sinh(F_from_vinf(vinf)) - F_from_vinf(vinf) );
    
    % Solve for v_inf that yields the desired time of flight T
    f_vinf = @(vinf) tof_from_vinf(vinf) - T;
    % Initial guess ~ average radial speed (works well for this problem)
    vinf_required_i = fzero(f_vinf, rf/T);   % [m/s]
    
    % Speeds at perihelion
    v_escape_rp_i  = sqrt(2*muSun/rp_i);                         % escape at rp
    v_p_required_i = sqrt(vinf_required_i^2 + v_escape_rp_i^2);  % needed at rp
    v_p_before_i   = sqrt(v_inf_jup_departure^2 + v_escape_rp_i^2);
    
    dv_sun_i = v_p_required_i - v_p_before_i;   % [m/s]
    
    dv_sun_kms(i)   = dv_sun_i/1e3;
    dv_total_kms(i) = dv_GEO_kms + dv_J_kms + dv_sun_kms(i);
end

figure;
plot(rp_vals_AU, dv_total_kms, '-o'); hold on;
plot(rp_vals_AU, dv_sun_kms, '--o');
xlabel('Perihelion distance r_p [AU]');
ylabel('\DeltaV [km/s]');
legend('Total mission \DeltaV','Solar Oberth \DeltaV','Location','northwest');
title(sprintf('\\DeltaV vs. perihelion radius (650 AU in %d yr)', T_yr));
grid on;