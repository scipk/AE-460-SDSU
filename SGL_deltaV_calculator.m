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

%% Planetary Data
R_earth   = 6378e3;         % Earth radius [m]
R_venus   = 6052e3;         % Venus radius [m]
R_jupiter = 71492e3;        % Jupiter radius [m]
r_GEO     = 42164e3;        % GEO radius [m]

%% Mission Inputs
rp_AU   = 0.046;                        % Solar perihelion [AU]
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
v_perihelion_required = sqrt(vinf_required^2 + v_escape_at_rp^2);

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

% We'll determine the required v_inf based on the trajectory
% For E-V-E-J sequence, we need approximately 10 km/s v_inf from Earth
v_inf_earth_escape = 10e3;              % Desired v_inf at Earth escape [m/s]

% Circular velocity at GEO
v_circular_GEO = sqrt(muEarth/r_GEO); % (Brown, Eq. 3.6)

% Velocity needed at GEO for hyperbolic escape with v_inf
v_escape_GEO = sqrt(v_inf_earth_escape^2 + 2*muEarth/r_GEO);    % (Brown, Eq. 3.15)

% Delta-V from GEO
dv_GEO = v_escape_GEO - v_circular_GEO;

fprintf('GEO altitude: %.0f km\n', (r_GEO - R_earth)/1e3);
fprintf('v_circular at GEO: %.2f km/s\n', v_circular_GEO/1e3);
fprintf('Desired v_inf from Earth: %.2f km/s\n', v_inf_earth_escape/1e3);
fprintf('v needed at GEO: %.2f km/s\n', v_escape_GEO/1e3);
fprintf('Delta-V (from GEO): %.2f km/s\n\n', dv_GEO/1e3);

%% Patched Conic 2: Venus Flyby
fprintf('========================================\n');
fprintf('PATCHED CONIC 2: VENUS FLYBY\n');
fprintf('========================================\n');

v_inf_venus_in = 5.5e3;                 % Arrival v_inf at Venus [m/s]
h_venus = 300e3;                        % Flyby altitude [m]
r_venus_periapsis = R_venus + h_venus;

% Unpowered gravity assist - magnitude preserved
v_inf_venus_out_helio = 9e3;            % Heliocentric v_inf after Venus [m/s]

fprintf('v_inf in (Venus-relative): %.2f km/s\n', v_inf_venus_in/1e3);
fprintf('Periapsis altitude: %.0f km\n', h_venus/1e3);
fprintf('Delta-V (unpowered GA): 0.00 km/s\n');
fprintf('v_inf out (heliocentric): %.2f km/s\n\n', v_inf_venus_out_helio/1e3);

%% Patched Conic 3: Earth Flyby
fprintf('========================================\n');
fprintf('PATCHED CONIC 3: EARTH FLYBY\n');
fprintf('========================================\n');

v_inf_earth_in = 6e3;                   % Arrival v_inf at Earth [m/s]
h_earth = 300e3;                        % Flyby altitude [m]
r_earth_periapsis = R_earth + h_earth;

% Unpowered gravity assist
v_inf_earth_out_helio = 14e3;           % Heliocentric v_inf after Earth [m/s]

fprintf('v_inf in (Earth-relative): %.2f km/s\n', v_inf_earth_in/1e3);
fprintf('Periapsis altitude: %.0f km\n', h_earth/1e3);
fprintf('Delta-V (unpowered GA): 0.00 km/s\n');
fprintf('v_inf out (heliocentric): %.2f km/s\n\n', v_inf_earth_out_helio/1e3);

%% Patched Conic 4: Jupiter Powered Flyby
fprintf('========================================\n');
fprintf('PATCHED CONIC 4: JUPITER POWERED FLYBY\n');
fprintf('========================================\n');

v_inf_jupiter_in = 10e3;                % Arrival v_inf at Jupiter [m/s]
r_jupiter_periapsis = 3 * R_jupiter;    % Periapsis at 3 R_J [m]
dv_jupiter = 2e3;                       % Powered burn [m/s]

% Velocity at periapsis before burn
v_p_jupiter_before = sqrt(v_inf_jupiter_in^2 + 2*muJupiter/r_jupiter_periapsis);    % (Hale, Eq. 5-4-1)

% Powered gravity assist formula (Oberth effect)
v_inf_jupiter_out = sqrt(v_inf_jupiter_in^2 + 2*v_p_jupiter_before*dv_jupiter + dv_jupiter^2);

fprintf('v_inf in (Jupiter-relative): %.2f km/s\n', v_inf_jupiter_in/1e3);
fprintf('Periapsis: %.1f R_J (%.0f km)\n', r_jupiter_periapsis/R_jupiter, r_jupiter_periapsis/1e3);
fprintf('v at periapsis (before): %.2f km/s\n', v_p_jupiter_before/1e3);
fprintf('Delta-V (at periapsis): %.2f km/s\n', dv_jupiter/1e3);
fprintf('v_inf out (heliocentric): %.2f km/s\n\n', v_inf_jupiter_out/1e3);

%% Patched Conic 5: Solar Oberth Maneuver
fprintf('========================================\n');
fprintf('PATCHED CONIC 5: SOLAR OBERTH\n');
fprintf('========================================\n');

% Velocity at perihelion before burn
v_perihelion_before = sqrt(v_inf_jupiter_out^2 + v_escape_at_rp^2); % (Brown, Eq. 3.57)

% Required delta-v at solar perihelion
dv_solar_oberth = v_perihelion_required - v_perihelion_before;

fprintf('v_inf from Jupiter: %.2f km/s\n', v_inf_jupiter_out/1e3);
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