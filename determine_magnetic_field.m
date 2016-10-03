clc, clear all
format bank
tic %Initiate runtime counter

% ***INPUTS***
sim_time = 283; %minutes in simulation
    %95.63 = 1 orbit
time_step = 10; %seconds between data sample
    %DO NOT SET BELOW 0.01 LEST YOU WISH TO MEET THE KRAKEN
orb_alt = 550; %mean altitude of orbit, km
init_pos = [0 0]; %initial position, [latitude longitude]
inclin = 0; %orbital inclination, degrees
    %0 < inclin < 89 means the satellite rotates with earth
    %91 < inclin < 179 means the satellite rotates against earth
eccent = 0; %eccentricity of orbit, 0 is circular
RA_ascend = 0; %right ascension of ascending node
arg_per   = 0; %argument of perigee
true_anom = 0;  %true anomaly of departure

% ***CONSTANTS***
g_parameter = 3.986E14; %graviational parameter of earth, m^3 / s^2
earth_radius = 6378; %radius of earth, km
semi_major = orb_alt + earth_radius; %semi_major axis height, km 
B_field_orb = []; %Define an empty matrix to store orbital components
secs_per_day = 86164; %seconds in a solar day

% ***CALCULATE***
orb_vel = sqrt(g_parameter/(semi_major * 1000)); %orbital velocity, m/s
init_quat = [semi_major 0 0 0]'; %initial measured quaternion
init_sat_ang_vel = orb_vel/sqrt(2); %initial satellite angular velocity
init_wheel_ang_vel = orb_vel/sqrt(2); %initial reaction wheel angular vel.
period = (2*pi*(semi_major * 1000)) / orb_vel;
num_orb = (sim_time * 60) / period;
fprintf('Simulating %2.2f Orbits \n', num_orb)
time = datenum(datetime);
time_of_day = mod(time, secs_per_day); % seconds in current day
gha = (time_of_day)/secs_per_day*360; %greenwich hour angle, degrees
% Find the latitude and longitude from orbit specs
[lat, long] = Orbit3D(RA_ascend, arg_per, true_anom, inclin, semi_major,...
    eccent, time_step, num_orb);

%{
%Calculate GMSTo
    JD = juliandate(datetime('today'));
    JD0 = NaN(size(JD)); %julian date of previous midnight
    JDmin = floor(JD)-.5;
    JDmax = floor(JD)+.5;
    JD0(JD > JDmin) = JDmin(JD > JDmin);
    JD0(JD > JDmax) = JDmax(JD > JDmax);
    H = (JD-JD0).*24;       %Time in hours past previous midnight
    D = JD - 2451545.0;     %Compute the number of days since J2000
    D0 = JD0 - 2451545.0;   %Compute the number of days since J2000
    T = D./36525;           %Compute the number of centuries since J2000
    %Calculate GMST in hours (0h to 24h), convert to degrees
GMSTo = mod(6.697374558 + 0.06570982441908.*D0  + 1.00273790935.*H + ...
    0.000026.*(T.^2),24).*15;
%}

% Call IGRF to create a matrix of magnetic field values in the NED frame
[Bx_NED, By_NED, Bz_NED] = igrf(time, lat, long, orb_alt);
B_field_NED = [Bx_NED, By_NED, Bz_NED];

% Convert magnetic field vector from NED to orbital frame, row by row
NED_matrix_row = 1; %initialize a counter for the following for loop

% Creates a matrix of mag field components in orbital frame
for i = 0:(time_step/60):sim_time
    % Pass the value of i in seconds to the ned_to_orb function
    ned_orb_matrix = ned_to_orb(inclin, period, i*60);
    % Gives the individual vectors for each timestep, and adds them to
    % larger matrix
    rotated_B_vec = ned_orb_matrix * B_field_NED(NED_matrix_row,:)';
    B_field_orb = [B_field_orb; rotated_B_vec'];
    NED_matrix_row = NED_matrix_row + 1;
end

% NED frame component graph
total_time = (0:(time_step / 60):sim_time)';
figure()
plot(total_time, B_field_NED)
axis([0 sim_time -3.5*10^4 3.5*10^4])
title('Magnetic Field (NED frame)')
legend('Bx', 'By', 'Bz')
xlabel('time, minutes')
ylabel('nano-Tesla')

% Orbital frame component graph
figure()
plot(total_time, B_field_orb)
axis([0 sim_time -3.5*10^4 3.5*10^4])
title('Magnetic Field (orbital frame)')
legend('Bx', 'By', 'Bz')
xlabel('time, minutes')
ylabel('nano-Tesla')

runtime = toc; %end runtime counter
fprintf('Total Simulation Time = %4.2f \n', runtime)

%NED to ECI dcm
%vern_ang = long + gha; %angle between vernal equinox and current point
%dcm = [-sind(lat)*cosd(vern_ang) -sind(lat)*sind(vern_ang)  cosd(lat);
%       -sind(vern_ang)            cosd(vern_ang)                    0;
%       -cosd(lat)*cos(vern_ang)  -cosd(lat)*sind(vern_ang) -sind(lat)];...
%direction cosine matrix