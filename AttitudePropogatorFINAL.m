


% Attitude Propagator

%% MAIN ATTITUDE PROPAGATOR
% PRIMARY FUNCTION OF THIS SCRIPT
function [] = AttitudePropagator_modified()
clc
tic
    addpath('igrf')
    addpath('Aero_Model')
    addpath('math')
    format longg
    
    % Satellite moment of inertia matrix for 6U satellite
    mass = 10; % kg
    dims = [.353 .2263 .1]; % x,y,z dimensions of  cube in meters
    moi_matrix = (1/12*mass)*[(dims(2)^2+dims(3)^2)  0  0; 0 (dims(1)^2+dims(3)^2) 0; 0 0 (dims(1)^2+dims(2)^2)];

    % Wheel moment of inertia matrix
    wheel_moi = [.00000523 0 0; 0 .00000523 0; 0 0 .00000523];
    %keeps track of how many time the function containing magmoment has ran
    global magmoment_function_counter
    magmoment_function_counter = 0;
    %{
    for some reason that function calculating mag moment gets called 
    27 times before the actual good data actually comes in, these first 27
    rows might be used by other functions since the mag moment was like
    over 10^30 or something ridiculous, so i made the function to only
    start adding data to the 'time_vs_mag_moment_matrix' after 27th call,
    the previous counter variable just keeps track of how many time it was
    called. hope this helps. not sure if accurate. see the comment in 
    'quat_bi_der_ang_accel' to see how to see the original output
    %}
    global not_magmoment_calls;
    not_magmoment_calls = 27;
    global time_vs_mag_moment_matrix
    time_vs_mag_moment_matrix = [];
    
    global sat_props
    circ_alt = 550; % (km) circular orbit altitude
    earth_radius = 6378; % km
    semi_major = circ_alt + earth_radius;
    grav_param = 398600;
    orbit_vel = sqrt(grav_param/semi_major);
    orbit_vel/sqrt(2)
    init_pos_vel = [semi_major 0 0 0 orbit_vel/sqrt(2) orbit_vel/sqrt(2)]';
    % [quat0 quat1 quat2 quat3 sat_ang_vel(wx wy wz) wheel_ang_vel(x,y,z)]
    % Angular velocities are in radians/sec
    % All quaternions transform from the orbital to the body frame
    % Note that these are frame quaternions. rotating the orbital frame
    % about the eigenaxis by angle theta aligns it with the body frame
    % The first component in the quaternion is the scalar component. The
    % last three are the vector components
    init_quat_orb_body = [1 0 0 0]';
    init_quat_eci_body = get_eci_body_q(init_quat_orb_body, init_pos_vel);
    init_attitude = [init_quat_eci_body' .01 .4 .2 0 0 0]';
    % [rx ry rz vx vy vz] in the eci frame
    
    bi_init = [0 0 0 0]'; % magnetic field in booms due to hysteresis
    
    sat_props = struct('sat_moi', moi_matrix, 'torque_vect', [0 0 0]', 'mass', 10, ...
    'des_quat_orbital', [1 0 0 0]', 'cur_quat_ang_vel', init_attitude, ...
    'cur_pos_vel', init_pos_vel, 'date', datenum([2015, 12, 12]), ...
    'att_diff_eq', @magDetumbleDiffEq, 'orb_diff_eq', @orbitDiffEq, ...
    'wheel_moi_mat', wheel_moi, 'cur_bi',bi_init);
    
    % approx_ang_vels = calibrate_system(sat_props)
    max_ang_vel = 25; % maximum angular velocity magnitude that you can start steady state pointing at
    detumble_step = 500; % seconds before checking if the satellite is sufficiently detumbled
    tspan = [0 10000]; % Propagation time (seconds)
    ang_vel_mag = norm(init_attitude(5:7));
    options = odeset('RelTol',1e-13, 'MaxStep', 100); % error bound
    all_times = 0;
    
    % Holds position and orientation data
    all_states = [init_pos_vel' init_attitude' bi_init'];
    diffEq = @orbit_attitude_diff_eq;
    
    % Keep integrating with the detumble controller until rates are low
    % enough
%     while ang_vel_mag > max_ang_vel && tspan(1) < tspan(2)
%          [time state_vect] = (diffEq, [tspan(1) tspan(1)+detumble_step],...
%             [sat_props.cur_pos_vel; sat_props.cur_quat_ang_vel; sat_props.cur_bi], options);
%          
%          tspan(1) = tspan(1) + detumble_step
%          result_size = size(state_vect);
%          num_vals = result_size(1);
%          ang_vel_mag = norm(state_vect(num_vals, 11:13));
%          all_times = [all_times; time];
%          
%          all_states = [all_states; get_state_orb_q(state_vect)];
%     end
    

    modes = {@magDetumbleDiffEq, @attitudeDiffEq}; % attitude control modes stored in a cell array
    times_for_modes = [600, 1e-10]; % Times for each mode used to be 8 hours, 15.8 hours (10 orbits)
    for i = 1:length(modes)
        sat_props.att_diff_eq = modes{i}; % Attitude control mode
        steady_state_time = times_for_modes(i); % seconds that the steady state controller will need
        tspan(2) = tspan(1) + steady_state_time; % integrate for specified time
        [time, state_vect] = ode45(diffEq, tspan, [sat_props.cur_pos_vel; sat_props.cur_quat_ang_vel; sat_props.cur_bi], options); % call to integrator
        all_times = [all_times; time];
        
        % all_states has the eci position and velocity as the first 6 elemets,
        % the orbital to body quaternion as the next
        % 4 elements, the satellite angular velocities as the next three,
        % and the wheel angular velocities as the last 4
        all_states = [all_states; get_state_orb_q(state_vect)];
        tspan(1) = tspan(2);
    end
    
%     % Simulate b-dot with aero drag and rotating orbital frame
%     sim_time = 60000; % seconds to simulate
%     [all_times state_vect] = sim_rotating_orb_frame(sim_time);
%     all_quat_ang_vels = state_vect(:, 7:16);
    % Figure out the magnetorquer powers and torques during the simulation
    [mag_torques, total_mag_power, mag_field_body_list] = get_mag_perf(all_times, all_states, sat_props.date);
    
    data_points = size(all_states);
    % quat_der = get_quat_der(all_quat_ang_vels(1,:)')'
    err_quat = zeros(data_points(1),4);
    for i = 1:data_points(1)
        err_quat(i,:) = get_rel_quat(all_states(i,7:10)', sat_props.des_quat_orbital)';
    end
    
    matrix = [all_states(:,7:10), all_states(:,14:16)];
    %quat_der_ang_accel = momDumpDiffEq(all_times, )

    time_label = 'Time (hours)';
    % Plot results
    
    % Convert time to hours
    %
    all_times = all_times/3600;
    figure(1)
    set(gcf,'color','w');
    plot(all_times, all_states(:,7:10))
    xlabel(time_label)
    ylabel('quaternion (dimensionless)')
    legend('q0','q1','q2','q3')
    title('Quaternion Over Time (orbital to body)')
    %}
    %subplot(2,2,2)
    figure(2)
    set(gcf,'color','w');
    plot(all_times / 3600, all_states(:,11:13)/(2*pi)*60)
    title('Satellite Angular Velocity over Time')
    xlabel('Time (Hours)')
    ylabel('Satellite Angular Velocity (rpm)')
    legend('wx','wy','wz')
    %{
    %subplot(2,2,3)
    figure(3)
    set(gcf,'color','w');
    plot(all_times, all_states(:,14:16)/(2*pi)*60)
    title('Reaction Wheel Angular Velocity Over Time')
    xlabel(time_label)
    ylabel('Reaction wheel angular velocity (rpm)')
    legend('x wheel', 'y wheel', 'z wheel')
    %
    figure(4)
    set(gcf,'color','w');
    plot(all_times, err_quat(:, 2:4))
    title('Error Quaternion Over Time')
    xlabel(time_label)
    ylabel('Error Quaternion')
    legend('err quat 1', 'err quat 2', 'err quat 3')
    %}
    % Calculates the orientation of the satellite long axis vector relative
    % to the orbital frame
    long_axis_orbital = zeros(data_points(1),3);
    for i = 1:data_points(1)
        long_axis_orbital(i,:) = v_rot_q([0 0 1]', all_states(i, 7:10)')';
    end
    pointing_err_ang = acosd(long_axis_orbital(:,1));
    %{
    figure(5)
    subplot(2,2,1)
    plot(all_times, long_axis_orbital)
    title('Satellite Long Axis Orientation Over Time')
    xlabel(time_label)
    ylabel('Satellite Long Axis Vector in Orbital Frame')
    legend('x component', 'y component', 'z component')
    subplot(2,2,2)
    plot(all_times, pointing_err_ang)
    title('Pointing Error Over time')
    xlabel(time_label)
    ylabel('Pointing Error (deg)')
    subplot(2,2,3);
    plot(all_times, total_mag_power);
    title('total magnetorquer power consumption over time');
    xlabel(time_label)
    ylabel('Power (watts)');
    set(gcf,'color','w');
    subplot(2,2,4);
    set(gcf,'color','w');
    plot(all_times, all_states(:,11:13)/(2*pi)*60)
    title('Satellite Angular Velocity Over Time')
    xlabel(time_label)
    ylabel('Satellite angular velocity (rpm)')
    legend('wx','wy','wz')
    %}
    figure(6)
    plot(time_vs_mag_moment_matrix(:,1)/3600, time_vs_mag_moment_matrix(:,2:4));
    xlabel('Time (Hours)');
    ylabel('Desired Magnetic Moment (Am^2)')
    legend('µx','µy','µz');
    title('Desired Magnetic Moment over time')
    %
    figure(7)
    plot([0; all_times]/3600, mag_field_body_list');
    xlabel('Time (Hours)')
    ylabel('Magnetic Field (µT)')
    legend('Bx','By','Bz')
    title('Magnetic Field over Time')
    %}
    %disp(all_states(:,14:16))
runtime = toc
end

%% CONVERT ECI/BODY TO ORB/BODY
function state_vect_orb_q = get_state_orb_q(state_vect_eci_q)
    % Given a state vector (20 dimensional) that contains eci to body
    % quaternions, convert the quaternions to orbital to body and return
    % the state vector
    state_vec_size = size(state_vect_eci_q); % dimensions of states vector
    % Quaternions are currently eci to body. Represent them in orbital
    % to body for easier interpretation
    orb_body_quat = zeros(state_vec_size(1) ,4); % preallocate for speed
    for j = 1:state_vec_size(1)
        orb_body_quat(j,:) = get_orb_body_q(state_vect_eci_q(j,7:10)', state_vect_eci_q(j, 1:6)')';
    end
    
    state_vect_orb_q = [state_vect_eci_q(:, 1:6) orb_body_quat state_vect_eci_q(:, 11:end)];
end


%% CREATE STATE VECTOR
%{
This function integrates position in orbit and attitude at the same time.
The input is a column vector containing, position, velocity, attitude
quaternion, and angular velocity. The output state vector will contain
velocity, acceleration, quaternion derivative, and angular acceleration
%}
function state_vect_der = orbit_attitude_diff_eq(time, state_vect)
    global sat_props;
    pos_vel_vect = state_vect(1:6);
    quat_ang_vel_bi = state_vect(7:20);
    vel_acc_vect = sat_props.orb_diff_eq(time, pos_vel_vect);
    quat_bi_der_ang_accel = sat_props.att_diff_eq(time, quat_ang_vel_bi);
    % rate of change of magnetic hysteresis field strength
    
    state_vect_der = [vel_acc_vect; quat_bi_der_ang_accel];
end


%% ORBITAL DIFFERENTIAL EQUATION
function vel_acc_vect_orb = orbitDiffEq (time, pos_vel_vector)
% Orbit Propagator. Calculates orbital velocity and acceleration based on
% input position and velocity
    global sat_props
    gravity_param = 398600;
    pos_vect = pos_vel_vector(1:3);
    vel_vect = pos_vel_vector(4:6);
    sat_props.cur_pos_vel = pos_vel_vector;
    grav_accel_vect = (-gravity_param/((norm(pos_vect))^3)).*pos_vect;
    vel_acc_vect_orb = [vel_vect; grav_accel_vect];
end


%% MAGNETORQUER DETUMBLE DIFFERENTIAL EQUATION
%{
This function simulates the outputs of the B-Dot magnetorquer detumble
algorithm using the IGRF magnetic field model. Magnetic hysteresis
effects of booms are also taken into account
need to have quaternion from eci frame to body frame and angular velocity
of body frame relative to eci frame. Also need the current magnetic field
induced in the booms along each boom length vector due to past hysteresis effects
%}
function quat_bi_der_ang_accel = magDetumbleDiffEq(time, quat_ang_vel_bi)
    global sat_props
    panel_ang = 110; % panel angle in degrees relative to the cubesat sides
    wheel_ang_vel = quat_ang_vel_bi(8:10);
    bi = quat_ang_vel_bi(11:14);
    % Frame quaternion from eci frame to body frame
    pos_quat_eci = quat_ang_vel_bi(1:4)/norm(quat_ang_vel_bi(1:4));
    ang_vel = quat_ang_vel_bi(5:7); % angular velocity relative to the eci frame
    % Get the orientation quaternion from the orbital frame to the body
    % frame
    orb_body_quat = get_orb_body_q(pos_quat_eci, sat_props.cur_pos_vel);
    % must provide orbital to body quaternion
    mag_field_body = get_mag_field_body(sat_props.cur_pos_vel(1:3), sat_props.cur_pos_vel(4:6),...
        orb_body_quat, time, sat_props.date);
    [mag_torque, des_mag_moment] = get_bdot_torques(mag_field_body, ang_vel);
    global time_vs_mag_moment_matrix
    global not_magmoment_calls
    global magmoment_function_counter
    if (magmoment_function_counter > not_magmoment_calls)
        %remove the if condition to add all data into the matrix if you
        %wanna the first 27 skewed data and figure out what they were.
        time_vs_mag_moment_matrix = [time_vs_mag_moment_matrix;time,des_mag_moment'];   
    else
        magmoment_function_counter = magmoment_function_counter + 1;
    end
    
      
    % get magnetic moment and rate of change of magnetic field in booms due
    % to hysteresis
    % [bi_der, mag_mom_hyst] = MagneticDamping3(ang_vel, mag_field_body, panel_ang-90, bi, boom_lens, boom_area);
    % hyst_torque = cross(mag_mom_hyst, mag_field_body);
    bi_der = [0 0 0 0]'; % ignoring hysteresis torques for now
    earth_radius = 6378; %km
    orbit_pos = sat_props.cur_pos_vel(1:3);
    density_alt = norm(orbit_pos) - earth_radius; % Altitude for density calculation purposes
    [aero_force, aero_torque] = aer(density_alt, orb_body_quat, panel_ang, norm(sat_props.cur_pos_vel(4:6)));
    roll_torque = [0 0 0]'; % Parasitic roll torques that may be caused by boom imperfections
    torque = aero_torque + mag_torque +  roll_torque;
    ang_mom = sat_props.sat_moi*ang_vel + sat_props.wheel_moi_mat*wheel_ang_vel;
    ang_accel = calc_ang_accel(torque, ang_vel, ang_mom, sat_props.sat_moi);
    
    quat_bi_der_ang_accel = zeros(7,1);
    quat_bi_der_ang_accel(1:4) = get_quat_der(quat_ang_vel_bi);
    quat_bi_der_ang_accel(5:7) = ang_accel;
    quat_bi_der_ang_accel(8:10) = [0 0 0]';
    quat_bi_der_ang_accel(11:14) = bi_der;
    sat_props.cur_quat_ang_vel_bi = quat_ang_vel_bi;
    sat_props.cur_bi = bi;
    
    
    
   
end


%% WHEEL DETUMBLE ACCELERATIONS
%{
This function outputs the satellite angular accelerations that result
when the satellite is in detumble mode
This is a reaction wheel detumble mode
%}
function quat_der_ang_accel = detumbleDiffEq(time, quat_ang_vel)
    global sat_props;
    der_gain = .1;
    quat_der_ang_accel(1:4) = get_quat_der(quat_ang_vel);
    quat_der_ang_accel = zeros(7,1);
    % Get the quaternion derivative
    quat_der_ang_accel(1:4) = get_quat_der(quat_ang_vel);
    % The desired angular acceleration will be opposite the direction of
    % the velocity vector
    ang_vel = quat_ang_vel(5:7);
    des_ang_accel = -der_gain*ang_vel;
    quat_der_ang_accel(5:7) = des_ang_accel;
    sat_props.cur_quat_ang_vel = quat_der_ang_accel;
end


%% STEADY STATE CONTROL LOOP
%{
This function simulates the steady state control loop and returns the angular
acceleration of the satellite after control is applied
Quaternions are body frame to orbital frame
Inputs are position quaternion (orbital to body), satellite angular
velocity, and wheel angular velocities
%}
function quat_der_ang_accel = attitudeDiffEq(time, quat_ang_vel)
    global sat_props;
    der_gain = .1;
    prop_gain = (der_gain^2)/2; % Critically damped system
    wheel_ang_vel = quat_ang_vel(8:10);
    % current orbital to body quaternion
    orb_body_quat = get_orb_body_q(quat_ang_vel(1:4), sat_props.cur_pos_vel);
    
    des_quat_orbital = sat_props.des_quat_orbital; % desired quaternion relative to the orbital frame
    
    % Position quaternion relative to the maneuvering frame
    errQuat = get_rel_quat(orb_body_quat, des_quat_orbital); % error quaternion (from current to desired orientation)
    % keep rotation angle between zero and 180 to avoid sign flip
    if errQuat(1) < 0
        errQuat = -errQuat;
    end
    
    ang_vel = quat_ang_vel(5:7); % Angular velocity vector about the body axes
    
    quat_der_ang_accel = zeros(7,1);
    quat_der_eci = get_quat_der(quat_ang_vel);
    
    derivative_term = ang_vel;
    
    des_ang_accel = prop_gain*errQuat(2:4) - der_gain*derivative_term; % assume inverse dynamics pd gets it right
    ang_momentum = sat_props.sat_moi*ang_vel + sat_props.wheel_moi_mat*wheel_ang_vel;
    req_torques = calc_req_torques(des_ang_accel, ang_vel, ang_momentum, sat_props.sat_moi);
    % disturbances = get_aero_disturbance();
    wheel_ang_accel = get_wheel_ang_accel(req_torques);
    ang_accel = calc_ang_accel(req_torques, ang_vel, ang_momentum, sat_props.sat_moi);
    % quaternion derivative in the eci frame. Rate of change of quaternion
    % from eci to body frame
    quat_der_ang_accel(1:4) = quat_der_eci; 
    quat_der_ang_accel(5:7) = ang_accel;
    quat_der_ang_accel(8:10) = wheel_ang_accel;
    quat_der_ang_accel(11:14) = [0 0 0 0]'; % forget hysteresis torque on this one
    sat_props.cur_quat_ang_vel = quat_ang_vel;
end
    

%% WHEEL ACCELERATION FOR REQUIRED TORQUE
%{
This function calculates the wheel angular accelerations requried (x,y,z)
to get a desired torque vector on the satellite body. This is based on
the desired torque vector and the wheel moment of inertias
%}
function wheel_ang_accels = get_wheel_ang_accel(des_torque)
    global sat_props
    wheel_ang_accels = inv(sat_props.wheel_moi_mat)*des_torque*(-1);
end
    

%% FREE SATELLITE ROTATION
%{
This function simulates the satellite rotating freely with no control
torques
%}
function quat_der_ang_accel = free_rotate(time, quat_ang_vel)
    global sat_props
    moi_matrix = sat_props.sat_moi;
    ang_vel = quat_ang_vel(5:7); % Angular velocity vector about the body axes
    quat_der_ang_accel = zeros(7,1);
    angular_momentum = moi_matrix*ang_vel;   
    ang_accel = calc_ang_accel([0 0 0]', ang_vel, angular_momentum, moi_matrix);
    % Calculate quaternion derivative
    quat_der_ang_accel(1:4) = get_quat_der(quat_ang_vel(1:7));
    % Include angular acceleration in the derivative matrix
    quat_der_ang_accel(5:7) = ang_accel;
    quat_der_ang_accel(8:10) = [0 0 0]';
end


%% WHEEL DESATURATION MODE
%{
This is the differential equation associated with momentum dumping mode.
In this mode, the magnetorquers constantly act to reduce reaction wheel
angular momentum without changing satellite attitude.
quat_ang_vel = [q0 q1 q2 q3 wx wy wz wheel_x wheel_y wheel_z]'
%}
function quat_der_ang_accel = momDumpDiffEq(time, quat_ang_vel)
    global sat_props;
    pos_quat = quat_ang_vel(1:4);
    sat_ang_vel = quat_ang_vel(5:7);
    wheel_ang_vel = quat_ang_vel(8:10);
    
    desat_gain = 1e6; % Gain for the desaturation controller
    
    % Calculate total angular momentum of the wheels
    wheel_ang_mom = sat_props.wheel_moi_mat*wheel_ang_vel;
    
    % Get the magnetic field in the satellite body frame in Tesla
    orbit_pos = sat_props.cur_pos_vel(1:3); % km
    orbit_vel = sat_props.cur_pos_vel(4:6); % km/s
    mag_field_body = get_mag_field_body(orbit_pos, orbit_vel, pos_quat, time, sat_props.date);
   
    des_mag_mom = desat_gain*cross(wheel_ang_mom, mag_field_body);
    magnetic_torque = cross(des_mag_mom, mag_field_body);
    wheel_torque = -magnetic_torque;
    total_torque = magnetic_torque + wheel_torque;
    
    ang_mom = sat_props.sat_moi*sat_ang_vel + sat_props.wheel_moi_mat*wheel_ang_vel;
    sat_ang_accel = calc_ang_accel(total_torque, sat_ang_vel, ang_mom, sat_props.sat_moi);
    wheel_ang_accel = get_wheel_ang_accel(wheel_torque);
    
    quat_der_ang_accel = zeros(10,1);
    quat_der_ang_accel(1:4) = get_quat_der(quat_ang_vel);
    quat_der_ang_accel(5:7) = sat_ang_accel;
    quat_der_ang_accel(8:10) = wheel_ang_accel;
    sat_props.cur_quat_ang_vel = quat_ang_vel;
end


%% DETERMINE DART TORQUE
function aero_torque = get_aero_disturbance()
    alt = 550; % km
    fin_length = .330
    zfin_width = .20931
    yfin_width = .083
    area = (.2263*.1)+(fin_length*sin(radians(panel_ang))*zfin_width)+(fin_length*sin(radians(panel_ang))*yfin_width); % m
    velocity = 7588; % m/s
    cd=4;
    force = 1/2*velocity^2*cd*area*getDensity(alt);
    cm_cp_offset = .02;
    t_mag = cm_cp_offset*force/sqrt(3);
    aero_torque = [t_mag t_mag t_mag]';
end


%% ESTIMATE ANGULAR VELOCITY
function [ang_vel_est] = calibrate_system(sat_props)
    sample_torque = [.01 .02 .015]'; % n*m
    sat_moi = sat_props.sat_moi;
    wheel_moi = sat_props.wheel_moi_mat;
    true_ang_vel = sat_props.cur_quat_ang_vel(5:7);
    wheel_ang_vel = sat_props.cur_quat_ang_vel(5:7);
    ang_momentum = sat_moi*true_ang_vel + wheel_moi*wheel_ang_vel
    ang_accel = calc_ang_accel(sample_torque, true_ang_vel, ang_momentum, sat_moi)
    wheel_ang_accel = get_wheel_ang_accel(sample_torque);
    
    % Calculate z angular velocity based on known x angular velocity
    function omega_z = get_omega_z(omega_x) 
        num = wheel_moi(1,1)*omega_x*wheel_ang_vel(3) - wheel_moi(1,1)*wheel_ang_accel(2) ...
            - sat_moi(2,2)*ang_accel(2);
        denom = (sat_moi(1,1) - sat_moi(3,3))*omega_x + wheel_moi(1,1)*wheel_ang_vel(1);
        omega_z = num/denom;
    end
    
    % Calculate y angular velocity based on known x angular velocity
    function omega_y = get_omega_y(omega_x)
        num = wheel_moi(1,1)*get_omega_z(omega_x)*wheel_ang_vel(2) - wheel_moi(1,1)*wheel_ang_accel(1) ...
            - sat_moi(1,1)*ang_accel(1);
        denom = (sat_moi(3,3) - sat_moi(2,2))*get_omega_z(omega_x) + wheel_moi(1,1)*wheel_ang_vel(3);
        omega_y = num/denom;
    end

    % Error should equal zero when the correct x angular velocity is
    % provided to this function. Use this function as an argument to fzero
    function error = root_func(omega_x)
        num = wheel_moi(1,1)*get_omega_y(omega_x)*wheel_ang_vel(1) - wheel_moi(1,1)*wheel_ang_accel(3) ...
            - sat_moi(3,3)*ang_accel(3);
        denom = (sat_moi(2,2) - sat_moi(1,1))*get_omega_y(omega_x) + wheel_moi(1,1)*wheel_ang_vel(3);
        error = omega_x - num/denom;
    end
        
    % Calculate the estimated angular velocity
    options = optimset('TolX', 1e-50)
    ang_vel_est = zeros(3,1);    
    ang_vel_est(1) = fzero(@root_func, true_ang_vel(1)*1.1, options);
    ang_vel_est(2) = get_omega_y(ang_vel_est(1));
    ang_vel_est(3) = get_omega_z(ang_vel_est(1));
     
end

%% end of program