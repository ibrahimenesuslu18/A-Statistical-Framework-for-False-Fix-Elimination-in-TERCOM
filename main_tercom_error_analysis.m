clear;clc;close all;
% --- Synthetic Simulation File for the Manuscript Titled "" --- %
% Initialize the terrain server

lla0 = [38.43, 41.37, 2500]; % initial LLA

terrain         = Terrain(5,0,".\\terrainServer", ...
    ".\\terrainServer\\server.dll"); % object constructor
terrain.init(lla0(1), lla0(2)); % Call init 

% Define system parameters
dt = 0.05; % cycle time

% Continious time point mass model with no angular moments
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0]; % System matrix
B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1]; % Input matrix

C = eye(6); % Dummy observation matrix

sys = ss(A, B, C, []); % Continous system model

sys_d = c2d(sys, dt, "zoh"); % Discritize the system

% Discrete model
A_d = sys_d.A; % Discrete system matrix 
B_d = sys_d.B; % Discrete input matrix

% Simulation parameters allocation and setup
time = 0:dt:100; % Time vector allocation
N = numel(time); % Number of iteratations

% Noise vector allocations
mu_INS              = zeros(6, 1); % Initial noise mean 
mu_IMU              = zeros(3, 1); % IMU drift noise mean
mu_motion           = zeros(3, 1); % Acceleration noise mean

sigma_INS           = blkdiag(100 * eye(3), 1e0 * eye(3)); % Initial noise covariance 
sigma_IMU           = blkdiag(1e-6 * eye(2), 1e-6); % IMU drift noise covariance
sigma_motion        = blkdiag(1e1 * eye(2), 1e-1); % Acceleration noise covariance

w_initial_noise     = mvnrnd(mu_INS, sigma_INS, 1)'; % Initial noise vector
w_IMU_drift_noise   = mvnrnd(mu_IMU, sigma_IMU, N)'; % IMU drift noise array
w_motion_noise      = mvnrnd(mu_motion, sigma_motion, N-1)'; % Acceleration noise array

% State and noise parameter allocations
x_truth             = zeros(6, N); % Ground truth state vector
x_INS               = zeros(6, N); % INS state vector

lla_truth           = zeros(N, 3); % Ground truth LLA(lat lon alt)
lla_INS             = zeros(N, 3); % INS LLA(lat lon alt)

u_IMU               = zeros(3, N-1); % IMU input
w_IMU_drift         = zeros(3, N); % IMU Drift input (Brown motion)

% initialize the states and construct input
x_truth(:, 1)       = [0; 0; 0;150; 155; 0]; % Ground truth state vector
x_INS(:, 1)         = x_truth(:, 1) + w_initial_noise; % initializa the INS state vector

lla_truth(1, :) = ned2lla(x_truth(1:3, 1)', lla0, "ellipsoid");
lla_INS(1, :) = ned2lla(x_INS(1:3, 1)', lla0, "ellipsoid");

u = [0.5; 0.5; 0.001]; % Input vector (constant acceleration)

radalt_meas = zeros(1, N); % Radalt measurement array

% Initial height above terrain
hat_hot_struct          = terrain.hatHot(lla_truth(1, 1), lla_truth(1, 2), lla_truth(1, 3), 1);
radalt_meas(1)          = hat_hot_struct.heightAboveTerrain;

% --- Simulation Loop ---%
for k = 1:N-1
    % if sqrt(sum(x_truth(4:5, k).^2)) >= 250 
    %     u = [0; 0; 0];
    % end
    x_truth(:, k + 1)       = A_d * x_truth(:, k) + B_d * (u + w_motion_noise(:, k)); % state transition equation
    lla_truth(k + 1, :)     = ned2lla(x_truth(1:3, k + 1)', lla0, "ellipsoid"); % ground truth LLA

    u_IMU(:, k)             = u + w_motion_noise(:, k) + w_IMU_drift(:, k) + w_IMU_drift_noise(:, k); % IMU input
    x_INS(:, k + 1)         = A_d * x_INS(:, k) + B_d * u_IMU(:, k); % INS state transition equation
    lla_INS(k + 1, :)       = ned2lla(x_INS(1:3, k + 1)', lla0, "ellipsoid"); % INS LLA

    % Synthetic radalt data generation
    hat_hot_struct          = terrain.hatHot(lla_truth(k + 1, 1), lla_truth(k + 1, 2), lla_truth(k + 1, 3), 1); % request terrain server
    radalt_meas(k)          = hat_hot_struct.heightAboveTerrain; % Pull the height above terrain response
    if mod(k, 40) == 0
        calculate_tercom_pos(pred_pos_INS, terrain, sigma_INS, radalt_meas);
    end
    
    w_IMU_drift(:, k + 1)   =  w_IMU_drift(:, k) + dt * w_IMU_drift_noise(:, k); % Integrate the IMU drift noise component
end

% --- Plotting ---
abs_err = abs(x_truth - x_INS);
L2_norm_err = sqrt(sum((x_truth(1:2, :) - x_INS(1:2, :)).^2));
figure();
plot(x_truth(2, :), x_truth(1, :), "LineWidth", 4, "LineStyle", "--");grid on;hold on;
plot(x_INS(2, :), x_INS(1, :), "LineWidth", 4);
xlabel("East Position (m)");ylabel("North Position (m)"); title("INS vs Ground truth");
legend("Ground Truth", "INS");

figure();
subplot(3, 1, 1);
plot(time, abs_err(1, :), "LineWidth", 2);grid on;title("North Error");
xlabel("time");ylabel("Northward Absolute Error (m)");
subplot(3, 1, 2);
plot(time, abs_err(2, :), "LineWidth", 2);grid on;title("East Error");
xlabel("time");ylabel("Northward Absolute Error (m)");
subplot(3, 1, 3);
plot(time, L2_norm_err, "LineWidth", 2);grid on;title("Euclidian Error");
xlabel("time");ylabel("L2 Norm Error (m)");
sgtitle("Error graphs over time");

figure();
subplot(3, 1, 1);
plot(time, w_IMU_drift(1, :), "LineWidth", 2);
subplot(3, 1, 2);
plot(time, w_IMU_drift(2, :), "LineWidth", 2);
subplot(3, 1, 3);
plot(time, w_IMU_drift(3, :), "LineWidth", 2);
