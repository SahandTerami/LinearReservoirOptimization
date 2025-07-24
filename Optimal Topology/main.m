clc
clear
close all

delete('*.txt');
delete('*.csv');
%% Parameters
K = 3;    % Total Number of Input Frequencies
N = 30;   % Total Number of Reservoir's Nodes

beta1 = 10^-8;   % Regularization parameter
beta2 = 10^-1;   % Harmonic mean Term weight
gamma = 6;       % REservoir's Constant Parameter

a = [1.1, 1.7, 2.1];    % Input signal coefficients
b = [2.2, 1.0, 1.6];    % Output signal coefficients
phi = [-0.5, 0.9, 1.1]; % Phase shifts
omega = [1, 3, 5];      % Frequencies
c = param(N)';          % Input weights(Decoupled Reservoir)

% Time vector for signal evaluation (for a specific time period)
T_total = 5000;                 % Total Time Steps
T_trans = 1000;                 % Transient Steps
T_train = 3000;                 % Training Time Steps
T_test = T_total - T_trans - T_train;

dt = 0.01;                      % Time Step Length
Final_t = T_total*dt;           % Final Time
len_trans = T_trans;            % Training Phase Length 
len_train = T_train;            % Test Phase Length 

%% Calculations

% Time Vector Creation
t_train = (T_trans+1)*dt:dt:(T_trans + T_train)*dt;
t_test = (T_trans + T_train + 1)*dt:dt:(T_total)*dt;
t = linspace(0, Final_t, T_total);              % Time Vector

% Initial Eigenvalues
[lambda_init] = initlambda(N);  

% Magnitude and Phase Matrices
[M_init,theta_init] = Mtheta(N,K,gamma,c,omega,lambda_init);

% Reservoir States Matrices Before The Optimization
Omega_init = generate_input(t, a, M_init, omega, theta_init);
Omega_train_init = Omega_init(len_trans+1:len_trans+len_train,:);
Omega_test_init = Omega_init(len_trans+len_train+1:end,:);

% Target Signal
[y_train,y_test] = signal_creation(omega,t, b, phi, len_trans, len_train);

% Initial Output Weights
[kappa_init] = initkappa(Omega_train_init, beta1, y_train); 

% Initial Normalized Training Error
rc_train_init = Omega_train_init*kappa_init';
Initial_Normal_Error_train = (1/sqrt(T_train))*norm(rc_train_init-y_train);
Initial_NRMSE_train = norm(rc_train_init-y_train)/norm(y_train);

% Initial Normalized Test Error
rc_test_init = Omega_test_init*kappa_init';
Initial_Normal_Error_test = (1/sqrt(T_test))*norm(rc_test_init-y_test);
Initial_NRMSE_test = norm(rc_test_init-y_test)/norm(y_test);

%%%%%%% Optimization %%%%%%%
[params] = optim(K, N, beta1, beta2, gamma, a, b, phi, omega, c, kappa_init, lambda_init);

% Run Julia script
% Make sure Julia is in your system Path, or provide full path
system('julia_Optimization.jl');

% Read The Optimal Values Resulted From Julia
optlambda = readmatrix('λ.csv');
optkappa = readmatrix('κ.csv');
optM = readmatrix('M.csv');
opttheta = readmatrix('θ.csv');
opterror_value = str2double(fileread('error.txt'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reservoir States Matrix After The Optimization
Omega = generate_input(t, a, optM, omega, opttheta);
Omega_train = Omega(len_trans+1:len_trans+len_train,:);
Omega_test = Omega(len_trans+len_train+1:end,:);

% Optimal Normalized Training Error
rc_train = Omega_train*optkappa;
Normal_Error_train = (1/sqrt(T_train))*norm(rc_train-y_train);
NRMSE_train = norm(rc_train-y_train)/norm(y_train);

% Optimal Normalized Test Error
rc_test = Omega_test*optkappa;
Normal_Error_test = (1/sqrt(T_test))*norm(rc_test-y_test);
NRMSE_test = norm(rc_test-y_test)/norm(y_test);
%% Save Data and Plot

data = [Initial_Normal_Error_train           Initial_NRMSE_train
        Initial_Normal_Error_test            Initial_NRMSE_test];

rowHeaders = {'Training Phase (Initial RC)', 'Testing Phase (Initial RC)'};          % Row descriptions
colHeaders = {'', 'Normalized Error (Eq.20-21)', 'NRMSE'};   % First empty, then column names
output = cell(size(data,1)+1, size(data,2)+1);
output(1, :) = colHeaders;
output(2:end, 1) = rowHeaders';
output(2:end, 2:end) = num2cell(data);

data2 = [Normal_Error_train           NRMSE_train
        Normal_Error_test             NRMSE_test];

rowHeaders = {'Training Phase (Optimal RC)', 'Testing Phase (Optimal RC)'};          % Row descriptions
colHeaders = {'', 'Normalized Error (Eq.20-21)', 'NRMSE'};   % First empty, then column names
output2 = cell(size(data2,1)+1, size(data2,2)+1);
output2(1, :) = colHeaders;
output2(2:end, 1) = rowHeaders';
output2(2:end, 2:end) = num2cell(data2);

% disp('Original Frequency Problem Before the Fourier Transform')
disp(output)
disp('')
% disp('Reduced Frequency Problem After the Fourier Transform')
disp(output2)
disp('')
% Plot
figure(1)
hold on;
plot(t_train, Omega_train * optkappa, '-', 'DisplayName', 'Optimal Reservoir (Training)', 'MarkerSize', 10, 'LineWidth', 4);
plot(t_train, (Omega_train_init * kappa_init'), '-.g', 'DisplayName', 'Non-Optimal Reservoir (Training)', 'MarkerSize', 10, 'LineWidth', 4);
xlabel('Time (t)');
ylabel('Amplitude');
plot(t_train, y_train, '-r', 'DisplayName', 'Output Data', 'LineWidth', 4);  % Actual output signal
box on;
set(gca, 'fontsize', 20);
legend('show', 'NumColumns', 3);

figure(2)
hold on;
plot(t_test, Omega_test * optkappa, '-', 'DisplayName', 'Optimal Reservoir (Training)', 'MarkerSize', 10, 'LineWidth', 4);
plot(t_test, (Omega_test_init * kappa_init'), '-.g', 'DisplayName', 'Non-Optimal Reservoir (Training)', 'MarkerSize', 10, 'LineWidth', 4);
xlabel('Time (t)');
ylabel('Amplitude');
plot(t_test, y_test, '-r', 'DisplayName', 'Output Data', 'LineWidth', 4);  % Actual output signal
box on;
set(gca, 'fontsize', 20);
legend('show', 'NumColumns', 3);