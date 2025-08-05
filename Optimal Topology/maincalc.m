function [overalldata, t_train, t_test, y_train, y_test, Omega_train, Omega_test, Omega_train_init, Omega_test_init, optkappa, kappa_init] = maincalc(K, N, beta1, beta2, gamma, T_total, T_trans, T_train, dt, a, b, omega, phi,  multi_num)
    delete('*.txt');
    c = param(N)';          % Input weights(Decoupled Reservoir)
    T_test = T_total - T_trans - T_train;
    Final_t = T_total*dt;           % Final Time

    % Calculations
    % Time Vector Creation
    t_train = (T_trans+1)*dt:dt:(T_trans + T_train)*dt;
    t_test = (T_trans + T_train + 1)*dt:dt:(T_total)*dt;
    t = linspace(0, Final_t, T_total);              % Time Vector
    
    len_trans = T_trans;            % Training Phase Length 
    len_train = T_train;            % Test Phase Length 
    old_error_value = 100;

    for ii = 1:multi_num
        delete('*.csv');
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
        oplambda = readmatrix('λ.csv');
        opkappa = readmatrix('κ.csv');
        opM = readmatrix('M.csv');
        optheta = readmatrix('θ.csv');
        operror_value = readmatrix('error.csv');
        
        if operror_value<old_error_value
            optlambda = oplambda;
            optkappa = opkappa;
            optM = opM;
            opttheta = optheta;
            opterror_value = operror_value;
            old_error_value = opterror_value;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
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
    
    % Save Data and Plot
    overalldata = [Initial_Normal_Error_train           Normal_Error_train          Initial_NRMSE_train         NRMSE_train
                   Initial_Normal_Error_test            Normal_Error_test           Initial_NRMSE_test          NRMSE_test];

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
end