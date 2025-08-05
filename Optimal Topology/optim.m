function [params] = optim(K, N, beta1, beta2, gamma, a, b, phi, omega, c, kappa_init, lambda_init)

        % Define parameters for Julia
        params = [K, N, beta1, beta2, gamma, a, b, phi, omega, c, kappa_init, lambda_init];
        
        % Save to file to being read in Julia
        writematrix(params, 'params.txt');
end