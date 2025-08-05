function q1 = generate_input(t, a, M, omega, theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  24/07/2025  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Developer: Sahand Tangerami %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function Generates The State Matrix of The Reservoir Computer. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                              %
%       t:          Total Simulation Time Vector                      %
%       a:          Amplitudes of Input Signals                       %
%       M:          Reservoir State Magnitude Matrix                  %
%       omega:      Input Signal Frequencies                          %
%       theta:      Reservoir State Phase Matrix                      %
% Output:                                                             %
%       q1:         State Matrix of The Reservoir Computer            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    [N, K] = size(M);            % Get number of neurons N and frequency components K
    T = length(t);              % Number of time steps
    u = zeros(N, T);            % Initialize output

    for i = 1:N
        for k = 1:K
            u(i, :) = u(i, :) + a(k) * M(i, k) * cos(omega(k) * t + theta(i, k));
        end
    end

    q1 = u';  % Transpose to get T x N
end