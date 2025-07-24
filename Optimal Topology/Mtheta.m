function [M,theta] = Mtheta(N,K,gamma,c,omega,lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  24/07/2025  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Developer: Sahand Tangerami %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function Generates The Reservoir State Magnitude and Phase     %
% Matrices.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                              %
%       N:          Total Number of Reservoir Nodes                   %
%       K:          Total Number of Input Frequencies                 %
%       gamma:      Reservoir Computer's Constant                     %
%       c:          Input Mask                                        %
%       omega:      Input Signal Frequencies                          %
%       lambda:     Eigenvalues of The Reservoir Computer             %
% Output:                                                             %
%       M:          Reservoir State Magnitude Matrix                  %
%       theta:      Reservoir State Phase Matrix                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:N
        for k = 1:K
            % Compute M_ik (magnitude)
            M(i, k) = abs(gamma * c(i) / (1j * omega(k) - gamma * (lambda(i) - 1)));

            % Compute theta_ik (phase)
            theta(i, k) = angle(gamma * c(i) / (1j * omega(k) - gamma * (lambda(i) - 1)));
            if theta(i,k)>0
                theta(i,k)=theta(i,k)-pi;
            end
        end
    end
end