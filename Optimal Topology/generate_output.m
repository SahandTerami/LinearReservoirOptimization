function y = generate_output(t, b, omega, phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  24/07/2025  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Developer: Sahand Tangerami %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function Generates The Output Signal as a Vector.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                              %
%       t:          Total Simulation Time Vector                      %
%       b:          Amplitudes of Output Signals                      %
%       omega:      Input Signal Frequencies                          %
%       phi:        Output Signal Phase Change                        %
% Output:                                                             %
%       y:         Vector of Output Signal                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    N = length(b);
    y = zeros(size(t));
    
    for i = 1:N
        y = y + b(i) * cos(omega(i)*t + phi(i));
    end
    
    y = y'; % Return as column vector
end