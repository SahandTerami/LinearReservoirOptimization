using DelimitedFiles
using JuMP, Ipopt
using LinearAlgebra
using Distributions

data = readdlm("params.txt", ',')
params = vec(data)

K = Int(params[1])       # Total Number of Input Frequencies
N = Int(params[2])       # Total Number of Reservoir's Nodes

β1 = params[3]      # Regularization parameter
β2 = params[4]      # Harmonic mean Term weight
γ  = params[5]

a = Vector{Float64}(undef, K)
b = Vector{Float64}(undef, K)
ϕ = Vector{Float64}(undef, K)
ω = Vector{Float64}(undef, K)
c_i = Vector{Float64}(undef, N)
c_t = Vector{Float64}(undef, N)
c_k = Vector{Float64}(undef, N)

a  = params[6:6+K-1]
b  = params[6+K:6+K+K-1]
ϕ  = params[6+K+K:6+K+K+K-1]
ω  = params[6+K+K+K:6+K+K+K+K-1]

c_i = params[6+K+K+K+K:6+K+K+K+K+N-1]
c_k = params[6+K+K+K+K+N:6+K+K+K+K+N+N-1]
c_t = params[6+K+K+K+K+N+N:6+K+K+K+K+N+N+N-1]

#=
println("K = ", K)
println("N = ", N)
println("β1 = ", β1)
println("β2 = ", β2)
println("γ = ", γ)
println("a = ", a)
println("b = ", b)
println("ϕ = ", ϕ)
println("ω = ", ω)
println("c = ", c_i)
=#

m = Model(Ipopt.Optimizer)
register(m, :sum, 1, sum; autodiff = true)
set_attribute(m, "print_info_string", "yes")
@variable(m,-20<= λ[1:N] <= 0) # Reservoir's eigenvalues
@variable(m, κ[1:N]) # Reservoir's output weights

@variable(m, 0 <=M[1:N, 1:K]) #Transfer function magnitude
@variable(m, -pi<= θ[1:N, 1:K] <= pi) #Transfer function phase

@variable(m, ϵ_cos[1:K]) #Estimation error
@variable(m, ϵ_sin[1:K]) #Estimation error

@variable(m, H>= 10^-3) #Harmonic mean
Ω = zeros(2 * K, N)  # Initialize Ω as a matrix of size 2K by N
ϵ_cos0_vals = zeros(K)
ϵ_sin0_vals = zeros(K)

    for j = 1:K
        cϵ_cos0 = 0
        cϵ_sin0 = 0

        for i = 1:N
            c = c_t[i]
            set_start_value(λ[i], c)
            
            ck = c_k[i]
            set_start_value(κ[i], ck)

            cM = cM = (γ * c_i[i]) / sqrt((γ^2)*(c - 1)^2 + ω[j]^2)
            set_start_value(M[i,j], cM)

            cθ = atan(ω[j] / (γ * (c - 1)))
            set_start_value(θ[i,j], cθ)

            cϵ_cos0 = cM*cos(cθ)*ck + cϵ_cos0
            cϵ_sin0 = cM*sin(cθ)*ck + cϵ_sin0
        end
        cϵ_cos = a[j]*cϵ_cos0
        set_start_value(ϵ_cos[j], cϵ_cos)

        cϵ_sin = a[j]*cϵ_sin0
        set_start_value(ϵ_sin[j], cϵ_sin)

        ϵ_cos0_vals[j] = cϵ_cos0
        ϵ_sin0_vals[j] = cϵ_sin0
    end

w = Vector{Float64}(undef, K)
for i in 1:K
    w[i] = 1/ω[i];
end
oversum1 = sqrt(sum(w[i] * (ϵ_cos0_vals[i]^2 + ϵ_sin0_vals[i]^2) for i in 1:K))
#println("Init error = ", oversum1)
###################################
# Set the parameters of the optimization
set_optimizer_attribute(m, "max_iter", 5000)            # Set the maximum number of iterations the solver is allowed to run  
set_optimizer_attribute(m, "tol", 1e-8)                 # Set the main convergence tolerance 
set_optimizer_attribute(m, "acceptable_tol", 1e-8)      # Set the acceptable convergence tolerance   
set_optimizer_attribute(m, "constr_viol_tol", 1e-8)     # Set the maximum allowed constraint violation   
set_optimizer_attribute(m, "dual_inf_tol", 1e-8)        # Set the tolerance for dual infeasibility   
set_optimizer_attribute(m, "compl_inf_tol", 1e-8)       # Set the tolerance for complementarity conditions  
set_optimizer_attribute(m, "acceptable_constr_viol_tol", 1e-8)      # Set the acceptable constraint violation tolerance 
###################################
# Evaluate the objective function before optimization
@NLobjective(m, Min, (sum(w[i]*(ϵ_cos[i]^2 + ϵ_sin[i]^2) for i in 1:K)) + sum(β1*(κ[j]*κ[j]) for j in 1:N)+ β2*(1/H)); # Weighted Objective
###################################
# Define the nonlinear constraints
@NLconstraint(m, con1[j in 1:K], ϵ_cos[j] == a[j]*sum(M[i,j]*cos(θ[i,j])*κ[i] for i in 1:N) - b[j]*cos(ϕ[j]));
@NLconstraint(m, con2[j in 1:K], ϵ_sin[j] == a[j]*sum(M[i,j]*sin(θ[i,j])*κ[i] for i in 1:N) - b[j]*sin(ϕ[j]));
@NLconstraint(m, con3[i in 1:N, j in 1:K], M[i,j]^2 * ((γ^2)*(λ[i]-1)^2 +     ω[j]^2) - (γ^2)*c_i[i]^2 == 0);       # Magnitude constraint for c / (jω - λ + 1)
@NLconstraint(m, con4[i in 1:N, j in 1:K], γ*tan(θ[i,j])*(λ[i]-1) - ω[j] == 0);     # Phase constraint for c / (jω - λ + 1)
@NLconstraint(m, con5[i in 1:N, j in i+1:N], H*sum(1/abs((λ[i]-λ[j]))) <= (N*(N-1)/2) + 10^-6);     # Harmonic mean Constraint
@NLconstraint(m, con7[i in 1:N], λ[i]*γ-γ+5<=10^-6);    
###################################
optimize!(m)
###################################

# Print output
oversum = sqrt(sum(w[i] * (value(ϵ_cos[i])^2 + value(ϵ_sin[i])^2) for i in 1:K))
println("\n")
println("Error in the Frequency Domain: ", oversum)     # Print error
println("\n")

println("Eigenvalues (λ) = ", value.(λ),"\n")       # Print optimal eigenvalues
println("Readout weights (κ) = ", value.(κ),"\n")   # Print optimal readout wieghts

println("Magnitude Matrix (M) = ")          # Print the magnitude matrix corresponding to the optimal eigenvalues
show(stdout, "text/plain", value.(M))
println("\n")

println("Phase Angle Matrix (θ) = ")        # Print the phase angle matrix corresponding to the optimal eigenvalues
show(stdout, "text/plain", value.(θ))
println("\n end")

# Save values
writedlm("λ.csv", value.(λ), ',')
writedlm("κ.csv", value.(κ), ',')
writedlm("M.csv", value.(M), ',')
writedlm("θ.csv", value.(θ), ',')

# Optional: save error too
open("error.txt", "w") do io
    write(io, string(oversum))
end


