
using JuMP, Ipopt
using LinearAlgebra
using Distributions
#########System Parameters

K = 3 # Total Number of Input Frequencies
N = 10 # Total Number of Reservoir's Nodes

a = Vector{Float64}(undef, K)
b = Vector{Float64}(undef, K)
ϕ = Vector{Float64}(undef, K)
ω = Vector{Float64}(undef, K)
c_i = Vector{Float64}(undef, N)
c_t = Vector{Float64}(undef, N)
c_k = Vector{Float64}(undef, N)

β1 = 10^-6      # Regularization parameter
β2 = 10^-1      # Harmonic mean Term weight
γ = 6

a[1], a[2], a[3] = 1.1, 1.7, 2.1    # Input signal coefficients
b[1], b[2], b[3] = 2.2, 1.0, 1.6    # Output signal coefficients
ϕ[1], ϕ[2], ϕ[3] = -0.5, 0.9, 1.1   # Phase shifts
ω[1], ω[2], ω[3] = 1, 3, 5          # Frequencies
mode = "manual"    # Set "manual" if you want to optimize your linear reservoir computer
                    # Set "auto" if you want to generate an optimized linear reservoir computer

# Decoupled reservoir's constant (c_i in eq. 12)
c_i[1], c_i[2], c_i[3], c_i[4], c_i[5], c_i[6], c_i[7], c_i[8], c_i[9], c_i[10]  = -1.99996560298497, -0.376357788431904, -2.86186350280611, 1.09105941750842, 1.36749055891047,
                                                                                   -1.02253798457437, 0.39469504442572 , 0.628226445668294, -1.65984541973612, -1.04932811954560 
# Initial guess of eigenvalues (Use your eigenvalues)
c_t[1], c_t[2], c_t[3], c_t[4], c_t[5], c_t[6], c_t[7], c_t[8], c_t[9], c_t[10]  = -1.9492762713584446e-6, -9.241228381769005, -0.13111971127125596, -1.6634201608952646,
                                                                                   -0.6426750348914327, -1.5313450034202654, -1.3961526699233624, -9.440358211745474,                                                                                   -0.5114473241311479, -1.795644654160857
# Initial guess of output weights (Use your output wieghts)
c_k[1], c_k[2], c_k[3], c_k[4], c_k[5], c_k[6], c_k[7], c_k[8], c_k[9], c_k[10]  = 118.29028024543335, 80.33561344451377, -146.6316735527552, -45.533238640181416, 80.05327262026651,
                                                                                   -41.40889456129366, -13.92479555692042, 134.11444974864088, 102.02275656816647, -41.7437310857018

m = Model(Ipopt.Optimizer)
set_attribute(m, "print_info_string", "yes")
@variable(m,-20<= λ[1:N] <= 0) # Reservoir's eigenvalues
@variable(m, κ[1:N]) # Reservoir's output weights

@variable(m, 0 <=M[1:N, 1:K]) #Transfer function magnitude
@variable(m, -pi<= θ[1:N, 1:K] <= pi) #Transfer function phase

@variable(m, ϵ_cos[1:K]) #Estimation error
@variable(m, ϵ_sin[1:K]) #Estimation error

@variable(m, H>= 10^-3) #Harmonic mean
Ω = zeros(2 * K, N)  # Initialize Ω as a matrix of size 2K by N

if mode == "manual"
    for j = 1:K
        cϵ_cos0 = 0
        cϵ_sin0 = 0

        for i = 1:N
            c = c_t[i]
            set_start_value(λ[i], c)
            
            ck = c_k[i]
            set_start_value(κ[i], ck)

            cM = c_i[i]*sqrt(((1-c)/(ω[j]^2+c^2))^2+(ω[j]/(ω[j]^2+c^2))^2)
            set_start_value(M[i,j], cM)

            cθ = atan((-ω[j])/(1-c))
            set_start_value(θ[i,j], cθ)

            cϵ_cos0 = cM*cos(cθ)*ck + cϵ_cos0
            cϵ_sin0 = cM*sin(cθ)*ck + cϵ_sin0
        end
        cϵ_cos = a[j]*cϵ_cos0
        set_start_value(ϵ_cos[j], cϵ_cos)

        cϵ_sin = a[j]*cϵ_sin0
        set_start_value(ϵ_sin[j], cϵ_sin)
    end
elseif mode == "auto"
    for j = 1:K
        cϵ_cos0 = 0
        cϵ_sin0 = 0

        for i = 1:N
            c = -1/(1*i)
            set_start_value(λ[i], c)
            
            ck = -5*i
            set_start_value(κ[i], ck)

            cM = c_i[i]*sqrt(((1-c)/(ω[j]^2+c^2))^2+(ω[j]/(ω[j]^2+c^2))^2)
            set_start_value(M[i,j], cM)

            cθ = atan((-ω[j])/(1-c))
            set_start_value(θ[i,j], cθ)

            cϵ_cos0 = cM*cos(cθ)*ck + cϵ_cos0
            cϵ_sin0 = cM*sin(cθ)*ck + cϵ_sin0
        end
        cϵ_cos = a[j]*cϵ_cos0
        set_start_value(ϵ_cos[j], cϵ_cos)

        cϵ_sin = a[j]*cϵ_sin0
        set_start_value(ϵ_sin[j], cϵ_sin)
    end
end

w = Vector{Float64}(undef, K)
for i in 1:K
    w[i] = 1/ω[i];
end

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

