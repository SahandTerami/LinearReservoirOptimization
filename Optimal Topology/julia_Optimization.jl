using DelimitedFiles, JuMP, Ipopt, LinearAlgebra, Distributions, MathOptInterface
const MOI = MathOptInterface

# Read params.txt
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



w = 1 ./ ω  # Optimization Weights

# Function of Building Model
function build_model()
    m = Model(Ipopt.Optimizer)

    # Set the parameters of the optimization
    set_optimizer_attribute(m, "max_iter", 10000)            # Set the maximum number of iterations the solver is allowed to run  
    set_optimizer_attribute(m, "tol", 1e-10)                 # Set the main convergence tolerance 
    set_optimizer_attribute(m, "acceptable_tol", 1e-10)      # Set the acceptable convergence tolerance   
    set_optimizer_attribute(m, "constr_viol_tol", 1e-10)     # Set the maximum allowed constraint violation 
    set_optimizer_attribute(m, "dual_inf_tol", 1e-10)        # Set the tolerance for dual infeasibility  
    set_optimizer_attribute(m, "compl_inf_tol", 1e-10)       # Set the tolerance for complementarity conditions
    set_optimizer_attribute(m, "acceptable_constr_viol_tol", 1e-10)      # Set the acceptable constraint violation tolerance 
    set_optimizer_attribute(m, "print_level", 5)

    @variable(m,-20 <= λ[1:N] <= 0)  # Eigenvalues
    @variable(m, κ[1:N])             # Readout weights

    @variable(m, 0 <= M[1:N, 1:K])   # Transfer magnitude
    @variable(m, -pi <= θ[1:N, 1:K] <= pi) # Transfer phase

    @variable(m, ϵ_cos[1:K]) # Error cosine part
    @variable(m, ϵ_sin[1:K]) # Error sine part

    @variable(m, H >= 1e-3)  # Harmonic mean

    # Evaluate the objective function before optimization
    @NLobjective(m, Min, sum(w[i]*(ϵ_cos[i]^2 + ϵ_sin[i]^2) for i=1:K) + sum(β1*κ[j]^2 for j=1:N) + β2*(1/H))   # Weighted Objective

    @NLconstraint(m, [j=1:K], ϵ_cos[j] == a[j]*sum(M[i,j]*cos(θ[i,j])*κ[i] for i=1:N) - b[j]*cos(ϕ[j]))
    @NLconstraint(m, [j=1:K], ϵ_sin[j] == a[j]*sum(M[i,j]*sin(θ[i,j])*κ[i] for i=1:N) - b[j]*sin(ϕ[j]))
    @NLconstraint(m, [i=1:N, j=1:K], M[i,j]^2 * ((γ^2)*(λ[i]-1)^2 + ω[j]^2) - (γ^2)*c_i[i]^2 == 0)         # Magnitude constraint for c / (jω - λ + 1)
    @NLconstraint(m, [i=1:N, j=1:K], γ*tan(θ[i,j])*(λ[i]-1) - ω[j] == 0)                                   # Phase constraint for c / (jω - λ + 1)
    @NLconstraint(m, sum(1/abs(λ[i]-λ[j]) for i=1:N-1, j=i+1:N) * H <= (N*(N-1)/2) + 1e-6)                 # Harmonic mean Constraint
    @NLconstraint(m, [i=1:N], λ[i]*γ - γ + 5 <= 1e-6)

    return m, λ, κ, M, θ, ϵ_cos, ϵ_sin, H
end

# Function to set initial values
function set_start_values!(m, λ, κ, M, θ, λ_start, κ_start)
    for i in 1:N
        set_start_value(λ[i], λ_start[i])
        set_start_value(κ[i], κ_start[i])
    end
    for j in 1:K
        for i in 1:N
            c = λ_start[i]
            cM = c_i[i]*sqrt(((1-c)/(ω[j]^2 + c^2))^2 + (ω[j]/(ω[j]^2 + c^2))^2)
            set_start_value(M[i,j], cM)
            cθ = atan((-ω[j])/(1-c))
            set_start_value(θ[i,j], cθ)
        end
    end
end

# Run optimization and return results
function run_optimization!(m, λ, κ, M, θ, ϵ_cos, ϵ_sin, H, λ_start, κ_start)
    set_start_values!(m, λ, κ, M, θ, λ_start, κ_start)
    optimize!(m)

    status = termination_status(m)
    if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
        oversum = sqrt(sum(w[i]*(value(ϵ_cos[i])^2 + value(ϵ_sin[i])^2) for i in 1:K))
        return oversum, value.(λ), value.(κ), M, θ
    else
        println("Optimization failed with status: ", status)
        return Inf, λ_start, κ_start, M, θ
    end
end

# Optimization loop with iterative initial values
function iterative_optimization(num_iter, λ_init, κ_init)
    m, λ, κ, M, θ, ϵ_cos, ϵ_sin, H = build_model()

    λ_current = copy(λ_init)
    κ_current = copy(κ_init)
    best_cost = Inf
    best_λ = λ_current
    best_κ = κ_current
    best_M = M
    best_θ = θ

    for iter in 1:num_iter
        println("\n=== Iteration $iter ===")
        obj, λ_new, κ_new, M_new, θ_new = run_optimization!(m, λ, κ, M, θ, ϵ_cos, ϵ_sin, H, λ_current, κ_current)
        println("Objective: ", obj)

        λ_current = λ_new
        κ_current = κ_new

        if obj < best_cost
            best_cost = obj
            best_λ = λ_new
            best_κ = κ_new
            best_M = value.(M_new)
            best_θ = value.(θ_new)
        end
    end

    return best_λ, best_κ, best_cost, best_M, best_θ
end

# Initial value from parameter file
num_iterations = 50
final_λ, final_κ, final_cost, final_M, final_θ = iterative_optimization(num_iterations, c_t, c_k)

println("\nBest optimized λ: ", final_λ)        # Print optimal eigenvalues
println("Best optimized κ: ", final_κ)          # Print optimal readout wieghts
println("\nBest Error: ", final_cost)           # Print error
println("\nM: ", final_M)                       # Print the magnitude matrix corresponding to the optimal eigenvalues
println("\nθ: ", final_θ)                       # Print the phase angle matrix corresponding to the optimal eigenvalues


# Save results to file
writedlm("λ.csv", final_λ, ',')
writedlm("κ.csv", final_κ, ',')
writedlm("error.csv", final_cost, ',')
writedlm("M.csv", final_M, ',')
writedlm("θ.csv", final_θ, ',')
