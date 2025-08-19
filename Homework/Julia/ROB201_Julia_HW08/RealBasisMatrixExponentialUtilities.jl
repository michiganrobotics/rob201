#=
Utility functions for solving linear ODE systems dx/dt = Ax, x(0) = x0.

Uses eigenvalue/eigenvector decomposition and construction via a real basis
to directly compute the time-domain solution x(t).

Final version: Removes aggressive final simplification to preserve solution structure.

Requires: SymPy, PyCall, LinearAlgebra

Current Location Context: Ann Arbor, MI, USA
Current Date/Time Context: Friday, March 28, 2025 at 8:06 AM
=#

using SymPy
using PyCall
using LinearAlgebra

# Import Python's SymPy via PyCall if not already done globally
# sympy = pyimport("sympy")

"""
    cleanUp(A::AbstractArray{<:Number}, tol::Real=1e-10) -> AbstractArray

Sets the real and imaginary parts of elements in array `A` to zero if their
absolute value is below the tolerance `tol`. Used for cleaning numerical results.
"""
function cleanUp(A::AbstractArray{<:Number}, tol::Real=1e-10)
    B = copy(A) # Work on a copy
    for i in eachindex(B)
        real_part = abs(real(B[i])) < tol ? 0.0 : real(B[i])
        imag_part = abs(imag(B[i])) < tol ? 0.0 : imag(B[i])
        B[i] = real_part + imag_part * im
    end
    return B
end


"""
    solve_ode_real_basis(A_in::AbstractMatrix{<:Real}, x0_in::AbstractVector{<:Real}, t::Sym; tol::Real=1e-8)

Solves the linear system dx/dt = Ax, x(0) = x0 using eigenvalue/eigenvector
decomposition and construction via a real basis. Directly computes x(t).

# Arguments
- `A_in`: The system matrix (n x n real matrix).
- `x0_in`: The initial condition vector (n x 1 real vector).
- `t`: SymPy symbol for the time variable (must be real=True).
- `tol`: Numerical tolerance for cleanup and checking if eigenvalues are real.

# Returns
- A `NamedTuple` `(x_t=...)` containing the symbolic time-domain solution vector `x_t`.
- Returns `nothing` if errors occur (e.g., basis matrix non-invertible).
"""
function solve_ode_real_basis(A_in::AbstractMatrix{<:Real}, x0_in::AbstractVector{<:Real}, t::Sym; tol::Real=1e-8)
    n = length(x0_in)
    if size(A_in) != (n, n)
        error("Matrix A dimension ($(size(A))) must match length of x0 ($(n))")
    end
    if !t.is_real
         @warn "Time symbol 't' should be created with real=True for correct symbolic functions (exp, cos, sin)."
    end

    println("Using real basis method derived from eigenvalues/eigenvectors...")

    A = Float64.(copy(A_in))
    x0 = Float64.(copy(x0_in))

    # --- Step 1: Eigenvalue Decomposition ---
    local eigen_decomp
    try
        eigen_decomp = eigen(A)
    catch e
        @error "Eigenvalue decomposition failed for matrix A: $A. Error: $e"
        return nothing
    end
    eigen_vals_raw = eigen_decomp.values
    eigen_vecs_raw = eigen_decomp.vectors

    eigen_vals = cleanUp(eigen_vals_raw, tol)
    eigen_vecs = cleanUp(eigen_vecs_raw, tol)
    #println("Cleaned Eigenvalues: ", eigen_vals, "\n")

    # --- Steps 2 & 3: Build Real Basis P and Store Eigenvalue Info ---
    basis_vecs = Vector{Float64}[]
    real_eig_info = [] # (lambda, vector_index_in_P)
    complex_eig_info = [] # (a, omega_positive, vector_R_index_in_P, vector_I_index_in_P)
    processed_indices = Set{Int}()
    current_basis_idx = 0

    for i in 1:n
        if i in processed_indices continue end

        λ = eigen_vals[i]
        v = eigen_vecs[:, i]

        if abs(imag(λ)) < tol # Real Case
            λ_r = real(λ)
            v_r = real.(v)
            push!(basis_vecs, v_r)
            current_basis_idx += 1
            push!(real_eig_info, (λ_r, current_basis_idx))
            push!(processed_indices, i)
        else # Complex Case
             conj_idx = -1
             for j in 1:n
                 if i != j && !(j in processed_indices) && abs(eigen_vals[j] - conj(λ)) < tol*10
                     conj_idx = j
                     break
                 end
             end

             if conj_idx == -1
                 @error "Could not find distinct conjugate partner for eigenvalue $λ (Index $i). Aborting."
                 return nothing
             end

             if imag(λ) < 0 # Ensure we use the one with positive imaginary part
                 λ = eigen_vals[conj_idx]
                 v = eigen_vecs[:, conj_idx]
             end

             a = real(λ)
             ω = imag(λ) # ω > 0
             v_R = real.(v)
             v_I = imag.(v)

             push!(basis_vecs, v_R); current_basis_idx += 1; r_idx = current_basis_idx
             push!(basis_vecs, v_I); current_basis_idx += 1; i_idx = current_basis_idx

             push!(complex_eig_info, (a, ω, r_idx, i_idx))
             push!(processed_indices, i)
             push!(processed_indices, conj_idx)
        end
    end

    if length(basis_vecs) != n
         @error "Failed to construct a full basis. Found $(length(basis_vecs)) basis vectors for dimension $n."
         return nothing
    end

    P = hcat(basis_vecs...)

    # --- Step 4: Find Initial Condition Coefficients α ---
    local α_vec
    try
        println("Solving V*α = x0 for coefficients α ...")
        α_vec = P \ x0
        α_vec = cleanUp(α_vec, tol)
        println("Coefficients α ≈ ")
        display(α_vec)
        println("\n")
    catch e
        @error "Failed to solve P*c = x0. Matrix P may be singular or ill-conditioned. Error: $e"
        return nothing
    end

    # --- Step 5: Construct Solution x(t) ---
    println("Constructing symbolic solution x(t)...")
    x_t_sym = Sym.(zeros(n))

    try
        # Process real eigenvalue components
        for (λ_r, basis_idx) in real_eig_info
            c_r = α_vec[basis_idx]
            if abs(c_r) < tol continue end
            v_r = basis_vecs[basis_idx]
            term = c_r * sympy.exp(λ_r * t) * Sym.(v_r)
            x_t_sym += term
        end

        # Process complex eigenvalue components
        for (a, ω, r_idx, i_idx) in complex_eig_info
            d = α_vec[r_idx]
            e = α_vec[i_idx]
            if abs(d) < tol && abs(e) < tol continue end # Skip if mode not excited

            v_R = basis_vecs[r_idx]
            v_I = basis_vecs[i_idx]
            v_R_sym = Sym.(v_R)
            v_I_sym = Sym.(v_I)

            term = sympy.exp(a*t) * (
                     (d*sympy.cos(ω*t) + e*sympy.sin(ω*t))*v_R_sym +
                     (-d*sympy.sin(ω*t) + e*sympy.cos(ω*t))*v_I_sym
                   )
            x_t_sym += term
        end
    catch e
        @error "Failed during symbolic construction of x(t). Error: $e"
        return nothing
    end

    # --- Final Steps ---
    println("Applying minimal cleanup (real)...")
    try
        # Ensure result is real
        x_t_sym = real.(x_t_sym) # Use Julia real broadcast

        # --- AGGRESSIVE SIMPLIFY REMOVED ---
        # x_t_sym = sympy.simplify.(x_t_sym) # REMOVED

        println("Solution construction complete.")
    catch e
        @warn "Minimal cleanup (real.) failed for x(t). Error: $e. Returning raw sum."
    end

    # Return the symbolic sum (without aggressive simplification)
    return (x_t=x_t_sym,)

end


# === Helper Functions for Creating Callable Solutions ===

"""
    generate_function(symbolic_expr::Sym, time_var::Sym)

Converts SymPy expression to callable Julia function using lambdify. Returns NaN on error.
"""
function generate_function(symbolic_expr::Sym, time_var::Sym)
    local lambda_func
    try
        lambda_func = sympy.lambdify(time_var, symbolic_expr, modules="numpy", cse=true)
    catch e1
        try
            @warn "Lambdify with numpy failed: $e1. Falling back to math."
            lambda_func = sympy.lambdify(time_var, symbolic_expr, modules="math", cse=true)
        catch e2
            @error "Lambdify failed: $e2. Returning NaN function."
            return t_val -> NaN
        end
    end

    function julia_func(t_val)
         input = isa(t_val, AbstractArray) ? t_val : [t_val]
         local result, numeric_result
         try
             result = lambda_func(input)
             numeric_result = PyAny(result)
             numeric_result = real.(numeric_result)
             return isa(t_val, AbstractArray) ? Float64.(numeric_result) : Float64(first(numeric_result))
         catch e
            @warn "Evaluation/Conversion failed in generated function: $e. Returning NaN."
            return isa(t_val, AbstractArray) ? fill(NaN, size(t_val)) : NaN
         end
     end
    return julia_func
end

"""
    mySolutionFactory(solution_vector_t::AbstractVector{<:Sym}, time_var::Sym)

Creates function evaluating symbolic solution vector at time points. Returns matrix (rows=vars, cols=times).
"""
function mySolutionFactory(solution_vector_t::AbstractVector{<:Sym}, time_var::Sym)
    solution_functions = map(solution_vector_t) do sol
        try generate_function(sol, time_var)
        catch e; @error "Failed generation for $sol: $e"; t_val -> NaN end
    end
    num_vars = length(solution_functions)

    function solutionMatrix(t_vec::AbstractVector{<:Real})
        t_vals = collect(t_vec)
        num_times = length(t_vals)
        result_matrix = Matrix{Float64}(undef, num_vars, num_times)
        for i in 1:num_vars
             try result_matrix[i, :] = solution_functions[i](t_vals)
             catch e; @warn "Eval component $i failed: $e"; result_matrix[i, :] .= NaN end
        end
        return result_matrix
    end
    return solutionMatrix
end