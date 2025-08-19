# Compute the Laplace tranform of dx/dt = Ax, x(0) = x0 and the inverse Laplace Transform
function Laplace(AA, x00, s, t, aTol=1e-8)
    x0 = copy(x00)
    A = copy(AA)
    
    # makes sure A and x0 are real valued
    x0 = real.(x0) 
    A = real.(A)
    
    # Obtain e-values and e-vectors   
    F = eigen(A)
    V = cleanUp(F.vectors, aTol)
    lambda = cleanUp(F.values, aTol)
        
    # Check e-values are distinct
    distinct = are_eigenvalues_distinct_within_tolerance(F.values, 100*aTol)
    if distinct == false
        println("Eigenvalues must be distinct for this method to work")
        return false
    end
        
    # Transform coordinates
    n = length(x0) 
    inV = V \ I(n)
    #display(inV)
    inV = cleanUp(inV, aTol)
    x0Bar = inV*x0; x0Bar = cleanUp(x0Bar, aTol)
    #Abar = inV*A*V same as Lambda below

    # Matrix in new coordinates is diagonal with e-values on the diagonal
    Lambda = diagm(lambda)

    # Define 's' as a complex symbol and 't' as a real symbol
    s = sympy.symbols("s")
    t = sympy.symbols("t", real=True)
    
    # Compute the Lapalce Transform in the diagonalizing coordinates
   X_sBar =  ( 1.0 ./(s .- lambda) ) .* x0Bar
    
    # Place in original coordinates
    V = cleanUp(V, aTol)
    X_s = V * X_sBar    
    X_s = sympy.simplify.(X_s)
  
    
    # Compute the inverse Laplace Transform of the transformed model
   sol_t =  X_sBar.applyfunc(decompose_and_inverse_laplace)
    
    # Solution back in the original coordinates    
    if false
        x_t = SymPy.simplify.(V * sol_t)
        x_t = SymPy.real.(x_t) 
    else
        x_t = SymPy.real.(V * sol_t)  
        x_t = SymPy.simplify.(x_t)    
    end
    x_t = SymPy.simplify.(x_t)
    
    
    # CleanUp answers
    X_s = cleanUpExpression.(X_s, aTol)
    x_t = cleanUpExpression.(x_t, aTol)
    

    return (x_t=x_t, X_s=X_s)
end


# removes terms that are less than a given tolerance
function cleanUpExpression(expression, aTol=1e-8)
    #println("Type of Expression = ", typeof(expression))
    if typeof(expression) == Float64
        return abs(expression) < aTol ? 0.0 : expression
    end

    if expression isa Sym
        terms = expression.as_ordered_terms()
        clean_terms = []
        for term in terms
            #println("Inspecting term: ", term)
            term_str = string(term)
            # Extract coefficients including those in scientific notation
            matches = eachmatch(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", term_str)
            coeffs = [m.match for m in matches]
            #println("Coefficients found: ", coeffs)
            keep_term = true
            for coeff in coeffs
                if occursin(r"e-", coeff) && abs(parse(Float64, coeff)) < aTol
                    #println("Coefficient ", coeff, " is below tolerance, zeroing out: ", coeff)
                    term_str = replace(term_str, coeff => "0")
                    updated_term = sympify(term_str)
                    if updated_term == 0
                        keep_term = false
                        #println("Resulting term after zeroing is zero.")
                        break
                    end
                end
            end
            if keep_term
                #println("Keeping term: ", sympify(term_str))
                push!(clean_terms, sympify(term_str))
            end
        end
        # Rebuild the expression from the cleaned terms
        clean_expr = isempty(clean_terms) ? 0.0 : sum(clean_terms)
        #println("Clean expression: ", clean_expr)
        return simplify(clean_expr)
    end

    return expression
end

# Function to clean up coefficients of a polynomial or a rational function
function cleanUpRationalFunction(rationalFunc, tol=1e-5)
    # Extract the numerator and denominator of the rational function
    num = sympy.numer(rationalFunc)
    denom = sympy.denom(rationalFunc)

    # Convert numerator to polynomial and get its coefficients
    num_poly = sympy.Poly(num, s)
    coeffs = num_poly.coeffs()

    
    # Clean up the coefficients
    cleaned_coeffs = cleanUp(coeffs, tol)


    # Rebuild the numerator from the cleaned coefficients
    new_num = sum([cleaned_coeffs[i] * s^(length(cleaned_coeffs)-i) for i in 1:length(cleaned_coeffs)])

    # Reassemble the rational function
    new_rational_func = new_num / denom
    new_rational_func = sympy.cancel.(new_rational_func)

    return sympy.simplify(new_rational_func)
end


function cleanUp(A, tol=1e-10)
    B = copy(A)
    for i in eachindex(B)
            real_part = abs(real(B[i])) < tol ? 0.0 : real(B[i])
            # Clean up the imaginary part
            imag_part = abs(imag(B[i])) < tol ? 0.0 : imag(B[i])
            # Reconstruct the complex number
            B[i] = real_part + imag_part * im
    end
    return B
end


# Function to check if all eigenvalues are distinct within a given tolerance
function are_eigenvalues_distinct_within_tolerance(Lambda, tol)
    # Sort the eigenvalues by magnitude to make comparison easier
    sorted_lambda = sort(Lambda, by=abs)
    # Check differences between consecutive eigenvalues
    for i in 1:length(sorted_lambda)-1
        if abs(sorted_lambda[i] - sorted_lambda[i+1]) < tol
            return false
        end
    end
    return true
end


# Function to handle decomposition and inverse Laplace transform of a symbolic expression
function decompose_and_inverse_laplace(Y)
    # Decompose Y if it is a rational function
    numerator = sympy.Poly(sympy.numer(Y), s)
    denominator = sympy.denom(Y)

    # Get coefficients and monomials
    coeffs = numerator.coeffs()
    monomials = numerator.monoms()

    # Calculate the inverse Laplace for each term and sum
    transforms = []
    for (coeff, monomial) in zip(coeffs, monomials)
        power = monomial[1]
        Yi = coeff * s^power / denominator
        transform = sympy.inverse_laplace_transform(Yi, s, t, noconds=true)
        push!(transforms, transform)
    end
    sum(transforms)
end

# Define a function to generate the callable solution from symbolic expressions
function mySolutionFactory(solution_t)
    # Convert each symbolic solution component into a callable Julia function
    solution_functions = [generate_function(sol, t) for sol in solution_t]

    # Return a new function that computes all components for any given t and returns a matrix
    function solutionVector(t_vec)
        # Apply each function to each time point and store results
        solutions = [func.(t_vec) for func in solution_functions]
        # Concatenate results into a matrix where each row is a time point and each column is a variable
        ColumnForm = hcat(solutions...)
        # return ColumnForm
        RowForm = transpose(ColumnForm)
        return RowForm
    end
    
    return solutionVector
end

# Helper function to convert a symbolic expression to a callable function
function generate_function(sol_t, t_sym)
    # Convert symbolic expression to a string
    sol_str = string(sol_t)
    # Replace SymPy's Heaviside function with Julia's step function for clarity
    sol_str = replace(sol_str, "Heaviside" => "(t -> t >= 0 ? 1 : 0)")
    # Evaluate the string to create a Julia function
    eval(Meta.parse("t -> $sol_str"))
end

# Not currently used
# Function to perform Laplace transform retaining coefficients
function myLaplaceTransform(f, t, s)
    # Extract the coefficient of the exponential term if present
    coeff, exp_term = f.args[1], f.args[2]
    
    # Perform Laplace transform on the exponential term
    L_exp_term = sympy.laplace_transform(f, t, s, noconds=True)
    
    # Return the product of the coefficient and the Laplace transform of the exponential term
    coeff * L_exp_term
end

