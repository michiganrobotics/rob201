using Random, LinearAlgebra, Statistics
Random.seed!(2024)

function test_limit_at_zero(f; target=Inf, threshold=1e5, tol=1e-3, verbose=false)
    # sample values of x approaching 0 (e.g., 1e-5, 1e-6, ..., 1e-8) and 
    # check if f(x) approaches "target" value
    #
    xs = [1e-5, 1e-6, 1e-7, 1e-8]
    xs = sort(vcat(-reverse(xs), xs))  # sample from both sides approaching zero
    ys = f.(xs)

    if verbose  # Set to true to see the function values
        println("Testing $(f)(x) as x → 0:")
        for (x, y) in zip(xs, ys)
            println("$(f)($x) = $y")
        end
    end

    result = false
    if target == Inf
        result = maximum(ys) > threshold
    elseif target == 0
        result = all(abs.(ys) .< tol)
    else
        result = all(abs.(ys .- target) .< tol)
    end

    if result
        println("✅ $(f)(x) appears to approach $target as x → 0")
    else
        println("❌ $(f)(x) does not appear to approach $target as x → 0")
    end

    return result
end

function test_all_non_negative(f; verbose=false)
    # sample values of x over the domain (e.g., -0.5, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.5)
    # and check if f(x) is non-negative
    #
    xs = [0.5, 0.1, 0.05, 0.01]
    xs = sort(vcat(-reverse(xs), xs))  # sample from both sides approaching zero
    ys = f.(xs)
    all_non_negative = all(y -> y >= 0, ys)  # Check if all values are greater or equal to zero

    if verbose  # Set to true to see the function values
        println("Testing if $(f)(x) >= 0 over the domain $(first(xs)) to $(last(xs)):")
        for (x, y) in zip(xs, ys)
            println("$(f)($x) = $y")
        end
    end

    if all_non_negative
        println("✅ $(f)(x) is non-negative over the domain.")
    else
        println("❌ $(f)(x) is not non-negative everywhere in the domain.")
    end

    return all_non_negative
end

function test_no_limit_at_zero(f; verbose=false, relative_threshold=0.2)
    # sample values of x approaching 0 (e.g. 1e-5, 1e-6, 1e-7, 1e-8)
    # and check if f(x) approaches a limit (Declaring "no limit" if function fluctuates wildly near 0 relative to its mean)
    #
    xs = [1e-5, 1e-6, 1e-7, 1e-8]
    xs = vcat(-reverse(xs), xs)
    ys = f.(xs)

    # Filter valid outputs
    valid = filter(y -> isfinite(y), ys)
    if length(valid) < 2
        println("❌ Not enough valid values to determine limit behavior.")
        return false
    end

    μ = mean(valid)
    σ = maximum(valid) - minimum(valid)

    # Normalize variation by average magnitude
    rel_variation = σ / (abs(μ) + eps())  # add eps() to avoid divide by zero

    if verbose
        println("Testing for NO LIMIT at x → 0:")
        for (x, y) in zip(xs, ys)
            println("$f($x) = $y")
        end
        println("Mean: $μ, Raw Variation: $σ, Relative Variation: $rel_variation")
    end

    if rel_variation > relative_threshold
        println("✅ $f(x) appears to have NO LIMIT at x → 0 (relative variation = $rel_variation)")
        return true
    else
        println("❌ $f(x) appears to approach a limit at x → 0 (relative variation = $rel_variation)")
        return false
    end
end

function test_limit_at_infinity(f; target=5.0, tol=1e-3, verbose=false)
    # Sample increasing values of x approaching ∞
    xs = [1e5, 1e6, 1e7, 1e8]
    ys = f.(xs)

    if verbose
        println("Testing $f(x) → $target as x → ∞:")
        for (x, y) in zip(xs, ys)
            println("$f($x) = $y")
        end
    end

    # Check whether all large-x values are close to target
    close_enough = all(abs.(ys .- target) .< tol)

    if close_enough
        println("✅ $f(x) appears to converge to $target as x → ∞")
        return true
    else
        println("❌ $f(x) does not appear to converge to $target as x → ∞")
        return false
    end
end

function test_strictly_increasing(f; domain=1:0.5:100, verbose=false)
    xs = collect(domain)
    ys = f.(xs)

    is_strictly_increasing = all(diff(ys) .> 0)

    if verbose
        println("Testing if $f(x) is strictly increasing over domain $(first(xs)) to $(last(xs))")
        for (x, y) in zip(xs, ys)
            println("$f($x) = $y")
        end
    end

    if is_strictly_increasing
        println("✅ $f(x) is strictly increasing.")
    else
        println("❌ $f(x) is NOT strictly increasing.")
    end

    return is_strictly_increasing
end