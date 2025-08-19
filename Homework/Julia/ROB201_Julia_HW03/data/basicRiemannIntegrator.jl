function myRiemannIntegralValueAndError(f, a, b, N, aTol=1e-6)
    # Works for b < a 
    if abs(b-a) < aTol
        println("a and b are too close together for most numerical integration methods.")
        return NaN
    end
    N = floor(Int64, N)
    if N < 2
        N = 2
    end
    Δx = (b-a)/N
#     # You can uncomment these messages if you want to see them
#     if b < a
#         println("Integrating backward with Δx = $Δx")
#     else
#         println("Integrating forward with Δx = $Δx")
#     end
    x = range(a, b, N+1)
    y = f.(x)
    upperSum = 0.0
    lowerSum = 0.0
    for i = 1:N
        lowerSum += Δx *minimum([y[i], y[i+1]])
        upperSum += Δx *maximum([y[i], y[i+1]])
    end
    if b < a # This ensures the upperSum is larger or equal to the lowerSum
             # Could remove and use instead
             # myPMerror =  abs(upperSum - lowerSum)/2.0
        temp = lowerSum
        lowerSum = upperSum
        upperSum = temp
    end    
    myEstimatedIntegral = (upperSum + lowerSum)/2.0
    myPMerror =  (upperSum - lowerSum)/2.0
    return myEstimatedIntegral, myPMerror
end