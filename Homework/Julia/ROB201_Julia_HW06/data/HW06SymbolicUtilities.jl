# Define the model parameters
function modelParametersSlinky()
    g, m, L, k, ell0 = (9.81, 1e-3, 7e-4, 50, 7e-4)
    return g, m, L, k, ell0
end

# Some model terms can be too long to parse effectively. This function breaks the terms into smaller 
# chunks
#
function fixlength(s1, s2, len, indent="    ")
#=
    FIXLENGTH Returns a string which has been divided up into < LEN
    character chunks with the '...' line divide appended at the end
    of each chunk.
    FIXLENGTH(S1, S2, L) is string S1 with string S2 used as break points into
    chunks less than length L.
=#
    
    tmp = s1
    K = "" 
    count = 0
    
    while length(tmp) > len
        I = Int[]
        
        for c in s2
            J = findall(occursin(string(c)), tmp)  # Convert character to string here
            
            if (c == '+' || c == '-') && !isempty(J)
                JJ = Int[]
                for cc in J
                    if tmp[cc] != 'e'
                        push!(JJ, cc)
                    end
                end
                I = vcat(I, JJ)
            else
                I = vcat(I, J)
            end
        end
        
        sort!(I)
        
        if isempty(I) && count == 0
            K = ""
            error("s2 does not exist in s1")
        end
        
        II = findall(x -> x <= len, I)
        
        if isempty(II)
            K = ""
            error("Cannot fixlength of s1 on basis of s2. Count = $(count)")
        end
        
        #K *= tmp[1:I[II[end]]] * "...\n" * indent
        K *= tmp[1:I[II[end]]] * "\n" * indent * indent * indent
        tmp = tmp[I[II[end]]+1:end]
        count += 1
    end
    
    K *= tmp
    return K
end


using Dates
using Symbolics

# Uses function fixlength

function  writeEOM(fcn_name, D, C, G, modelParamString,variableNamesString)

s2 = "+-*" # favorable places to break lines

println("[creating ", uppercase(fcn_name), ".jl]")

open("$(fcn_name).jl", "w") do fid
    n = length(q)
    println(fid, "function $fcn_name(q, dq)")
    println(fid, "# ", uppercase(fcn_name))
    println(fid, "# ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(fid, "#")
    println(fid, "# Author: Grizzle")
    println(fid, "#")
    println(fid, "# Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau ")
    println(fid, "# The Robot Equations: From Lagrange's Equations of Motion")
    println(fid, "#")
    println(fid, modelParamString)
    println(fid, "#")
    println(fid, "# Variable names for the model")
    println(fid, variableNamesString) 
    println(fid, "#")
    println(fid, "D = zeros($n, $n)")
    for i in 1:n
        for j in 1:n
            Temp0 = D[i,j]
            if !isequal(Temp0, 0)
                Temp1 = string(Temp0)
                Temp1 = replace(Temp1, "//" => "/")
                Temp1=fixlength(Temp1,s2,65)
                println(fid, "  D[$i, $j] = $Temp1")
            end
        end
    end

    println(fid, "#")
    println(fid, "C = zeros($n, $n)")
    for i in 1:n
        for j in 1:n
            Temp0 = C[i,j]
            if !isequal(Temp0, 0)
                Temp1 = string(Temp0)
                Temp1 = replace(Temp1, "//" => "/")
                Temp1=fixlength(Temp1,s2,65)
                println(fid, "  C[$i, $j] = $Temp1")
            end
        end
    end

    println(fid, "#")
    println(fid, "G = zeros($n)")
    for i in 1:n
        Temp1 = string(G[i])
        Temp1 = replace(Temp1, "//" => "/")
        Temp1=fixlength(Temp1,s2,65)
        println(fid, "  G[$i] = $Temp1")
    end

    
    println(fid, "#")
    println(fid, "  return (D=D, C=C, G=G)")
    println(fid, "end")
end

println("File $(fcn_name).jl created successfully!")
end