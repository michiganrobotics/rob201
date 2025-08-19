using LinearAlgebra

function jointTorques(q, vx)
    # Link Parameters
    g, L1, L2, L3, m1, m2, m3 = (9.81, 1.0, 0.9, 2.5, 8.0, 5.0, 8.0)
    # Drag resistance 
    Kw = 3e-1
    # Variable names for the model
    q1, q2, q3 = q 
    # Gravity vector
    G = zeros(eltype(q), 3)
    G[1] = g*m2*(L1*cos(q1) + L2*cos(q1 + q2)) + g*m3*(L1*cos(q1) + L2*cos(q1 + q2) + L3*cos(q1 + q2 + q3)) + L1*g*m1*cos(q1)
    G[2] = g*m3*(L2*cos(q1 + q2) + L3*cos(q1 + q2 + q3)) + L2*g*m2*cos(q1 + q2)
    G[3] = L3*g*m3*cos(q1 + q2 + q3)
    # Identity matrix B
    B = Matrix{eltype(q)}(I, 3, 3)
    # Jacobian of the third link's center of mass
    J1transpose = [-0.5L3*sin(q1 + q2 + q3) - L1*sin(q1) - L2*sin(q1 + q2) - 0.5L3*sin(q1 + q2 + q3);  
        -L2*sin(q1 + q2) - 0.5L3*sin(q1 + q2 + q3); -0.5L3*sin(q1 + q2 + q3)]
    J2transpose = [L1*cos(q1) + L2*cos(q1 + q2) + 0.5L3*cos(q1 + q2 + q3) ; 
        L2*cos(q1 + q2) + 0.5L3*cos(q1 + q2 + q3) ;  0.5L3*cos(q1 + q2 + q3)]
    Jacp3cTranspose = [J1transpose  J2transpose]
    # Wind forces
    Fxwind = Kw * abs(sin(q1 + q2 + q3)) * vx * abs(vx)
    Fywind = 0.0 * Kw * abs(cos(q1 + q2 + q3)) * vx * abs(vx)
    Fwind = [Fxwind; Fywind]
    # Calculate torques
    tau = (-G + Jacp3cTranspose * Fwind)   
    return tau
end


# Define the model parameters
function modelParameters()
    g, L1, L2, L3, m1, m2, m3 = (9.81, 1.0, 0.9, 2.5, 8.0, 5.0, 8.0)
    return (L1=L1, L2=L2, L3=L3)
end




function plot_hemisphere(x_center, y_center, radius; color=:grey)
    θ = LinRange(0, π, 100)  # 100 points from 0 to π
    x = x_center .+ radius .* cos.(θ)
    y = y_center .+ radius .* sin.(θ)
    plot!(x, y, fill=(0, color), label=nothing)
end

# Function to create the vertices of a rectangle
function create_flag(x, y, width, height, phi)
    # Define the corners of the rectangle
    corners = [
        (0, 0),
        (width, 0),
        (width, height),
        (0, height)
    ]
    
    # Rotation matrix
    R = [
        cos(phi) -sin(phi);
        sin(phi)  cos(phi)
    ]
    
    # Rotate and translate the corners
    rotated_corners = [(R * [cx, cy]) .+ [x, y] for (cx, cy) in corners]
    
    # Separate the x and y coordinates
    x_coords = [point[1] for point in rotated_corners]
    y_coords = [point[2] for point in rotated_corners]
    
    # Close the rectangle by adding the first point again
    push!(x_coords, x_coords[1])
    push!(y_coords, y_coords[1])
    
    return (x=x_coords, y=y_coords)
end

function linkPostions(q)
    th1,th2,th3 = q
    params = modelParameters()
    p0 = [0.0; 0.0]
    p1 = p0 + [params.L1 * cos(th1); params.L1 * sin(th1)]
    p2 = p1 + [params.L2 * cos(th1+th2); params.L2 * sin(th1+th2)]
    p3 = p2 + [params.L3 * cos(th1+th2+th3); params.L3 * sin(th1+th2+th3)]
    p3c =  p2 + (2.0/3.0)*[params.L3 * cos(th1+th2+th3); params.L3 * sin(th1+th2+th3)]
    return (p0=p0, p1=p1, p2=p2, p3=p3, p3c=p3c, phi = (th1+th2+th3))
end


function plot_armFlagBearer(q; line_thickness=5, ball_size=10)
    positions = linkPostions(q)
#
    p0 = positions.p0
    p1 = positions.p1
    p2 = positions.p2
    p3 = positions.p3  
    p3c = positions.p3c
    phi = positions.phi   
   
    
    # Colors for the lines for each call
    line_colors = [:blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :black]
    
    # Increment the call count and determine the color for this call
    call_count = 0
    current_color = line_colors[mod(call_count, length(line_colors)) + 1]        
             
    # Plot a rectangle at the end of the last link as a flag    
    # Parameters
    width = 1.2  # Width of the rectangle
    height = 0.6  # Height of the rectangle
    # Create the flag's vertices
    Flag = create_flag(p3[1], p3[2], width, height, phi + pi/2)
    # Plot the flag
    plot!(Flag.x, Flag.y, seriestype = :shape, fillalpha=0.8, fillcolor=:yellow, 
         linecolor=:blue, linewidth=10, legend=false)
    #scatter!([x], [y], legend=false)

    
    # Plot the lines
    plot!([p0[1], p1[1], p2[1], p3[1]], [p0[2], p1[2], p2[2], p3[2]], 
        linewidth=line_thickness, color=current_color, label="Configuration $(call_count)")        
  
    
    # Plot the balls
    scatter!([p0[1]], [p0[2]], color=:black, markersize=ball_size, label=nothing)
    scatter!([p1[1], p2[1]], [p1[2], p2[2]], color=:red, markersize=ball_size, label=nothing)
    scatter!([p3[1]], [p3[2]], color=:yellow, markersize=1.5*ball_size, label=nothing)

end



