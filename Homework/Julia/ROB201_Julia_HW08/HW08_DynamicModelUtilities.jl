function modelParameters3LinkWalker(th3Des=pi/8)
#
    r=1;   # length of each leg
    m=5;   # mass of a leg
    Mh=15; # mass of hips
    Mt=10; #mass of torso
    L=0.5; # distance between hips and torso mass
    g=9.8; # acceleration due to gravity
#    
    return  g, r, L, m, Mh, Mt, th3Des
end

function dyn_mod_3LinkWalkerWithController(q, dq)
# DYN_MOD_3LINKWALKERWITHCONTROLLER
# 2024-06-26 09:34:40
#
# Author: Grizzle
#
# Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau 
# The Robot Equations: From Lagrange's Equations of Motion
#
g, r, L, m, Mh, Mt, th3Des = modelParameters3LinkWalker()
#
# Variable names for the model
th1, th2, th3 = q 
dth1, dth2, dth3 = dq
#
D = zeros(3, 3)
  D[1, 1] = Mh*(r^2) + 2Mt*(r^2)*(0.5(cos(th1)^2) + 0.5(sin(th1)^2)) + 2m*
            (r^2)*(0.625(cos(th1)^2) + 0.625(sin(th1)^2))
  D[1, 2] = -0.5m*(r^2)*cos(th1)*cos(th2) - 0.5m*(r^2)*sin(th1)*sin(th2)
  D[1, 3] = L*Mt*r*cos(th1)*cos(th3) + L*Mt*r*sin(th1)*sin(th3)
  D[2, 1] = -0.5m*(r^2)*cos(th1)*cos(th2) - 0.5m*(r^2)*sin(th1)*sin(th2)
  D[2, 2] = 2m*(r^2)*(0.125(cos(th2)^2) + 0.125(sin(th2)^2))
  D[3, 1] = L*Mt*r*cos(th1)*cos(th3) + L*Mt*r*sin(th1)*sin(th3)
  D[3, 3] = 2Mt*(L^2)*(0.5(cos(th3)^2) + 0.5(sin(th3)^2))
#
C = zeros(3, 3)
  C[1, 2] = dth2*(0.5m*(r^2)*cos(th1)*sin(th2) - 0.5m*(r^2)*sin(th1)*
            cos(th2))
  C[1, 3] = dth3*(L*Mt*r*sin(th1)*cos(th3) - L*Mt*r*cos(th1)*sin(th3))
  C[2, 1] = dth1*(0.5m*(r^2)*sin(th1)*cos(th2) - 0.5m*(r^2)*cos(th1)*
            sin(th2))
  C[3, 1] = dth1*(L*Mt*r*cos(th1)*sin(th3) - L*Mt*r*sin(th1)*cos(th3))
#
G = zeros(3)
  G[1] = -Mh*g*r*sin(th1) - Mt*g*r*sin(th1) - (3/2)*g*m*r*sin(th1)
  G[2] = (1/2)*g*m*r*sin(th2)
  G[3] = -L*Mt*g*sin(th3)
#
B = zeros(3, 2)
  B[1, 1] = -1
  B[2, 2] = -1
  B[3, 1] = 1
  B[3, 2] = 1
#
JacG = zeros(3, 3)
  JacG[1, 1] = -Mh*g*r*cos(th1) - Mt*g*r*cos(th1) - (3/2)*g*m*r*cos(th1)
  JacG[2, 2] = (1/2)*g*m*r*cos(th2)
  JacG[3, 3] = -L*Mt*g*cos(th3)
#
h = zeros(2)
  h[1] = th3 - th3Des
  h[2] = th1 + th2
#
Jac_h = zeros(2, 3)
  Jac_h[1, 3] = 1
  Jac_h[2, 1] = 1
  Jac_h[2, 2] = 1
#
Lfh = zeros(2)
  Lfh[1] = dth3
  Lfh[2] = dth1 + dth2
#
JacqLfh = zeros(2, 3)
#
  return (D=D, C=C, G=G, B=B, JacG=JacG, h=h, Jac_h=Jac_h, Lfh=Lfh, JacqLfh=JacqLfh)
end

function impact_terms_3LinkWalker(q, dq)
# IMPACT_TERMS_3LINKWALKER
# 2024-06-26 10:04:05
#
# Author: Grizzle
#
#
g, r, L, m, Mh, Mt, th3Des = modelParameters3LinkWalker()
#
# Variable names for the model
th1, th2, th3, x, z = q 
dth1, dth2, dth3, dx, dz = dq
#
De = zeros(5, 5)
  De[1, 1] = (r^2)*(2Mh + 2Mt)*(0.5(cos(th1)^2) + 0.5(sin(th1)^2)) + 2m*(r^2)*
            (0.625(cos(th1)^2) + 0.625(sin(th1)^2))
  De[1, 2] = -0.5m*(r^2)*cos(th1)*cos(th2) - 0.5m*(r^2)*sin(th1)*sin(th2)
  De[1, 3] = L*Mt*r*cos(th1)*cos(th3) + L*Mt*r*sin(th1)*sin(th3)
  De[1, 4] = r*(Mh + Mt)*cos(th1) + 1.5m*r*cos(th1)
  De[1, 5] = -Mh*r*sin(th1) - Mt*r*sin(th1) - 1.5m*r*sin(th1)
  De[2, 1] = -0.5m*(r^2)*cos(th1)*cos(th2) - 0.5m*(r^2)*sin(th1)*sin(th2)
  De[2, 2] = 2m*(r^2)*(0.125(cos(th2)^2) + 0.125(sin(th2)^2))
  De[2, 4] = -0.5m*r*cos(th2)
  De[2, 5] = 0.5m*r*sin(th2)
  De[3, 1] = L*Mt*r*cos(th1)*cos(th3) + L*Mt*r*sin(th1)*sin(th3)
  De[3, 3] = 2Mt*(L^2)*(0.5(cos(th3)^2) + 0.5(sin(th3)^2))
  De[3, 4] = L*Mt*cos(th3)
  De[3, 5] = -L*Mt*sin(th3)
  De[4, 1] = r*(Mh + Mt)*cos(th1) + 1.5m*r*cos(th1)
  De[4, 2] = -0.5m*r*cos(th2)
  De[4, 3] = L*Mt*cos(th3)
  De[4, 4] = Mh + Mt + 2m
  De[5, 1] = -Mh*r*sin(th1) - Mt*r*sin(th1) - 1.5m*r*sin(th1)
  De[5, 2] = 0.5m*r*sin(th2)
  De[5, 3] = -L*Mt*sin(th3)
  De[5, 5] = Mh + Mt + 2m
#
Jac_footSwing = zeros(2, 5)
  Jac_footSwing[1, 1] = r*cos(th1)
  Jac_footSwing[1, 2] = -r*cos(th2)
  Jac_footSwing[1, 4] = 1
  Jac_footSwing[2, 1] = -r*sin(th1)
  Jac_footSwing[2, 2] = r*sin(th2)
  Jac_footSwing[2, 5] = 1
#
  return (De=De, Jac_footSwing=Jac_footSwing)
end

# For use in the ODE simulator
function threeLinkDynamics(dx, x, p, t)
    n = floor(Int, length(x) / 2) # In Julia n/2 is a Float64
    q = x[1:n]
    dq = x[n+1:end]
    model = dyn_mod_3LinkWalkerWithController(q, dq)
    # Control Params and Terms
    wn = 5
    zeta = 1.1
    Kp = wn^2
    Kd = 2*zeta*wn
    L2fh = model.JacqLfh*dq + model.Jac_h*(model.D \ (-model.C*dq - model.G))
    LgLfh =  model.Jac_h*(model.D \ model.B)
    # Control Signal
    u = LgLfh \ (- L2fh - Kp*model.h - Kd*model.Lfh)
    # ODE
    dx1 = dq
    dx2 = (model.D) \ (-model.C*dq-model.G + model.B*u) # note the use of backslash
    dx[1:n] = dx1
    dx[n+1:end] = dx2
    return dx
end

# For use in the ODE simulator
function threeLinkImpactModel(q, dq)
    q_minus = [q;0.0;0.0]
    dq_minus = [dq;0.0;0.0]
    model = impact_terms_3LinkWalker(q_minus, dq_minus)
    De = model.De
    Jac_footSwing = model.Jac_footSwing
    A = [De -Jac_footSwing'; Jac_footSwing zeros(2,2)]
    ans = A \ [De * dq_minus; zeros(2,1)]
    
    # Initialize Swap matrix
    R = zeros(3, 5)
    R[1, 2] = 1
    R[2, 1] = 1
    R[3, 3] = 1
 
    # Compute q_plus and dq_plus
    q_plus = R * q_minus
    dq_plus = R * ans[1:5]   
    return q_plus, dq_plus
end

# construct collision detection function. If output is 0.0, the 
# solver knows there has been an impact with the ground
function condition(x, t, integrator)
    pi/8 - x[1]
end

# construct the Jump map. It indicates what happens after the end of a step
function affect!(integrator)
    q = integrator.u[1:3]
    dq = integrator.u[4:6]
    q_plus, dq_plus = threeLinkImpactModel(q, dq)
    integrator.u[1:3] =  q_plus
    integrator.u[4:6] = dq_plus
end

function plot_biped(th1, th2, th3, r, L, p_footStance)
    # Calculate positions based on the given angles
    p_m1 = p_footStance + [r/2*sin(th1), r/2*cos(th1)]
    p_Mh = p_footStance + [r*sin(th1), r*cos(th1)]
    p_Mt = p_Mh + [L*sin(th3), L*cos(th3)]
    p_m2 = p_Mh - [r/2*sin(th2), r/2*cos(th2)]
    p_footSwing = p_Mh - [r*sin(th2), r*cos(th2)]

    # Create a plot with predefined limits
    plt = plot(xlim=(-0.5, 9), ylim=(-0.1, 2), aspect_ratio=:equal, legend=false)

    # Plot each leg and the torso with different colors
    plot!(plt, [p_footStance[1], p_Mh[1]], [p_footStance[2], p_Mh[2]], line=(:red, 2), label="Stance Leg")
    plot!(plt, [p_Mh[1], p_Mt[1]], [p_Mh[2], p_Mt[2]], line=(:black, 2), label="Torso")
    plot!(plt, [p_Mh[1], p_footSwing[1]], [p_Mh[2], p_footSwing[2]], line=(:blue, 2), label="Swing Leg")

    plt
    return p_footSwing
end 

function animate_biped(th1, th2, th3, r, L)
    p_footStance = [0.0; 0.0]
    p_footSwing = p_footStance
    th1_last = 0.0
    anim = @animate for i in 1:length(th1)
       if abs(th1[i] - th1_last)>pi/6
           p_footStance[1] = p_footSwing[1]
        end
       p_footSwing = plot_biped(th1[i], th2[i], th3[i], r, L, p_footStance)
       th1_last = th1[i]
    end
   return anim
end

function evaluate_answers(answer_code)
    # This is an obfuscating mechanism using a pseudo-hash function with an additional obfuscation key
    function pseudo_hash(c, i, obfuscation_key)
        ascii_val = Int(c)  # Convert character to its ASCII integer value
        # Apply a more complex hashing function that uses both the position and a key
        return (ascii_val % 16 + i * (ascii_val % 3) + obfuscation_key[i] % 7) % 11
    end
    
    # Define a constant obfuscation key (could be randomly generated for each session or exam)
    obfuscation_key = [3, 8, 5, 6]

    # Compute the pseudo-hash for each student answer based on its position and the obfuscation key
    student_hashes = [pseudo_hash(Char(ans), idx, obfuscation_key) for (idx, ans) in enumerate(answer_code)]

    # Compare hashes to determine the number of correct answers
    # Obfuscated problem hashes, precomputed using the same obfuscation key
    problem_hashes = [6, 6, 4, 8] # Example hashes that you would compute beforehand
    num_correct = sum([student_hashes[i] == problem_hashes[i] for i in 1:length(answer_code)])
    
    return num_correct
end

function evaluate_answers02(answer_code)
    # This is an obfuscating mechanism using a pseudo-hash function with an additional obfuscation key
    function pseudo_hash(c, i, obfuscation_key)
        ascii_val = Int(c)  # Convert character to its ASCII integer value
        # Apply a more complex hashing function that uses both the position and a key
        return (ascii_val % 16 + i * (ascii_val % 3) + obfuscation_key[i] % 7) % 11
    end
    
    # Define a constant obfuscation key (could be randomly generated for each session or exam)
    obfuscation_key = [3, 8, 5, 6]

    # Compute the pseudo-hash for each student answer based on its position and the obfuscation key
    student_hashes = [pseudo_hash(Char(ans), idx, obfuscation_key) for (idx, ans) in enumerate(answer_code)]

    # Compare hashes to determine the number of correct answers
    # Obfuscated problem hashes, precomputed using the same obfuscation key
    problem_hashes = [5, 6, 7, 2] # Example hashes that you would compute beforehand
    num_correct = sum([student_hashes[i] == problem_hashes[i] for i in 1:length(answer_code)])
    
    return num_correct
end

