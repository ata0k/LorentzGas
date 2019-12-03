
using LinearAlgebra
using Plots
using StaticArrays
function RATIONAL_APPROXIMATION(α, ε)
    T = typeof(ε)
    h1,h2 = 1,0
    k1, k2 = 0,1
    b = α
    while abs(k1*α - h1) > ε 
        a = floor(Int,b)
        h1,h2 = a*h1 +h2,h1
        k1, k2 = a*k1 +k2, k1
        b = T(1)/(b-a)
    end
    return k1,h1
end 
function EFFICIENT_DISC_COLLISION(m, b, r)
    T = typeof(r)
    corte = 1e-15
    if T == BigFloat
        corte = 1e-30
    end
    kn = 0
    b1 = b
    ε = r*√(m^2 +T(1))
    if b < ε || (T(1)-b) < ε
        if b < T(1)/T(2)
            q, p = RATIONAL_APPROXIMATION(m, T(2)*ε)
        else
            q, p = RATIONAL_APPROXIMATION(m, T(2)*(T(1)-ε))
        end 
        b = mod(m*q+b,T(1))
        kn = kn+q
    end
    contador = 1
    while (b > ε && T(1)-b > ε) && contador<1e6
        contador += 1
        if b < T(1)/T(2)
            q, p = RATIONAL_APPROXIMATION(m, T(2)*b)
        else
            q, p = RATIONAL_APPROXIMATION(m, T(2)*(T(1)-b))
        end
        b = mod(m*q+b,T(1))
        kn += q
        if abs(b-b1) < corte
            @show "está dentro de un canal"
            return (Inf, Inf)
        end 
    end
    q = kn
    p = floor(Int,m*q)+1  
    return (q, p)
end

function LOCAL_STEP(x,v,r)
    T = typeof(r)
    e1 = [T(1),T(0)]
    e2 = [T(0),T(1)]
    n = floor.(x)
    x = x-n
    t1 = (T(1)-x[1])/v[1]
    b1 = x[2]+t1*v[2]
    t2 = -x[1]/v[1]
    b2 = x[2]+t2*v[2]
    ε = r/v[1]
    if (x-e1)⋅v< T(0)
        if abs(b1)<ε
            return e1+n, 0
        end
    end
    if (x-e2)⋅v< T(0)
        if abs(1-b2)<ε
            return e2+n, 0
        end
    end
    if (x-e1-e2)⋅v< T(0)
        if abs(1-b1)<ε
            return e1+e2+n, 0
        end
    end
    if (x-2*e2-e1)⋅v< T(0)
        if abs(2-b1)<ε
            return e1+2*e2+n, 0
        end
    end
    m = v[2]/v[1]
    ϵ = r*√(m^2 +T(1))
    b = b1-floor(b1)
    if ε>b || ε> 1-b
        s1, collider = LOCAL_STEP([1,b1]+n,v,r)
        return s1, collider
    end
    return [1,b1]+n, 1
end
function FIND_NEXT_DISC_FIRST_OCTANT(x,v,r)
    T = typeof(r)
    if norm(round.(x)-x)<r
        return round.(x)
    end
    xp, collided = LOCAL_STEP(x,v,r)
    if collided == 0
        return xp
    end
    m = v[2]/v[1]
    b = xp[2]
    b = b-floor(b)
    ϵ = r*√(m^2 +T(1))
    p,q = EFFICIENT_DISC_COLLISION(m,b,r)
    xpp = round.(xp)+[p,q]-[T(0),round(b)]
    return xpp
end
function OCTANT(v)
    θ = atan(v[2],v[1])
    if θ<0
        θ += 2π
    end
    return ceil(Int, 4θ/π)
end
function TRANSFORMATION(n, T)
    R1 = [T(0) T(1);
         -T(1) T(0)]
    R2 = [T(0) T(1);
         T(1) T(0)]
    if mod(n,2) == 1
        return R1^((n-1)/2)
    else
        return R2*(R1^((n-2)/2))
    end
end
function FIND_NEXT_DISC(x,v,r)
    n = OCTANT(v)
    T = typeof(r)
    TT = TRANSFORMATION(n, T)
    x2 = TT*x
    v2 = TT*v
    c = FIND_NEXT_DISC_FIRST_OCTANT(x2,v2,r) 
    if c[1] == Inf
        return [Inf, Inf]
    end
    return inv(TT)*c
end
function COLLISION(x,c,v,r)
    x = SVector(x[1],x[2])
    c = SVector(c[1],c[2])
    v = SVector(v[1],v[2])
    v2 = v⋅v
    B = ((x.-c)⋅v)/v2
    C = ((x.-c)⋅(x.-c)-r^2)/v2
    if B^2 -C<0
    #    @show (B)
#        return Inf
        return -B
    end
    t = -B -√(B^2 -C)
    return t
end
function POST_COLLISION_VELOCITY(x,c,v)
    n = (x-c)/norm(x-c)
    v = v-2*(v⋅n)*n
    v = v/norm(v)
    return v
end
function LORENTZ_GAS(x,v,r,steps)
    δt = 1e-10
    if typeof(r) == BigFloat
        δt = 1e-20
    end
    for i in 1:steps
        c = FIND_NEXT_DISC(x,v,r)
        if c[1] == Inf
            return "is in a chanel"
        end
        t = COLLISION(x,c,v,r)
        x = x+v*t
        v = POST_COLLISION_VELOCITY(x,c,v)
        x = x+v*δt
    end
    return x,v
end
function LORENTZ_GAS_T(x,v,r,tiempo; test = false)
    if test
        X = Float64[]
        Y = Float64[]
        δt = 1e-10
        if typeof(r) == BigFloat
            δt = 1e-20
        end
        t = 0
        v0 = copy(v)
        while t<tiempo
            v0 = copy(v)
            c = FIND_NEXT_DISC(x,v,r)
            if c[1] == Inf
                return "is in a chanel"
            end
            tt = COLLISION(x,c,v,r)
            t += tt
            x = x+v*tt
            v = POST_COLLISION_VELOCITY(x,c,v)
            x = x+v*δt
            push!(X,x[1])
            push!(Y,x[2])
            t += δt
        end
        Δt = t-tiempo
        if δt<Δt
            Δt -= δt
            x = x-v*δt
        else
            x = x-v*Δt
            return x,v, X,Y
        end
        x = x-v0*Δt
        return x,v0, X, Y
    else
        δt = 1e-10
        if typeof(r) == BigFloat
            δt = 1e-20
        end
        t = 0
        v0 = copy(v)
        while t<tiempo
            v0 = copy(v)
            c = FIND_NEXT_DISC(x,v,r)
            if c[1] == Inf
                return "is in a chanel"
            end
            tt = COLLISION(x,c,v,r)
            t += tt
            x = x+v*tt
            v = POST_COLLISION_VELOCITY(x,c,v)
            x = x+v*δt
            t += δt
        end
        Δt = t-tiempo
        if δt<Δt
            Δt -= δt
            x = x-v*δt
        else
            x = x-v*Δt
            return x,v
        end
        x = x-v0*Δt
        return x,v0
    end
end
function INITIAL_POSITION(r)
    T = typeof(r)
    e1 = [T(1), T(0)]
    e2 = [T(1), T(1)]
    e3 = [T(0), T(1)]    
    x = rand(2)
    test = true
    while test
        if norm(x)>r && norm(x-e1)>r && norm(x-e2)>r && norm(x-e3)>r
            test = false
        else
            x = rand(2)
        end
    end
    return x
end
function DRAW_CIRCLE(x,r; color = :blue)
    X = zeros(101)
    Y = zeros(101)
    for i in 0:100
        X[i+1] = r*cos(2π/100*i)+x[1]
        Y[i+1] = r*sin(2π/100*i)+x[2]
    end
    plot!(X,Y, color = color, linewidth = 0.5)
end  
