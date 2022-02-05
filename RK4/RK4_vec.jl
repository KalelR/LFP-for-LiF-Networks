using Plots
function f(t, y)
    # [-y[1], -5*y[2]]
    [-t, -5*t]
end

function integraStep(t, y, h, f)
    k1 = h * f(t, y);
    k2 = h * f(t+h/2, y+k1/2 );
    k3 = h * f(t+h/2, y+k2/2 );
    k4 = h * f(t+h,   y+k3 );

    yn = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    # y = yn;
    t += h;
    return yn, t
end

function solve()
    y = [10.,20.];
    t_step = 0.1
    t = 0.
    m_y, v_t = [Float64[] for i=1:2], Float64[]
    for aux=0:1
        push!(m_y[1], y[1])
        push!(m_y[2], y[2])
        y, t = integraStep(t,y,t_step,f)
        println(t, " ", y)
        push!(v_t, t)
    end
    return v_t, m_y
end

v_t, m_y = solve()

plot(v_t, m_y[1])
plot!(v_t, m_y[2])
