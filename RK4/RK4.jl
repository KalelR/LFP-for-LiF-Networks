using Plots
function f(t, y)
    -t
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
    y = 10.;
    t_step = 0.1
    t = 0.
    v_y, v_t = Float64[], Float64[]
    for aux=0:100
        y, t = integraStep(t,y,t_step,f)
        println(t, " ", y)
        push!(v_y, y)
        push!(v_t, t)
    end
    return v_t, v_y
end

v_t, v_y = solve()

plot(v_t, v_y)
