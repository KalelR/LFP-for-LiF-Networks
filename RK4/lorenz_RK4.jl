
using Plots
pyplot()

# function lorenz!(du,u,p,t)
#   x,y,z = u
#   σ,ρ,β = p
#   du[1] =  σ*(y-x)
#   du[2] =  x*(ρ-z) - y
#   du[3] =  x*y - β*z
#   return du
# end
#
function lorenz!(t, y)
    p = [10.0,28.0,8/3]
    x,y,z = y
    σ,ρ,β = p
    du = zeros(3)
    du[1] =  σ*(y-x)
    du[2] =  x*(ρ-z) - y
    du[3] =  x*y - β*z
    return du
end


#--- f = dy/dx
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

function solver(y0, t_span, dt, p)
    y = y0;
    t, t_end = t_span
    n_vars = length(y0)
    aux_max = Int64(floor((t_end-t)/dt))
    v_t = zeros(Float64, aux_max)
    m_y = zeros(Float64, (n_vars, aux_max))
    for aux=1:aux_max
        m_y[:,aux] = y
        v_t[aux] = t
        y, t = integraStep(t,y,dt,lorenz!)
    end
    return v_t, m_y
end

p = [10.0,28.0,8/3]
y0 = [1.0, 0.0, 0.0]
tspan = (0.0,1000.0)
dt = 0.01
v_t, m_y = @time solver(y0, tspan, dt, p)

# plot(v_t, m_y[1,:])
# plot!(v_t, m_y[2])
# scatter3d(m_y[1,:], m_y[2,:], m_y[3,:],   markersize=2, markerstrokewidth=0)
# plot3d(m_y[1,:], m_y[2,:], m_y[3,:],   markersize=2, markerstrokewidth=0)
# plot(m_y, vars=(1,2,3))

savefig("lorenz_rk4_x.png")
