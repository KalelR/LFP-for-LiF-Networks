using Plots
pyplot()

using BenchmarkTools

function LIF!(t, u, du, p)
    gl::Float64, vr::Float64, C::Float64, I::Float64 = p
    v = u[1]::Float64
    # println(I, " ", typeof(I))

    du[1] = (-gl*(v - vr) + I)/C
    return du
end


#--- f = dy/dx
function integraStep(t::Float64, y, dy, h, f, p)
    k1 = h * f(t, y, dy, p);
    # println(k1, " ", y)
    k2 = h * f(t+h/2, y+k1/2, dy, p );
    k3 = h * f(t+h/2, y+k2/2, dy, p );
    k4 = h * f(t+h,   y+k3, dy, p );

    yn = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    # y = yn;
    t += h;
    return yn, t
end
#----

# t = 0.
# y = 0.
# dy = 0.
# h = 0.01
# @btime integraStep(t,y,dy,h,LIF!,p)

#----solver function, receives ...
function solver(y0, t_span, dt, p,f)
    y = y0;
    t, t_end = t_span
    n_vars = length(y0)
    aux_max = Int64(floor((t_end-t)/dt))
    v_t = zeros(Float64, aux_max)
    m_y = zeros(Float64, (n_vars, aux_max))
    dy = similar(y)
    vr = -70.;
    v_th = 0.;
    for aux=1:aux_max
        # println(m_y[:,aux], " ", y)
        m_y[:,aux] = y
        v_t[aux] = t
        # y, t = @inbounds @fastmath integraStep(t,y, dy, dt, lorenz!, p)
        y, t = integraStep(t,y, dy, dt, f, p)
        if(y[1] >= v_th)
            y = [vr];
        end

    end
    return v_t, m_y
end
#--------

#------- main code
#model params
gl = 4.5 #nS
vr = -70.0 #mV
I = 500.0 #pA
C = 100.0 #pF
vth = 0.0 #mV

#integrator params
u0 = collect(range(-60, -65, length=N))
# u0 = [-60.0]
tspan = [0.0, 1000.0]
dt = 0.1
p = [gl, vr, C, I]

v_t, m_y = @btime @inbounds @fastmath solver(u0, tspan, dt, p, LIF!)

plot(v_t, m_y[1,:], marker=:circle, markerstrokewidth=0, markersize=1)
# plot!(xlim=(0,2))
# plot!(v_t, m_y[2])
# scatter3d(m_y[1,:], m_y[2,:], m_y[3,:],   markersize=2, markerstrokewidth=0)
# plot3d(m_y[1,:], m_y[2,:], m_y[3,:],  lw=1, markersize=2, markerstrokewidth=0)
# plot(m_y, vars=(1,2,3))

# savefig("lorenz_rk4_x.png")



##--tests
