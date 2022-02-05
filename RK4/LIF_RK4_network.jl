using Plots
pyplot()

using BenchmarkTools

function LIF!(t, u, du, p)
    gl::Float64, vr::Float64,C::Float64, I::Float64 = p
    # v = u[1]::Float64
    v = u
    # println(I, " ", typeof(I))

    du = @. (-gl*(v - vr) + I)/C
    # println(v, " ", du)
    return du
end


#--- f = dy/dx
function integraStep(t::Float64, y, dy, h, f, p)
    # for i=1:N
        k1 = h * f(t, y, dy, p);
        # println(k1, " ", y)
        k2 = h * f(t+h/2, y+k1/2, dy, p );
        k3 = h * f(t+h/2, y+k2/2, dy, p );
        k4 = h * f(t+h,   y+k3, dy, p );

        yn = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    # end
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
function solver(y0, t_span, dt, p,f, v_reset, v_th, N)
    v_y = y0;
    v_dy = similar(v_y)
    t, t_end = t_span
    n_vars = length(y0)
    aux_max = Int64(floor((t_end-t)/dt))
    v_t = zeros(Float64, aux_max)
    m_y = zeros(Float64, (n_vars, aux_max))
    for aux=1:aux_max
        m_y[:,aux] = v_y
        v_t[aux] = t
        # v_y, t = @inbounds @fastmath integraStep(t, v_y, v_dy, dt, f, p)
        v_y, t = integraStep(t, v_y, v_dy, dt, f, p)
        for i=1:N if(v_y[1*i] >= v_th) v_y[i] = v_reset; end end #reset potential

    end
    return v_t, m_y
end
#--------

#------- main code
#model params
gl = 4.5 #nS
v_reset = -70.0 #mV
I = 500.0 #pA
C = 100.0 #pF
vth = -10. #mV
N = 100;
Er = 0.0;
#integrator params
u0 = collect(range(-60, -65, length=N))
# u0 = [-60.0]
tspan = [0.0, 1000.0]
dt = 0.05
p = [gl, Er, C, I]

# v_t, m_y = @btime @inbounds @fastmath solver(u0, tspan, dt, p, LIF!, v_reset, vth, N)
v_t, m_y =  @btime solver(u0, tspan, dt, p, LIF!, v_reset, vth, N)

plot(v_t, m_y[1,:], marker=:circle, markerstrokewidth=0, markersize=1)
# plot!(xlim=(0,2))
# plot!(v_t, m_y[2])
# scatter3d(m_y[1,:], m_y[2,:], m_y[3,:],   markersize=2, markerstrokewidth=0)
# plot3d(m_y[1,:], m_y[2,:], m_y[3,:],  lw=1, markersize=2, markerstrokewidth=0)
# plot(m_y, vars=(1,2,3))

# savefig("lorenz_rk4_x.png")



##--tests
