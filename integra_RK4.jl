 

#--- f = dy/dx
function integraStep(t::Float64, y, dy, h, f, p, aux)
    # for i=1:N
        k1 = h * f(t, y, dy, p, aux, 1);
        # println(k1, " ", y)
        k2 = h * f(t+h/2, y+k1/2, dy, p, aux, 2 );
        k3 = h * f(t+h/2, y+k2/2, dy, p, aux, 3 );
        k4 = h * f(t+h,   y+k3, dy, p, aux,4 );

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
function solver(y0, t_span, dt, p,f, v_reset, v_th, N, m_tS, v_numS, m_tS_ext, v_numS_ext)
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
        v_y, t = integraStep(t, v_y, v_dy, dt, f, p, aux)
        # length()
        for i=1:N
             if(v_y[1*i] >= v_th)
                 v_y[i] = v_reset;
                 v_numS[i] +=1; m_tS[i][v_numS[i]] = t;
            end
            if(v_numS_ext[i] < length(m_tS_ext[i]))
                if (t >= m_tS_ext[i][v_numS_ext[i]+1])
                v_numS_ext[i] += 1;
                end 
            end
        end #reset potential

    end
    return v_t, m_y
end
#