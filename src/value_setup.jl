# ------------------------------------
# set volume
# ------------------------------------
function set_volume(nodes, cellxmax, cellymax, volume)
    for j in 1:cellymax
        for i in 1:cellxmax
            vec_r1x = nodes[i+1,j+1,1] - nodes[i,j,1]
            vec_r1y = nodes[i+1,j+1,2] - nodes[i,j,2]
            vec_r2x = nodes[i,j+1,1] - nodes[i+1,j,1]
            vec_r2y = nodes[i,j+1,2] - nodes[i+1,j,2]

            volume[i,j] = abs(vec_r1x*vec_r2y - vec_r1y*vec_r2x) /2
        end
    end
    return volume
end 

# ------------------------------------
# set dx and dy for local time stepping
# ------------------------------------
function set_dx_lts(dx, dy, nodes, cellxmax, cellymax)
    # define at cell boundaries
    # dx = zeros(cellxmax+1, cellymax)
    # dy = zeros(cellxmax, cellymax+1)
    
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            x1 = nodes[i,j,1]
            y1 = nodes[i,j,2]
            x2 = nodes[i,j+1,1]
            y2 = nodes[i,j+1,2]
            
            if (x2 - x1) == 0.0
                a = 0.0
                b  = y1 - a*x1
            else
                a = (y2-y1) / (x2-x1)
                b = y1 - a*x1
            end

            ccx = 0.125 * (nodes[i,j,1] + nodes[i+1,j,1] + nodes[i,j+1,1] + nodes[i+1,j+1,1])
            ccy = 0.125 * (nodes[i,j,2] + nodes[i+1,j,2] + nodes[i,j+1,2] + nodes[i+1,j+1,2])
            
            dx1 = 0.0
            if a == 0.0
                dx1 = abs(ccx - b)
            else
                dx1 = abs(-a*ccx + ccy -b)/abs(a)
            end

            ccx = 0.125 * (nodes[i-1,j,1] + nodes[i,j,1] + nodes[i-1,j+1,1] + nodes[i,j+1,1])
            ccy = 0.125 * (nodes[i-1,j,2] + nodes[i,j,2] + nodes[i-1,j+1,2] + nodes[i,j+1,2])
            
            dx2 = 0.0
            if a == 0.0
                dx2 = abs(ccx - b)
            else
                dx2 = abs(-a*ccx + ccy -b)/abs(a)
            end
            dx[i,j] = 0.5 * (dx1+dx2)
        end
    end

    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1        
            x1 = nodes[i,j,1]
            y1 = nodes[i,j,2]
            x2 = nodes[i+1,j,1]
            y2 = nodes[i+1,j,2]
            
            if (x2 - x1) == 0.0
                a = 0.0
                b  = y1 - a*x1
            else
                a  = (y2-y1) / (x2-x1)
                b  = y1 - a*x1
            end

            ccx = 0.125 * (nodes[i,j,1] + nodes[i+1,j,1] + nodes[i,j+1,1] + nodes[i+1,j+1,1])
            ccy = 0.125 * (nodes[i,j,2] + nodes[i+1,j,2] + nodes[i,j+1,2] + nodes[i+1,j+1,2])
            
            dy1 = 0.0
            if a == 0.0
                dy1 = abs(ccy - b)
            else
                dy1 = abs(-a*ccx + ccy -b)/abs(a)
            end

            ccx = 0.125 * (nodes[i,j,1] + nodes[i+1,j,1] + nodes[i,j-1,1] + nodes[i+1,j-1,1])
            ccy = 0.125 * (nodes[i,j,2] + nodes[i+1,j,2] + nodes[i,j-1,2] + nodes[i+1,j-1,2])
            
            dy2 = 0.0
            if a == 0.0
                dy2 = abs(ccy - b)
            else
                dy2 = abs(-a*ccx + ccy -b)/abs(a)
            end
            dy[i,j] = 0.5 * (dy1+dy2)
        end
    end

    return dx, dy
end

# ------------------------------------
# set dtau for local time stepping
# ------------------------------------
function set_lts(dtau, lambda_facex, lambda_facey, Qbase, cellxmax, cellymax, mu, dx, dy,
                vecAx, vecAy, volume, specific_heat_ratio, cfl)
    g = specific_heat_ratio
    
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            rho_av = 0.5 * (Qbase[i,j,1] + Qbase[i-1,j,1])
            u_av   = 0.5 * (Qbase[i,j,2] + Qbase[i-1,j,2])
            v_av   = 0.5 * (Qbase[i,j,3] + Qbase[i-1,j,3])
            mu_av  = 0.5 * (   mu[i,j]   +    mu[i-1,j]  )
            
            ap = (g * Qbase[  i,j,4] / Qbase[  i,j,1])^0.5
            am = (g * Qbase[i-1,j,4] / Qbase[i-1,j,1])^0.5
            a_av = 0.5 * (ap + am)

            U   = u_av*vecAx[i,j,1] + v_av*vecAx[i,j,2]
            lambda_facex[i,j] = abs(U) + a_av + 2*mu_av/(rho_av*dx[i,j])
        end
    end
    
    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1
            rho_av = 0.5 * (Qbase[i,j,1] + Qbase[i,j-1,1])
            u_av   = 0.5 * (Qbase[i,j,2] + Qbase[i,j-1,2])
            v_av   = 0.5 * (Qbase[i,j,3] + Qbase[i,j-1,3])
            mu_av  = 0.5 * (   mu[i,j]   +    mu[i,j-1]  )
            
            ap = (g * Qbase[  i,j,4] / Qbase[  i,j,1])^0.5
            am = (g * Qbase[i,j-1,4] / Qbase[i,j-1,1])^0.5
            a_av = 0.5 * (ap + am)

            V   = u_av*vecAy[i,j,1] + v_av*vecAy[i,j,2]
            lambda_facey[i,j] = abs(V) + a_av + 2*mu_av/(rho_av*dy[i,j])
        end
    end
    for j in 2:cellymax-1
        for i in 2:cellxmax-1
            a1 = lambda_facex[  i,  j]
            a2 = lambda_facex[i+1,  j]
            a3 = lambda_facey[  i,  j]
            a4 = lambda_facey[  i,j+1]
            lmax = maximum([a1,a2,a3,a4])

            dtau[i,j] = cfl * volume[i,j] / lmax
        end
    end

    return dtau
end

# ------------------------------------
# set gas constant 
# ------------------------------------
function set_gasconst(Qbase, cellxmax, cellymax, Rs)
    for s in 1:nch
        for i in 1:cellxmax
            for j in 1:cellymax
                Rhat[i,j] += Qbase[i,j,npr+s] * Rs[s]
            end
        end
    end
    return Rhat
end

# ------------------------------------
# set Minf for AUSM+up (don't use)
# ------------------------------------
function set_Minf(bdcon, specific_heat_ratio, Rd)
    rho = 0 
    u   = 0 
    v   = 0 
    p   = 0 

    for i in 1:4
        if Int(bdcon[i][1]) == 0 || Int(bdcon[i][1]) == 5
            rho = bdcon[i][2]
            u   = bdcon[i][3]
            v   = bdcon[i][4]
            p   = bdcon[i][5]
        elseif Int(bdcon[i][1]) == 6
            rho = bdcon[i][2]
            u   = bdcon[i][3]
            v   = bdcon[i][4]
            T   = bdcon[i][8]
            p   = rho*Rd*T
        end
    end

    a = (specific_heat_ratio * p / rho)^0.5
    u = (u^2 + v^2)^0.5
    M = u/a

    return M
end

# ------------------------------------
# Conversion from primitive variables to conserved variables
# ------------------------------------
function base_to_conservative(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """ 
    
    for j in 1:cellymax
        for i in 1:cellxmax
            Qcon[i,j,1] = Qbase[i,j,1]
            Qcon[i,j,2] = Qbase[i,j,1]*Qbase[i,j,2]
            Qcon[i,j,3] = Qbase[i,j,1]*Qbase[i,j,3]
            Qcon[i,j,4] = Qbase[i,j,4]/(specific_heat_ratio-1)+Qbase[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2
            for l in 1:nch
                Qcon[i,j,npr+l] = Qbase[i,j,1]*Qbase[i,j,npre+l]
            end
        end
    end
    return Qcon
end

# ------------------------------------
# Conversion from conserved variables to primitive variables
# ------------------------------------
function conservative_to_base(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """

    for j in 1:cellymax
        for i in 1:cellxmax
            Qbase[i,j,1] = Qcon[i,j,1]
            Qbase[i,j,2] = Qcon[i,j,2]/Qcon[i,j,1]
            Qbase[i,j,3] = Qcon[i,j,3]/Qcon[i,j,1]
            Qbase[i,j,4] = (Qcon[i,j,4]-Qcon[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2)*(specific_heat_ratio-1)
            for l in 1:nch
                Qbase[i,j,npr+l] = Qcon[i,j,npre+l]/Qcon[i,j,1]
            end
        end
    end
    return Qbase
end

# ------------------------------------
# Conversion from Cartesian coordinate system to general coordinate system
# ------------------------------------
function setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume)
    for l in 1:nval
        for j in 1:cellymax
            for i in 1:cellxmax
                Qcon_hat[i,j,l] = Qcon[i,j,l] * volume[i,j]
            end
        end
    end
    return Qcon_hat
end

# ------------------------------------
# Conversion from general coordinate system to Cartesian coordinate system
# ------------------------------------
function Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume)
    for l in 1:nval
        for j in 1:cellymax
            for i in 1:cellxmax
                Qcon[i,j,l] = Qcon_hat[i,j,l] / volume[i,j]
            end
        end
    end
    return Qcon
end