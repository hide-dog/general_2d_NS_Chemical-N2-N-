function set_volume(nodes, cellxmax, cellymax)
    volume = zeros(cellxmax, cellymax)
    for i in 1:cellxmax
        for j in 1:cellymax
            vec_r1x = nodes[i+1,j+1,1] - nodes[i,j,1]
            vec_r1y = nodes[i+1,j+1,2] - nodes[i,j,2]
            vec_r2x = nodes[i,j+1,1] - nodes[i+1,j,1]
            vec_r2y = nodes[i,j+1,2] - nodes[i+1,j,2]

            volume[i,j] = abs(vec_r1x*vec_r2y - vec_r1y*vec_r2x) /2
        end
    end
    return volume
end 

function set_mu(Qbase,cellxmax,cellymax,specific_heat_ratio,Rd)
    mu = zeros(cellxmax,cellymax)

    # サザーランドの式
    # https://www.jstage.jst.go.jp/article/jsam1937/37/4/37_4_694/_pdf/-char/ja
    mu0 = 1.82e-5 # 基準粘度[Pa s]
    T0 = 293.15       # 基準温度[K]
    C = 117       # サザーランド定数[K]

    for i in 1:cellxmax
        for j in 1:cellymax
            # T = p/(rho Rd )
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)

            if T <0
                println("T is munus ( 'A' ) !!")
                println("i : "*string(i)*"  j : "*string(j))
                println(Qbase[i,j,:])
            end
            
            mu[i,j] = mu0 * (T/T0)^1.5 * (T0+C)/(T+C)
        end
    end

    
    return mu
end


function set_lambda(Qbase,cellxmax,cellymax,mu,specific_heat_ratio,Rd)
    
    lambda = zeros(cellxmax,cellymax)
    
    # サザーランドの式
    lam0 = 22.3*10^(-3)  # 基準熱伝導率　[W/mK]
    T0 = 273.15       # 基準温度[K]
    C = 125       # サザーランド定数[K]

    for i in 1:cellxmax
        for j in 1:cellymax
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)
            lambda[i,j] = lam0*((T0+C)/(T+C))*(T/T0)^1.5
        end
    end

    return lambda
end

function set_Minf(bdcon, specific_heat_ratio)
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
        end
    end

    a = (specific_heat_ratio * p / rho)^0.5
    u = (u^2 + v^2)^0.5
    M = u/a

    return M
end

function set_gasconst(Qbase,cellxmax,cellymax,nval,nch,R)
    Rhat = zeros(cellxmax, cellymax)
    mw   = [28e-3, 14e-3]
    npre = nval - nch

    for i in 1:cellxmax
        for j in 1:cellymax
            for s in 1:nch
                Rs = R/mw[s]
                Rhat[i,j] += Qbase[i,j,npre+s] * Rs
            end
        end
    end
    return Rhat
end