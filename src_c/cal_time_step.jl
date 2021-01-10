function time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:nval
                Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + dt*RHS[i,j,l]
            end
        end
    end
    return Qcon_hat
end 

function one_wave(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio,volume)
    # セル中心のonewave
    A_adv_hat_p = zeros(cellxmax,cellymax,4,4)
    A_adv_hat_m = zeros(cellxmax,cellymax,4,4)
    B_adv_hat_p = zeros(cellxmax,cellymax,4,4)
    B_adv_hat_m = zeros(cellxmax,cellymax,4,4)
    A_beta_shig = zeros(cellxmax,cellymax)
    B_beta_shig = zeros(cellxmax,cellymax)
    beta = 1.1
    for i in 2:cellxmax
        for j in 2:cellymax
            for k in 1:2 #A,B
                jacob_temp = zeros(4,4)
                if k == 1
                    kx_av = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
                    ky_av = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
                elseif k == 2
                    kx_av = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
                    ky_av = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
                end

                rho = Qbase[i,j,1]
                u = Qbase[i,j,2]
                v = Qbase[i,j,3]
                e = Qcon[i,j,4]
                p = Qbase[i,j,4]

                g = specific_heat_ratio
                
                Z = kx_av*u+ky_av*v
                q2   = 0.5*(u^2+v^2)
                b1c2 = 0.5*q2*(g-1)
                gebyrho = g*e/rho

                jacob_temp[1,1] = 0.0
                jacob_temp[1,2] = kx_av
                jacob_temp[1,3] = ky_av
                jacob_temp[1,4] = 0.0
            
                jacob_temp[2,1] = -u*Z + kx_av*b1c2
                jacob_temp[2,2] = Z - (g-2)*kx_av*u
                jacob_temp[2,3] = ky_av*u - kx_av*(g-1)*v
                jacob_temp[2,4] = (g-1)*kx_av
            
                jacob_temp[3,1] = -v*Z + ky_av*b1c2
                jacob_temp[3,2] = kx_av*v - ky_av*(g-1)*u
                jacob_temp[3,3] = Z - ky_av*(g-2)*v
                jacob_temp[3,4] = (g-1)*ky_av
            
                jacob_temp[4,1] = Z*(-gebyrho + 2*b1c2)
                jacob_temp[4,2] = kx_av*(gebyrho-b1c2) - (g-1)*u*Z
                jacob_temp[4,3] = ky_av*(gebyrho-b1c2) - (g-1)*v*Z
                jacob_temp[4,4] = g*Z

                c = (g*rho/p)^0.5
                shigma = abs(Z) + c*(kx_av^2+ky_av^2)^0.5

                I_temp = zeros(4,4)
                for l in 1:4
                    I_temp[l,l] = beta * shigma
                end

                if k == 1
                    for l in 1:4
                        for m in 1:4
                            A_adv_hat_p[i,j,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                            A_adv_hat_m[i,j,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                        end
                    end
                    A_beta_shig[i,j] = beta * shigma
                elseif k ==2
                    for l in 1:4
                        for m in 1:4
                            B_adv_hat_p[i,j,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                            B_adv_hat_m[i,j,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                        end
                    end
                    B_beta_shig[i,j] = beta * shigma
                end

            end
        end
    end

    return A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig
end

function central_diff_jacobian(Qbase,Qcon,cellxmax,cellymax,mu,lambda,vecAx,vecAy,specific_heat_ratio,volume)
    jalphaP = zeros(cellxmax, cellymax)
    jbetaP = zeros(cellxmax, cellymax)

    # 角セルの処理
    Qbase[1,1,1] = 1.0
    Qbase[cellxmax,1,1] = 1.0
    Qbase[1,cellymax,1] = 1.0
    Qbase[cellxmax,cellymax,1] = 1.0
    volume[1,1,1] = 1.0
    volume[cellxmax,1,1] = 1.0
    volume[1,cellymax,1] = 1.0
    volume[cellxmax,cellymax,1] = 1.0
    
    for i in 1:cellxmax
        for j in 1:cellymax
            rho = Qbase[i,j,1]
            xi_x_av = 0.5*(vecAx[i,j,1] + vecAx[i+1,j,1]) / volume[i,j]
            xi_y_av = 0.5*(vecAx[i,j,2] + vecAx[i+1,j,2]) / volume[i,j]
            alpha = (xi_x_av^2+xi_y_av^2)^0.5
            #alpha = xi_x_av^2+xi_y_av^2
            jalphaP[i,j] = alpha / volume[i,j] * (2*mu[i,j]/rho)

            eta_x_av = 0.5*(vecAy[i,j,1] + vecAy[i,j+1,1]) / volume[i,j]
            eta_y_av = 0.5*(vecAy[i,j,2] + vecAy[i,j+1,2]) / volume[i,j]
            beta = (xi_x_av^2+xi_y_av^2)^0.5
            #beta = xi_x_av^2+xi_y_av^2
            jbetaP[i,j] = beta / volume[i,j] * (2*mu[i,j]/rho)
        end
    end
    return jalphaP, jbetaP
end

function lusgs(dt,Qcon_hat,Qconn_hat,delta_Q,A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig,jalphaP, jbetaP,RHS,cellxmax,cellymax,volume)
    D = zeros(cellxmax,cellymax)
    Lx = zeros(cellxmax,cellymax,4,4)
    Ly = zeros(cellxmax,cellymax,4,4)
    Ux = zeros(cellxmax,cellymax,4,4)
    Uy = zeros(cellxmax,cellymax,4,4)
    
    I = zeros(4,4)
    for i in 1:4
        I[i,i] = 1.0
    end
    
    # L,U
    for i in 1:cellxmax
        for j in 1:cellymax
            for l in 1:4
                for m in 1:4
                    Lx[i,j,l,m] = dt*(A_adv_hat_p[i,j,l,m] + jalphaP[i,j]*I[l,m])
                    Ly[i,j,l,m] = dt*(B_adv_hat_p[i,j,l,m] + jbetaP[i,j]*I[l,m])
                    Ux[i,j,l,m] = dt*(A_adv_hat_m[i,j,l,m] - jalphaP[i,j]*I[l,m])
                    Uy[i,j,l,m] = dt*(B_adv_hat_m[i,j,l,m] - jbetaP[i,j]*I[l,m])
                end
            end
        end
    end

    LdQ = zeros(cellxmax,cellymax,4)
    UdQ = zeros(cellxmax,cellymax,4)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                for m in 1:4
                    LdQ[i,j,l] += Lx[i-1,j,l,m]*delta_Q[i-1,j,m] + Ly[i,j-1,l,m]*delta_Q[i,j-1,m]
                    UdQ[i,j,l] += Ux[i+1,j,l,m]*delta_Q[i+1,j,m] + Uy[i,j+1,l,m]*delta_Q[i,j+1,m]
                end
            end
        end
    end

    # diagonal
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            D[i,j] = 2.0 + dt*(A_beta_shig[i,j]+2*jalphaP[i,j] + B_beta_shig[i,j]+2*jbetaP[i,j])
            #D[i,j] = volume[i,j]*2 + dt*(A_beta_shig[i,j]+2*jalphaP[i,j] + B_beta_shig[i,j]+2*jbetaP[i,j])
        end
    end    
    
    #println(A_beta_shig[50,50])
    #println(B_beta_shig[50,50])
    #println(jalphaP[50,50])
    #println(jbetaP[50,50])
    
    
    RHS_temp = zeros(cellxmax,cellymax,4)
    # RHS
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                RHS_temp[i,j,l] = - (Qcon_hat[i,j,l]-Qconn_hat[i,j,l])*1.0 + dt*RHS[i,j,l]
            end
        end
    end

    #println(RHS_temp[:,50,:])
    #throw(UndefVarError(:x))

    # lower sweep
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                delta_Q[i,j,l] = D[i,j]^(-1) * (RHS_temp[i,j,l]+LdQ[i,j,l])
            end
        end
    end             

    # upepr sweep
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:4
                delta_Q[i,j,l] = delta_Q[i,j,l] - D[i,j]^(-1) * UdQ[i,j,l]
            end
        end
    end
    
    #println(Qcon_hat[:,50,:])
    #println(delta_Q[:,50,:])
    #throw(UndefVarError(:x))

    return delta_Q
end

