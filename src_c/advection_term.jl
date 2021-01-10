function AUSM_plusup(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, Minf, volume, nval, nch)
    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)
    g         = specific_heat_ratio

    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            # i-1セル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            
            rhoL = Qbase[i-1,j,1]
            UL = Qcon[i-1,j,2] / Qcon[i-1,j,1]*cell_vecAx[1] + Qcon[i-1,j,3] / Qcon[i-1,j,1]*cell_vecAx[2]
            pL = Qbase[i-1,j,4]
            
            # iセル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            
            rhoR = Qbase[i,j,1]
            #UR = Qbase[i,j,2]*cell_vecAx[1] + Qbase[i,j,3]*cell_vecAx[2]
            UR = Qcon[i,j,2] / Qcon[i,j,1]*cell_vecAx[1] + Qcon[i,j,3] / Qcon[i,j,1]*cell_vecAx[2]
            pR = Qbase[i,j,4]
            
            # AUSM+up
            mdot, ph = AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g)
#=
            if rhoL == 0
                println(" L  ")
                println(i)
                println(j)
            end
            if rhoR == 0
                println("   ")
                println(i)
                println(j)
            end
=#
            # flux half
            temp_vecX = zeros(nval)
            temp_vecX[2] = vecAx[i,j,1]
            temp_vecX[3] = vecAx[i,j,2]

            Lpsi = zeros(nval)
            Rpsi = zeros(nval)

            for l in 1:nval
                Lpsi[l] = Qcon[i-1,j,l] / Qcon[i-1,j,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i-1,j,4] + Qbase[i-1,j,4]) / Qcon[i-1,j,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]

            if mdot > 0
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecX[l]
                end
            else
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecX[l]
                end
            end            
        end
    end

    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
        
            # j-1セル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            
            rhoL = Qbase[i,j-1,1]
            VL = Qbase[i,j-1,2]*cell_vecAy[1] + Qbase[i,j-1,3]*cell_vecAy[2]
            pL = Qbase[i,j-1,4]
            
            # jセル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            
            rhoR = Qbase[i,j,1]
            VR = Qbase[i,j,2]*cell_vecAy[1] + Qbase[i,j,3]*cell_vecAy[2]
            pR = Qbase[i,j,4]

            # AUSM+up
            mdot, ph = AUSM_plusup_half(rhoL, rhoR, VL, VR, pL, pR, Minf, g)
            #=
            if i == 100 && j == 50
                println(" y ")
                println(mdot)
                println(ph)
            end
            =#
            
            # flux half
            temp_vecY = zeros(nval)
            temp_vecY[2] = vecAy[i,j,1]
            temp_vecY[3] = vecAy[i,j,2]

            Lpsi = zeros(nval)
            Rpsi = zeros(nval)

            for l in 1:nval
                Lpsi[l] = Qcon[i,j-1,l] / Qcon[i,j-1,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i,j-1,4] + Qbase[i,j-1,4]) / Qcon[i,j-1,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]
            
            if mdot > 0
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecY[l]
                end
            else
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecY[l]
                end
            end
            
            #=
            if i == 35 && j ==2
                println(" y ")
                println(mdot)
                println(temp_vecY)
                println(Lpsi)
                println(Rpsi)
            end
            =#

        end
    end
    
    return E_adv_hat, F_adv_hat
end

function AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g)
    # param
    beta  = 1/8
    kp    = 0.25
    sigma = 1.0
    ku    = 0.75
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    Mbar = (( UL^2 + UR^2 ) / ( 2 * ah^2 ))^0.5
    Mo   = (min(1,max( Mbar^2, Minf^2 )))^0.5
    fa   = Mo * (2-Mo)

    alpha = 3/16 * ( -4 + 5*fa )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        Mtp  = 0.25*(ML + 1)^2
        Mtm  = -0.25*(ML - 1)^2
        #M2m1 = (ML^2 - 1)^2
        M_p4 = Mtp * (1 - 16*beta*Mtm)
        #M_p4 = (1-ML)* (Mtp + beta*M2m1) + ML*0.5*(ML + abs(ML))
        p_p5 = Mtp * ((2-ML) - 16*alpha*ML*Mtm)
    end
    

    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        Mtp  = 0.25*(MR + 1)^2
        Mtm  = -0.25*(MR - 1)^2
        M_m4 = Mtm * (1 + 16*beta*Mtp)
        p_m5 = Mtm * ((-2-MR) + 16*alpha*MR*Mtp)
    end
    
    #=
    Mtp  = 0.25*(ML + 1)^2
    Mtm  = -0.25*(ML - 1)^2
    M2m1 = (ML^2 - 1)^2
    M_p4 = (1-ML)* (Mtp + beta*M2m1) + ML*0.5*(ML + abs(ML))
    
    Mtp  = 0.25*(MR + 1)^2
    Mtm  = -0.25*(MR - 1)^2
    M2m1 = (MR^2 - 1)^2
    M_m4 = (1-MR)* (Mtm - beta*M2m1) + MR*0.5*(MR - abs(MR))
    =#

    # M half
    Mp = -kp * max( 1 - sigma*Mbar^2, 0) * (pR - pL)/(rhoh*ah^2)

    # Mh = M_p4 + M_m4 + Mp/fa
    Mh = M_p4 + M_m4
    
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    pu = -ku * p_p5 * p_m5 * (rhoL + rhoR) * ah *(UR - UL)

    ph = p_p5*pL + p_m5*pR + fa * pu
    return mdot, ph, ah, Mh
end

function AUSM(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio)
    E_adv_hat = zeros(cellxmax+1,cellymax,4)
    F_adv_hat = zeros(cellxmax,cellymax+1,4)
    
    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            # i-1セル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            
            rho = Qbase[i-1,j,1]
            U = Qbase[i-1,j,2]*cell_vecAx[1] + Qbase[i-1,j,3]*cell_vecAx[2]
            p = Qbase[i-1,j,4]
            a_im1 = (specific_heat_ratio * p / rho)^0.5
            M = U/a_im1

            if abs(M) <= 1
                M_p = 0.25 * (M+1)^2
                p_p = 0.25 *p * (1+M)^2 *(2-M)
            elseif abs(M) > 1
                M_p = 0.5 * (M+abs(M))
                p_p = 0.5 *p * (M+abs(M))/M
            end

            # iセル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            

            rho = Qbase[i,j,1]
            U = Qbase[i,j,2]*cell_vecAx[1] + Qbase[i,j,3]*cell_vecAx[2]
            p = Qbase[i,j,4]
            a_i = (specific_heat_ratio * p / rho)^0.5
            M = U/a_i

            if abs(M) <= 1
                M_m = -0.25 * (M-1)^2
                p_m = 0.25 *p * (1-M)^2 *(2+M)
            elseif abs(M) > 1
                M_m = 0.5 * (M-abs(M))
                p_m = 0.5 *p * (M-abs(M))/M
            end

            # bound
            M_half = M_p + M_m
            p_half = p_p + p_m

            temp_vecX = zeros(4)
            temp_vecX[2] = vecAx[i,j,1]
            temp_vecX[3] = vecAx[i,j,2]

            LaQc = zeros(4)
            RaQc = zeros(4)
            for l in 1:3
                LaQc[l] = a_im1*Qcon[i-1,j,l]
                RaQc[l] = a_i*Qcon[i,j,l]
            end
            # hat{E}+1/J*p
            LaQc[4] = a_im1*(Qcon[i-1,j,4]+Qbase[i-1,j,4])
            RaQc[4] = a_i*(Qcon[i,j,4]+Qbase[i,j,4])

            for l in 1:4
                E_adv_hat[i,j,l] = 0.5 * M_half*(LaQc[l] + RaQc[l]) - 0.5*abs(M_half)*(RaQc[l] - LaQc[l]) + p_half*temp_vecX[l]
            end
        end
    end

    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
        
            # j-1セル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            

            rho = Qbase[i,j-1,1]
            V = Qbase[i,j-1,2]*cell_vecAy[1] + Qbase[i,j-1,3]*cell_vecAy[2]
            p = Qbase[i,j-1,4]
            a_jm1 = (specific_heat_ratio * p / rho)^0.5
            M = V/a_jm1

            if abs(M) <= 1
                M_p = 0.25 * (M+1)^2
                p_p = 0.25 *p * (1+M)^2 *(2-M)
            elseif abs(M) > 1
                M_p = 0.5 * (M+abs(M))
                p_p = 0.5 *p * (M+abs(M))/M
            end

            # jセル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            
            rho = Qbase[i,j,1]
            V = Qbase[i,j,2]*cell_vecAy[1] + Qbase[i,j,3]*cell_vecAy[2]
            p = Qbase[i,j,4]
            a_j = (specific_heat_ratio * p / rho)^0.5
            M = V/a_j

            if abs(M) <= 1
                M_m = -0.25 * (M-1)^2
                p_m = 0.25 *p * (1-M)^2 *(2+M)
            elseif abs(M) > 1
                M_m = 0.5 * (M-abs(M))
                p_m = 0.5 *p * (M-abs(M))/M
            end

            # bound
            M_half = M_p + M_m
            p_half = p_p + p_m

            temp_vecY = zeros(4)
            temp_vecY[2] = vecAy[i,j,1]
            temp_vecY[3] = vecAy[i,j,2]
            
            LaQc = zeros(4)
            RaQc = zeros(4)
            for l in 1:3
                LaQc[l] = a_jm1*Qcon[i-1,j,l]
                RaQc[l] = a_j*Qcon[i,j,l]
            end
            # hat{E}+1/J*p
            LaQc[4] = a_jm1*(Qcon[i-1,j,4]+Qbase[i-1,j,4])
            RaQc[4] = a_j*(Qcon[i,j,4]+Qbase[i,j,4])


            for l in 1:4
                F_adv_hat[i,j,l] = 0.5 * M_half*(LaQc[l] + RaQc[l]) - 0.5*abs(M_half)*(RaQc[l] - LaQc[l]) + p_half*temp_vecY[l]
            end
        end
    end

    return E_adv_hat, F_adv_hat
end

function setup_cell_flux_hat(Qcon,cellxmax,cellymax,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus)    
    cell_E_hat_plas  = zeros(cellxmax,cellymax,4)
    cell_E_hat_minus = zeros(cellxmax,cellymax,4)
    cell_F_hat_plas  = zeros(cellxmax,cellymax,4)
    cell_F_hat_minus = zeros(cellxmax,cellymax,4)

    for i in 1:cellxmax
        for j in 1:cellymax
            for l in 1:4
                for m in 1:4
                    # Ehat = Ahat * Qhat = Ahat * Q/J = 1/J*Ahat * Q
                    # 1/J*Ahatの計算をvecAを使ったため，ここではQconを掛けてる
                    cell_E_hat_plas[i,j,l]  += cell_Ahat_plas[i,j,l,m]*Qcon[i,j,m]
                    cell_E_hat_minus[i,j,l] += cell_Ahat_minus[i,j,l,m]*Qcon[i,j,m]
                    cell_F_hat_plas[i,j,l]  += cell_Bhat_plas[i,j,l,m]*Qcon[i,j,m]
                    cell_F_hat_minus[i,j,l] += cell_Bhat_minus[i,j,l,m]*Qcon[i,j,m]
                end
            end
        end
    end

    return cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus
end

function cal_jacobi(Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy,volume)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    
    q=u^2+v^2
    c=(rp/rho)^0.5
    theta=Z=U= A^xi_x * u   +   A^xi_y * v
    H=(e+p)/rho
    b1=q^2/2 * (r-1)/(c^2)
    b2=(r-1)/(c^2)

    xi_x_bar= xi_x / (xi_x^2+xi_y^2)^0.5
    xi_y_bar= xi_y / (xi_x^2+xi_y^2)^0.5
    U_bar=Z_bar=theta_bar=


    A = Rleft * Lambda * Rright
    """

    cell_Ahat_plas  = zeros(cellxmax,cellymax,4,4)
    cell_Ahat_minus = zeros(cellxmax,cellymax,4,4)
    cell_Bhat_plas  = zeros(cellxmax,cellymax,4,4)
    cell_Bhat_minus = zeros(cellxmax,cellymax,4,4)
    r=specific_heat_ratio

    for i in 1:cellxmax
        for j in 1:cellymax
            q=(Qbase[i,j,2]^2+Qbase[i,j,3]^2)^0.5
            c=0
            try c=(r*Qbase[i,j,4]/Qbase[i,j,1])^0.5
            catch
                println("\n"*string(i)*","*string(j)*" Qbase error")
                println(Qbase[i,j,:])
                println("\n")
                throw(UndefVarError(:x))
            end

            H  = (Qcon[i,j,4]+Qbase[i,j,4])/Qbase[i,j,1]
            b1 = q^2/2 * (r-1)/(c^2)
            b2 = (r-1)/(c^2)
            u  = Qbase[i,j,2]
            v  = Qbase[i,j,3]

            k_x = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            k_y = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            
            Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
            A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

            for l in 1:4
                for m in 1:4
                    cell_Ahat_plas[i,j,l,m]  = A_plas[l,m]
                    cell_Ahat_minus[i,j,l,m] = A_minus[l,m]
                end
            end

            k_x = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            k_y = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])

            Rleft_E,Lambda_E,Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
            A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

            for l in 1:4
                for m in 1:4
                    cell_Bhat_plas[i,j,l,m]  = A_plas[l,m]
                    cell_Bhat_minus[i,j,l,m] = A_minus[l,m]
                end
            end
        end
    end

    return cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus
end


function Eigenvalue_vector(q,c,H,b1,b2,u,v,k_x,k_y)
    Z= k_x*u + k_y*v
    k_x_bar=k_x / (k_x^2+k_y^2)^0.5
    k_y_bar=k_y / (k_x^2+k_y^2)^0.5
    Z_bar=Z / (k_x^2+k_y^2)^0.5
    
    Lambda = [Z-c*(k_x^2+k_y^2)^0.5 
              Z
              Z+c*(k_x^2+k_y^2)^0.5 
              Z]

    Rleft = [[1 1 1 0]
             [u-k_x_bar*c u u+k_x_bar*c -k_y_bar]
             [v-k_y_bar*c v v+k_y_bar*c k_x_bar]
             [H-c*Z_bar q^2/2 H+c*Z_bar -(k_y_bar*u-k_x_bar*v)]]

    Rright = [[(b1+Z_bar/c)/2 -(k_x_bar/c+b2*u)/2 -(k_y_bar/c+b2*v)/2 b2/2]
              [1-b1 b2*u b2*v -b2]
              [(b1-Z_bar/c)/2 (k_x_bar/c-b2*u)/2 (k_y_bar/c-b2*v)/2 b2/2]
              [k_y_bar*u-k_x_bar*v -k_y_bar k_x_bar 0]]
    
    return Rleft,Lambda,Rright
end

function RLmbdaR(R,Lam,Rm)
    Lamp=zeros(4,4)
    Lamm=zeros(4,4)
    for i in 1:length(Lam)
        Lamp[i,i]=(Lam[i]+abs(Lam[i]))/2
        Lamm[i,i]=(Lam[i]-abs(Lam[i]))/2
    end
    Ap=nn_inner_product(nn_inner_product(R,Lamp),Rm)
    Am=nn_inner_product(nn_inner_product(R,Lamm),Rm)
    return Ap,Am
end

function nn_inner_product(a,b)
    temp=zeros(size(a)[1],size(a)[1])
    for i in 1:size(a)[1]
        for j in 1:size(a)[1]
            for k in 1:size(a)[1]
                temp[i,j] += a[i,k]*b[k,j]
            end
        end
    end
    return temp    
end


function FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,cellxmax,cellymax)
    E_hat = zeros(cellxmax+1,cellymax,4)
    F_hat = zeros(cellxmax,cellymax+1,4)

    Threads.@threads for i in 2:cellxmax+1-1
        for j in 2:cellymax-1
            for k in 1:4
                E_hat[i,j,k] = cell_E_hat_plas[i-1,j,k] + cell_E_hat_minus[i,j,k]
            end
        end
    end
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax+1-1
            for k in 1:4
                F_hat[i,j,k] = cell_F_hat_plas[i,j-1,k] + cell_F_hat_minus[i,j,k]
            end
        end
    end

    return E_hat,F_hat
end