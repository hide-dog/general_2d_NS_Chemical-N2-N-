function set_wdot(Qbase, cellxmax, cellymax, Rd, nch, nre, nval)
    
    wdot = zeros(cellxmax, cellymax, nval)
    npre = nval - nch

    mw_N2_N = [28e-3, 14e-3]
    avgdro  = 6.022*10^23

    # N2 + N2 --> N + N + N2
    # N2 + N  --> N + N + N
    fmol = [[2, 0],
            [1, 1]]
    # N2 + N2 <-- N + N + N2
    # N2 + N  <-- N + N + N
    bmol = [[1, 2],
            [0, 3]]

    for i in 2:cellxmax-1
        for j in 2:cellymax-1
        
            # Initial setting
            total = 0.0
            rho_mw = zeros(nch)
            for s in 1:nch
                rho_mw[s] = Qbase[i,j,npre+s] / mw_N2_N[s]   # [mol/m3] = [kg/m3] / [kg/mol] 
                total     = total + rho_mw[s]
            end
            total = log10(total*avgdro)

            if total <= 20.0
                nsb = 1
            elseif total <= 21.0
                nsb = 2
            elseif total <= 22.0
                nsb = 3
            elseif total <= 23.0
                nsb = 4
            elseif total <= 24.0
                nsb = 5
            else
                nsb = 6
            end
            
            T   = Qbase[i,j,npre]/(Qbase[i,j,1]*Rd)
            
            keq = k_eq(T, nsb, nre)                 # [mol/(cm^3 s)]
            kf  = freaction_rate(T, nre)            # [mol/(cm^3 s)]
            kb  = breaction_rate(kf, keq, nre)      # [mol/(cm^3 s)]

            # Forward and backward rate for reaction ir
            lf = ones(nre)
            lb = ones(nre)
            for r in 1:nre
                for s in 1:nch
                    # 並行定数の単位に合わせるため，モル濃度を[mol/cm^3]の単位で計算している
                    lf[r] = lf[r]*1.0e-6*rho_mw[s]^fmol[r][s]
                    lb[r] = lb[r]*1.0e-6*rho_mw[s]^bmol[r][s]
                end
                lf[r] = lf[r]*kf[r]
                lb[r] = lb[r]*kb[r]
            end

            for s in 1:nch
                wdot[i,j,npre+s] = 0.0
                for r in 1:nre
                    temp = (fmol[r][s] - bmol[r][s])*(lf[r] - lb[r])
                    wdot[i,j,npre+s] += temp 
                end
                # モル濃度[mol/cm^3]から[mol/m^3]に変換するため，1.0e6を乗じている，
                wdot[i,j,npre+s] = 1.0e6 * wdot[i,j,npre+s] * mw_N2_N[s]
            end
        end
    end
    return wdot
end

function chemical_value(Qbase, cellxmax, cellymax, Rd, nval, nch)

    chmu        = zeros(cellxmax, cellymax)
    chlambda_tr = zeros(cellxmax, cellymax)
    chD         = zeros(cellxmax, cellymax, nch)
    chi         = zeros(cellxmax, cellymax, nch)
    # N2 - N2
    # N2 - N
    # N  - N
    num_nncross = binomial(BigInt(2*nch-1), BigInt(nch)) # 重複組み合わせ nHr = n+r-1Cr
    npre = nval - nch

    mw = [28.0e-3, 14.0e-3]  # 原子量[kg/mol]
    avgdro  = 6.022*10^23    # アボガドロ定数
    kb = 1.380/10^(23)      # ボルツマン定数
    pi = 3.1415              # pi

    ms = zeros(nch)          # 原子一個当たりの質量[kg]
    for s in 1:nch
        ms[s] = mw[s] / avgdro
    end

    for i in 1:cellxmax
        for j in 1:cellymax
        
            total = 0.0
            for s in 1:nch
                total += Qbase[i,j,npre+s]/mw[s]
            end
            for s in 1:nch
                chi[i,j,s] = (Qbase[i,j,npre+s]/mw[s]) /total
            end

            T = Qbase[i,j,npre]/(Qbase[i,j,1]*Rd)
            p = Qbase[i,j,npre]
            
            # 対象行列に格納
            # [ N2-N2   N2-N ]
            # [   0.0    N-N ]
            #
            ps11 = collision_cross_section_pisigma11(T, nch)
            ps22 = collision_cross_section_pisigma22(T, nch)
        
            deltaij1 = zeros(nch, nch)
            deltaij2 = zeros(nch, nch)
            
            # 以下では行列のすべてを計算しているが、下三角は0なのでキャンセルされる
            for si in 1:nch
                for sj in 1:nch
                    mi = ms[si]
                    mj = ms[sj]
                        
                    temp = 2*mi*mj / (pi*kb*T*(mi+mj))
                    deltaij1[si,sj] = 8/3 * temp^0.5 * ps11[si,sj]
                    deltaij2[si,sj] = 16/5 * temp^0.5 * ps22[si,sj]
                end
            end

            chmu[i,j] = 0.0
            temp      = 0.0
            for si in 1:nch
                for sj in 1:nch
                    temp += chi[i,j,sj] * deltaij2[si,sj]
                end
                chmu[i,j] += ms[si]*chi[i,j,si]/temp
            end


            chlambda_tr[i,j] = 0.0
            temp =0.0
            for si in 1:nch
                for sj in 1:nch
                    alphaij = 1 + (1-ms[si]/ms[sj]) * (0.45-2.54*ms[si]/ms[sj]) / (1+ms[si]/ms[sj])^2
                    temp += alphaij * chi[i,j,sj] * deltaij2[si, sj]
                end
                chlambda_tr[i,j] += chi[i,j,si]/temp
            end
            chlambda_tr[i,j] = 15/4 * kb * chlambda_tr[i,j]

            # 対象の行と列を足して(si,si)を引く
            # [ N2-N2   N2-N ]
            # [   0.0    N-N ]
            #
            #=
            (例) N2
                sum[1,si] + sum[si,1] - [si,si]*2 = [ N2-N2  N2-N ] + [ N2-N2 0.0 ] - [N2-N2]*2
                                                  = [ N2-N ]
            
            Dijを導入すると0割りが発生するため、inverse Dijを計算している
            =#
            inv_Dij = zeros(nch,nch)
            for si in 1:nch
                for sj in 1:nch
                    inv_Dij[si,sj] = (p*deltaij1[si,sj])/kb*T
                end
            end

            for si in 1:nch
                tempsum = 0
                for sj in 1:nch
                    tempsum += chi[i,j,sj]*inv_Dij[si,sj]
                    tempsum += chi[i,j,sj]*inv_Dij[sj,si]
                end
                tempsum -= chi[i,j,si]*inv_Dij[si,si]*2
                rho  = Qbase[i,j,1]
                rhos = Qbase[i,j,npre+si]
                Cs   = rhos / rho
                if chi[i,j,si] == 0.0 || tempsum == 0.0
                    chD[i,j,si] = 0.0
                else
                    chD[i,j,si] = (1-Cs)/tempsum * Cs/chi[i,j,si]
                end
            end
        end
    end

    return chmu, chlambda_tr, chD, chi
end
