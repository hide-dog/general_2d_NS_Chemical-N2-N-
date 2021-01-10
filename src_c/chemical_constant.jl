function collision_cross_section_pisigma11(T, nch)

    # pisigma
    pisigma11 = zeros(nch, nch)
    
    # [A, B, C, D]
    A = 0.0
    B = 0.0
    C = 0.0
    D = 0.0

    # A review of reaction rates and thermodynamic and transport properties 
    # for an 11-species air model for chemical and 
    # thermal nosequilibrium calculations to 30,000 K

    # N2-N2
    A = 0.0
    B = -0.0112
    C = -0.1182
    D = 4.8464
    
    pisigma11[1,1] = pisigma(T, A, B, C, D)

    # N2-N
    A = 0.0
    B = -0.0194
    C = 0.0119
    D = 4.1055

    pisigma11[1,2] = pisigma(T, A, B, C, D)
    pisigma11[2,1] = 0.0

    # N-N
    A = 0.0
    B = -0.0033
    C = -0.0572
    D = 5.0452

    pisigma11[2,2] = pisigma(T, A, B, C, D)
    
    return pisigma11
end

function collision_cross_section_pisigma22(T, nch)

    # pisigma
    pisigma22 = zeros(nch, nch)
    
    # [A, B, C, D]
    A = 0.0
    B = 0.0
    C = 0.0
    D = 0.0

    # A review of reaction rates and thermodynamic and transport properties 
    # for an 11-species air model for chemical and 
    # thermal nosequilibrium calculations to 30,000 K

    # N2-N2
    A = 0.0
    B = -0.0203
    C = 0.0683
    D = 4.0900
    
    pisigma22[1,1] = pisigma(T, A, B, C, D)

    # N2-N
    A = 0.0
    B = -0.0190
    C = 0.0239
    D = 4.1782

    pisigma22[1,2] = pisigma(T, A, B, C, D)
    pisigma22[2,1] = 0.0

    # N-N
    A = 0.0
    B = -0.0118
    C = -0.0960
    D = 4.3252

    pisigma22[2,2] = pisigma(T, A, B, C, D)
    
    return pisigma22
end

function pisigma(T, A, B, C, D)
    ans = exp(D) * T^(A^2*log(T) + B*log(T) + C)
    return ans
end

function freaction_rate(Tfr, nreaction)
    kf = zeros(nreaction)

    # N2+N2 -> 2N+N2
    Cr = 7.0*10^(21)
    sr = -1.60
    tr = 113200

    kf[1] = reaction_rate_arrhcnius(Tfr, Cr, sr, tr)

    # N2+N  -> 2N+N
    Cr = 3.0*10^(22)
    sr = -1.60
    tr = 113200

    kf[2] = reaction_rate_arrhcnius(Tfr, Cr, sr, tr)

    return kf
end

function reaction_rate_arrhcnius(Tfr, Cr, sr, tr)
    k = Cr * Tfr^(sr) * exp(-tr/Tfr)
    return k
end

function k_eq(T, nsb, nreaction)
    keq = zeros(nreaction)
    
    # Equilibrium constants curve-fits 
    # C. Park, “Nonequilibrium Hypersonic Aerothermodynamics”, Wiley, New York, 1990.

    # N2+N2 -> 2N+N2
    ar1 = [3.490700,
            2.072300,
            1.606000,
            1.535100,
            1.476600,
            1.476600]
            
    ar2 = [8.313300,
            1.389700,
            1.573200,
            1.606100,
            1.629100,
            1.629100]
    
    ar3 = [4.097800,
            2.061700,
            1.392300,
            1.299300,
            1.215300,
            1.215300]

    ar4 = [-1.272800,
            -1.182800,
            -1.153300,
            -1.149400,
            -1.145700,
            -1.145700]
    
    ar5 = [7.487000,
            1.510500,
            -4.543000,
            -6.980000,
            -9.444000,
            -9.444000]


    keq[1] = keq_curve_fitting(T, ar1[nsb], ar2[nsb], ar3[nsb], ar4[nsb], ar5[nsb])

    # N2+N  -> 2N+N
    # 上式と一緒
    
    keq[2] = keq_curve_fitting(T, ar1[nsb], ar2[nsb], ar3[nsb], ar4[nsb], ar5[nsb])
        
    return keq
end

function keq_curve_fitting(T, A1, A2, A3, A4, A5)
    Z = 10^4/T
    keq = exp(A1/Z + A2 + A3*log(Z) + A4*Z + A5*Z^2)
    keq = min(keq, 700.0)
    return keq
end

function breaction_rate(kf, keq, nreaction)
    kb = zeros(nreaction)
    for i in 1:nreaction
        kb[i] = kf[i] / keq[i]
    end
    return kb
end

function chemical_enthalpy(T)
    # 複数の温度モードを考慮した方がよい
    delta_hs_0 = [0.0, 3.364e7]  # N2,N[J/kg]
    mw  = [28e-3, 14e-3]          # [kg/mol]
    R   = 8.314                   # [J/K/mol]
    nch = length(delta_hs_0)

    num_molecule = 1             # 二原子分子の数
    theta_vib = [3.353]

    Rs = zeros(nch)              # 化学種気体定数
    for s in 1:nch
        Rs[s] = R / mw[s]
    end

    hs = zeros(nch)
    for s in 1:num_molecule
        etr   = 1.5 * Rs[s] * T
        erot  = Rs[s] * T
        evib  = Rs[s] * theta_vib[s] / (exp(theta_vib[s]/T)-1)
        hs[s] = etr + erot + evib + delta_hs_0[s] + Rs[s]*T
    end
    for s in num_molecule+1:nch
        etr   = 1.5 * Rs[s] * T
        hs[s] = etr + delta_hs_0[s] + Rs[s]*T
    end

    return hs
end