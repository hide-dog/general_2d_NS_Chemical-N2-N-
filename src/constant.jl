# ------------------------------------
# define constant 
# ------------------------------------
const nval = 6          # number of conserved variables
const nch  = 2          # number of chemical species
const npr  = nval - nch # number of primitive variables (rho, u, v, p)
const nre  = 2          # number of chemical reactions
const nmol = 1          # 二原子分子の数
const R    = 8.314      # gas constant, J/(K mol)

const avgdro  = 6.022*10^23    # アボガドロ定数
const kb_const = 1.380/10^(23) # ボルツマン定数
const pi = 3.141592            # pi

function array_const(Rs)
    # N2,N
    delta_hs_0 = [0.0, 3.364e7]   # [J/kg]
    mw         = [28e-3, 14e-3]   # [kg/mol]

    # N2
    theta_vib    = [3.353]
    
    # 化学種の気体定数
    for s in 1:nch
        Rs[s] = R / mw[s]
    end
    
    return delta_hs_0, mw, theta_vib, Rs
end

# ------------------------------------
# collision_cross_section
#
# A review of reaction rates and thermodynamic and transport properties 
# for an 11-species air model for chemical and 
# thermal nosequilibrium calculations to 30,000 K
#
# ------------------------------------
function ccs_const()
    # N2-N2
    # N2-N
    # N-N
    ccs11_A = [0.0, 0.0, 0.0]
    ccs11_B = [-0.0112, -0.0194, -0.0033]
    ccs11_C = [-0.1182, 0.0119, -0.0572]
    ccs11_D = [4.8464, 4.1055, 5.0452]

    ccs22_A = [0.0, 0.0, 0.0]
    ccs22_B = [-0.0203, -0.0190, -0.0118]
    ccs22_C = [-0.0683, 0.0239, -0.0960]
    ccs22_D = [4.0900, 4.1782, 4.3252]

    ccs11_A = Array{Float16}(ccs11_A)
    ccs11_B = Array{Float16}(ccs11_B)
    ccs11_C = Array{Float16}(ccs11_C)
    ccs11_D = Array{Float16}(ccs11_D)
    ccs22_A = Array{Float16}(ccs22_A)
    ccs22_B = Array{Float16}(ccs22_B)
    ccs22_C = Array{Float16}(ccs22_C)
    ccs22_D = Array{Float16}(ccs22_D)
    return ccs11_A, ccs11_B, ccs11_C, ccs11_D, ccs22_A, ccs22_B, ccs22_C, ccs22_D
end

function rr()
    # N2+N2 -> 2N+N2
    # N2+N  -> 2N+N
    rr_Cr = [7.0*10^(21), 3.0*10^(22)]
    rr_sr = [-1.60, -1.60]
    rr_tr = [113200, 113200]
    
    rr_Cr = Array{Float16}(rr_Cr)
    rr_sr = Array{Float16}(rr_sr)
    rr_tr = Array{Float16}(rr_tr)
    return rr_Cr, rr_sr, rr_tr
end

# ------------------------------------
# Equilibrium constants curve-fits 
#
# C. Park, “Nonequilibrium Hypersonic Aerothermodynamics”, Wiley, New York, 1990.
#
# ------------------------------------
function eq_const()
    
    # N2+N2 -> 2N+N2
    ar1 = [[3.490700, 2.072300, 1.606000, 1.535100, 1.476600, 1.476600],
        [3.490700, 2.072300, 1.606000, 1.535100, 1.476600, 1.476600]]
            
    ar2 = [[8.313300, 1.389700, 1.573200, 1.606100, 1.629100, 1.629100],
        [8.313300, 1.389700, 1.573200, 1.606100, 1.629100, 1.629100]]

    ar3 = [[4.097800, 2.061700, 1.392300, 1.299300, 1.215300, 1.215300],
        [4.097800, 2.061700, 1.392300, 1.299300, 1.215300, 1.215300]]

    ar4 = [[-1.272800, -1.182800, -1.153300, -1.149400, -1.145700, -1.145700],
        [-1.272800, -1.182800, -1.153300, -1.149400, -1.145700, -1.145700]]

    ar5 = [[7.487000, 1.510500, -4.543000, -6.980000, -9.444000, -9.444000],
        [7.487000, 1.510500, -4.543000, -6.980000, -9.444000, -9.444000]]

    ar1 = Array{Float16}(ar1)
    ar2 = Array{Float16}(ar2)
    ar3 = Array{Float16}(ar3)
    ar4 = Array{Float16}(ar4)
    ar5 = Array{Float16}(ar5)

    return ar1, ar2, ar3, ar4, ar5 
end
