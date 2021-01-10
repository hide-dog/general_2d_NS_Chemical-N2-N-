# ------------------------------------
# Secure the array for explicit and implicit
# ------------------------------------
function common_allocation(cellxmax, cellymax)
    # define at cell center
    Qbase  = zeros(cellxmax, cellymax, nval)           # primitive variables
    volume = zeros(cellxmax, cellymax)                 # volume

    Qcon     = zeros(cellxmax, cellymax, nval)         # conserved variables
    Qcon_hat = zeros(cellxmax, cellymax, nval)         # conserved variables in the general coordinate system

    Rhat = zeros(cellxmax, cellymax)                   # gas const Considering chemical species

    chmu        = zeros(cellxmax, cellymax)            # viscosity
    chlambda_tr = zeros(cellxmax, cellymax)            # 熱伝導係数
    chD         = zeros(cellxmax, cellymax, nch)       # 拡散係数
    chi         = zeros(cellxmax, cellymax, nch)       # 化学種濃度
    
    wdot = zeros(cellxmax, cellymax, nval)             # 化学生成量
    RHS  = zeros(cellxmax, cellymax, nval)             # right hand side
    
    # define at cell boundaries
    dx = zeros(cellxmax+1, cellymax)                   # distance from cell center
    dy = zeros(cellxmax, cellymax+1)                   # distance from cell center

    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)    # flux of advection in the x-direction
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)    # flux of advection in the y-direction

    E_vis_hat = zeros(cellxmax+1,   cellymax, nval)    # flux of viscosity in the x-direction
    F_vis_hat = zeros(  cellxmax, cellymax+1, nval)    # flux of viscosity in the y-direction

    # for chemical reaction
    pisigma11 = zeros(nch, nch)                        # value of collision cross section
    pisigma22 = zeros(nch, nch)                        # value of collision cross section
    kf  = zeros(nre)                                    # forward reaction rate constant
    keq = zeros(nre)                                    # Equilibrium constant
    kb  = zeros(nre)                                    # backward reaction rate constant
    hs = zeros(nch)                                    # enthalpy of each chemical species
    Rs = zeros(nch)                                    # 化学種気体定数

    return Qbase, volume, dx, dy, Qcon, Qcon_hat, Rhat, chmu, chlambda_tr, chD, chi, 
            E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, wdot, RHS,
            pisigma11, pisigma22, kf, keq, kb, hs, Rs
end 

function allocation_implicit(cellxmax, cellymax)
    # define at cell center
    Qbasen     = zeros(cellxmax, cellymax, nval)        # primitive variables for inner iteration
    Qconn      = zeros(cellxmax, cellymax, nval)        # conserved variables for inner iteration
    Qconn_hat  = zeros(cellxmax, cellymax, nval)        # conserved variables in general coordinate system for inner iteration
    Qbasem     = zeros(cellxmax, cellymax, nval)        # primitive variables for inner iteration

    dtau         = zeros(cellxmax, cellymax)            # computational time steps

    A_adv_hat_p = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix A+ for one-wave approximation
    A_adv_hat_m = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix A- for one-wave approximation
    B_adv_hat_p = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix B+ for one-wave approximation
    B_adv_hat_m = zeros(cellxmax, cellymax, nval, nval) # Jacobian matrix B- for one-wave approximation
    A_beta_shig = zeros(cellxmax, cellymax)             # A sigma for one-wave approximation
    B_beta_shig = zeros(cellxmax, cellymax)             # B sigma for one-wave approximation

    jalphaP = zeros(cellxmax, cellymax)                 # Jacobian matrix for viscosity
    jbetaP  = zeros(cellxmax, cellymax)                 # Jacobian matrix for viscosity

    delta_Q      = zeros(cellxmax, cellymax, nval)      # delta Q for lusgs
    delta_Q_temp = zeros(cellxmax, cellymax, nval)      # temporary delta Q for lusgs

    D  = zeros(cellxmax, cellymax)                      # Diagonal of LHS
    Lx = zeros(cellxmax, cellymax, nval, nval)          # lower of LHS
    Ly = zeros(cellxmax, cellymax, nval, nval)          # lower of LHS
    Ux = zeros(cellxmax, cellymax, nval, nval)          # upper of LHS
    Uy = zeros(cellxmax, cellymax, nval, nval)          # upper of LHS

    LdQ = zeros(cellxmax, cellymax, nval)               # lower of LHS
    UdQ = zeros(cellxmax, cellymax, nval)               # upper of LHS

    RHS_temp = zeros(cellxmax, cellymax, nval)          # temporary RHS

    # define at cell boundaries
    lambda_facex = zeros(cellxmax+1, cellymax)          # lambda for computational time steps
    lambda_facey = zeros(cellxmax, cellymax+1)          # lambda for computational time steps

    # misc
    norm2 = zeros(nval)                                 # Residuals by norm-2
    I = zeros(nval, nval)                               # identity matrix
    for l in 1:nval
        I[l,l] = 1.0
    end

    return Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
            A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
            jalphaP, jbetaP, delta_Q, delta_Q_temp, D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp,
            norm2, I
end