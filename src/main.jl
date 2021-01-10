using ProgressMeter
using Dates

function main()
    
    # Start of time measurement
    start_t = now()

    out_dir  = "result"         # output dir
    PARAMDAT = "PARAMDAT.json"  # directory
    fwrite   = "write"          # write file 
    
    # read grids and parameter
    xmax, ymax, nodes, vecAx, vecAy = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl,
    init_rho, init_u, init_v, init_p, init_N2, init_N, init_T,
    specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)
    
    # number of cells
    cellxmax = xmax - 1
    cellymax = ymax - 1
    
    # allocation
    Qbase, volume, dx, dy, Qcon, Qcon_hat, Rhat, chmu, chlambda_tr, chD, chi, 
    E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, wdot, RHS,
    pisigma11, pisigma22, kf, keq, kb, hs, Rs    = common_allocation(cellxmax, cellymax)

    delta_hs_0, mw, theta_vib, Rs = array_const(Rs)
    
    # set initial condition
    Qbase, restartnum = set_initQbase(Qbase, cellxmax, cellymax, restart_file, init_rho, init_u, init_v, init_p, init_T,
                                init_N2, init_N, specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rs)
    
    # set volume, dx and dy
    volume = set_volume(nodes, cellxmax, cellymax, volume)
    dx, dy = set_dx_lts(dx, dy, nodes, cellxmax, cellymax)
    reset_write(fwrite)

    # write number of threads
    print("threads num : ")
    println(Threads.nthreads())

    # check boundary condition
    check_bd(bdcon)
        
    # main loop
    if time_integ == "1"
        # exlicit scheme
        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            # step number
            evalnum = t + restartnum

            Rhat = set_gasconst(Qbase, cellxmax, cellymax, Rs)
            
            # set conserved variables in the general coordinate system
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, specific_heat_ratio, Rhat)
            Qcon     = base_to_conservative(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume)
            
            # set viscosity and thermal Conductivity
            chmu, chlambda_tr, chD, chi = chemical_value(chmu, chlambda_tr, chD, chi, Qbase, cellxmax, cellymax)
            
            # advection_term
            E_adv_hat, F_adv_hat = AUSM_plus(E_adv_hat, F_adv_hat, Qbase, Qcon, cellxmax, cellymax, 
                                            vecAx, vecAy, specific_heat_ratio, volume, nval)
                        
            # viscos_term
            E_vis_hat, F_vis_hat = central_diff(E_vis_hat, F_vis_hat, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
            
            # sorce term
            wdot = set_wdot(Qbase, cellxmax, cellymax, Rd, nch, nre, nval)
            
            # RHS
            RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)
            
            # calculate primitive variables
            Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            Qbase = conservative_to_base(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            
            # output
            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
            end
            
            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
        end
    elseif time_integ == "2"
        # implicit

        # allocation for implicit
        Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
        A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
        jalphaP, jbetaP, delta_Q, delta_Q_temp, D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp,
        norm2, I = allocation_implicit(cellxmax, cellymax, nval)

        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            # step number
            evalnum = t + restartnum

            Rhat = set_gasconst(Qbase, cellxmax, cellymax, Rs)
            
            # write physical time 
            output_physicaltime(fwrite, t, dt)
            
            # copy
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbasen[i,j,l] = Qbase[i,j,l]
                        Qbasem[i,j,l] = Qbase[i,j,l]
                    end
                end
            end

            # set conserved variables in the general coordinate system for inner iteration
            Qconn     = base_to_conservative(Qbasen, Qconn, cellxmax, cellymax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, Qconn_hat, cellxmax, cellymax, volume, nval)     

            # start inner iteration
            for tau in 1:in_nt
                # set conserved variables in the general coordinate system
                Qbasem = set_boundary(Qbasem, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval)
                Qcon     = base_to_conservative(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)

                # set viscosity and thermal Conductivity
                mu     = set_mu(mu, Qbasem, cellxmax, cellymax, specific_heat_ratio, Rd)
                lambda = set_lambda(lambda, Qbasem, cellxmax, cellymax, mu, specific_heat_ratio, Rd)

                # set inner time step by local time stepping
                dtau   = set_lts(dtau, lambda_facex, lambda_facey, Qbase, cellxmax, cellymax, mu, dx, dy,
                                vecAx, vecAy, volume, specific_heat_ratio, cfl)
                                
                #advection_term
                E_adv_hat, F_adv_hat = AUSM_plus(E_adv_hat, F_adv_hat, Qbasem, Qcon, cellxmax, cellymax, 
                                                vecAx, vecAy, specific_heat_ratio, volume, nval)

                # viscos_term
                E_vis_hat, F_vis_hat = central_diff(E_vis_hat, F_vis_hat, Qbasem, Qcon, cellxmax, cellymax, mu, lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
                
                # RHS
                RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume)
            
                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, 
                A_beta_shig, B_beta_shig = one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig, I,
                                                    Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
                # lusgs_viscos_term
                jalphaP, jbetaP = central_diff_jacobian(jalphaP, jbetaP, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                                                        vecAx, vecAy, specific_heat_ratio, volume, nval)
                
                # LUSGS
                ite = 0
                while true
                    # copy delta_Q
                    for l in 1:nval
                        for j in 1:cellymax
                            for i in 1:cellxmax
                                delta_Q_temp[i,j,l] = delta_Q[i,j,l]
                            end
                        end
                    end
                    
                    # Reversing the left-hand side by lusgs
                    delta_Q = lusgs(D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp, I, dt, dtau, Qcon_hat, Qconn_hat, delta_Q,
                                    A_adv_hat_p,  A_adv_hat_m,  B_adv_hat_p,  B_adv_hat_m,  A_beta_shig,  B_beta_shig,
                                     jalphaP,  jbetaP, RHS, cellxmax, cellymax, volume, nval)
                    
                    # cal Residuals by norm-2
                    res   = set_res(delta_Q, delta_Q_temp, cellxmax, cellymax, nval)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax, init_small, nval)
                    
                    # Find out if the results converged
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok
                        break
                    end
                    if ite % 100 ==0
                        println(" now cal norm2 ")
                        println(norm2)
                    end

                    ite += 1
                end

                # output inner time
                output_innertime(fwrite, tau, norm2, nval)
                
                # Updating
                for l in 1:nval
                    for j in 2:cellymax-1
                        for i in 2:cellxmax-1
                            Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + delta_Q[i,j,l]
                        end
                    end
                end
                
                # calculate primitive variables
                Qcon = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
                Qbasem = conservative_to_base(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)
            end
            # End of the inner iteration

            # Updating in physical time
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbase[i,j,l] = Qbasem[i,j,l]
                    end
                end
            end
            
            # output
            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                output_result(stepnum, Qbase, cellxmax, cellymax, specific_heat_ratio, Rhat, out_file_front, out_ext, out_dir)
            end

            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
        end
    end
    
    # end of time measurement
    end_t = now()

    # output of calculation time
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax)
end


# -- main --
main()
#throw(UndefVarError(:x))