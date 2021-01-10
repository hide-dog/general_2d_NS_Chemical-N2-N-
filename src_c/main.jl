    using ProgressMeter

function main()
    out_dir  = "result"
    PARAMDAT = "PARAMDAT.json"
    fwrite   = "write"
    
    nch  = 2    # rho,u,v,p
    nval = 6    # rho,u,v,p + N2, N
    nre  = 2    # N2+N2 -> 2N+N2
                # N2+N  -> 2N+N
    Rd   = 287.0
    R    = 8.314

    xmax, ymax, nodes, vecAx, vecAy = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau,
    init_rho, init_u, init_v, init_p, init_N2, init_N, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)

    Qbase, cellxmax, cellymax, restartnum = set_initQbase(xmax, ymax, restart_file, init_rho, init_u, init_v, init_p, init_N2, init_N, 
                                                        specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval)

    # init Delta_Qcon_hat
    volume = set_volume(nodes, cellxmax, cellymax)
    Minf   = set_Minf(bdcon, specific_heat_ratio)
    reset_write(fwrite)
    
    # main loop
    print("threads num : ")
    println(Threads.nthreads())

    # check bd
    check_bd(bdcon)
    
    #throw(UndefVarError(:x))
    # main loop
    prog = Progress(nt,1)
    for t in 1:nt
        next!(prog)
                
        evalnum = t+restartnum
        if time_integ == "1"
            # exlicit scheme

            #println(Qbase[1,3,:])
            
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, nval)
            Qcon     = base_to_conservative(Qbase, cellxmax, cellymax, specific_heat_ratio, nval, nch)
            Qcon_hat = setup_Qcon_hat(Qcon, cellxmax, cellymax, volume, nval, nch)

            # initial_setup
            chmu, chlambda_tr, chD, chi = chemical_value(Qbase, cellxmax, cellymax, Rd, nval, nch)

            # mu     = set_mu(Qbase,cellxmax,cellymax,specific_heat_ratio,Rd)
            # lambda = set_lambda(Qbase,cellxmax,cellymax,mu,specific_heat_ratio,Rd)
            
            # RHS
            # sorce term
            wdot = set_wdot(Qbase, cellxmax, cellymax, Rd, nch, nre, nval)

            #println(Qbase[1,3,:])
            #println(Qcon[1,3,:])
            # advection_term
            #E_adv_hat, F_adv_hat = AUSM(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio)
            E_adv_hat, F_adv_hat = AUSM_plusup(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, Minf, volume, nval, nch)
            
            #println("AUSM")
            
            #println(Qcon[40,101,:])
            #println(Qcon[40,102,:])
            #println(E_adv_hat[2,10,:])
            #println(E_adv_hat[35,3,:])
            #println(F_adv_hat[2,100,:])
            #println(wdot[2,2,:])
            #println(F_adv_hat[35,3,:])
            
            #=
            cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus = cal_jacobi(Qbase,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy,volume)
            cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus = setup_cell_flux_hat(Qcon,cellxmax,cellymax,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus)
            E_adv_hat,F_adv_hat = FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,cellxmax,cellymax)
            =#

            #println("FVS")
            #println(E_adv_hat[52,50,:])
            #println(E_adv_hat[35,3,:])
            #println(F_adv_hat[52,50,:])
            #println(F_adv_hat[35,3,:])
            
            
            # viscos_term
            E_vis_hat, F_vis_hat = central_diff(Qbase, Qcon, cellxmax, cellymax, chmu, chlambda_tr, chD, chi,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval, nch)
            #E_vis_hat = zeros(cellxmax+1, cellymax+1, nval)
            #F_vis_hat = zeros(cellxmax+1, cellymax+1, nval)

            RHS = setup_RHS(cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, wdot, nval, volume)
            #=
            println(E_adv_hat[2,50,:])
            println(F_adv_hat[2,50,:])
            println(E_vis_hat[2,50,:])
            println(F_vis_hat[2,50,:])
            println(wdot[2,50,:])
            println(RHS[2,50,:])
            =#
            

            #=
            println("          ")
            println(E_adv_hat[100,2,:])
            println(F_adv_hat[100,2,:])
            println(E_adv_hat[100,3,:])
            println(F_adv_hat[100,3,:])
            println(Qcon_hat[100,2,:])
            println(RHS[100,2,:])
            =#

            #println(RHS[2,4,:])
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)

            Qcon  = Qhat_to_Q(Qcon_hat, cellxmax, cellymax, volume, nval, nch)
            Qbase = conservative_to_base(Qcon, cellxmax, cellymax, specific_heat_ratio, nval, nch)

        elseif time_integ == "2"
            # Qcn
            Qbasen = copy(Qbase)
            Qconn = base_to_conservative(Qbasen,cellxmax,cellymax,specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn,cellxmax,cellymax,volume)

            Qbasem = copy(Qbase)

            for tau in 1:in_nt
                # LHS (A_adv_hat=jacobian)
                Qbasem = set_boundary(Qbasem,cellxmax,cellymax,vecAx,vecAy,bdcon)
                
                Qcon = base_to_conservative(Qbasem,cellxmax,cellymax,specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon,cellxmax,cellymax,volume)

                # initial_setup
                mu = set_mu(Qbasem,cellxmax,cellymax,specific_heat_ratio,Rd)
                lambda = set_lambda(Qbasem,cellxmax,cellymax,mu,specific_heat_ratio,Rd)
                                
                # RHS
                #advection_term
                #E_adv_hat, F_adv_hat = AUSM(Qbasem,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio)
                cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus = cal_jacobi(Qbasem,Qcon,cellxmax,cellymax,specific_heat_ratio,vecAx,vecAy,volume)
                cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus = setup_cell_flux_hat(Qcon,cellxmax,cellymax,cell_Ahat_plas,cell_Ahat_minus,cell_Bhat_plas,cell_Bhat_minus)
                E_adv_hat,F_adv_hat = FVS(cell_E_hat_plas,cell_E_hat_minus,cell_F_hat_plas,cell_F_hat_minus,cellxmax,cellymax)

                # viscos_term
                E_vis_hat, F_vis_hat = central_diff(Qbasem,Qcon,cellxmax,cellymax,mu,lambda,vecAx,vecAy,specific_heat_ratio,volume,Rd)
                
                RHS = setup_RHS(cellxmax,cellymax,E_adv_hat, F_adv_hat,E_vis_hat, F_vis_hat)

                #println(E_adv_hat[:,50,:])
                #println(E_adv_hat[:,50,:])

                #println(RHS[:,50,:])

                #throw(UndefVarError(:x))


                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig = one_wave(Qbasem,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio,volume)
                # lusgs_viscos_term
                jalphaP, jbetaP = central_diff_jacobian(Qbasem,Qcon,cellxmax,cellymax,mu,lambda,vecAx,vecAy,specific_heat_ratio,volume)
                #jalphaP = zeros(cellxmax, cellymax)
                #jbetaP = zeros(cellxmax, cellymax)

                # LUSGS
                delta_Q = zeros(cellxmax,cellymax,4)

                ite=0
                while true
                    delta_Q_temp = copy(delta_Q)
                    delta_Q = lusgs(dt,Qcon_hat,Qconn_hat,delta_Q,A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig,jalphaP, jbetaP,RHS,cellxmax,cellymax,volume)
                    
                    # println(RHS[:,50,:])
                    #println(delta_Q[:,50,:])

                    res = set_res(delta_Q,delta_Q_temp,cellxmax,cellymax)
                    norm2 = check_converge(res,RHS,cellxmax,cellymax,init_small)

                    ite += 1
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok
                        break
                    end
                    if ite % 100 ==0
                        println(norm2)
                    end
                end
                #println(ite)
                #println(delta_Q[:,50,:])
                
                for i in 2:cellxmax-1
                    for j in 2:cellymax-1
                        for k in 1:4
                            Qcon_hat[i,j,k] = Qcon_hat[i,j,k] + delta_Q[i,j,k]
                        end
                    end
                end

                Qcon = Qhat_to_Q(Qcon_hat,cellxmax,cellymax,volume)
                Qbasem = conservative_to_base(Qcon,cellxmax,cellymax,specific_heat_ratio)

                #println("xxx")
            end            
            Qbase = copy(Qbasem)
        end
        
        if round(evalnum) % every_outnum == 0
            println("\n")
            println("nt_______________________________"*string(round(evalnum)))
            output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
        end

        check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)

    end
end

# -- main --
main()