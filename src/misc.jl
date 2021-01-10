# ------------------------------------
# set initial conditions
# ------------------------------------
function set_initQbase(Qbase, cellxmax, cellymax, restart_file, init_rho, init_u, init_v, init_p, init_T,
                        init_N2, init_N, specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rs)

    restart_check = 0
    try Qbase = setup_restart_value(Qbase, cellxmax, cellymax, out_dir, restart_file)
        println("Restart "*restart_file)
        restart_check = 2
    catch 
        restart_check = 1
    end

    if restart_check == 1
        Qbase = setup_init_value(Qbase, cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N)
        println("Start Initial condition")
        restartnum = 0
        Rhat = set_gasconst(Qbase, cellxmax, cellymax, Rs)
        output_result(0, Qbase, cellxmax, cellymax, specific_heat_ratio, Rhat, out_file_front, out_ext, out_dir)
    end

    return Qbase, restartnum
end

function setup_init_value(Qbase, cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N)
    for j in 1:cellymax
        for i in 1:cellxmax
            Qbase[i,j,1] = init_rho
            Qbase[i,j,2] = init_u
            Qbase[i,j,3] = init_v
            Qbase[i,j,4] = init_p
            Qbase[i,j,5] = init_N2
            Qbase[i,j,6] = init_N
        end
    end
    return Qbase
end

function setup_restart_value(Qbase, cellxmax, cellymax, out_dir, restart_file)
    skipnum = 1
    fff = []
    open("result/"*restart_file, "r") do f
        fff = read(f, String)
    end 
    fff = split(fff,"\n", keepempty = false)
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end
    
    k = 1
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            temp = split(fff[k+skipnum]," ")
            for l in 1:nval
                Qbase[i,j,l] = parse(Float64,temp[l]) 
            end
            k = k+1
        end
    end
    return Qbase
end

# ------------------------------------
# find out if results were diverge
# ------------------------------------
function check_divrege(Qbase, cellxmax, cellymax, Rhat, fwrite)
    ite = 0
    for j in 2:cellymax-1
        for i in 2:cellxmax-1        
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rhat[i,j])
            if T < 0 || isequal(T, NaN) == true
                open( fwrite, "a" ) do f
                    ai = @sprintf("%4.0f", i)
                    aj = @sprintf("%4.0f", j)

                    write(f, "\n")
                    write(f, " diverge ")
                    write(f, "\n")
                    write(f, " i = "*ai)
                    write(f, "\n")
                    write(f, " j = "*aj)
                    write(f, "\n")
                end
                ite = 1
            end
        end
    end

    if ite == 1
        println("\n")
        println("\n")
        println(" T<0 ")
        println(" diverge ")
        println("\n")
        throw(UndefVarError(:x))
    end
end