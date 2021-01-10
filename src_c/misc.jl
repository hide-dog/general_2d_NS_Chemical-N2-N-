function set_initQbase(xmax, ymax, restart_file, init_rho, init_u, init_v, init_p, init_N2, init_N, 
                        specific_heat_ratio, out_file_front, out_ext, out_dir, restart_num, Rd, nval)
    Qbase=[]
    cellxmax = xmax - 1
    cellymax = ymax - 1

    restart_check = 0
    try Qbase = setup_restart_value(cellxmax, cellymax, out_dir, restart_file, nval)
        println("Restart "*restart_file)
        restart_check = 2
    catch 
        restart_check = 1
    end

    if restart_check == 1
        Qbase = setup_init_value(cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N, nval)
        println("Start Initial condition")
        restart_num = 0
        output_result(0, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
    end

    return Qbase, cellxmax, cellymax, restart_num
end

function setup_init_value(cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N, nval)
    Qbase=zeros(cellxmax, cellymax, nval)
    for i in 1:cellxmax
        for j in 1:cellymax
            Qbase[i,j,1] = init_rho
            Qbase[i,j,2] = init_u
            Qbase[i,j,3] = init_v
            Qbase[i,j,4] = init_p
            Qbase[i,j,5] = init_rho * init_N2
            Qbase[i,j,6] = init_rho * init_N
        end
    end
    return Qbase
end

function setup_restart_value(cellxmax, cellymax, out_dir, restart_file, nval)
    Qbase = zeros(cellxmax, cellymax, nval)

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
            for k in 1:nval
                Qbase[i,j,k] = parse(Float64,temp[k]) 
            end
            k = k+1
        end
    end
    return Qbase
end

function check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite)
    ite = 0
    for i in 1:cellxmax
        for j in 1:cellymax
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)
            if T < 0
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