using Printf

function output_result(stepnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval)
    
    stepnum = string(stepnum)
    while length(stepnum) < 6
        stepnum = "0"*stepnum
    end
    
    fff = out_dir*"/"*out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], p[Pa], massfrac_N2[kg/m^3], massfrac_N[kg/m^3], T[K]\n")
        for i in 2:cellxmax-1
            for j in 2:cellymax-1
                for l in 1:nval
                    a = @sprintf("%8.8e", Qbase[i,j,l])
                    write(f, a*" ")
                end
                T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)
                a = @sprintf("%8.8e", T)
                write(f, a*"\n")
            end
        end
    end
    println("\nwrite "*fff)
end

function reset_write(fwrite)
    open( fwrite, "w" ) do f 
    end
end