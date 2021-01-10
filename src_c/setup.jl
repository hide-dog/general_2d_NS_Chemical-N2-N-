
function base_to_conservative(Qbase, cellxmax, cellymax, specific_heat_ratio, nval, nch)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    npre = nval - nch
    Qcon = zeros(cellxmax, cellymax, nval)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            Qcon[i,j,1] = Qbase[i,j,1]
            Qcon[i,j,2] = Qbase[i,j,1]*Qbase[i,j,2]
            Qcon[i,j,3] = Qbase[i,j,1]*Qbase[i,j,3]
            Qcon[i,j,4] = Qbase[i,j,4]/(specific_heat_ratio-1)+Qbase[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2
            for l in 1:nch
                Qcon[i,j,npre+l] = Qbase[i,j,1]*Qbase[i,j,npre+l]
            end
        end
    end
    return Qcon
end

function conservative_to_base(Qcon, cellxmax, cellymax, specific_heat_ratio, nval, nch)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    npre = nval - nch
    Qbase = zeros(cellxmax, cellymax, nval)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            Qbase[i,j,1]=Qcon[i,j,1]
            Qbase[i,j,2]=Qcon[i,j,2]/Qcon[i,j,1]
            Qbase[i,j,3]=Qcon[i,j,3]/Qcon[i,j,1]
            Qbase[i,j,4]=(Qcon[i,j,4]-Qcon[i,j,1]*(Qbase[i,j,2]^2+Qbase[i,j,3]^2)/2)*(specific_heat_ratio-1)
            for l in 1:nch
                Qbase[i,j,npre+l] = Qcon[i,j,npre+l]/Qcon[i,j,1]
            end
        end
    end
    return Qbase
end

function setup_Qcon_hat(Qcon, cellxmax, cellymax, volume, nval, nch)
    Qcon_hat = zeros(cellxmax, cellymax, nval)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for l in 1:nval
                Qcon_hat[i,j,l] = Qcon[i,j,l]*volume[i,j]
            end
        end
    end
    return Qcon_hat
end

function Qhat_to_Q(Qcon_hat, cellxmax, cellymax, volume, nval, nch)
    Qcon = zeros(cellxmax, cellymax, nval)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for l in 1:nval
                Qcon[i,j,l] = Qcon_hat[i,j,l]/volume[i,j]
            end
        end
    end
    return Qcon
end