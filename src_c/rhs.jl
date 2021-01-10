function setup_RHS(cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, wdot, nval, volume)
    RHS = zeros(cellxmax, cellymax, nval)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for l in 1:nval
                RHS[i,j,l] = - (E_adv_hat[i+1,  j, l] - E_vis_hat[i+1,  j, l]) +
                             (E_adv_hat[  i,  j, l] - E_vis_hat[  i,  j, l]) -
                             (F_adv_hat[  i,j+1, l] - F_vis_hat[  i,j+1, l]) +
                             (F_adv_hat[  i,  j, l] - F_vis_hat[  i,  j, l]) -
                             wdot[i,j,l] * volume[i,j]
            end
        end
    end
    return RHS
end