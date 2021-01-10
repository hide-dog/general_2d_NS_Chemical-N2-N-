# ------------------------------------
# calcucation of the viscosity term by central difference
# ------------------------------------
function central_diff(E_vis_hat, F_vis_hat, Qbase, Qcon, cellxmax, cellymax, mu, lambda,
                        vecAx, vecAy, specific_heat_ratio, volume, Rd, nval)
    
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            rho_av = 0.5*(Qbase[i,j,1] + Qbase[i-1,j,1]) 
            u_av   = 0.5*(Qbase[i,j,2] + Qbase[i-1,j,2])
            v_av   = 0.5*(Qbase[i,j,3] + Qbase[i-1,j,3])
            T_av   = Qbase[i,j,4]/(Qbase[i,j,1]*Rd) - Qbase[i-1,j,4]/(Qbase[i-1,j,1]*Rd)
            
            mu_av     = 0.5*(mu[i-1,j] + mu[i,j]) 
            lambda_av = (0.5*(1/lambda[i-1,j] + 1/lambda[i,j]))^(-1)

            volume_av = 0.5*(volume[i,j] + volume[i-1,j])
    
            dudxi = Qbase[i,j,2] - Qbase[i-1,j,2]
            dvdxi = Qbase[i,j,3] - Qbase[i-1,j,3]
            dTdxi = Qbase[i,j,4]/(Qbase[i,j,1]*Rd) - Qbase[i-1,j,4]/(Qbase[i-1,j,1]*Rd)
    
            dudeta = 0.25 * (Qbase[i,j+1,2] - Qbase[i,j-1,2] + Qbase[i-1,j+1,2] - Qbase[i-1,j-1,2])
            dvdeta = 0.25 * (Qbase[i,j+1,3] - Qbase[i,j-1,3] + Qbase[i-1,j+1,3] - Qbase[i-1,j-1,3])
            dTdeta = 0.25* (Qbase[i,j+1,4]/(Qbase[i,j+1,1]*Rd) - Qbase[i,j-1,4]/(Qbase[i,j-1,1]*Rd) +
                            Qbase[i-1,j+1,4]/(Qbase[i-1,j+1,1]*Rd) -  Qbase[i-1,j-1,4]/(Qbase[i-1,j-1,1]*Rd))           

            vecAy_xav    = 0.25*( vecAy[i,j,1] + vecAy[i,j+1,1] + vecAy[i-1,j,1] + vecAy[i-1,j+1,1] )
            vecAy_yav    = 0.25*( vecAy[i,j,2] + vecAy[i,j+1,2] + vecAy[i-1,j,2] + vecAy[i-1,j+1,2] )
            
            tau_xx = 2/3*mu_av*(2*( vecAx[i,j,1] * dudxi + vecAy_xav * dudeta )
                    - (vecAx[i,j,2]*dvdxi + vecAy_yav*dvdeta ))          
                    
            tau_yy = 2/3*mu_av*(( vecAx[i,j,1]*dudxi -vecAy_xav*dudeta )
                    + 2*(vecAx[i,j,2]*dvdxi + vecAy_yav*dvdeta ))          
                    
            tau_xy = mu_av*(                
                    + (vecAx[i,j,2]*dudxi + vecAy_yav*dudeta ) 
                    + ( vecAx[i,j,1]*dvdxi + vecAy_xav*dvdeta ) )
            
            dTdx = vecAx[i,j,1]*dTdxi + vecAy_xav*dTdeta
            dTdy = vecAx[i,j,2]*dTdxi + vecAy_yav*dTdeta
       
            betax = u_av*tau_xx + v_av*tau_xy + lambda_av * dTdx
            betay = u_av*tau_xy + v_av*tau_yy + lambda_av * dTdy
            
            E_vis_hat[i,j,1] = 0.0
            E_vis_hat[i,j,2] = (vecAx[i,j,1]*tau_xx + vecAx[i,j,2]*tau_xy) / volume_av
            E_vis_hat[i,j,3] = (vecAx[i,j,1]*tau_xy + vecAx[i,j,2]*tau_yy) / volume_av
            E_vis_hat[i,j,4] = (vecAx[i,j,1]*betax + vecAx[i,j,2]*betay)   / volume_av
        end
    end
    
    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1
            rho_av = 0.5*(Qbase[i,j,1] + Qbase[i,j-1,1]) 
            u_av = 0.5*(Qbase[i,j,2] + Qbase[i,j-1,2])
            v_av = 0.5*(Qbase[i,j,3] + Qbase[i,j-1,3])
            T_av   = Qbase[i,j,4]/(Qbase[i,j,1]*Rd) - Qbase[i,j-1,4]/(Qbase[i,j-1,1]*Rd)

            mu_av     = 0.5*(mu[i,j-1] + mu[i,j])
            lambda_av = (0.5*(1/lambda[i,j-1] + 1/lambda[i,j]))^(-1)

            volume_av = 0.5*(volume[i,j] + volume[i,j-1])
    
            dudxi = 0.25 * (Qbase[i+1,j,2] - Qbase[i-1,j,2] + Qbase[i+1,j-1,2] - Qbase[i-1,j-1,2])
            dvdxi = 0.25 * (Qbase[i+1,j,3] - Qbase[i-1,j,3] + Qbase[i+1,j-1,3] - Qbase[i-1,j-1,3])
            dTdxi = 0.25 * (Qbase[i+1,j,4]/(Qbase[i+1,j,1]*Rd) - Qbase[i-1,j,4]/(Qbase[i-1,j,1]*Rd) +
                            Qbase[i+1,j-1,4]/(Qbase[i+1,j-1,1]*Rd) - Qbase[i-1,j-1,4]/(Qbase[i-1,j-1,1]*Rd))

            dudeta = Qbase[i,j,2] - Qbase[i,j-1,2]
            dvdeta = Qbase[i,j,3] - Qbase[i,j-1,3]
            dTdeta = (Qbase[i,j,4]/(Qbase[i,j,1]*Rd) - Qbase[i,j-1,4]/(Qbase[i,j-1,1]*Rd))
                    
        
            vecAx_xav    = 0.25*( vecAx[i,j,1] + vecAx[i+1,j,1] + vecAx[i,j-1,1] + vecAx[i+1,j-1,1] )
            vecAx_yav    = 0.25*( vecAx[i,j,2] + vecAx[i+1,j,2] + vecAx[i,j-1,2] + vecAx[i+1,j-1,2] )
        
            tau_xx = 2/3*mu_av*(2*( vecAx_xav * dudxi + vecAy[i,j,1] * dudeta )
                    - (vecAx_yav*dvdxi + vecAy[i,j,2]*dvdeta ))          
                    
            tau_yy = 2/3*mu_av*(( vecAx_xav*dudxi -vecAy[i,j,1]*dudeta )
                    + 2*(vecAx_yav*dvdxi + vecAy[i,j,2]*dvdeta ))          
                    
            tau_xy = mu_av*(                
                    + (vecAx_yav*dudxi + vecAy[i,j,2]*dudeta ) 
                    + ( vecAx_xav*dvdxi + vecAy[i,j,1]*dvdeta ) )
            
            dTdx = vecAx_xav*dTdxi + vecAy[i,j,1]*dTdeta
            dTdy = vecAx_yav*dTdxi + vecAy[i,j,2]*dTdeta
                   
            betax = u_av*tau_xx + v_av*tau_xy + lambda_av * dTdx
            betay = u_av*tau_xy + v_av*tau_yy + lambda_av * dTdy

            F_vis_hat[i,j,1] = 0.0
            F_vis_hat[i,j,2] = (vecAy[i,j,1]*tau_xx + vecAy[i,j,2]*tau_xy) / volume_av
            F_vis_hat[i,j,3] = (vecAy[i,j,1]*tau_xy + vecAy[i,j,2]*tau_yy) / volume_av
            F_vis_hat[i,j,4] = (vecAy[i,j,1]*betax + vecAy[i,j,2]*betay)   / volume_av
        end
    end
    return E_vis_hat, F_vis_hat
end
