function [Dxustar,Dyvstar] = DxDyRC(A_U_U,A_V_V,Dx,Dy,u_vec,v_vec,p_vec,space,FV,ib,Area,Length)

        Dxustar = zeros(space.dimY*space.dimX,1);
        Dyvstar = zeros(space.dimY*space.dimX,1);

        Dxp = Dx*p_vec;
        Dyp = Dy*p_vec;

        % Rhie Chow for south node
        u_bar_s = 0.5*(u_vec(ib.P)+u_vec(ib.s));
        A_bar_U_P_s = 0.5*(A_U_U(ib.P,ib.P) + A_U_U(ib.s,ib.s));
        dpx_bar_s = 0.5*(Dxp(ib.P)+Dxp(ib.s));
        uRC_s = u_bar_s + Area.S_s(ib.s).*(1./diag(A_bar_U_P_s)).*(dpx_bar_s-Dxp(ib.s));
        
        v_bar_s = 0.5*(v_vec(ib.P)+v_vec(ib.s));
        A_bar_V_P_s = 0.5*(A_V_V(ib.P,ib.P) + A_V_V(ib.s,ib.s));
        dpy_bar_s = 0.5*(Dyp(ib.P)+Dyp(ib.s));
        vRC_s = v_bar_s + Area.S_s(ib.s).*(1./diag(A_bar_V_P_s)).*(dpy_bar_s-Dyp(ib.s));
        
        % Rhie Chow for east node
        u_bar_e = 0.5*(u_vec(ib.P)+u_vec(ib.e));
        A_bar_U_P_e = 0.5*(A_U_U(ib.P,ib.P) + A_U_U(ib.e,ib.e));
        dpx_bar_e = 0.5*(Dxp(ib.P)+Dxp(ib.e));
        uRC_e = u_bar_e + Area.S_e(ib.e).*(1./diag(A_bar_U_P_e)).*(dpx_bar_e-Dxp(ib.e));
        
        v_bar_e = 0.5*(v_vec(ib.P)+v_vec(ib.e));
        A_bar_V_P_e = 0.5*(A_V_V(ib.P,ib.P) + A_V_V(ib.e,ib.e));
        dpy_bar_e = 0.5*(Dyp(ib.P)+Dyp(ib.e));
        vRC_e = v_bar_e + Area.S_e(ib.e).*(1./diag(A_bar_V_P_e)).*(dpy_bar_e-Dyp(ib.e));
        
        % Rhie Chow for west node
        u_bar_w = 0.5*(u_vec(ib.P)+u_vec(ib.w));
        A_bar_U_P_w = 0.5*(A_U_U(ib.P,ib.P) + A_U_U(ib.w,ib.w));
        dpx_bar_w = 0.5*(Dxp(ib.P)+Dxp(ib.w));
        uRC_w = u_bar_w + Area.S_w(ib.w).*(1./diag(A_bar_U_P_w)).*(dpx_bar_w-Dxp(ib.w));
        
        v_bar_w = 0.5*(v_vec(ib.P)+v_vec(ib.w));
        A_bar_V_P_w = 0.5*(A_V_V(ib.P,ib.P) + A_V_V(ib.w,ib.w));
        dpy_bar_w = 0.5*(Dyp(ib.P)+Dyp(ib.w));
        vRC_w = v_bar_w + Area.S_w(ib.w).*(1./diag(A_bar_V_P_w)).*(dpy_bar_w-Dyp(ib.w));
        
        % Rhie Chow for north node
        u_bar_n = 0.5*(u_vec(ib.P)+u_vec(ib.n));
        A_bar_U_P_n = 0.5*(A_U_U(ib.P,ib.P) + A_U_U(ib.n,ib.n));
        dpx_bar_n = 0.5*(Dxp(ib.P)+Dxp(ib.n));
        uRC_n = u_bar_n + Area.S_n(ib.n).*(1./diag(A_bar_U_P_n)).*(dpx_bar_n-Dxp(ib.n));
        
        v_bar_n = 0.5*(v_vec(ib.P)+v_vec(ib.n));
        A_bar_V_P_n = 0.5*(A_V_V(ib.P,ib.P) + A_V_V(ib.n,ib.n));
        dpy_bar_n = 0.5*(Dyp(ib.P)+Dyp(ib.n));
        vRC_n = v_bar_n + Area.S_n(ib.n).*(1./diag(A_bar_V_P_n)).*(dpy_bar_n-Dyp(ib.n));
        
        Dxustar(ib.P) = (Length.dy_sw_se(ib.P).*uRC_s + Length.dy_se_ne(ib.P).*uRC_e + Length.dy_ne_nw(ib.P).*uRC_n...
            + Length.dy_nw_sw(ib.P).*uRC_w)./Area.S_P(ib.P);
        Dyvstar(ib.P) = -1*(Length.dx_sw_se(ib.P).*vRC_s + Length.dx_se_ne(ib.P).*vRC_e + Length.dx_ne_nw(ib.P).*vRC_n...
            + Length.dx_nw_sw(ib.P).*vRC_w)./Area.S_P(ib.P);
end