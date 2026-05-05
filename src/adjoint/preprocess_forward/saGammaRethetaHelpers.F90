










module saGammaRethetaHelpers
    !
    !  Helper functions for the SA-sLM2015 transition model
    !  (Piotrowski & Zingg, AIAA J., 2020, doi:10.2514/1.J059784).
    !
    !  All smooth approximations and correlation functions are collected
    !  here so they can be reused in residual, Jacobian, and debug-output
    !  code without duplication.
    !
    !  Source terms are NOT activated in this module — these are pure
    !  helper utilities.
    !
    use constants, only: realType, intType, zero, one, two, half, third, fourth, eps
    implicit none
    save

    ! ===================================================================
    !  Model constants (from paper)
    ! ===================================================================
    real(kind=realType), parameter :: slm_cthetat = 0.03_realType
    real(kind=realType), parameter :: slm_sigmathetat = 2.0_realType
    real(kind=realType), parameter :: slm_sigmaf = 1.0_realType
    real(kind=realType), parameter :: slm_ca1 = 2.0_realType
    real(kind=realType), parameter :: slm_ca2 = 0.06_realType
    real(kind=realType), parameter :: slm_ce1 = 1.0_realType
    real(kind=realType), parameter :: slm_ce2 = 50.0_realType
    real(kind=realType), parameter :: slm_ccrossflow = 0.6_realType

    ! Smooth min/max parameters
    real(kind=realType), parameter :: slm_p_smooth = 300.0_realType

contains

    ! ===================================================================
    !  φ_p(g1, g2, p): overflow-safe smooth max/min (Algorithm 1)
    !  p > 0 → smooth max;  p < 0 → smooth min
    ! ===================================================================
    function phi_p(g1, g2, p) result(phi)
        implicit none
        real(kind=realType), intent(in) :: g1, g2, p
        real(kind=realType) :: phi
        real(kind=realType) :: diff

        diff = p * (g2 - g1)

        if (abs(diff) < 1.0e-15_realType) then
            phi = (g1 + g2) * half + log(two) / abs(p) - log(two) / p
        else if (diff > 20.0_realType) then
            phi = g2 - log(two) / p
        else if (diff < -20.0_realType) then
            phi = g1 - log(two) / p
        else
            phi = (g1 + g2) * half + log(one + exp(diff)) / p - log(two) / p
        end if
    end function phi_p

    ! ===================================================================
    !  Smooth max:  smooth_max(g1, g2) ≈ max(g1, g2)
    !  Uses p = +300
    ! ===================================================================
    function smooth_max(g1, g2) result(val)
        implicit none
        real(kind=realType), intent(in) :: g1, g2
        real(kind=realType) :: val

        val = phi_p(g1, g2, slm_p_smooth)
    end function smooth_max

    ! ===================================================================
    !  Smooth min:  smooth_min(g1, g2) ≈ min(g1, g2)
    !  Uses p = -300
    ! ===================================================================
    function smooth_min(g1, g2) result(val)
        implicit none
        real(kind=realType), intent(in) :: g1, g2
        real(kind=realType) :: val

        val = phi_p(g1, g2, -slm_p_smooth)
    end function smooth_min

    ! ===================================================================
    !  Smooth Fonset (Eqs. 46-47)
    !  Fonset1 = sqrt( (Re_S/(2.6*Re_theta_c))^2 + R_T^2 )
    !  Fonset  = (tanh(6*(Fonset1 - 1.35)) + 1) / 2
    ! ===================================================================
    function smooth_Fonset(Re_S, Re_theta_c, R_T) result(val)
        implicit none
        real(kind=realType), intent(in) :: Re_S, Re_theta_c, R_T
        real(kind=realType) :: val
        real(kind=realType) :: Fonset1, ratio, Re_theta_c_safe

        Re_theta_c_safe = max(Re_theta_c, eps)
        ratio = Re_S / (2.6_realType * Re_theta_c_safe)
        Fonset1 = sqrt(ratio * ratio + R_T * R_T)
        val = (tanh(6.0_realType * (Fonset1 - 1.35_realType)) + one) / two
    end function smooth_Fonset

    ! ===================================================================
    !  Smooth Fturb (Eq. 48)
    !  Fturb = (1 - Fonset) * exp(-R_T)
    ! ===================================================================
    function smooth_Fturb(Fonset, R_T) result(val)
        implicit none
        real(kind=realType), intent(in) :: Fonset, R_T
        real(kind=realType) :: val

        val = (one - Fonset) * exp(-R_T)
    end function smooth_Fturb

    ! ===================================================================
    !  Smooth Flength (Eqs. 49-50)
    !  Flength1 = exp(-3e-2 * (ReThetaTilde - 460))
    !  Flength  = 44 - (44 - (0.50 - 3e-4*(ReThetaTilde-596))) /
    !             (1 + Flength1)^(1/6)
    ! ===================================================================
    function smooth_Flength(ReThetaTilde) result(val)
        implicit none
        real(kind=realType), intent(in) :: ReThetaTilde
        real(kind=realType) :: val
        real(kind=realType) :: Flength1, base

        Flength1 = exp(-3.0e-2_realType * (ReThetaTilde - 460.0_realType))
        base = one + Flength1
        ! (1 + Flength1)^(1/6):
        val = 44.0_realType - (44.0_realType - &
              (0.50_realType - 3.0e-4_realType * (ReThetaTilde - 596.0_realType))) / &
              (base**(one / 6.0_realType))
    end function smooth_Flength

    ! ===================================================================
    !  Smooth Re_theta_c (Eq. 51)
    !  Re_theta_c = 0.67*ReThetaTilde + 24*sin(ReThetaTilde/240 + 0.5) + 14
    ! ===================================================================
    function smooth_ReThetaC(ReThetaTilde) result(val)
        implicit none
        real(kind=realType), intent(in) :: ReThetaTilde
        real(kind=realType) :: val

        val = 0.67_realType * ReThetaTilde &
            + 24.0_realType * sin(ReThetaTilde / 240.0_realType + 0.5_realType) &
            + 14.0_realType
    end function smooth_ReThetaC

    ! ===================================================================
    !  Smooth F(lambda_theta) (Eqs. 54-57)
    !  F1 = 1 + 0.275*(1 - exp(-35*lam))*exp(-Tu/0.5)
    !  F2 = smooth_max(F1, 1)
    !  F3 = 1 - (-12.986*lam - 123.66*lam^2 - 405.689*lam^3)*exp(-(Tu/1.5)^1.5)
    !  F  = smooth_min(F2, F3)
    ! ===================================================================
    function smooth_FlambdaTheta(lambda_theta, Tu) result(val)
        implicit none
        real(kind=realType), intent(in) :: lambda_theta, Tu
        real(kind=realType) :: val
        real(kind=realType) :: F1val, F2val, F3val, lam

        lam = lambda_theta

        F1val = one + 0.275_realType * (one - exp(-35.0_realType * lam)) &
                * exp(-Tu / 0.5_realType)

        F2val = smooth_max(F1val, one)

        F3val = one - (-12.986_realType * lam &
                       - 123.66_realType * lam * lam &
                       - 405.689_realType * lam * lam * lam) &
                * exp(-(Tu / 1.5_realType)**1.5_realType)

        val = smooth_min(F2val, F3val)
    end function smooth_FlambdaTheta

    ! ===================================================================
    !  Re_theta correlation (Eqs. 8-9)
    !  Uses Tu_inf (freestream turbulence intensity, %)
    !  and lambda_theta (pressure gradient parameter)
    ! ===================================================================
    function ReTheta_correlation(Tu, lambda_theta) result(val)
        implicit none
        real(kind=realType), intent(in) :: Tu, lambda_theta
        real(kind=realType) :: val
        real(kind=realType) :: Flam, Tu_safe

        Tu_safe = max(Tu, 0.027_realType)  ! Protect 1/Tu^2

        Flam = smooth_FlambdaTheta(lambda_theta, Tu_safe)

        if (Tu_safe <= 1.3_realType) then
            val = (1173.51_realType &
                   - 589.428_realType * Tu_safe &
                   + 0.2196_realType / (Tu_safe * Tu_safe)) * Flam
        else
            val = 331.50_realType * (Tu_safe - 0.5658_realType)**(-0.671_realType) * Flam
        end if
    end function ReTheta_correlation

    ! ===================================================================
    !  Boundary layer proxies
    ! ===================================================================
    subroutine BL_proxies(ReThetaTilde, mu, rho, U, Re_ref, &
                          theta_BL, delta_BL)
        implicit none
        real(kind=realType), intent(in)  :: ReThetaTilde, mu, rho, U, Re_ref
        real(kind=realType), intent(out) :: theta_BL, delta_BL
        real(kind=realType) :: U_safe

        U_safe = max(U, eps)
        theta_BL = ReThetaTilde * mu / (rho * U_safe) / Re_ref
        delta_BL = 7.5_realType * theta_BL   ! 15/2
    end subroutine BL_proxies

    ! ===================================================================
    !  F_wake (Eq. 5)
    !  F_wake = exp(-Re_S / 1e6)
    ! ===================================================================
    function F_wake(Re_S) result(val)
        implicit none
        real(kind=realType), intent(in) :: Re_S
        real(kind=realType) :: val

        val = exp(-Re_S / 1.0e6_realType)
    end function F_wake

    ! ===================================================================
    !  F_theta_t (Eq. 3)
    !  F_theta_t = F_wake * exp(-(d/delta)^4)
    ! ===================================================================
    function F_theta_t(Fwake_val, d, delta) result(val)
        implicit none
        real(kind=realType), intent(in) :: Fwake_val, d, delta
        real(kind=realType) :: val
        real(kind=realType) :: delta_safe, ratio

        delta_safe = max(delta, eps)
        ratio = d / delta_safe
        val = Fwake_val * exp(-(ratio**4))
    end function F_theta_t

    ! ===================================================================
    !  Strain-rate Reynolds number  Re_S = rho*d^2*S/mu * Re
    ! ===================================================================
    function strain_Re_S(rho, d, S, mu, Re_ref) result(val)
        implicit none
        real(kind=realType), intent(in) :: rho, d, S, mu, Re_ref
        real(kind=realType) :: val
        real(kind=realType) :: mu_safe

        mu_safe = max(mu, eps)
        val = rho * d * d * S / mu_safe * Re_ref
    end function strain_Re_S

    ! ===================================================================
    !  Velocity magnitude with floor
    ! ===================================================================
    function velocity_mag(u, v, w_vel) result(val)
        implicit none
        real(kind=realType), intent(in) :: u, v, w_vel
        real(kind=realType) :: val

        val = sqrt(u * u + v * v + w_vel * w_vel)
        val = max(val, eps)
    end function velocity_mag

    ! ===================================================================
    !  Streamwise vorticity and crossflow helicity
    !  Omega_streamwise = |U_hat . omega_vec|
    !  H_crossflow = d * Omega_streamwise / U
    ! ===================================================================
    subroutine crossflow_helicity(u_vel, v_vel, w_vel, U_mag, &
                                  vortx, vorty, vortz, d, &
                                  Omega_streamwise, H_crossflow)
        implicit none
        real(kind=realType), intent(in)  :: u_vel, v_vel, w_vel, U_mag
        real(kind=realType), intent(in)  :: vortx, vorty, vortz, d
        real(kind=realType), intent(out) :: Omega_streamwise, H_crossflow
        real(kind=realType) :: U_safe, ux, uy, uz

        U_safe = max(U_mag, eps)
        ux = u_vel / U_safe
        uy = v_vel / U_safe
        uz = w_vel / U_safe

        Omega_streamwise = abs(ux * vortx + uy * vorty + uz * vortz)
        H_crossflow = d * Omega_streamwise / U_safe
    end subroutine crossflow_helicity

    ! ===================================================================
    !  Pressure gradient parameter lambda_theta
    !  lambda_theta = (rho * theta_BL^2 / mu) * dU/ds * Re
    !
    !  dU/ds = (u/U)*dU/dx + (v/U)*dU/dy + (w/U)*dU/dz
    !  dU/dx = (u*du/dx + v*dv/dx + w*dw/dx) / U
    ! ===================================================================
    function lambda_theta_calc(rho, theta_BL, mu, Re_ref, &
                               u_vel, v_vel, w_vel, U_mag, &
                               dudx, dudy, dudz, &
                               dvdx, dvdy, dvdz, &
                               dwdx, dwdy, dwdz) result(val)
        implicit none
        real(kind=realType), intent(in) :: rho, theta_BL, mu, Re_ref
        real(kind=realType), intent(in) :: u_vel, v_vel, w_vel, U_mag
        real(kind=realType), intent(in) :: dudx, dudy, dudz
        real(kind=realType), intent(in) :: dvdx, dvdy, dvdz
        real(kind=realType), intent(in) :: dwdx, dwdy, dwdz
        real(kind=realType) :: val
        real(kind=realType) :: U_safe, gradU_x, gradU_y, gradU_z, gradU_s, mu_safe

        U_safe = max(U_mag, eps)
        mu_safe = max(mu, eps)

        gradU_x = (u_vel * dudx + v_vel * dvdx + w_vel * dwdx) / U_safe
        gradU_y = (u_vel * dudy + v_vel * dvdy + w_vel * dwdy) / U_safe
        gradU_z = (u_vel * dudz + v_vel * dvdz + w_vel * dwdz) / U_safe

        gradU_s = (u_vel / U_safe) * gradU_x + (v_vel / U_safe) * gradU_y &
             + (w_vel / U_safe) * gradU_z

        val = (rho * theta_BL * theta_BL / mu_safe) * gradU_s * Re_ref

        ! Clamp to prevent correlation blowup
        val = max(-0.1_realType, min(0.1_realType, val))
    end function lambda_theta_calc

    ! ===================================================================
    !  Crossflow sink helper (Eqs. 15-22)
    !  Returns Re_scf
    ! ===================================================================
    function crossflow_Re_scf(H_crossflow, R_T, h, theta_t) result(val)
        implicit none
        real(kind=realType), intent(in) :: H_crossflow, R_T, h, theta_t
        real(kind=realType) :: val
        real(kind=realType) :: DeltaHcf, DeltaH_plus, DeltaH_minus
        real(kind=realType) :: f_plus, f_minus, h_safe, theta_t_safe, ratio

        DeltaHcf = H_crossflow * (one + min(R_T, 0.4_realType))

        DeltaH_plus  = max(0.1066_realType - DeltaHcf, zero)
        DeltaH_minus = max(-(0.1066_realType - DeltaHcf), zero)

        f_plus  = 6200.0_realType * DeltaH_plus &
                + 50000.0_realType * DeltaH_plus * DeltaH_plus

        f_minus = 75.0_realType * tanh(DeltaH_minus / 0.0125_realType)

        h_safe = max(h, eps)
        theta_t_safe = max(theta_t, eps)
        ratio = h_safe / theta_t_safe
        ratio = max(ratio, eps)

        val = -35.088_realType * log(ratio) + 319.51_realType &
            + f_plus - f_minus
    end function crossflow_Re_scf

    ! ===================================================================
    !  Time scale t (Eq. 7)
    !  t = 500*mu / (rho*U^2) * (1/Re)
    ! ===================================================================
    function time_scale_t(mu, rho, U, Re_ref) result(val)
        implicit none
        real(kind=realType), intent(in) :: mu, rho, U, Re_ref
        real(kind=realType) :: val
        real(kind=realType) :: U_safe, Re_safe

        U_safe = max(U, eps)
        Re_safe = max(Re_ref, eps)
        val = 500.0_realType * mu / (rho * U_safe * U_safe) / Re_safe
    end function time_scale_t

end module saGammaRethetaHelpers
