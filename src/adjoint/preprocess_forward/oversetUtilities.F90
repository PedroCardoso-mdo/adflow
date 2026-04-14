










module oversetUtilities
contains

    ! --------------------------------------------------
    !           Tapenade Routine BELOW this point
    ! --------------------------------------------------

    subroutine fracToWeights(frac, weights)
        use constants
        implicit none
        real(kind=realType), intent(in), dimension(3) :: frac
        real(kind=realType), intent(out), dimension(8) :: weights

        weights(1) = (one - frac(1)) * (one - frac(2)) * (one - frac(3))
        weights(2) = (frac(1)) * (one - frac(2)) * (one - frac(3))
        weights(3) = (one - frac(1)) * (frac(2)) * (one - frac(3))
        weights(4) = (frac(1)) * (frac(2)) * (one - frac(3))
        weights(5) = (one - frac(1)) * (one - frac(2)) * (frac(3))
        weights(6) = (frac(1)) * (one - frac(2)) * (frac(3))
        weights(7) = (one - frac(1)) * (frac(2)) * (frac(3))
        weights(8) = (frac(1)) * (frac(2)) * (frac(3))
    end subroutine fracToWeights

    subroutine fracToWeights2(frac, weights)
        use constants
        implicit none
        real(kind=realType), intent(in), dimension(3) :: frac
        real(kind=realType), intent(out), dimension(8) :: weights

        weights(1) = (one - frac(1)) * (one - frac(2)) * (one - frac(3))
        weights(2) = (frac(1)) * (one - frac(2)) * (one - frac(3))
        weights(3) = (frac(1)) * (frac(2)) * (one - frac(3))
        weights(4) = (one - frac(1)) * (frac(2)) * (one - frac(3))

        weights(5) = (one - frac(1)) * (one - frac(2)) * (frac(3))
        weights(6) = (frac(1)) * (one - frac(2)) * (frac(3))
        weights(7) = (frac(1)) * (frac(2)) * (frac(3))
        weights(8) = (one - frac(1)) * (frac(2)) * (frac(3))

    end subroutine fracToWeights2

    subroutine newtonUpdate(xCen, blk, frac0, frac)

        ! This routine performs the newton update to recompute the new
        ! "frac" (u,v,w) for the point xCen. The actual search is performed
        ! on the the dual cell formed by the cell centers of the 3x3x3 block
        ! of primal nodes. This routine is AD'd with tapenade in both
        ! forward and reverse.

        use constants
        implicit none

        ! Input
        real(kind=realType), dimension(3), intent(in) :: xCen
        real(kind=realType), dimension(3, 3, 3, 3), intent(in) :: blk
        real(kind=realType), dimension(3), intent(in) :: frac0
        ! Output
        real(kind=realType), dimension(3), intent(out) :: frac

        ! Working
        real(kind=realType), dimension(3, 1:8) :: xn
        real(kind=realType) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
        real(kind=realType) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, val
        real(kind=realType) :: f(3), x(3)
        integer(kind=intType), dimension(8), parameter :: indices = [1, 2, 4, 3, 5, 6, 8, 7]
        integer(kind=intType) :: i, j, k, ii, ll
        real(kind=realType), parameter :: adtEps = 1.e-25_realType
        real(kind=realType), parameter :: thresConv = 1.e-10_realType

        ! Compute the cell center locations for the 8 nodes describing the
        ! dual cell. Note that this must be counter-clockwise ordering.

        ii = 0
        do k = 1, 2
            do j = 1, 2
                do i = 1, 2
                    ii = ii + 1
                    xn(:, indices(ii)) = eighth * ( &
                                         blk(i, j, k, :) + &
                                         blk(i + 1, j, k, :) + &
                                         blk(i, j + 1, k, :) + &
                                         blk(i + 1, j + 1, k, :) + &
                                         blk(i, j, k + 1, :) + &
                                         blk(i + 1, j, k + 1, :) + &
                                         blk(i, j + 1, k + 1, :) + &
                                         blk(i + 1, j + 1, k + 1, :))
                end do
            end do
        end do

        ! Compute the coordinates relative to node 1.

        do i = 2, 8
            xn(:, i) = xn(:, i) - xn(:, 1)
        end do

        ! Compute the location of our seach point relative to the first node.
        x = xCen - xn(:, 1)

        ! Modify the coordinates of node 3, 6, 8 and 7 such that
        ! they correspond to the weights of the u*v, u*w, v*w and
        ! u*v*w term in the transformation respectively.

        xn(1, 7) = xn(1, 7) + xn(1, 2) + xn(1, 4) + xn(1, 5) &
                   - xn(1, 3) - xn(1, 6) - xn(1, 8)
        xn(2, 7) = xn(2, 7) + xn(2, 2) + xn(2, 4) + xn(2, 5) &
                   - xn(2, 3) - xn(2, 6) - xn(2, 8)
        xn(3, 7) = xn(3, 7) + xn(3, 2) + xn(3, 4) + xn(3, 5) &
                   - xn(3, 3) - xn(3, 6) - xn(3, 8)

        xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
        xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
        xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)

        xn(1, 6) = xn(1, 6) - xn(1, 2) - xn(1, 5)
        xn(2, 6) = xn(2, 6) - xn(2, 2) - xn(2, 5)
        xn(3, 6) = xn(3, 6) - xn(3, 2) - xn(3, 5)

        xn(1, 8) = xn(1, 8) - xn(1, 4) - xn(1, 5)
        xn(2, 8) = xn(2, 8) - xn(2, 4) - xn(2, 5)
        xn(3, 8) = xn(3, 8) - xn(3, 4) - xn(3, 5)

        ! Set the starting values of u, v and w based on our previous values

        u = frac0(1); v = frac0(2); w = frac0(3); 
        ! The Newton algorithm to determine the parametric
        ! weights u, v and w for the given coordinate.

        NewtonHexa: do ll = 1, 15

            ! Compute the RHS.

            uv = u * v; uw = u * w; vw = v * w; wvu = u * v * w

            f(1) = xn(1, 2) * u + xn(1, 4) * v + xn(1, 5) * w &
                   + xn(1, 3) * uv + xn(1, 6) * uw + xn(1, 8) * vw &
                   + xn(1, 7) * wvu - x(1)
            f(2) = xn(2, 2) * u + xn(2, 4) * v + xn(2, 5) * w &
                   + xn(2, 3) * uv + xn(2, 6) * uw + xn(2, 8) * vw &
                   + xn(2, 7) * wvu - x(2)
            f(3) = xn(3, 2) * u + xn(3, 4) * v + xn(3, 5) * w &
                   + xn(3, 3) * uv + xn(3, 6) * uw + xn(3, 8) * vw &
                   + xn(3, 7) * wvu - x(3)

            ! Compute the Jacobian.

            a11 = xn(1, 2) + xn(1, 3) * v + xn(1, 6) * w + xn(1, 7) * vw
            a12 = xn(1, 4) + xn(1, 3) * u + xn(1, 8) * w + xn(1, 7) * uw
            a13 = xn(1, 5) + xn(1, 6) * u + xn(1, 8) * v + xn(1, 7) * uv

            a21 = xn(2, 2) + xn(2, 3) * v + xn(2, 6) * w + xn(2, 7) * vw
            a22 = xn(2, 4) + xn(2, 3) * u + xn(2, 8) * w + xn(2, 7) * uw
            a23 = xn(2, 5) + xn(2, 6) * u + xn(2, 8) * v + xn(2, 7) * uv

            a31 = xn(3, 2) + xn(3, 3) * v + xn(3, 6) * w + xn(3, 7) * vw
            a32 = xn(3, 4) + xn(3, 3) * u + xn(3, 8) * w + xn(3, 7) * uw
            a33 = xn(3, 5) + xn(3, 6) * u + xn(3, 8) * v + xn(3, 7) * uv

            ! Compute the determinant. Make sure that it is not zero
            ! and invert the value. The cut off is needed to be able
            ! to handle exceptional cases for degenerate elements.

            val = a11 * (a22 * a33 - a32 * a23) + a21 * (a13 * a32 - a12 * a33) &
                  + a31 * (a12 * a23 - a13 * a22)
            val = sign(one, val) / max(abs(val), adtEps)

            ! Compute the new values of u, v and w.

            du = val * ((a22 * a33 - a23 * a32) * f(1) &
                        + (a13 * a32 - a12 * a33) * f(2) &
                        + (a12 * a23 - a13 * a22) * f(3))
            dv = val * ((a23 * a31 - a21 * a33) * f(1) &
                        + (a11 * a33 - a13 * a31) * f(2) &
                        + (a13 * a21 - a11 * a23) * f(3))
            dw = val * ((a21 * a32 - a22 * a31) * f(1) &
                        + (a12 * a31 - a11 * a32) * f(2) &
                        + (a11 * a22 - a12 * a21) * f(3))

            u = u - du; v = v - dv; w = w - dw

            ! Exit the loop if the update of the parametric
            ! weights is below the threshold

            val = sqrt(du * du + dv * dv + dw * dw)
            if (val <= thresConv) then
                exit NewtonHexa
            end if

        end do NewtonHexa

        ! We would *like* that all solutions fall inside the hexa, but we
        ! can't be picky here since we are not changing the donors. So
        ! whatever the u,v,w is we have to accept. Even if it is greater than
        ! 1 or less than zero, it shouldn't be by much.

        frac(1) = u
        frac(2) = v
        frac(3) = w

    end subroutine newtonUpdate

end module oversetUtilities
