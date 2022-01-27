module algorithm
  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  public :: calc_e_kin, calc_force
contains
  subroutine calc_force(natom, xyz, f, l, e_pot_lj)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of particles
    real(wp), dimension(3, natom), intent(in) :: xyz
    !> distance vector between two particles
    real(wp), dimension(3) :: xyz_r
    !> squared distance between two particles
    real(wp) :: xyzij2

    !> forces
    real(wp), dimension(3, natom), intent(out) :: f
    !> force: derivate of LJ potential
    real(wp) :: ff

    !> length of cubic box
    real(wp), intent(in) :: l

    !> parameter in LJ potential
    real(wp), parameter :: sigma = 3.405_wp
    real(wp), parameter :: sigma2 = sigma*sigma
    real(wp) :: sigmar2
    real(wp) :: sigmar6
    real(wp) :: sigmar12

    !> parameter in LJ potential
    real(wp), parameter :: epsilon = 120.0_wp

    !> cutoff parameters for interaction radius
    real(wp) :: rcutoff
    real(wp) :: rcutoff2

    !> potential energy
    real(wp), intent(out) :: e_pot_lj

    !> iteration variable
    integer :: i, j

    rcutoff = 0.50_wp*l
    rcutoff2 = rcutoff*rcutoff

    f = 0.0_wp
    e_pot_lj = 0.0_wp
    do i = 1, natom - 1
      do j = i + 1, natom
        ! periodic boundary conditions
        xyz_r = xyz(:, i) - xyz(:, j)
        xyz_r = xyz_r - l*anint(xyz_r/l)

        xyzij2 = sum(xyz_r**2)
        if (xyzij2 .lt. rcutoff2) then
          sigmar2 = sigma2/xyzij2
          sigmar6 = sigmar2*sigmar2*sigmar2
          sigmar12 = sigmar6*sigmar6

          ! update force
          ff = (48.0_wp*epsilon*(sigmar12 - (0.5_wp*sigmar6)))/xyzij2
          f(:, i) = f(:, i) + ff*xyz_r
          f(:, j) = f(:, j) + ff*xyz_r

          ! energy evaluation
          e_pot_lj = e_pot_lj + sigmar12 - sigmar6 + 3.83738839608178386E-3_wp
        end if
      end do
    end do
    e_pot_lj = 4.0_wp*epsilon*e_pot_lj
  end subroutine calc_force

  subroutine calc_e_kin(natom, v, mass, e_kin)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) :: v

    !> mass of particle
    real(wp), intent(in) :: mass

    !> potential energy
    real(wp), intent(out) :: e_kin

    !> iteration variable
    integer :: i

    e_kin = 0.0_wp
    do i = 1, natom
      e_kin = e_kin + 0.5_wp*mass*sum(v(:, i)**2)
    end do

  end subroutine calc_e_kin

  !*********************************************************
  !****************** HARMONIC OSCILLATOR ******************
  !*********************************************************

  subroutine calc_force_harm(natom, xyz, f, a)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) :: xyz

    !> forces
    real(wp), dimension(3, natom), intent(out) :: f

    !> whatever
    real(wp), dimension(3, natom), intent(in) :: a

    !> force constant
    real(wp), parameter :: k = 1.0_wp

    !> iteration
    integer :: i

    f = 0.0_wp
    do i = 1, natom
      f(:, i) = -k*(xyz(:, i) - a(:, i))
    end do

  end subroutine calc_force_harm

  subroutine calc_pot_harm(natom, xyz, a, k, e_pot_harm)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) :: xyz

    !> whatever
    real(wp), dimension(3, natom), intent(in) :: a

    !> force constant
    real(wp), intent(in) :: k

    !> potential energy
    real(wp), intent(out) :: e_pot_harm

    !> iteration variable
    integer :: i

    e_pot_harm = 0.0_wp
    do i = 1, natom
      e_pot_harm = e_pot_harm + 0.5_wp*k*sum((xyz(:, i) - a(:, i))**2)
    end do

  end subroutine calc_pot_harm

end module algorithm
