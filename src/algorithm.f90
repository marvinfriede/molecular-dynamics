module algorithm
  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  public :: calc_force, calc_pot_harm, calc_e_kin
contains
  subroutine calc_force(natom, xyz, f, a, l, k)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) :: xyz

    !> forces
    real(wp), dimension(3, natom), intent(out) :: f

    !> whatever
    real(wp), dimension(3, natom), intent(in) :: a

    !> length of cubic box
    real(wp), intent(in) :: l

    !> force constant
    real(wp), intent(in) :: k

    !> iteration
    integer :: i

    f = 0.0_wp
    do i = 1, natom
      f(:, i) = -k*(xyz(:, i) - a(:, i))
    end do

  end subroutine calc_force

  subroutine calc_pot_harm(natom, xyz, a, l, k, e_pot_harm)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(in) :: xyz

    !> whatever
    real(wp), dimension(3, natom), intent(in) :: a

    !> length of cubic box
    real(wp), intent(in) :: l

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
end module algorithm
