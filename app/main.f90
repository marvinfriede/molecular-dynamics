program main
  !> prints a matrix quantity to screen
  !  examples:
  !  call write_vector(vec, name='vector')
  !  call write_matrix(mat, name='matrix')
  !  call write_matrix(mat, name='packed matrix')
  use print_matrix, only: write_matrix

  use build_box
  use algorithm
  implicit none

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  !> number of atoms
  integer, parameter :: natom = 108

  !> length of cubic box
  real(wp), parameter :: l = 17.158_wp

  !> atom coordinates of the system
  real(wp), dimension(3, natom) :: xyz

  !> forces
  real(wp), dimension(3, natom) :: f

  !> velocities
  real(wp), dimension(3, natom) :: v

  !> whatever
  real(wp), dimension(3, natom) :: a

  !> force constant
  real(wp), parameter :: k = 1.0_wp

  !> number of iterations
  integer, parameter :: itime = 200

  !> length of time step
  real(wp), parameter :: delta = 0.005_wp

  !> mass of particle
  real(wp), parameter :: mass = 39.948_wp

  !> temperature to rescale velocities
  integer, parameter :: t_rescale = 8

  !> potential energy
  real(wp) :: e_pot_harm
  !> kinetic energy
  real(wp) :: e_kin
  !> total energy
  real(wp) :: e_tot

  !> iteration variables
  integer :: i, j, w

  a = 0.0_wp
  v = 0.0_wp

  !> write energies
  open (unit=15, file="energy.csv")
  write (15, "(A27)") "time,e_kin,e_pot_harm,e_tot"

  !> pass 1 for sc and 4 for fcc
  call build_grid(xyz, natom, l, 4)

  !> calc for initial
  call calc_e_kin(natom, v, mass, e_kin)
  call calc_pot_harm(natom, xyz, a, l, k, e_pot_harm)
  e_tot = e_kin + e_pot_harm

  !> write trajectory for starting positions
  open (unit=16, file="trj.xyz")
  write (16, "(I3)") natom
  write (16, "(F15.10)") e_tot
  do w = 1, natom
    write (16, *) "Ar", xyz(1, w), xyz(2, w), xyz(3, w)
  end do

  !> evaluate forces and hence f/m from positions
  call calc_force(natom, xyz, f, a, l, k)

  time_prop: do i = 1, itime
    do j = 1, natom
      ! propagate velocities (half step)
      v(:, j) = v(:, j) + 0.5_wp*f(:, j)/mass*delta

      ! propagate positions
      xyz(:, j) = xyz(:, j) + v(:, j)*delta
    end do

    !> evaluate forces and hence f/m from positions
    call calc_force(natom, xyz, f, a, l, k)

    do j = 1, natom
      ! propagate velocities again (half step)
      v(:, j) = v(:, j) + 0.5_wp*f(:, j)/mass*delta
    end do

    call calc_e_kin(natom, v, mass, e_kin)
    call calc_pot_harm(natom, xyz, a, l, k, e_pot_harm)
    e_tot = e_kin + e_pot_harm

    ! write energies to console
    write (*, "(A16,2X,F15.8)") "Kinetic Energy =", e_kin
    write (*, "(A18,F15.8)") "Potential Energy =", e_pot_harm
    write (*, "(A14,4X,F15.8)") "Total Energy =", e_tot

    ! write energies to file
    write (15, "(F15.10,A1,F15.10,A1,F15.10,A1,F15.10)") delta*i, ",", e_kin, ",", e_pot_harm, ",", e_tot

    ! write trajectory
    write (16, "(I3)") natom
    write (16, "(F15.10)") e_tot
    do w = 1, natom
      write (16, *) "Ar", xyz(1, w), xyz(2, w), xyz(3, w)
    end do

  end do time_prop

end program main
