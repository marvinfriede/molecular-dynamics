program main
  !> prints a matrix quantity to screen for debug
  !  examples:
  !  call write_vector(vec, name='vector')
  !  call write_matrix(mat, name='matrix')
  !  call write_matrix(mat, name='packed matrix')
  !use print_matrix, only: write_matrix

  use build_box
  use algorithm
  implicit none

  intrinsic :: sqrt

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

  !> number of iterations
  integer, parameter :: itime = 200

  !> length of time step
  real(wp), parameter :: delta = 0.005_wp

  !> mass of particle
  real(wp), parameter :: mass = 39.948_wp

  !> temperature
  real(wp) :: T
  !> temperature to rescale velocities
  integer, parameter :: T_req = 8
  !> factor to rescale velocities
  real(wp) :: T_scaler

  !> potential energy
  real(wp) :: e_pot
  !> kinetic energy
  real(wp) :: e_kin
  !> total energy
  real(wp) :: e_tot

  !> iteration variables
  integer :: i, j, w

  ! initial temperature scaling (1 means no scaling)
  T_scaler = 1.0_wp

  ! initial velocities
  v = 0.0_wp

  ! pass 1 for sc and 4 for fcc
  call build_grid(xyz, natom, l, 4)

  ! evaluate forces and hence f/m from positions
  call calc_force(natom, xyz, f, l, e_pot)

  ! write energies
  open (unit=15, file="out/energy.csv")
  write (15, "(A22)") "time,e_kin,e_pot,e_tot"

  ! write trajectory for starting positions
  open (unit=16, file="out/trj.xyz")
  write (16, "(I3)") natom
  write (16, "(F15.8)") e_pot ! e_kin = 0 in first step (v = 0)
  do w = 1, natom
    write (16, *) "Ar", xyz(1, w), xyz(2, w), xyz(3, w)
  end do

  !*************************************************************
  !************************* MAIN LOOP *************************
  !*************************************************************

  main_loop: do i = 1, itime
    ! scale factor for velocities (skip first iteration because v = 0)
    if (i .gt. 1) then
      T = e_kin/3.0_wp/natom
      T_scaler = sqrt(T_req/T)
    end if

    do j = 1, natom
      ! scale and propagate velocities (half step)
      v(:, j) = T_scaler*(v(:, j) + 0.5_wp*f(:, j)/mass*delta)

      ! propagate positions
      xyz(:, j) = xyz(:, j) + v(:, j)*delta
    end do

    ! evaluate forces and hence f/m from positions
    call calc_force(natom, xyz, f, l, e_pot)

    do j = 1, natom
      ! propagate velocities again (half step) (no scaling?)
      v(:, j) = v(:, j) + 0.5_wp*f(:, j)/mass*delta
    end do

    call calc_e_kin(natom, v, mass, e_kin)
    e_tot = e_kin + e_pot

    ! write energies to console
    write (*, "(A16,2X,F15.8)") "Kinetic Energy =", e_kin
    write (*, "(A18,F15.8)") "Potential Energy =", e_pot
    write (*, "(A14,4X,F15.8)") "Total Energy =", e_tot

    ! write energies to file
    write (15, "(F15.8,A1,F15.8,A1,F15.8,A1,F15.8)") delta*i, ",", e_kin, ",", e_pot, ",", e_tot

    ! write trajectory
    write (16, "(I3)") natom
    write (16, "(F15.8)") e_tot
    do w = 1, natom
      write (16, *) "Ar", xyz(1, w), xyz(2, w), xyz(3, w)
    end do
  end do main_loop
end program main
