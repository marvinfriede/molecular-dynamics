module build_box
  implicit none
  private

  !> floating point accuracy
  integer, parameter :: wp = selected_real_kind(15)

  public :: build_grid
contains
  subroutine build_grid(xyz, natom, l, formula_units)
    implicit none

    !> number of atoms
    integer, intent(in) :: natom

    !> atom coordinates of the system
    real(wp), dimension(3, natom), intent(out) :: xyz

    !> number of formula units per unit cell
    integer, intent(in) :: formula_units

    !> number of lattice points along side of cube
    integer :: nl

    !> iteration
    integer :: n, i, j, k

    !> coordinates
    real(wp) :: x, y, z

    !> length of cubic box
    real(wp), intent(in) :: l
    !> half length of cubic box
    real(wp) :: hl
    !> distance of atoms along side of cube
    real(wp) :: dl

    ! output file: fort.14
    open (unit=14, file="out/box.xyz")
    write (14, "(1X,I3)") natom
    write (14, *) ""

    nl = int((natom/formula_units)**(1.0_wp/3.0_wp))

    if ((nl*formula_units)**3 .lt. natom) then
      nl = nl + 1
    end if

    hl = l/2.0_wp
    dl = l/nl

    n = 1
    outer: do i = 0, nl - 1
      do j = 0, nl - 1
        do k = 0, nl - 1
          x = i*dl - hl
          y = j*dl - hl
          z = k*dl - hl
          xyz(1, n) = x
          xyz(2, n) = y
          xyz(3, n) = z
          write (14, *) "Ar", x, y, z

          if (formula_units == 4) then
            x = i*dl - hl
            y = (j + 0.5_wp)*dl - hl
            z = (k + 0.5_wp)*dl - hl
            xyz(1, n + 1) = x
            xyz(2, n + 1) = y
            xyz(3, n + 1) = z
            write (14, *) "Ar", x, y, z

            x = (i + 0.5_wp)*dl - hl
            y = j*dl - hl
            z = (k + 0.5_wp)*dl - hl
            xyz(1, n + 2) = x
            xyz(2, n + 2) = y
            xyz(3, n + 2) = z
            write (14, *) "Ar", x, y, z

            x = (i + 0.5_wp)*dl - hl
            y = (j + 0.5_wp)*dl - hl
            z = k*dl - hl
            xyz(1, n + 3) = x
            xyz(2, n + 3) = y
            xyz(3, n + 3) = z
            write (14, *) "Ar", x, y, z
          end if

          n = n + formula_units
          if ((natom + 1) == n) exit outer
        end do
      end do
    end do outer

  end subroutine build_grid
end module build_box
