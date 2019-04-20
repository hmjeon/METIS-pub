!
! =============================================================================
!
! Module - Exam_METIS
! Last Updated : 01/09/2019, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of METIS, which allows scientists to build and solve
! the sequence design of complex DNA nanostructures.
! Copyright 2019 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! METIS is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! METIS is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
module Exam_METIS

    use Data_Prob
    use Data_Geom

    use Para
    use Math
    use Mani

    implicit none

    ! Triangular-mesh objects
    public Exam_METIS_Hexagon_Tri       ! 01. Hexagon
    public Exam_METIS_Hexagon_Hole      ! 02. Hexagon with a hole
    public Exam_METIS_Square_Tri        ! 03. Square
    public Exam_METIS_Circle            ! 04. Circle

    ! Quadrilateral-mesh objects
    public Exam_METIS_Six_Parallelogram ! 05. Six Parallelogram
    public Exam_METIS_Curved_Beam       ! 06. Curved Beam
    public Exam_METIS_Quarter_Circle    ! 07. Quarter Circle
    public Exam_METIS_Annulus           ! 08. Annulus

    ! Without the internal mesh
    public Exam_METIS_Triangle          ! 09. Triangle
    public Exam_METIS_Square            ! 10. Square
    public Exam_METIS_Hexagon           ! 11. Hexagon
    public Exam_METIS_Octagon           ! 12. Octagon

    ! Letters
    public Exam_METIS_Letter_A          ! 13. Letter, A
    public Exam_METIS_Letter_T          ! 14. Letter, T
    public Exam_METIS_Letter_G          ! 15. Letter, G
    public Exam_METIS_Letter_C          ! 16. Letter, C

    ! Without the internal mesh
    public Exam_METIS_Heart             ! Heart
    public Exam_METIS_Map_Korea         ! Korea
    public Exam_METIS_Map_China         ! China
    public Exam_METIS_Map_US            ! US
    public Exam_METIS_Map_Finland       ! Finland

    ! Variable vertex-number
    public Exam_METIS_N_Sided_Poly      ! N-sided polygon

    public Exam_METIS_Control_Plate

contains

! -----------------------------------------------------------------------------

! Example of the hexagon with the triangular mesh
subroutine Exam_METIS_Hexagon_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: x_width, y_width, del_x, del_y
    integer :: i, j, index, n_i_poi, n_j_poi, n, nx, ny
    character :: pn

    ! Set problem
    prob.name_prob = "01_Hexagon_Tri"
    call Mani_Set_Prob(prob, 'blue')

    ! Set options
    n = 6

    ! The number of points and faces
    geom.n_iniP = n + 1
    geom.n_face = n

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    geom.iniP(n+1).pos(1) = 0.0d0
    geom.iniP(n+1).pos(2) = 0.0d0
    geom.iniP(n+1).pos(3) = 0.0d0

    ! Set position vector
    do i = 1, n
        geom.iniP(i).pos(1) = 10.0d0 * dcos(2.0d0 * pi * dble(i) / dble(n))
        geom.iniP(i).pos(2) = 10.0d0 * dsin(2.0d0 * pi * dble(i) / dble(n))
        geom.iniP(i).pos(3) = 0.0d0
    end do

    ! Set connectivity
    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    do i = 1, n
        if(i == n) then
            geom.face(i).poi(1) = n + 1
            geom.face(i).poi(2) = i
            geom.face(i).poi(3) = 1
        else
            geom.face(i).poi(1) = n + 1
            geom.face(i).poi(2) = i
            geom.face(i).poi(3) = i + 1
        end if
    end do
end subroutine Exam_METIS_Hexagon_Tri

! -----------------------------------------------------------------------------

! Example of the hexagon with a hole
subroutine Exam_METIS_Hexagon_Hole(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "02_Hexagon_Hole"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 19

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   86.6328d0, -150.0566d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [    0.0000d0, -200.1111d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [    0.0000d0, -100.0021d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  -86.6328d0, -150.0566d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [ -173.2656d0, -100.0021d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  -86.6328d0,  -50.0545d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -173.2656d0,   -0.0000d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [ -173.2656d0,  100.0021d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [  -86.6328d0,   50.0545d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [  -86.6328d0,  150.0566d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [    0.0000d0,  200.1111d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [    0.0000d0,  100.0021d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [   86.6328d0,  150.0566d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  173.2656d0,  100.0021d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [   86.6328d0,   50.0545d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  173.2656d0,   -0.0000d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [  173.2656d0, -100.0021d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [   86.6328d0,  -50.0545d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  2,  1 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  4,  2 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  6,  5,  4 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  7,  5 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  9,  8,  7 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [  9, 10,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 12, 11, 10 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 12, 13, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 15, 14, 13 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 15, 16, 14 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 18, 17, 16 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 18,  1, 17 ]
    geom.face(13).n_poi = 6; allocate(geom.face(13).poi(6)); geom.face(13).poi(1:6) = [  9,  6,  3, 18, 15, 12 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 13, 12, 15 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 16, 15, 18 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  1, 18,  3 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  4,  3,  6 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  7,  6,  9 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [ 10,  9, 12 ]
end subroutine Exam_METIS_Hexagon_Hole

! -----------------------------------------------------------------------------

! Example of the square with the triangular mesh
subroutine Exam_METIS_Square_Tri(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "03_Sqaure_Tri"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 18

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  150.0000d0,  -50.0000d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  150.0000d0, -150.0000d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [   50.0000d0, -150.0000d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  -50.0000d0, -150.0000d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  -50.0000d0,  -50.0000d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [ -150.0000d0, -150.0000d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -150.0000d0,  -50.0000d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [ -150.0000d0,   50.0000d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -150.0000d0,  150.0000d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [  -50.0000d0,  150.0000d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [   50.0000d0,  150.0000d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [   50.0000d0,   50.0000d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  150.0000d0,  150.0000d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  150.0000d0,   50.0000d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [  -50.0000d0,   50.0000d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [   50.0000d0,  -50.0000d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  2,  1 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  5,  4,  3 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  5,  6,  4 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  5,  7,  6 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  5,  8,  7 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 10,  9,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 12, 11, 10 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 12, 13, 11 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [ 12, 14, 13 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 12,  1, 14 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 15, 12, 10 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [  5, 12, 15 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  8, 15, 10 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [  5, 15,  8 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 16,  1, 12 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [  3,  1, 16 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [  5, 16, 12 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  3, 16,  5 ]
end subroutine Exam_METIS_Square_Tri

! -----------------------------------------------------------------------------

! Example of the circle
subroutine Exam_METIS_Circle(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "04_Circle"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = 19
    geom.n_face = 24

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -42.8419d0,   0.3155d0, 0.0000d0 ]; geom.iniP(10).pos(1:3) = [   0.0000d0,  -0.2165d0, 0.0000d0 ]
    geom.iniP( 2).pos(1:3) = [ -37.2731d0, -21.2527d0, 0.0000d0 ]; geom.iniP(11).pos(1:3) = [   0.0000d0, -42.9701d0, 0.0000d0 ]
    geom.iniP( 3).pos(1:3) = [ -36.9329d0,  21.5902d0, 0.0000d0 ]; geom.iniP(12).pos(1:3) = [  11.3495d0,  19.5828d0, 0.0000d0 ]
    geom.iniP( 4).pos(1:3) = [ -22.7947d0,  -0.0741d0, 0.0000d0 ]; geom.iniP(13).pos(1:3) = [  11.4451d0, -19.8967d0, 0.0000d0 ]
    geom.iniP( 5).pos(1:3) = [ -21.6133d0, -37.1190d0, 0.0000d0 ]; geom.iniP(14).pos(1:3) = [  21.2310d0,  37.0880d0, 0.0000d0 ]
    geom.iniP( 6).pos(1:3) = [ -21.2310d0,  37.0880d0, 0.0000d0 ]; geom.iniP(15).pos(1:3) = [  21.6133d0, -37.1190d0, 0.0000d0 ]
    geom.iniP( 7).pos(1:3) = [ -11.4451d0, -19.8967d0, 0.0000d0 ]; geom.iniP(16).pos(1:3) = [  22.7947d0,  -0.0741d0, 0.0000d0 ]
    geom.iniP( 8).pos(1:3) = [ -11.3495d0,  19.5828d0, 0.0000d0 ]; geom.iniP(17).pos(1:3) = [  36.9329d0,  21.5902d0, 0.0000d0 ]
    geom.iniP( 9).pos(1:3) = [   0.0000d0,  42.7183d0, 0.0000d0 ]; geom.iniP(18).pos(1:3) = [  37.2731d0, -21.2527d0, 0.0000d0 ]
    geom.iniP(19).pos(1:3) = [  42.8419d0,   0.3155d0, 0.0000d0 ]

    ! Set point position vectors
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [ 12,  8, 10 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  9,  8, 12 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [ 11,  7,  5 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [ 10,  7, 13 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [ 15, 18, 13 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 11, 15, 13 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 13,  7, 11 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 10,  8,  4 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  4,  7, 10 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  1,  4,  3 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  3,  4,  8 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 16, 12, 10 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [ 10, 13, 16 ]
    geom.face(14).n_poi = 3; allocate(geom.face(14).poi(3)); geom.face(14).poi(1:3) = [ 16, 18, 19 ]
    geom.face(15).n_poi = 3; allocate(geom.face(15).poi(3)); geom.face(15).poi(1:3) = [ 16, 13, 18 ]
    geom.face(16).n_poi = 3; allocate(geom.face(16).poi(3)); geom.face(16).poi(1:3) = [ 17, 16, 19 ]
    geom.face(17).n_poi = 3; allocate(geom.face(17).poi(3)); geom.face(17).poi(1:3) = [ 12, 16, 17 ]
    geom.face(18).n_poi = 3; allocate(geom.face(18).poi(3)); geom.face(18).poi(1:3) = [  2,  5,  7 ]
    geom.face(19).n_poi = 3; allocate(geom.face(19).poi(3)); geom.face(19).poi(1:3) = [  2,  4,  1 ]
    geom.face(20).n_poi = 3; allocate(geom.face(20).poi(3)); geom.face(20).poi(1:3) = [  7,  4,  2 ]
    geom.face(21).n_poi = 3; allocate(geom.face(21).poi(3)); geom.face(21).poi(1:3) = [  6,  8,  9 ]
    geom.face(22).n_poi = 3; allocate(geom.face(22).poi(3)); geom.face(22).poi(1:3) = [  6,  3,  8 ]
    geom.face(23).n_poi = 3; allocate(geom.face(23).poi(3)); geom.face(23).poi(1:3) = [  9, 12, 14 ]
    geom.face(24).n_poi = 3; allocate(geom.face(24).poi(3)); geom.face(24).poi(1:3) = [ 12, 17, 14 ]
end subroutine Exam_METIS_Circle

! -----------------------------------------------------------------------------

! Example of the 6 parallelogram
subroutine Exam_METIS_Six_Parallelogram(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "05_Six_Parallelogram"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 13
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   86.6212d0,  -49.9800d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [   86.6212d0, -150.0756d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [    0.0000d0, -100.0278d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [    0.0000d0,   -0.0123d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  -86.6212d0, -150.0756d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  -86.6212d0,  -49.9800d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -173.2425d0,   -0.0123d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  -86.6212d0,   50.0355d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [  -86.6212d0,  150.0509d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [    0.0000d0,  100.0031d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [   86.6212d0,  150.0509d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [   86.6212d0,   50.0355d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  173.2425d0,   -0.0123d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [ 4,  3,  2,  1 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [ 4,  6,  5,  3 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [ 4,  8,  7,  6 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [ 4, 10,  9,  8 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 4, 12, 11, 10 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 4,  1, 13, 12 ]
end subroutine Exam_METIS_Six_Parallelogram

! -----------------------------------------------------------------------------

! Example of the curved beam
subroutine Exam_METIS_Curved_Beam(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "06_Curved_Beam"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 6

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  -99.6791d0, -171.5371d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [ -112.7442d0, -270.7541d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [ -227.7428d0, -270.7541d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [ -210.7970d0, -141.7849d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [ -160.9944d0,  -21.4826d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  -61.3893d0,  -79.0466d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  -81.6984d0,   81.7446d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [   -0.3326d0,    0.3788d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [   21.5288d0,  161.0406d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [   78.9634d0,   61.3061d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [  141.7017d0,  210.8432d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [  171.5833d0,   99.5959d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  270.8003d0,  227.7890d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  270.8003d0,  112.6610d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [  4,  3,  2,  1 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [  1,  6,  5,  4 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [  6,  8,  7,  5 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [  8, 10,  9,  7 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [ 10, 12, 11,  9 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 12, 14, 13, 11 ]
end subroutine Exam_METIS_Curved_Beam

! -----------------------------------------------------------------------------

! Example of the quarter_circle
subroutine Exam_METIS_Quarter_Circle(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "07_Quarter_Circle"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 8

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  111.6174d0, -156.9291d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  -22.6599d0, -156.9291d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  -44.9640d0,  -44.9559d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [   66.8960d0,  -67.3732d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [ -156.9372d0, -156.9291d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [ -156.9372d0,  -22.6518d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -156.9372d0,  111.6255d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  -67.3813d0,   66.9041d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -156.9372d0,  245.9028d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [   -9.4133d0,  203.1061d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [  122.4864d0,  122.4945d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [   22.1747d0,   22.0695d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  203.0981d0,   -9.4052d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  245.8947d0, -156.9291d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = 4; allocate(geom.face(1).poi(4)); geom.face(1).poi(1:4) = [  4,  3,  2,  1 ]
    geom.face(2).n_poi = 4; allocate(geom.face(2).poi(4)); geom.face(2).poi(1:4) = [  3,  6,  5,  2 ]
    geom.face(3).n_poi = 4; allocate(geom.face(3).poi(4)); geom.face(3).poi(1:4) = [  3,  8,  7,  6 ]
    geom.face(4).n_poi = 4; allocate(geom.face(4).poi(4)); geom.face(4).poi(1:4) = [  8, 10,  9,  7 ]
    geom.face(5).n_poi = 4; allocate(geom.face(5).poi(4)); geom.face(5).poi(1:4) = [  8, 12, 11, 10 ]
    geom.face(6).n_poi = 4; allocate(geom.face(6).poi(4)); geom.face(6).poi(1:4) = [ 12,  4, 13, 11 ]
    geom.face(7).n_poi = 4; allocate(geom.face(7).poi(4)); geom.face(7).poi(1:4) = [  4,  1, 14, 13 ]
    geom.face(8).n_poi = 4; allocate(geom.face(8).poi(4)); geom.face(8).poi(1:4) = [  8,  3,  4, 12 ]
end subroutine Exam_METIS_Quarter_Circle

! -----------------------------------------------------------------------------

! Example of the annulus
subroutine Exam_METIS_Annulus(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "08_Annulus"
    call Mani_Set_Prob(prob, 'red')

    ! The number of points and faces
    geom.n_iniP = 20
    geom.n_face = 11

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  263.7559d0,   -0.0000d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  213.3808d0, -154.9888d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  132.3943d0,  -96.2922d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  163.7487d0,   -0.0000d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [   81.4249d0, -250.8352d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [   50.5163d0, -155.7318d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  -81.4397d0, -250.8352d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  -50.5311d0, -155.7318d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -213.3957d0, -154.9888d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -132.4092d0,  -96.2922d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -263.7708d0,   -0.0000d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -163.6150d0,   -0.0000d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [ -213.3957d0,  154.9888d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [ -132.4092d0,   96.2922d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [  -81.4397d0,  250.8352d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  -50.5311d0,  155.7318d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [   81.4249d0,  250.8352d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [   50.5163d0,  155.7318d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [  213.3808d0,  154.9888d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  132.3943d0,   96.2922d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi =  4; allocate(geom.face(  1).poi(  4)); geom.face(  1).poi(1: 4) = [  4,  3,  2,  1 ]
    geom.face( 2).n_poi =  4; allocate(geom.face(  2).poi(  4)); geom.face(  2).poi(1: 4) = [  3,  6,  5,  2 ]
    geom.face( 3).n_poi =  4; allocate(geom.face(  3).poi(  4)); geom.face(  3).poi(1: 4) = [  6,  8,  7,  5 ]
    geom.face( 4).n_poi =  4; allocate(geom.face(  4).poi(  4)); geom.face(  4).poi(1: 4) = [  8, 10,  9,  7 ]
    geom.face( 5).n_poi =  4; allocate(geom.face(  5).poi(  4)); geom.face(  5).poi(1: 4) = [ 10, 12, 11,  9 ]
    geom.face( 6).n_poi =  4; allocate(geom.face(  6).poi(  4)); geom.face(  6).poi(1: 4) = [ 12, 14, 13, 11 ]
    geom.face( 7).n_poi =  4; allocate(geom.face(  7).poi(  4)); geom.face(  7).poi(1: 4) = [ 14, 16, 15, 13 ]
    geom.face( 8).n_poi =  4; allocate(geom.face(  8).poi(  4)); geom.face(  8).poi(1: 4) = [ 16, 18, 17, 15 ]
    geom.face( 9).n_poi =  4; allocate(geom.face(  9).poi(  4)); geom.face(  9).poi(1: 4) = [ 18, 20, 19, 17 ]
    geom.face(10).n_poi =  4; allocate(geom.face( 10).poi(  4)); geom.face( 10).poi(1: 4) = [ 20,  4,  1, 19 ]
    geom.face(11).n_poi = 10; allocate(geom.face( 11).poi( 10)); geom.face( 11).poi(1:10) = [ 20, 18, 16, 14, 12, 10,  8,  6,  3,  4 ]
end subroutine Exam_METIS_Annulus

! -----------------------------------------------------------------------------

! Example of the letter, A
subroutine Exam_METIS_Letter_A(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "13_Letter_A"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 11

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  307.7862d0, -138.2976d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  352.9592d0, -272.2034d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  148.0672d0, -265.7501d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  109.3474d0, -131.8443d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  -84.2513d0, -125.3910d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  -51.9849d0,  -18.9117d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [   80.3076d0,  -18.9117d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [ -127.8110d0, -260.9102d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -310.1165d0, -257.6835d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -271.3968d0, -128.6177d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -237.5170d0,  -17.2984d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -206.8639d0,   77.8876d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [ -163.3041d0,  206.9535d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [    7.7081d0,  197.2735d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [  -90.7046d0,  413.4588d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  122.2540d0,  413.4588d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [  175.4937d0,  258.5798d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [  240.0266d0,   68.2077d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  3,  2,  1 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  7,  6,  5,  4 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [ 10,  9,  8,  5 ]
    geom.face( 4).n_poi = 5; allocate(geom.face( 4).poi(5)); geom.face( 4).poi(1:5) = [  5,  6, 12, 11, 10 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  6, 14, 13, 12 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 14, 15, 13 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 14, 16, 15 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 14, 17, 16 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 14,  7, 18, 17 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [  7,  4,  1, 18 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [  7, 14,  6 ]
end subroutine Exam_METIS_Letter_A

! -----------------------------------------------------------------------------

! Example of the letter, T
subroutine Exam_METIS_Letter_T(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "14_Letter_T"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 16
    geom.n_face = 11

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [   94.4446d0, -249.3668d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [   94.4446d0, -405.4856d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  -51.3646d0, -396.6487d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  -48.4189d0, -247.8939d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  -54.3102d0,  -88.8294d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [   92.9717d0,  -97.6663d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  -48.4189d0,   71.7079d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [ -164.7717d0,   80.5448d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -132.3697d0,  220.4627d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [  -11.5985d0,  216.0442d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -263.4506d0,   96.7458d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -256.0865d0,  213.0986d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  109.1727d0,  221.9355d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  103.2815d0,   76.1264d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [  269.7101d0,  216.0442d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  266.7644d0,   73.1807d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi = 3; allocate(geom.face( 1).poi(3)); geom.face( 1).poi(1:3) = [  3,  2,  1 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  1,  4,  3 ]
    geom.face( 3).n_poi = 3; allocate(geom.face( 3).poi(3)); geom.face( 3).poi(1:3) = [  6,  5,  4 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  6,  7,  5 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [ 10,  9,  8,  7 ]
    geom.face( 6).n_poi = 4; allocate(geom.face( 6).poi(4)); geom.face( 6).poi(1:4) = [  9, 12, 11,  8 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 14, 13, 10 ]
    geom.face( 8).n_poi = 4; allocate(geom.face( 8).poi(4)); geom.face( 8).poi(1:4) = [ 14, 16, 15, 13 ]
    geom.face( 9).n_poi = 3; allocate(geom.face( 9).poi(3)); geom.face( 9).poi(1:3) = [  7,  6, 14 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [  4,  1,  6 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 14, 10,  7 ]
end subroutine Exam_METIS_Letter_T

! -----------------------------------------------------------------------------

! Example of the letter, G
subroutine Exam_METIS_Letter_G(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "15_Letter_G"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 24
    geom.n_face = 13

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  113.1343d0, -324.8157d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  -28.8381d0, -334.4956d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  -25.6115d0, -173.1633d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  105.0677d0, -178.0033d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [ -169.1972d0, -279.6426d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [ -283.7432d0, -171.5500d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -122.4109d0,  -94.1105d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [ -312.7830d0,  -16.6710d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -140.1574d0,    4.3022d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -311.1697d0,  128.5281d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -117.5709d0,  146.2746d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -230.5035d0,  262.4339d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  -62.7179d0,  346.3267d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [   38.9214d0,  199.5143d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [  134.1075d0,  351.1666d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  238.9735d0,  243.0740d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [  161.5340d0,  141.4346d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [  119.5876d0,  -57.0041d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [  268.0133d0,  -87.6572d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  259.9467d0, -218.3364d0, 0.0d0 ]
    geom.iniP(21).pos(1:3) = [    3.4283d0,  -61.8440d0, 0.0d0 ]
    geom.iniP(22).pos(1:3) = [   -6.2516d0,   70.4484d0, 0.0d0 ]
    geom.iniP(23).pos(1:3) = [  117.9742d0,   72.0618d0, 0.0d0 ]
    geom.iniP(24).pos(1:3) = [  250.2667d0,   31.7287d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  3,  2,  1 ]
    geom.face( 2).n_poi = 3; allocate(geom.face( 2).poi(3)); geom.face( 2).poi(1:3) = [  3,  5,  2 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  3,  7,  6,  5 ]
    geom.face( 4).n_poi = 4; allocate(geom.face( 4).poi(4)); geom.face( 4).poi(1:4) = [  7,  9,  8,  6 ]
    geom.face( 5).n_poi = 4; allocate(geom.face( 5).poi(4)); geom.face( 5).poi(1:4) = [  9, 11, 10,  8 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 11, 12, 10 ]
    geom.face( 7).n_poi = 4; allocate(geom.face( 7).poi(4)); geom.face( 7).poi(1:4) = [ 11, 14, 13, 12 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 14, 15, 13 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 14, 17, 16, 15 ]
    geom.face(10).n_poi = 4; allocate(geom.face(10).poi(4)); geom.face(10).poi(1:4) = [ 20, 19, 18,  4 ]
    geom.face(11).n_poi = 4; allocate(geom.face(11).poi(4)); geom.face(11).poi(1:4) = [ 23, 22, 21, 18 ]
    geom.face(12).n_poi = 4; allocate(geom.face(12).poi(4)); geom.face(12).poi(1:4) = [ 18, 19, 24, 23 ]
    geom.face(13).n_poi = 3; allocate(geom.face(13).poi(3)); geom.face(13).poi(1:3) = [  4,  1, 20 ]
end subroutine Exam_METIS_Letter_G

! -----------------------------------------------------------------------------

! Example of the letter, C
subroutine Exam_METIS_Letter_C(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    prob.name_prob = "16_Letter_C"
    call Mani_Set_Prob(prob, 'green')

    ! The number of points and faces
    geom.n_iniP = 18
    geom.n_face = 12

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [    4.1990d0, -164.0915d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  103.4937d0, -175.9475d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  102.0117d0, -313.7745d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  -62.4914d0, -291.5443d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  187.9683d0, -115.1851d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  248.7307d0, -235.2279d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -182.5342d0, -230.7819d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  -78.7935d0,  -84.0629d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -253.6707d0, -107.7751d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -253.6707d0,   34.4979d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -231.4406d0,  181.2169d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [  -89.1676d0,   70.0661d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [ -106.9517d0,  283.4756d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [   23.4652d0,  317.5618d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [   26.4292d0,  175.2889d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  150.9180d0,  299.7777d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [  156.8461d0,  135.2746d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [  254.6587d0,  221.2312d0, 0.0d0 ]

    ! Set face connectivity
    geom.face( 1).n_poi = 4; allocate(geom.face( 1).poi(4)); geom.face( 1).poi(1:4) = [  4,  3,  2,  1 ]
    geom.face( 2).n_poi = 4; allocate(geom.face( 2).poi(4)); geom.face( 2).poi(1:4) = [  3,  6,  5,  2 ]
    geom.face( 3).n_poi = 4; allocate(geom.face( 3).poi(4)); geom.face( 3).poi(1:4) = [  1,  8,  7,  4 ]
    geom.face( 4).n_poi = 3; allocate(geom.face( 4).poi(3)); geom.face( 4).poi(1:3) = [  8,  9,  7 ]
    geom.face( 5).n_poi = 3; allocate(geom.face( 5).poi(3)); geom.face( 5).poi(1:3) = [  8, 10,  9 ]
    geom.face( 6).n_poi = 3; allocate(geom.face( 6).poi(3)); geom.face( 6).poi(1:3) = [ 12, 11, 10 ]
    geom.face( 7).n_poi = 3; allocate(geom.face( 7).poi(3)); geom.face( 7).poi(1:3) = [ 12, 13, 11 ]
    geom.face( 8).n_poi = 3; allocate(geom.face( 8).poi(3)); geom.face( 8).poi(1:3) = [ 15, 14, 13 ]
    geom.face( 9).n_poi = 4; allocate(geom.face( 9).poi(4)); geom.face( 9).poi(1:4) = [ 15, 17, 16, 14 ]
    geom.face(10).n_poi = 3; allocate(geom.face(10).poi(3)); geom.face(10).poi(1:3) = [ 17, 18, 16 ]
    geom.face(11).n_poi = 3; allocate(geom.face(11).poi(3)); geom.face(11).poi(1:3) = [ 13, 12, 15 ]
    geom.face(12).n_poi = 3; allocate(geom.face(12).poi(3)); geom.face(12).poi(1:3) = [ 10,  8, 12 ]
end subroutine Exam_METIS_Letter_C

! -----------------------------------------------------------------------------

! Example of N-sided polygon
subroutine Exam_METIS_N_Sided_Poly(prob, geom, n)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom
    integer,        intent(in)    :: n

    integer :: i

    ! Set problem
    if(n == 4) prob.name_prob = "4_Sided_Poly"
    if(n == 5) prob.name_prob = "5_Sided_Poly"
    if(n == 6) prob.name_prob = "6_Sided_Poly"
    call Mani_Set_Prob(prob, 'blue')

    ! The number of points and faces
    geom.n_iniP = n + 1
    geom.n_face = n

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    geom.iniP(n+1).pos(1) = 0.0d0
    geom.iniP(n+1).pos(2) = 0.0d0
    geom.iniP(n+1).pos(3) = 0.0d0

    ! Set the position vector
    do i = 1, n
        geom.iniP(i).pos(1) = 10.0d0 * dcos(2.0d0 * pi * dble(i) / dble(n))
        geom.iniP(i).pos(2) = 10.0d0 * dsin(2.0d0 * pi * dble(i) / dble(n))
        geom.iniP(i).pos(3) = 0.0d0
    end do

    ! Set face connectivity
    do i = 1, geom.n_face
        geom.face(i).n_poi = 3
        allocate(geom.face(i).poi(3))
    end do

    do i = 1, n
        if(i == n) then
            geom.face(i).poi(1) = n + 1
            geom.face(i).poi(2) = i
            geom.face(i).poi(3) = 1
        else
            geom.face(i).poi(1) = n + 1
            geom.face(i).poi(2) = i
            geom.face(i).poi(3) = i + 1
        end if
    end do

    ! Set the orientation of the geometry
    if(n == 5) call Mani_Set_2DGeo_Ori(geom, [0.0d0, 0.0d0, 1.0d0], 18.0d0)
end subroutine Exam_METIS_N_Sided_Poly

! -----------------------------------------------------------------------------

! Example of the triangle
! https://www.calculatorsoup.com/calculators/geometry-plane/polygon.php
subroutine Exam_METIS_Triangle(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "09_Triangle"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 3
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [   0.0000d0,  57.7393d0,  0.0d0 ]
    geom.iniP(2).pos(1:3) = [  50.0000d0, -28.8697d0,  0.0d0 ]
    geom.iniP(3).pos(1:3) = [ -50.0000d0, -28.8697d0,  0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Triangle

! -----------------------------------------------------------------------------

! Example of the square
subroutine Exam_METIS_Square(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "10_Square"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 4
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [  50.0000d0,  50.0098d0, 0.0d0 ]
    geom.iniP(2).pos(1:3) = [  50.0000d0, -50.0098d0, 0.0d0 ]
    geom.iniP(3).pos(1:3) = [ -50.0000d0, -50.0098d0, 0.0d0 ]
    geom.iniP(4).pos(1:3) = [ -50.0000d0,  50.0098d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Square

! -----------------------------------------------------------------------------

! Example of the hexagon
subroutine Exam_METIS_Hexagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "11_Hexagon"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 6
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [   50.0046d0,  86.6112d0, 0.0d0 ]
    geom.iniP(2).pos(1:3) = [  100.0091d0,   0.0113d0, 0.0d0 ]
    geom.iniP(3).pos(1:3) = [   50.0046d0, -86.6225d0, 0.0d0 ]
    geom.iniP(4).pos(1:3) = [  -50.0046d0, -86.6225d0, 0.0d0 ]
    geom.iniP(5).pos(1:3) = [ -100.0091d0,   0.0113d0, 0.0d0 ]
    geom.iniP(6).pos(1:3) = [  -50.0046d0,  86.6112d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Hexagon

! -----------------------------------------------------------------------------

! Example of the octagon
subroutine Exam_METIS_Octagon(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "12_Octagon"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 8
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP(1).pos(1:3) = [   50.0000d0,  120.7409d0, 0.0d0 ]
    geom.iniP(2).pos(1:3) = [  120.7631d0,   50.0222d0, 0.0d0 ]
    geom.iniP(3).pos(1:3) = [  120.7631d0,  -50.0222d0, 0.0d0 ]
    geom.iniP(4).pos(1:3) = [   50.0000d0, -120.7409d0, 0.0d0 ]
    geom.iniP(5).pos(1:3) = [  -50.0000d0, -120.7409d0, 0.0d0 ]
    geom.iniP(6).pos(1:3) = [ -120.7631d0,  -50.0222d0, 0.0d0 ]
    geom.iniP(7).pos(1:3) = [ -120.7631d0,   50.0222d0, 0.0d0 ]
    geom.iniP(8).pos(1:3) = [  -50.0000d0,  120.7409d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Octagon

! -----------------------------------------------------------------------------

! Example of the heart
subroutine Exam_METIS_Heart(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "Heart"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 14
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [    0.0000d0,  109.5660d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [   75.9962d0,  174.5627d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  181.1909d0,  196.5616d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  286.3857d0,  134.9647d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  310.9845d0,    6.7711d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  255.5872d0, -136.4218d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  146.3927d0, -256.6157d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [    0.0000d0, -349.2111d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -146.3927d0, -256.6157d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -255.5872d0, -136.4218d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -310.9845d0,    6.7711d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -286.3857d0,  134.9647d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [ -181.1909d0,  196.5616d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  -75.9962d0,  174.5627d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Heart

! -----------------------------------------------------------------------------

! Example of the map of Korea
subroutine Exam_METIS_Map_Korea(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "Korea"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 27
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  314.6327d0,  358.8431d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [   14.2873d0,  154.4122d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [    0.9853d0,   50.7966d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  138.9062d0,  -47.2182d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  276.1269d0, -287.3545d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  302.0308d0, -500.1866d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  254.4236d0, -673.8129d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [   75.8966d0, -706.0177d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -119.4329d0, -847.4391d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -190.1435d0, -776.7284d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -111.7317d0, -554.0948d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -190.1435d0, -370.6671d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  -59.2238d0, -337.7621d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [ -154.4381d0, -201.2415d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [ -311.9619d0, -211.7431d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [ -377.7719d0, -107.4273d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [ -293.0591d0,  102.6044d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [ -440.7814d0,  188.0173d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [ -344.8669d0,  307.7354d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [ -199.2449d0,  375.6456d0, 0.0d0 ]
    geom.iniP(21).pos(1:3) = [  -96.3294d0,  518.4672d0, 0.0d0 ]
    geom.iniP(22).pos(1:3) = [   74.4964d0,  457.5580d0, 0.0d0 ]
    geom.iniP(23).pos(1:3) = [   92.6992d0,  569.5749d0, 0.0d0 ]
    geom.iniP(24).pos(1:3) = [  278.2272d0,  647.9867d0, 0.0d0 ]
    geom.iniP(25).pos(1:3) = [  328.6348d0,  755.8030d0, 0.0d0 ]
    geom.iniP(26).pos(1:3) = [  434.3508d0,  624.8833d0, 0.0d0 ]
    geom.iniP(27).pos(1:3) = [  303.4310d0,  509.3658d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Map_Korea

! -----------------------------------------------------------------------------

! Example of the map of China
subroutine Exam_METIS_Map_China(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "China"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 27
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  344.8121d0,   23.7372d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  445.6048d0, -137.0207d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  486.8575d0, -308.8361d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  378.8349d0, -495.1111d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  127.4912d0, -593.3521d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  -37.9449d0, -535.0880d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [ -197.0017d0, -593.3521d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [ -256.9670d0, -450.0309d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [ -253.5647d0, -350.0888d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [ -356.4838d0, -302.0315d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [ -574.6553d0, -333.5026d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [ -777.5165d0, -187.2044d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [ -831.1025d0,    9.2775d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [ -895.7459d0,  182.3687d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [ -692.8847d0,  228.7249d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [ -493.0005d0,  423.5057d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [ -283.3347d0,  185.7710d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [  -37.9449d0,  130.9092d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [  165.7669d0,  220.6444d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  322.6972d0,  340.1497d0, 0.0d0 ]
    geom.iniP(21).pos(1:3) = [  237.2148d0,  446.8964d0, 0.0d0 ]
    geom.iniP(22).pos(1:3) = [  373.7315d0,  606.8037d0, 0.0d0 ]
    geom.iniP(23).pos(1:3) = [  566.8111d0,  480.9192d0, 0.0d0 ]
    geom.iniP(24).pos(1:3) = [  707.1554d0,  472.4135d0, 0.0d0 ]
    geom.iniP(25).pos(1:3) = [  645.0637d0,  270.4029d0, 0.0d0 ]
    geom.iniP(26).pos(1:3) = [  498.7654d0,  142.3919d0, 0.0d0 ]
    geom.iniP(27).pos(1:3) = [  387.3406d0,  120.7023d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Map_China

! -----------------------------------------------------------------------------

! Example of the map of US
subroutine Exam_METIS_Map_US(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "Unite_States"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 28
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [ -724.2248d0,  434.4027d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [ -384.6789d0,  362.9877d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [  -74.6728d0,  356.4954d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [    7.7791d0,  248.3991d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  157.4261d0,  286.7035d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  238.9041d0,  101.3490d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  411.2740d0,  199.3824d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  484.3121d0,  271.4466d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [  544.6903d0,  401.6168d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [  663.8235d0,  327.2802d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [  562.8686d0,  231.1945d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [  501.1920d0,  110.4382d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  438.5416d0,   20.8448d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  473.9245d0, -116.4668d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [  333.3667d0, -229.1077d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [  385.6295d0, -325.5180d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [  406.4048d0, -516.0663d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [  281.1039d0, -346.6179d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [  181.1229d0, -348.5656d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  -66.2328d0, -368.6916d0, 0.0d0 ]
    geom.iniP(21).pos(1:3) = [ -199.3244d0, -519.9616d0, 0.0d0 ]
    geom.iniP(22).pos(1:3) = [ -272.3625d0, -390.4408d0, 0.0d0 ]
    geom.iniP(23).pos(1:3) = [ -385.3281d0, -389.1423d0, 0.0d0 ]
    geom.iniP(24).pos(1:3) = [ -479.7907d0, -295.6536d0, 0.0d0 ]
    geom.iniP(25).pos(1:3) = [ -798.5614d0, -198.5940d0, 0.0d0 ]
    geom.iniP(26).pos(1:3) = [ -913.1500d0,   22.4679d0, 0.0d0 ]
    geom.iniP(27).pos(1:3) = [ -918.0192d0,  230.8699d0, 0.0d0 ]
    geom.iniP(28).pos(1:3) = [ -856.0180d0,  438.9473d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Map_US

! -----------------------------------------------------------------------------

! Example of the map of Finland
subroutine Exam_METIS_Map_Finland(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! Set problem
    prob.name_prob = "Finland"
    call Mani_Set_Prob(prob, 'violet')

    ! The number of points and faces
    geom.n_iniP = 25
    geom.n_face = 1

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set point position vectors
    geom.iniP( 1).pos(1:3) = [  -89.9544d0,   518.0807d0, 0.0d0 ]
    geom.iniP( 2).pos(1:3) = [  -25.7360d0,   748.7797d0, 0.0d0 ]
    geom.iniP( 3).pos(1:3) = [   65.5697d0,   789.5629d0, 0.0d0 ]
    geom.iniP( 4).pos(1:3) = [  174.8321d0,   708.3009d0, 0.0d0 ]
    geom.iniP( 5).pos(1:3) = [  134.6576d0,   504.6892d0, 0.0d0 ]
    geom.iniP( 6).pos(1:3) = [  252.1376d0,   357.3827d0, 0.0d0 ]
    geom.iniP( 7).pos(1:3) = [  194.9194d0,   202.1630d0, 0.0d0 ]
    geom.iniP( 8).pos(1:3) = [  284.3990d0,    14.0733d0, 0.0d0 ]
    geom.iniP( 9).pos(1:3) = [  278.3119d0,  -139.6247d0, 0.0d0 ]
    geom.iniP(10).pos(1:3) = [  343.7477d0,  -262.2786d0, 0.0d0 ]
    geom.iniP(11).pos(1:3) = [  303.2688d0,  -359.9757d0, 0.0d0 ]
    geom.iniP(12).pos(1:3) = [  442.3578d0,  -460.1077d0, 0.0d0 ]
    geom.iniP(13).pos(1:3) = [  352.2695d0,  -638.7625d0, 0.0d0 ]
    geom.iniP(14).pos(1:3) = [  165.7016d0,  -902.9403d0, 0.0d0 ]
    geom.iniP(15).pos(1:3) = [ -257.6525d0, -1032.8987d0, 0.0d0 ]
    geom.iniP(16).pos(1:3) = [ -411.3504d0,  -870.9833d0, 0.0d0 ]
    geom.iniP(17).pos(1:3) = [ -379.0891d0,  -679.8500d0, 0.0d0 ]
    geom.iniP(18).pos(1:3) = [ -395.2197d0,  -495.7169d0, 0.0d0 ]
    geom.iniP(19).pos(1:3) = [ -112.4764d0,  -188.6254d0, 0.0d0 ]
    geom.iniP(20).pos(1:3) = [  -68.6497d0,   -43.1450d0, 0.0d0 ]
    geom.iniP(21).pos(1:3) = [ -182.4775d0,    79.5090d0, 0.0d0 ]
    geom.iniP(22).pos(1:3) = [ -186.1297d0,   392.0788d0, 0.0d0 ]
    geom.iniP(23).pos(1:3) = [ -387.0022d0,   565.2553d0, 0.0d0 ]
    geom.iniP(24).pos(1:3) = [ -286.5660d0,   669.0394d0, 0.0d0 ]
    geom.iniP(25).pos(1:3) = [ -209.8692d0,   525.9938d0, 0.0d0 ]

    ! Set face connectivity
    geom.face(1).n_poi = geom.n_iniP
    allocate(geom.face(1).poi(geom.n_iniP))
    do i = 1, geom.n_iniP
        geom.face(1).poi(i) = geom.n_iniP - i + 1
    end do
end subroutine Exam_METIS_Map_Finland

! -----------------------------------------------------------------------------

! Example of the controllable plate
subroutine Exam_METIS_Control_Plate(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    double precision :: ang_t = 0.0d0, ang_b = 0.0d0
    double precision :: del_x, del_y, angle, point_x, point_y, point_z
    integer :: i, j, index, x_num, y_num, n_i_poi, n_j_poi, n_i_elem, n_j_elem

    ! Set problem
    prob.name_prob = "Control_Plate"
    call Mani_Set_Prob(prob, 'blue')

    ! Geoemtry options
    x_num = 4
    y_num = 2
    ang_t = 90.0d0
    !ang_b = 45.0d0

    n_i_poi  = x_num+1
    n_j_poi  = y_num+1
    n_i_elem = n_i_poi-1
    n_j_elem = n_j_poi-1
    del_x    = dble(x_num)/(n_i_poi-1)
    del_y    = dble(y_num)/(n_j_poi-1)

    ! The number of points and faces
    geom.n_iniP = n_i_poi*n_j_poi
    geom.n_face = n_i_elem*n_j_elem

    allocate(geom.iniP(geom.n_iniP))
    allocate(geom.face(geom.n_face))

    ! Set position vector
    do j = 1, n_j_poi
        do i = 1, n_i_poi
            index = n_i_poi * (j - 1) + i

            geom.iniP(index).pos(1) = del_x * (i - 1)
            geom.inip(index).pos(2) = del_y * (j - y_num/2 - 1)
            geom.inip(index).pos(3) = 0.0d0
        end do
    end do

    ! Set connectivity
    do i = 1, geom.n_face
        geom.face(i).n_poi = 4
        allocate(geom.face(i).poi(4))
    end do

    do j = 1, n_j_elem
        do i = 1, n_i_elem
            index = n_i_elem * (j - 1) + i

            geom.face(index).poi(1) = n_i_poi * (j - 1) + i
            geom.face(index).poi(2) = n_i_poi * (j - 1) + i + 1
            geom.face(index).poi(3) = n_i_poi * j + i + 1
            geom.face(index).poi(4) = n_i_poi * j + i
        end do
    end do

    ! Make the beam twist
    do i = 1, geom.n_iniP

        ! x_num : ang_t = position : angle
        angle = Deg2Rad(ang_t)/dble(x_num)*geom.inip(i).pos(1)

        point_x = geom.inip(i).pos(1)
        point_y = dcos(angle) * geom.inip(i).pos(2) - dsin(angle) * geom.inip(i).pos(3)
        point_z = dsin(angle) * geom.inip(i).pos(2) + dcos(angle) * geom.inip(i).pos(3)

        geom.inip(i).pos(1) = point_x
        geom.inip(i).pos(2) = point_y
        geom.inip(i).pos(3) = point_z
    end do

    ! Make the beam bent
    do i = 1, geom.n_iniP

        ! x_num : ang_t = position : angle
        angle = Deg2Rad(ang_b)/dble(x_num)*geom.inip(i).pos(1)

        point_y = geom.inip(i).pos(2)
        point_x = dcos(angle) * geom.inip(i).pos(1) - dsin(angle) * geom.inip(i).pos(3)
        point_z = dsin(angle) * geom.inip(i).pos(1) + dcos(angle) * geom.inip(i).pos(3)

        geom.inip(i).pos(1) = -point_x
        geom.inip(i).pos(2) = +point_y
        geom.inip(i).pos(3) = -point_z
    end do
end subroutine Exam_METIS_Control_Plate

! -----------------------------------------------------------------------------

end module Exam_METIS