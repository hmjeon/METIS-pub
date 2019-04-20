!
! =============================================================================
!
! Module - Data_Geom
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
module Data_Geom

! -----------------------------------------------------------------------------

    ! Section data for arbitrary cross-section
    type :: SecType
        character(len=10) :: types          ! Lattice type, squre or honeycomb

        integer :: dir                      ! Section convection of caDNAnano
        integer :: maxR, minR               ! Maximum and minimum row
        integer :: maxC, minC               ! Maximum and minimum column
        integer :: n_row, n_col             ! Size of row and column
        integer :: ref_row                  ! Reference row to set t-axis
        integer :: ref_maxC, ref_minC       ! Maximum and minimum column in reference row

        integer, allocatable :: id(:)       ! Section ID
        integer, allocatable :: posR(:)     ! Row position number
        integer, allocatable :: posC(:)     ! Column position number

        ! Connetivity for self connection route
        ! -1 : neighbor connection, num : section number to be connected for self-connection
        integer, allocatable :: conn(:)
    end type SecType

! -----------------------------------------------------------------------------

    ! Point type data structure
    type :: PointType
        double precision :: pos(3)      ! Position vector
        double precision :: ori_pos(3)  ! Original position vector
    end type PointType

! -----------------------------------------------------------------------------

    ! Line type data structure
    type :: LineType
        integer :: iniL             ! Inital line number
        integer :: sec              ! Cross-section ID
        integer :: iniP(2)          ! Initial point connectivity
        integer :: poi(2)           ! Point connectivity
        integer :: neiF(2)          ! Neighbor faces  (direction)
        integer :: neiP(2, 2)       ! Neighbor points (point number, direction)
        integer :: neiL(2, 2)       ! Neighbor lines  (point number, direction)
                                    ! Point number : 1 - starting, 2 - ending
                                    ! Direction    : 1 - left, 2 - right
        integer :: n_xover          ! The number of scaffold crossovers
        double precision :: t(3, 3) ! Local vector at the center
    end type LineType

! -----------------------------------------------------------------------------

    ! Face type data structure
    type :: FaceType
        integer :: n_poi                ! The number of points
        integer, allocatable :: poi(:)  ! Connectivity
    end type FaceType
    
! -----------------------------------------------------------------------------

    ! Junction data structure
    type :: JuncType
        integer :: n_arm                ! The number of arms
        integer :: poi_c                ! Center position

        integer :: n_un_scaf = 0        ! # of unpaired nucleotide in scaf
        integer :: n_un_stap = 0        ! # of unpaired nucleotide in stap

        double precision :: ref_ang     ! Reference angle between two neighboring edges
        double precision :: tot_ang     ! Total angle at the junction
        double precision :: gap         ! Gap distance between junction and end of edges

        integer, allocatable :: iniL(:)         ! Initial line
        integer, allocatable :: modP(:)         ! Modified point
        integer, allocatable :: croP(:,:)       ! Sectional point (# of arms, # of sections)
        integer, allocatable :: node(:,:)       ! Nodes (# of arms, # of sections)
        integer, allocatable :: conn(:,:)       ! Node connectivity (node #, node # to be connected))
        integer, allocatable :: type_conn(:)    ! Section connection type, negihbor(1) and self(2)
    end type JuncType

! -----------------------------------------------------------------------------

    ! Geometry data type to manage section, point, line and face data
    type :: GeomType
        integer :: n_sec                    ! The number of sections
        integer :: n_iniP, n_modP, n_croP   ! The number of initial, modified and sectional points
        integer :: n_iniL, n_croL           ! The number of initial, sectional lines
        integer :: n_face                   ! The number of faces
        integer :: min_edge_length
        integer :: max_edge_length
        integer :: n_junc                   ! The number of juncs

        type(SecType) :: sec                        ! Section
        type(PointType), allocatable :: iniP(:)     ! Initial points
        type(PointType), allocatable :: modP(:)     ! Modified points
        type(PointType), allocatable :: croP(:)     ! Cross-sectional points
        type(LineType),  allocatable :: iniL(:)     ! Initial lines
        type(LineType),  allocatable :: croL(:)     ! Cross-sectional lines
        type(FaceType),  allocatable :: face(:)     ! Face
        type(JuncType),  allocatable :: junc(:)
    end type GeomType

! -----------------------------------------------------------------------------

end module Data_Geom