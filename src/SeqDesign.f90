!
! =============================================================================
!
! Module - SeqDesign
! Last Updated : 03/13/2019, by Hyungmin Jun (hyungminjun@outlook.com)
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
module SeqDesign

    use Data_Prob
    use Data_Mesh
    use Data_DNA

    use Section

    use Para
    use Mani
    use List
    use Math

    implicit none

    public  SeqDesign_Design

    private SeqDesign_Build_dnaTop
    private SeqDesign_Get_Rand_Sequence
    private SeqDesign_Get_Comp_Sequence
    private SeqDesign_Build_Strand
    private SeqDesign_Make_Linear_Stap_Single_Xover
    private SeqDesign_Make_Linear_Stap_Nick
    private SeqDesign_Make_Linear_Scaf_Nick_Outside
    private SeqDesign_Make_Linear_Scaf_Nick
    private SeqDesign_Make_Short_Scaf
    private SeqDesign_Break_Staple_Max_Length

    private SeqDesign_Move_Nick
    
    private SeqDesign_Build_Sequence_Design_Mix
    private SeqDesign_Build_Sequence_Design
    private SeqDesign_Avoid_Barrier
    private SeqDesign_Make_Short_Strand
    private SeqDesign_Count_Remainder
    private SeqDesign_Rebuild_Strand
    private SeqDesign_Build_Region_Staple
    private SeqDesign_Build_Region_Staple_1
    private SeqDesign_Order_Staple
    private SeqDesign_Print_14nt_Region_Simple
    private SeqDesign_CirGraph_Count_Edge
    private SeqDesign_CirGraph_Init_Variable
    private SeqDesign_Assign_Sequence
    private SeqDesign_Set_M13mp18
    private SeqDesign_Get_M13mp18
    private SeqDesign_Get_Lamda
    private SeqDesign_Import_Sequence
    private SeqDesign_Set_Rand_Sequence

    private SeqDesign_Chimera_Atom
    private SeqDesign_Chimera_Route
    private SeqDesign_Chimera_Sequence_Design
    private SeqDesign_Chimera_Strand

    type :: RegionType
        integer :: types,    length                 ! 1-vertex, 2-edge, region length
        integer :: sta_pos,  cen_pos,  end_pos      ! Position
        integer :: sta_base, cen_base, end_base     ! Base ID
        integer :: n_14nt,   n_4nt                  ! # of 14nt seeds
    end type RegionType

    ! GraphType data strucutre is corresponding to node
    type GraphType
        ! node2base(i) : From base to index
        ! base2node(i) : From index to base
        ! node(i)      : 0-nomarl, 1-first 14nt, 2-secondary 14nt, 3-4nt
        ! edge(i,1)    : [1]->2 : source of connectivity
        ! edge(i,2)    : 1->[2] : target of connectivity
        ! edge(i,3)    : Weight factor
        integer, allocatable :: node2base(:)
        integer, allocatable :: base2node(:)
        integer, allocatable :: node(:)
        integer, allocatable :: edge(:,:)
    end type GraphType
contains

! -----------------------------------------------------------------------------

! Design topology
subroutine SeqDesign_Design(prob, geom, mesh, dna)
    type(ProbType),  intent(inout) :: prob
    type(GeomType),  intent(inout) :: geom
    type(MeshType),  intent(inout) :: mesh
    type(DNAType),   intent(inout) :: dna

    ! Build dnaTop data from dna base data
    call SeqDesign_Build_dnaTop(dna)

    ! Build strand data from dnaTop
    call SeqDesign_Build_Strand(dna)

    ! Make linear staples by the nick
    call SeqDesign_Make_Linear_Stap_Nick(mesh, dna)

    ! Staple break rule ranging from 20 to 60-bp
    call SeqDesign_Break_Staple_Max_Length(prob, mesh, dna)

    ! Make linear scaffold by the nick
    call SeqDesign_Make_Linear_Scaf_Nick_Outside(prob, geom, mesh, dna)

    ! Make short scaffold strand
    call SeqDesign_Make_Short_Scaf(mesh, dna)

    ! Rebuild strand data from dnaTop
    call SeqDesign_Rebuild_Strand(dna)

    ! List in long length order of the staple
    call SeqDesign_Order_Staple(dna)

    ! Print 14nt region with various representations
    if(para_max_break_scaf == 0) then
        call SeqDesign_Print_14nt_Region_Simple(prob, geom, mesh, dna)
    end if

    ! Assign DNA sequence according to para_set_seq_scaf
    call SeqDesign_Assign_Sequence(prob, dna)

    ! Write atom model by dnaTop and strand data
    call SeqDesign_Chimera_Atom(prob, dna)

    ! Chimera topology route
    call SeqDesign_Chimera_Route(prob, geom, mesh, dna)

    ! Chimera sequence design
    call SeqDesign_Chimera_Sequence_Design(prob, geom, mesh, dna)

    ! Write Chimera file for strand and sequence
    call SeqDesign_Chimera_Strand(prob, dna)
end subroutine SeqDesign_Design

! -----------------------------------------------------------------------------

! Reset possible staple crossovers
subroutine SeqDesign_Reset_Possible_Stap_Xover(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Crossover based on cross-sectional edges
    type :: CroLType
        integer :: max_bp, min_bp
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, k, croL_cur, croL_com, sec_cur, sec_com, id_bp
    integer :: up_scaf1, dn_scaf1, up_scaf2, dn_scaf2
    integer :: up_cur, up_com, dn_cur, dn_com
    logical :: b_nei_up, b_nei_dn, b_scaf

    ! Print progress
    write(p_redir, "(a)"), "   * Reset possible staple Xovers"

    ! Reset connectivity for bases in staple at crossovers
    dna.n_xover_stap  = 0
    dna.n_sxover_stap = 0
    do i = 1, mesh.n_node

        ! Reset connection at crossovers
        if(dna.base_stap(i).xover /= -1) then
            if(dna.base_scaf(i).up == -1) then
                dna.base_stap(i).up = dna.base_scaf(i).dn
                dna.base_stap(i).dn = -1
            else if(dna.base_scaf(i).dn == -1) then
                dna.base_stap(i).up = -1
                dna.base_stap(i).dn = dna.base_scaf(i).up
            else
                dna.base_stap(i).up = dna.base_scaf(i).dn
                dna.base_stap(i).dn = dna.base_scaf(i).up
            end if
        end if

        dna.base_stap(i).xover = -1
    end do

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    croL(1:geom.n_croL).max_bp = -999999
    croL(1:geom.n_croL).min_bp =  999999

    ! Find maximum and minimum basepair ID in cross-sectional edges
    do i = 1, mesh.n_node
        croL_cur = mesh.node(i).croL
        id_bp    = mesh.node(i).bp

        ! Set maximum and minimum base ID
        if(croL(croL_cur).max_bp < id_bp) croL(croL_cur).max_bp = id_bp
        if(croL(croL_cur).min_bp > id_bp) croL(croL_cur).min_bp = id_bp
    end do

    ! Find the possible staple double crossovers
    dna.n_xover_stap  = 0
    dna.n_sxover_stap = 0
    do i = 1, mesh.n_node       ! Loop for current node

        ! Print progress bar
        call Mani_Progress_Bar(i, mesh.n_node)

        ! Loop for comparing node
        do j = i + 1, mesh.n_node

            ! Exception for the pre-constructed crossovers (due to double crossover)
            if(dna.base_stap(i).xover /= -1 .and. dna.base_stap(j).xover /= -1) then
                cycle
            end if

            ! It should be skipped when condition below
            ! Basepair ID and iniL shoud be the same and croL and section ID should be different
            if(mesh.node(i).bp   /= mesh.node(j).bp  ) cycle
            if(mesh.node(i).iniL /= mesh.node(j).iniL) cycle
            if(mesh.node(i).croL == mesh.node(j).croL) cycle
            if(mesh.node(i).sec  == mesh.node(j).sec ) cycle

            ! Find section ID
            sec_cur  = mesh.node(i).sec
            sec_com  = mesh.node(j).sec
            croL_cur = mesh.node(i).croL
            croL_com = mesh.node(j).croL
            id_bp    = mesh.node(i).bp

            ! To eliminate boundary staple crossovers
            if(croL(croL_cur).min_bp + para_gap_xover_bound_stap > id_bp) cycle
            if(croL(croL_cur).max_bp - para_gap_xover_bound_stap < id_bp) cycle
            if(croL(croL_com).min_bp + para_gap_xover_bound_stap > id_bp) cycle
            if(croL(croL_com).max_bp - para_gap_xover_bound_stap < id_bp) cycle

            ! Determine whether the node has crossover or not
            if(Section_Connection_Stap(geom, sec_cur, sec_com, id_bp) == .true.) then

                ! To eliminate crossover if there is neighboring scaffold crossover
                !
                !     dn_scaf1   i   up_scaf1
                !    *---*---*---*---*---*---*-->  : node i
                !            |       |
                ! <--*---*---*---*---*---*---*     : node j
                !     up_scaf2   j   dn_scaf2
                b_scaf   = .false.
                up_scaf1 = mesh.node(i).id
                dn_scaf1 = mesh.node(i).id
                up_scaf2 = mesh.node(j).id
                dn_scaf2 = mesh.node(j).id

                ! Check neighbor scaffold crossovers
                do k = 1, para_gap_xover_two
                    if( (dna.base_scaf(mesh.node(up_scaf1).id).xover == dna.base_scaf(mesh.node(dn_scaf2).id).id   .and. &
                         dna.base_scaf(mesh.node(dn_scaf2).id).xover == dna.base_scaf(mesh.node(up_scaf1).id).id ) .or.  &
                        (dna.base_scaf(mesh.node(dn_scaf1).id).xover == dna.base_scaf(mesh.node(up_scaf2).id).id   .and. &
                         dna.base_scaf(mesh.node(up_scaf2).id).xover == dna.base_scaf(mesh.node(dn_scaf1).id).id ) ) then
                        b_scaf = .true.
                        exit
                    else
                        up_scaf1 = mesh.node(up_scaf1).up
                        dn_scaf1 = mesh.node(dn_scaf1).dn
                        up_scaf2 = mesh.node(up_scaf2).up
                        dn_scaf2 = mesh.node(dn_scaf2).dn
                    end if
                end do

                if(b_scaf == .true.) cycle

                ! Find upper or downward neighboring crossovers
                ! Node numbering is opposite to staple ID
                up_cur   = mesh.node(i).dn
                up_com   = mesh.node(j).up
                id_bp    = mesh.node(up_cur).bp
                b_nei_up = Section_Connection_Stap(geom, sec_cur, sec_com, id_bp)

                dn_cur   = mesh.node(i).up
                dn_com   = mesh.node(j).dn
                id_bp    = mesh.node(dn_cur).bp
                b_nei_dn = Section_Connection_Stap(geom, sec_cur, sec_com, id_bp)

                ! Set current and previous or next crossovers (double crossover)
                ! For current crossover
                dna.n_xover_stap       = dna.n_xover_stap + 2
                dna.base_stap(i).xover = dna.base_stap(j).id
                dna.base_stap(j).xover = dna.base_stap(i).id

                ! For neighboring crossover
                if(b_nei_up == .true.) then
                    dna.base_stap(up_cur).xover = dna.base_stap(up_com).id
                    dna.base_stap(up_com).xover = dna.base_stap(up_cur).id

                    ! Set connectivity
                    dna.base_stap(i).up      = dna.base_stap(j).id
                    dna.base_stap(j).dn      = dna.base_stap(i).id
                    dna.base_stap(up_cur).dn = dna.base_stap(up_com).id
                    dna.base_stap(up_com).up = dna.base_stap(up_cur).id

                else if(b_nei_dn == .true.) then
                    dna.base_stap(dn_cur).xover = dna.base_stap(dn_com).id
                    dna.base_stap(dn_com).xover = dna.base_stap(dn_cur).id

                    ! Set connectivity
                    dna.base_stap(i).dn      = dna.base_stap(j).id
                    dna.base_stap(j).up      = dna.base_stap(i).id
                    dna.base_stap(dn_cur).up = dna.base_stap(dn_com).id
                    dna.base_stap(dn_com).dn = dna.base_stap(dn_cur).id
                end if
            end if
        end do
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * The number of possible staple Xovers : "//trim(adjustl(Int2Str(dna.n_xover_stap)))
    write(P_redir, "(a)")

    ! Deallocate memory
    deallocate(croL)
end subroutine SeqDesign_Reset_Possible_Stap_Xover

! -----------------------------------------------------------------------------

! Build dnaTop data from dna base data
subroutine SeqDesign_Build_dnaTop(dna)
    type(DNAType), intent(inout) :: dna

    integer   :: i, j, n_jump, id_up, id_down
    character :: seq_across

    ! Set the number of dna top
    dna.n_top = dna.n_base_scaf + dna.n_base_stap

    ! Allocate dna top data
    allocate(dna.top(dna.n_top))

    ! Copy data from bases in scaffold strand
    do i = 1, dna.n_base_scaf

        ! ID and node ID
        dna.top(i).id   = dna.base_scaf(i).id
        dna.top(i).node = dna.base_scaf(i).node

        ! Connectivity
        dna.top(i).up = dna.base_scaf(i).up
        dna.top(i).dn = dna.base_scaf(i).dn

        ! Crossover and across ID
        dna.top(i).xover = dna.base_scaf(i).xover

        if(dna.base_scaf(i).across == -1) then
            dna.top(i).across = -1
        else
            dna.top(i).across = dna.base_scaf(i).across + dna.n_base_scaf
        end if

        ! Strand, residue ID and sequence
        dna.top(i).strand  = -1
        dna.top(i).address = -1
        dna.top(i).b_14nt  = .false.
        dna.top(i).seq     = "N"

        ! Assign position vector
        dna.top(i).pos(1:3) = dna.base_scaf(i).pos(1:3)
    end do

    ! Copy data from bases in staple strand
    n_jump = dna.n_base_scaf
    do i = 1, dna.n_base_stap

        ! ID and node ID
        dna.top(i+n_jump).id   = dna.base_stap(i).id + n_jump
        dna.top(i+n_jump).node = dna.base_stap(i).node

        ! Connectivity
        if(dna.base_stap(i).up == -1) then
            dna.top(i+n_jump).up = -1
        else
            dna.top(i+n_jump).up = dna.base_stap(i).up + n_jump
        end if

        if(dna.base_stap(i).dn == -1) then
            dna.top(i+n_jump).dn = -1
        else
            dna.top(i+n_jump).dn = dna.base_stap(i).dn + n_jump
        end if

        ! Crossover and across ID
        if(dna.base_stap(i).xover == -1) then
            dna.top(i+n_jump).xover = -1
        else
            dna.top(i+n_jump).xover = dna.base_stap(i).xover + n_jump
        end if
        dna.top(i+n_jump).across = dna.base_stap(i).across

        ! Strand, residue ID and sequence
        dna.top(i+n_jump).strand  = -1
        dna.top(i+n_jump).address = -1
        dna.top(i+n_jump).b_14nt  = .false.
        dna.top(i+n_jump).seq     = "N"

        ! Sequence for poly Tn loop
        if(dna.base_stap(i).across == -1) dna.top(i+n_jump).seq = "T"

        ! Assign position vector
        dna.top(i+n_jump).pos(1:3) = dna.base_stap(i).pos(1:3)
    end do

    ! Print progress
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " +===================================================+"
    write(p_redir, "(a)"), " | 6. Sequence design                                |"
    write(p_redir, "(a)"), " +===================================================+"
    write(p_redir, "(a)")
    write(p_redir, "(a)"), "  6.1. Build dnaTop data"
    write(p_redir, "(a)"), "   * total # of nucleotides   : "//trim(adjustl(Int2Str(dna.n_top)))
    write(p_redir, "(a)"), "   * # of nts in the scaffold : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
    write(p_redir, "(a)"), "   * # of nts in staples      : "//trim(adjustl(Int2Str(dna.n_base_stap)))
    write(p_redir, "(a)")
end subroutine SeqDesign_Build_dnaTop

! -----------------------------------------------------------------------------

! Get random sequence
function SeqDesign_Get_Rand_Sequence result(seq)
    character :: seq
    integer :: random
    real :: drandom

    ! Random number generation (range from 1 to 4)
    call random_number(drandom)

    random = int(4.0*drandom) + 1

    if(random == 1 .or. random == 5) then
        seq = "A"
    else if(random == 2) then
        seq = "T"
    else if(random == 3) then
        seq = "G"
    else if(random == 4) then
        seq = "C"
    end if
end function SeqDesign_Get_Rand_Sequence

! -----------------------------------------------------------------------------

! Get complementary sequence
function SeqDesign_Get_Comp_Sequence(seq) result(com_seq)
    character, intent(in) :: seq

    character :: com_seq

    if(seq == "A") then
        com_seq = "T"
    else if(seq == "T") then
        com_seq = "A"
    else if(seq == "C") then
        com_seq = "G"
    else if(seq == "G") then
        com_seq = "C"
    else
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | The sequences are not assigned.                  |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if
end function SeqDesign_Get_Comp_Sequence

! -----------------------------------------------------------------------------

! Build strand data from dnaTop
subroutine SeqDesign_Build_Strand(dna)
    type(DNAType), intent(inout) :: dna

    type(StrandType), allocatable :: strand(:)
    logical,          allocatable :: b_visit(:)
    type(ListBase),   pointer :: list_base, ptr_base
    type(TopType) :: cur_base

    integer, parameter :: max_strand = 2000
    integer :: i, j, n_strand, n_base, int_base
    logical :: b_end

    ! Allocate and initialize strand data
    allocate(strand(max_strand))
    call Mani_Init_StrandType(strand, max_strand)

    ! Nullify the linked lists
    nullify(list_base, ptr_base)

    ! Allocate and initilize the b_visit data
    allocate(b_visit(dna.n_top))
    b_visit(1:dna.n_top) = .false.

    n_strand = 0
    do
        ! Check b_visit if there is 0 (0 means not visiting base)
        b_end = .true.
        do i = 1, dna.n_top
            if(b_visit(dna.top(i).id) == .false.) then
                b_end = .false.
                exit
            end if
        end do
        if(b_end == .true.) exit

        ! Increase the number of strands
        n_strand = n_strand + 1
        cur_base = dna.top(i)
        int_base = cur_base.id

        if(n_strand > max_strand) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | Exceed maximum number of strands                 |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if

        ! Find the first base in the current strand
        do
            if(cur_base.dn == -1) then
                ! cur_base is at the 3'-end of the strand
                strand(n_strand).b_circular = .false.
                exit
            else if(cur_base.dn == int_base) then
                ! cur_base goes back to the starting point
                strand(n_strand).b_circular = .true.
                cur_base = dna.top(int_base)
                exit
            end if

            ! Update current base
            cur_base = dna.top(cur_base.dn)

            if(b_visit(cur_base.id) == .true.) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Reached a visited base                           |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if
        end do

        ! Walk through the current strand
        n_base = 1

        ! Insert data to linked list
        allocate(ptr_base)
        ptr_base%id = cur_base.id
        list_base => List_Insert_Base(list_base, ptr_base)

        dna.top(cur_base.id).strand  = n_strand
        dna.top(cur_base.id).address = n_base
        dna.top(cur_base.id).b_14nt  = .false.
        b_visit(cur_base.id)         = .true.

        ! Loop to add a new base into current strand
        do
            ! Check for going out loop
            if(strand(n_strand).b_circular == .false.) then
                ! If non-circular and upper ID is equal to -1
                if(cur_base.up == -1) exit
            else
                ! If circular and base ID is equal to init base ID
                if(cur_base.up == int_base) exit
            end if

            cur_base = dna.top(cur_base.up)

            if(b_visit(cur_base.id) == .true.) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Reached a visited base.                          |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            n_base = n_base + 1

            ! Insert data to linked list
            allocate(ptr_base)
            ptr_base%id = cur_base.id
            list_base => List_Insert_Base(list_base, ptr_base)

            dna.top(cur_base.id).strand  = n_strand
            dna.top(cur_base.id).address = n_base
            dna.top(cur_base.id).b_14nt  = .false.
            b_visit(cur_base.id)         = .true.
        end do

        strand(n_strand).n_base = n_base

        ! Allocate array using the linked list
        allocate(strand(n_strand).base(n_base))

        ! Put strand data from linked list
        ptr_base => list_base
        do i = 1, n_base
            strand(n_strand).base(n_base+1-i) = ptr_base%id
            ptr_base => ptr_base%next
        end do
    end do

    ! Deallocate b_visit data
    deallocate(b_visit)

    ! Set the number of strands and copy data
    dna.n_strand = n_strand
    dna.n_scaf   = 1
    dna.n_stap   = dna.n_strand - 1

    allocate(dna.strand(dna.n_strand))
    do i = 1, dna.n_strand

        dna.strand(i).n_base     = strand(i).n_base
        dna.strand(i).b_circular = strand(i).b_circular

        ! Copy data from base to dna.base
        allocate(dna.strand(i).base(dna.strand(i).n_base))
        do j = 1, dna.strand(i).n_base
            dna.strand(i).base(j) = strand(i).base(j)
        end do

        ! Set strand type
        dna.strand(i).type1 = "stap"
        if(i == 1) dna.strand(i).type1 = "scaf"

        ! Deallocate strand base data
        deallocate(strand(i).base)
    end do

    ! Deallocate strand 
    deallocate(strand)

    ! Delete linked list allocated
    call List_Delete_Base(list_base)
    !call List_Delete_Base(ptr_base)

    ! Print progress
    write(p_redir, "(a)"), "  6.2. Build strand data"
    write(p_redir, "(a)"), "   * Total # of strands            : "//trim(adjustl(Int2Str(dna.n_strand)))
    write(p_redir, "(a)"), "   * # of strands for the scaffold : "//trim(adjustl(Int2Str(1)))
    write(p_redir, "(a)"), "   * # of strands for staples      : "//trim(adjustl(Int2Str(dna.n_strand - 1)))
    write(p_redir, "(a)"), "   * Detailed info. on DNA strand data"

    ! Print progress in detail
    if(p_detail == .true.) then
        do i = 1, dna.n_strand
            write(p_redir, "(i10, a$)"), i, " strand: "
            if(dna.strand(i).type1 == "scaf") then
                write(p_redir, "(a$)"), "scaffold, "
            else if(dna.strand(i).type1 == "stap") then
                write(p_redir, "(a$)"), "staple, "
            end if
            write(p_redir, "(a, i6$)"), "# of bases : ", dna.strand(i).n_base
            write(p_redir, "(a, l  )"), ", circular : ", dna.strand(i).b_circular

            ! Print bases numbering
            write(p_redir, "(a15, i6, a$)"), "Bases:", dna.strand(i).base(1), "->"
            do j = 2, dna.strand(i).n_base - 1
                if(mod(j, 20)== 0) then
                    write(p_redir, "(i6, a)"), dna.strand(i).base(j), "->"
                else if(mod(j, 20)== 1) then
                    write(p_redir, "(i21, a$)"), dna.strand(i).base(j), "->"
                else
                    write(p_redir, "(i6, a$)"), dna.strand(i).base(j), "->"
                end if
            end do
            write(p_redir, "(i6)"), dna.strand(i).base(dna.strand(i).n_base)
            write(p_redir, "(a )")
        end do
    end if
    write(p_redir, "(a)")
end subroutine SeqDesign_Build_Strand

! -----------------------------------------------------------------------------

! Make linear staples by the single crossover
subroutine SeqDesign_Make_Linear_Stap_Single_Xover(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, base, init_base, across, node, dn_base

    ! Loop to make non-circular staple strand
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).type1 == "scaf") cycle

        ! Find non-circular strand
        if(dna.strand(i).b_circular == .true.) then

            ! Find first base that has double crossover
            base = dna.strand(i).base(1)
            do j = 1, dna.strand(i).n_base
                if( dna.top(base).xover /= -1 .and. &
                    dna.top(dna.top(base).dn).xover /= -1 ) then
                    exit
                end if
                base = dna.top(base).dn
            end do

            ! Find downward base
            dn_base = dna.top(base).dn

            ! Disconnect between current and downward bases
            dna.top(base).dn    = -1
            dna.top(dn_base).up = -1

            ! If this point is crossover
            if( dna.top(base).xover    == dna.top(dn_base).id .and. &
                dna.top(dn_base).xover == dna.top(base).id ) then

                ! Delete single crossover
                dna.n_xover_stap  = dna.n_xover_stap  - 1
                dna.n_sxover_stap = dna.n_sxover_stap + 1

                dna.top(base).xover    = -1
                dna.top(dn_base).xover = -1
            end if
        end if
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.3. Make linear staples by the single Xover"
    write(p_redir, "(a)"), "   * # of staples : "//trim(adjustl(Int2Str(dna.n_base_stap)))
    write(p_redir, "(a)")
end subroutine SeqDesign_Make_Linear_Stap_Single_Xover

! -----------------------------------------------------------------------------

! Make linear staples by the nick
subroutine SeqDesign_Make_Linear_Stap_Nick(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, base, init_base, across, node, dn_base
    integer :: count, max_count, start_base, max_start_base

    ! Loop for all staples
    do i = 1, dna.n_strand

        ! Skip when strand is the scaffold or staple is linear
        if(dna.strand(i).type1      == "scaf" ) cycle
        if(dna.strand(i).b_circular == .false.) cycle

        ! Find the position where the double crossover exists
        base = dna.strand(i).base(1)
        do j = 1, dna.strand(i).n_base
            if( dna.top(base).xover /= -1 .and. &
                dna.top(dna.top(base).dn).xover /= -1 ) then
                exit
            end if
            base = dna.top(base).dn
        end do

        ! To make the single crossover if there are no proper regions
        init_base = base

        ! Find maximum length region
        count      = 1
        max_count  = 0
        start_base = base
        do j = 1, dna.strand(i).n_base

            base   = dna.top(base).id
            across = dna.top(base).across

            if(across == -1) then

                ! For base in Tn-loop
                if(max_count < count) then
                    max_count      = count
                    max_start_base = start_base
                end if
                count      = 0
                start_base = base
            else

                ! Check scaffold and staple crossovers
                if(dna.top(base).xover /= -1 .or. dna.top(across).xover /= -1) then

                    if(max_count < count) then
                        max_count      = count
                        max_start_base = start_base
                    end if

                    count      = 0
                    start_base = base
                end if
            end if

            count = count + 1
            base  = dna.top(base).up
        end do

        ! Exclude starting and ending bases
        max_count = max_count - 1

        ! The minimum gap b/w xover/Tn and first nick, para_gap_min_region = 8
        if( max_count > para_gap_min_region ) then

            ! Set base at the center of the maximum length region
            base = dna.top(max_start_base).id
            do j = 1, max_count / 2 + 1
                base = dna.top(base).up
            end do
        else

            ! Make the single crossover if there are proper regions
            base = init_base
        end if

        dn_base = dna.top(base).dn

        ! Disconnect between current and downward bases
        dna.top(base).dn    = -1
        dna.top(dn_base).up = -1

        ! If this is the crossover position
        ! In case of "max_count > para_gap_min_region", there are bug in making single xover
        if( dna.top(base).xover    == dna.top(dn_base).id .and. &
            dna.top(dn_base).xover == dna.top(base).id ) then

            ! Delete single crossover
            dna.n_xover_stap  = dna.n_xover_stap  - 1
            dna.n_sxover_stap = dna.n_sxover_stap + 1

            dna.top(base).xover    = -1
            dna.top(dn_base).xover = -1
        end if
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.3. Make linear staples by nick"
    write(p_redir, "(a)"), "   * # of staples : "//trim(adjustl(Int2Str(dna.n_base_stap)))
    write(p_redir, "(a)")
end subroutine SeqDesign_Make_Linear_Stap_Nick

! -----------------------------------------------------------------------------

! Build sequence design combined maximum and optimal cutting
subroutine SeqDesign_Build_Sequence_Design_Mix(prob, mesh, dna)
    type(ProbType), intent(inout) :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, cen_base, pre_base, up_base, pre_reg
    integer :: n_region, cn_tn, length, final_length, cen_pos, pre_pos, bgn_pos
    integer :: cng_para, up, xover, upxover
    integer :: cng_gap_nick
    logical :: b_ext, b_cut, b_14nt

    cng_gap_nick = 8

    ! Loop to make short staple based on mixed method
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).type1 == "scaf") cycle

        ! For short staple, para_max_break_stap <= 60
        if(dna.strand(i).n_base <= para_max_break_stap) cycle

        ! Count unpaired staple nucleotides if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! For vertex in case of DX tile design
        !if(cn_tn > 0 .and. prob.sel_edge_sec == 1 .and. dna.strand(i).n_base < 80) cycle

        ! --------------------------------------------------
        ! Build region of staple strands
        ! --------------------------------------------------
        ! Tn loop |<----->| : 7 nt poly T loop
        !     ~~~~~=======*===========*=====*~~~~~~~
        !                 |<--------->|<--->|
        !           11 and 5 nt region not including crossover and nick
        !
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Print information on region of staple strands
        write(p_redir, "(i6    )"), i
        write(p_redir, "(a, i5$)"), "-strd, # nts:", dna.strand(i).n_base
        write(p_redir, "(a, i5 )"), ", # nts@Tn:", cn_tn

        do j = 1, n_region
            write(p_redir, "(i11, a$)"), j, "-rgn"
            write(p_redir, "(a,  i3$)"), ", T:",        region(j).types
            write(p_redir, "(a,  i4$)"), ", regn len:", region(j).length
            write(p_redir, "(a,  i4$)"), ", str P:",    region(j).sta_pos
            write(p_redir, "(a,  i4$)"), ", cen P:",    region(j).cen_pos
            write(p_redir, "(a,  i4 )"), ", end P:",    region(j).end_pos
        end do
        write(p_redir, "(a)")
        write(p_redir, "(a)"), "      ----- Make nick -----"

        ! --------------------------------------------------
        ! Cut staples to make short multi staples with 14nt seeds
        ! --------------------------------------------------
        bgn_pos = 0
        pre_reg = 0
        b_cut   = .false.
        b_14nt  = .false.
        j       = 0

        do
            ! Check loop
            j = j + 1
            if(j == n_region + 1) exit

            ! Skip when the region length is smaller than 5
            if(j /= n_region .and. region(j).length < 5) cycle

            ! Cutted length from center to beginning
            length = region(j).cen_pos - bgn_pos
            if(j == n_region) length = dna.strand(i).n_base - bgn_pos + 1

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_break_stap) then
                b_cut = .true.
            end if

            ! Skip this current loop to prevent cutting this region
            if(j /= n_region .and. length <= para_max_break_stap) then
                if( (b_14nt == .false. .and. region(j).length + 1 >= 14 .and. region(j).types == 1) .or. &
                    (b_14nt == .false. .and. region(j).length + 2 >= 14 .and. region(j).types == 2) ) then
                    b_14nt = .true.
                    cycle
                end if
            end if

            if(b_cut == .true. .and. b_14nt == .true. .and. length <= para_max_break_stap) then

                ! --------------------------------------------------
                ! Optimal cutting to have 14nt seeds
                ! --------------------------------------------------
                ! Centered base ID and position
                cen_base = region(j).cen_base
                cen_pos  = region(j).cen_pos

                ! To avoid short region cutting
                if(region(j).length < cng_gap_nick) cycle

                ! If the final length is smaller than minimum length
                if(dna.strand(i).n_base - cen_pos < para_min_break_stap) cycle

                ! Print progress
                write(p_redir, "(i11, a$)"), j, "->cut rgn"
                write(p_redir, "(a,  i4$)"), ", len:",     cen_pos - bgn_pos
                write(p_redir, "(a,  i4$)"), ", cut pos:", cen_pos
                write(p_redir, "(a,  i4$)"), ", rmn len:", dna.strand(i).n_base - cen_pos
                write(p_redir, "(a      )"), "->14nt cut"

                ! Increase the number of staples
                dna.n_stap = dna.n_stap + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(cen_base).up
                dna.top(cen_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update beginning position and flags
                bgn_pos = cen_pos
                b_cut   = .false.
                b_14nt  = .false.
            else if(b_cut == .true. .and. length > para_max_break_stap) then

                ! --------------------------------------------------
                ! If the length exceeds the maximum cutting length
                ! --------------------------------------------------
                ! To avoid small region in previous region
                jj       = j
                b_ext    = .false.
                cng_para = para_gap_xover_nick

                ! Find previous region
                do
                    ! If the region length is longer than paramter value (default : 8)
                    if( (region(jj).types == 1 .and. region(jj).length >= cng_para * 2 + 2 + 1) .or. &
                        (region(jj).types == 2 .and. region(jj).length >= cng_para * 2 + 2) ) then
                        pre_base = region(jj).cen_base
                        pre_pos  = region(jj).cen_pos

                        ! Check cutted and remained length
                        if( pre_pos - bgn_pos >= para_min_break_stap .and. &
                            pre_pos - bgn_pos <= para_max_break_stap .and. &
                            dna.strand(i).n_base - pre_pos >= para_min_break_stap ) exit
                    end if

                    ! Go back previous region
                    jj = jj - 1

                    ! Exception that the region cutted already
                    if(jj == pre_reg .or. jj == 0) then

                        ! Reset parameter and region index
                        jj       = j
                        cng_para = cng_para - 1

                        if(cng_para == 0) then
                            b_ext    = .true.
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos
                            exit
                        end if

                        if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                        if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                        if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                        write(p_redir, "(a)"), "   ** Adjusting parameter, para_gap_xover_nick : From "//&
                            trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                            trim(adjustl(Int2Str(cng_para)))//" - "//&
                            trim(adjustl(Int2Str(cng_para * 2 + 2)))
                    end if
                end do

                ! No cutting due to exception, which may exceeds 60
                if(b_ext == .true.) then
                    write(p_redir, "(a$)"), "   |-last region with exception, remaining length : "
                    write(p_redir, "(i4)"), dna.strand(i).n_base - bgn_pos
                    if(j == n_region) cycle
                end if

                ! If the cutted length is smaller than para_min_break_stap
                if(pre_pos - bgn_pos < para_min_break_stap) then
                    cycle
                end if

                ! Check 14nt seed cutting
                up      = dna.top(region(jj).end_base).up
                xover   = dna.top(up).xover
                upxover = dna.top(dna.top(up).up).xover

                if( region(jj).end_pos + 1 - bgn_pos >= para_min_break_stap .and. &
                    region(jj).end_pos + 1 - bgn_pos <= para_max_break_stap .and. &
                    dna.strand(i).n_base - region(jj).end_pos + 1 >= para_min_break_stap .and. &
                    region(jj).length >= 12 .and. xover /= -1 .and. upxover /= -1 .and. para_set_stap_sxover == "on") then

                    ! Print progress
                    write(p_redir, "(i7,  a$)"), jj, "->cut rgn"
                    write(p_redir, "(a, i4$ )"), ", len:",     region(jj).end_pos + 1 - bgn_pos
                    write(p_redir, "(a, i4$ )"), ", cut pos:", region(jj).end_pos + 1
                    write(p_redir, "(a, i4$ )"), ", rmn len:", dna.strand(i).n_base - region(jj).end_pos + 1
                    write(p_redir, "(a      )"), "->max cut-single Xover"

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    if(dna.top(up).up == xover) dna.top(up).up = -1
                    if(dna.top(up).dn == xover) dna.top(up).dn = -1

                    if(dna.top(xover).up == up) dna.top(xover).up = -1
                    if(dna.top(xover).dn == up) dna.top(xover).dn = -1

                    dna.top(up).xover    = -1
                    dna.top(xover).xover = -1

                    dna.n_xover_stap  = dna.n_xover_stap  - 1
                    dna.n_sxover_stap = dna.n_sxover_stap + 1

                    ! Update starting position and flag
                    pre_reg = jj
                    bgn_pos = region(jj).end_pos + 1
                    b_cut   = .false.
                    b_14nt  = .false.
                    j       = jj
                else

                    ! Print progress
                    write(p_redir, "(i7,  a$)"), jj, "->cut rgn"
                    write(p_redir, "(a,  i4$)"), ", len:",     pre_pos - bgn_pos
                    write(p_redir, "(a,  i4$)"), ", cut pos:", pre_pos
                    write(p_redir, "(a,  i4$)"), ", rmn len:", dna.strand(i).n_base - pre_pos
                    write(p_redir, "(a      )"), "->max cut-nick"

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    ! Update starting position and flag
                    pre_reg = jj
                    bgn_pos = pre_pos
                    b_cut   = .false.
                    b_14nt  = .false.
                    j       = jj
                end if
            end if

            ! If remained length is smaller than para_max_break_stap, exit this loop
            if(dna.strand(i).n_base - bgn_pos < para_max_break_stap) then
                write(p_redir, "(a$ )"), "   |-->last regn, len:"
                write(p_redir, "(i4$)"), dna.strand(i).n_base - bgn_pos
                write(p_redir, "(a$ )"), ", cut pos: ++"
                write(p_redir, "(a$ )"), ", rmn len: ++"
                write(p_redir, "(a  )"), "->14nt cut"
                exit
            end if

            ! Check final staple and print information
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - bgn_pos

                ! If the last staple exceeds maximum length
                if(final_length > para_max_break_stap) then

                    jj       = j
                    cng_para = para_gap_xover_nick

                    ! Find previous region
                    do
                        ! If the region length is longer than paramter value (default : 8)
                        if(region(jj).length >= cng_para * 2 + 2) then
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos

                            ! Check cutted and remained length
                            if( pre_pos - bgn_pos >= para_min_break_stap .and. &
                                pre_pos - bgn_pos <= para_max_break_stap .and. &
                                dna.strand(i).n_base - pre_pos >= para_min_break_stap ) exit
                        end if

                        ! Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_reg .or. jj == 0) then

                            jj       = j
                            cng_para = cng_para - 1

                            if(cng_para == 0) then
                                deallocate(region)
                                return
                            end if
                        end if
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(p_redir, "(a)"), "   WARNING : This strand length exceeds para_max_break_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_break_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_break_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_break_stap) ) then
                        do k = 0, 11, 11
                            write(p_redir, "(i7, a, i4)"), j, "-final rgn3, len:", &
                                dna.strand(i).n_base - bgn_pos
                        end do
                        cycle
                    end if

                    ! Add # of staple
                    up_base    = dna.top(pre_base).up
                    dna.n_stap = dna.n_stap + 1
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    write(p_redir, "(i7, a, i4)"), jj, "-final rgn4, len:", &
                        final_length - (dna.strand(i).n_base - pre_pos)
                    write(p_redir, "(i7, a, i4)"), jj, "-final rgn5, len:", &
                        dna.strand(i).n_base - pre_pos
                else
                    write(p_redir, "(i7, a, i4)"), j, "-final rgn6, len:", &
                        dna.strand(i).n_base - bgn_pos
                end if
            end if
        end do
        write(p_redir, "(a)")

        ! Deallocate memory
        deallocate(region)
    end do
end subroutine SeqDesign_Build_Sequence_Design_Mix

! -----------------------------------------------------------------------------

! Staple break rule, Maximized staple length
subroutine SeqDesign_Break_Staple_Max_Length(prob, mesh, dna)
    type(ProbType), intent(inout) :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, cen_base, pre_base, up_base, pre_region, cng_para
    integer :: n_region, cn_tn, length, final_length, cen_pos, pre_pos, bgn_pos, base1, base2
    logical :: b_cut, b_ext

    ! Loop for all strands
    do i = 1, dna.n_strand

        ! Only staples
        if(dna.strand(i).type1 == "scaf") cycle

        ! For short staple, para_max_break_stap <= 60
        if(dna.strand(i).n_base <= para_max_break_stap) cycle

        ! Count unpaired nucleotides for the current staple
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! --------------------------------------------------
        ! Build the region of the current staple strand
        ! --------------------------------------------------
        ! Tn loop |<----->| : 7 nt poly T loop
        !     ~~~~~=======*===========*=====*~~~~~~~
        !                 |<--------->|<--->|
        !           11 and 5 nt region not including crossover and nick
        !
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Print information on the current staple region
        write(p_redir, "(i6$   )"), i
        write(p_redir, "(a, i5$)"), "-strd, # nts:", dna.strand(i).n_base
        write(p_redir, "(a, i5 )"), ", # nts@Tn:",  cn_tn
        do j = 1, n_region
            write(p_redir, "(i8, a$ )"), j, "-rgn"
            write(p_redir, "(a,  i2$)"), ", T:",     region(j).types
            write(p_redir, "(a,  i4$)"), ", len:",   region(j).length
            write(p_redir, "(a,  i4$)"), ", str P:", region(j).sta_pos
            write(p_redir, "(a,  i4$)"), ", cen P:", region(j).cen_pos
            write(p_redir, "(a,  i4 )"), ", end P:", region(j).end_pos
        end do
        write(p_redir, "(a)")
        write(p_redir, "(a)"), "      ----- Make nick -----"

        ! break the staple strand to make multiple compartments
        bgn_pos    = 0
        pre_region = 0
        b_cut      = .false.

        do j = 1, n_region

            ! Find the centered base and position
            cen_base = region(j).cen_base
            cen_pos  = region(j).cen_pos
            length   = cen_pos - bgn_pos

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_break_stap) b_cut = .true.

            ! Staple cutting with maximum length
            if(b_cut == .true. .and. length >= para_max_break_stap) then

                ! To avoid small region that will be 4nt region
                jj       = j - 1
                b_ext    = .false.
                cng_para = para_gap_xover_nick

                ! Find previous large enough region
                do
                    ! If the region length exceed, para_gap_xover_nick = 3 (default : 8)
                    if( (region(jj).types == 1 .and. region(jj).length >= cng_para * 2 + 2 + 1) .or. &
                        (region(jj).types == 2 .and. region(jj).length >= cng_para * 2 + 2) ) then
                        pre_base = region(jj).cen_base
                        pre_pos  = region(jj).cen_pos

                        ! Check remained staple length
                        if( pre_pos - bgn_pos >= para_min_break_stap .and. &
                            pre_pos - bgn_pos <= para_max_break_stap .and. &
                            dna.strand(i).n_base - pre_pos >= para_min_break_stap) exit
                    end if

                    !Go back previous region
                    jj = jj - 1

                    ! Exception that the region cutted already
                    if(jj == pre_region .or. jj == 0) then

                        ! Reset parameter and region index
                        jj       = j - 1
                        cng_para = cng_para - 1
                        if(cng_para == 0) then
                            b_ext    = .true.
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos
                            exit
                        end if

                        ! Print information on changing parameter
                        if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                        if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                        if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                        write(p_redir, "(a)"), "   ** Adjusting parameter, para_gap_xover_nick : From "//&
                            trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                            trim(adjustl(Int2Str(cng_para)))//" - "//&
                            trim(adjustl(Int2Str(cng_para * 2 + 2)))
                    end if
                end do

                ! No cutting due to exception, which might exceed 60 edge length
                if(b_ext == .true.) then
                    write(p_redir, "(a$)"), "   |-last region with exception, remaining length : "
                    write(p_redir, "(i4)"), dna.strand(i).n_base - bgn_pos
                    if(j == n_region) cycle
                end if

                ! Print progress
                write(p_redir, "(i7, a$)"), jj, "->cut rgn"
                write(p_redir, "(a, i4$)"), ", len:",     pre_pos - bgn_pos
                write(p_redir, "(a, i4$)"), ", cut pos:", pre_pos
                write(p_redir, "(a, i4$)"), ", rmn len:", dna.strand(i).n_base - pre_pos
                write(p_redir, "(a     )"), "->max cut-nick"

                ! 1 : with single crossover
                if(0 .and. region(jj).types == 1 .and. region(jj).length <= 9) then

                    ! Make single crossover
                    if(dna.top(dna.top(region(jj).end_base).up).across == -1) then
                        ! Back crossover
                        base1      = dna.top(region(jj).sta_base).dn
                        base2      = dna.top(base1).dn
                        pre_pos    = region(jj-1).end_pos + 1
                        pre_region = jj - 1
                    else
                        ! Front crossover
                        base1      = dna.top(region(jj).end_base).up
                        base2      = dna.top(base1).up
                        pre_pos    = region(jj).end_pos + 1
                        pre_region = jj
                    end if

                    ! Make single crossover
                    write(p_redir, "(a$     )"), "***** Single Xover"
                    write(p_redir, "(a,  i4$)"), ", len:",     pre_pos - bgn_pos
                    write(p_redir, "(a,  i4$)"), ", cut pos:", pre_pos
                    write(p_redir, "(a,  i4 )"), ", rmn len:", dna.strand(i).n_base - pre_pos

                    ! Cut crossover and make new connectivity
                    dna.top(base1).xover = -1
                    dna.top(base2).xover = -1
                    dna.top(base1).dn    = -1
                    dna.top(base2).up    = -1
                    dna.n_stap           = dna.n_stap + 1
                    dna.n_xover_stap     = dna.n_xover_stap - 1
                    dna.n_sxover_stap    = dna.n_sxover_stap + 1

                    ! Update flag
                    bgn_pos = pre_pos
                    b_cut   = .false.
                else

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    ! Update starting position and flag
                    pre_region = jj
                    bgn_pos    = pre_pos
                    b_cut      = .false.
                end if
            end if

            ! Check remained length whether it is smaller than para_max_break_stap
            if(dna.strand(i).n_base - bgn_pos <= para_max_break_stap) then
                write(p_redir, "(i7, a, i4)"), j, "-final rgn2, len:", &
                    dna.strand(i).n_base - bgn_pos
                exit
            end if

            ! Check if it is the final region
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - bgn_pos

                ! If the last staple exceeds maximum length
                if(final_length >= para_max_break_stap) then

                    jj       = j
                    cng_para = para_gap_xover_nick
                    do
                        ! Exception
                        if(jj == 0) then
                            deallocate(region)
                            return
                        end if

                        if(region(jj).length >= para_gap_xover_nick * 2 + 2) then
                            pre_base = region(jj).cen_base
                            pre_pos  = region(jj).cen_pos

                            ! Check remained staple length
                            if( pre_pos - bgn_pos >= para_min_break_stap .and. &
                                pre_pos - bgn_pos <= para_max_break_stap .and. &
                                dna.strand(i).n_base - pre_pos >= para_min_break_stap) exit
                        end if

                        !Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_region .or. jj == 0) then

                            ! Reset parameter and region index
                            jj       = j - 1
                            cng_para = cng_para - 1
                            if(cng_para == 0) then
                                b_ext    = .true.
                                pre_base = region(jj).cen_base
                                pre_pos  = region(jj).cen_pos
                                exit
                            end if

                            ! Print information on changing parameter
                            if(para_gap_xover_nick - cng_para == 1 .and. prob.n_cng_max_stap == 0) prob.n_cng_max_stap = 1
                            if(para_gap_xover_nick - cng_para == 2 .and. prob.n_cng_max_stap == 1) prob.n_cng_max_stap = 2
                            if(para_gap_xover_nick - cng_para == 3 .and. prob.n_cng_max_stap == 2) prob.n_cng_max_stap = 3

                            write(p_redir, "(a)"), "   ** Adjusting parameter, para_gap_xover_nick : From "//&
                                trim(adjustl(Int2Str(cng_para + 1)))//" to "//&
                                trim(adjustl(Int2Str(cng_para)))//" - "//&
                                trim(adjustl(Int2Str(cng_para * 2 + 2)))
                        end if
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(p_redir, "(a)"), "   WARNING : This strand length exceeds para_max_break_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_break_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_break_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_break_stap) ) then
                        write(p_redir, "(i7, a, i5)"), j, "-final rgn3, len:", &
                            dna.strand(i).n_base - bgn_pos
                        cycle
                    end if

                    ! Increase the number of staples
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    write(p_redir, "(i7, a, i5)"), jj, "-final rgn4, len:", &
                        final_length - (dna.strand(i).n_base - pre_pos)
                    write(p_redir, "(i7, a, i5)"), jj, "-final rgn5, len:", &
                        dna.strand(i).n_base - pre_pos
                else
                    write(p_redir, "(i7, a, i5)"), j, "-final rgn6, len:", &
                        dna.strand(i).n_base - bgn_pos
                end if
            end if
        end do
        write(p_redir, "(a)")

        deallocate(region)
    end do
end subroutine SeqDesign_Break_Staple_Max_Length

! -----------------------------------------------------------------------------

! Move the nick
subroutine SeqDesign_Move_Nick(dna, ref_up, ref_dn, dir, step)
    type(DNAType), intent(inout) :: dna
    integer,       intent(in)    :: ref_up
    integer,       intent(in)    :: ref_dn
    character(2),  intent(in)    :: dir
    integer,       intent(in)    :: step

    integer :: i, cur, up, dn

    up = ref_up
    dn = ref_dn

    ! Reconnect current staple
    dna.top(dn).dn = up
    dna.top(up).up = dn

    ! Move current staple
    if(trim(dir) == "up") then

        cur = up
        do i = 1, step
        cur = dna.top(cur).up
        end do
        up = cur
        dn = dna.top(cur).up
    else if(trim(dir) == "dn") then

        cur = dn
        do i = 1, step
        cur = dna.top(cur).dn
        end do
        dn = cur
        up = dna.top(cur).dn
    end if

    dna.top(dn).dn = -1
    dna.top(up).up = -1
end subroutine SeqDesign_Move_Nick

! -----------------------------------------------------------------------------

! Build staple region without crossovers and unpaired nucleotides
subroutine SeqDesign_Build_Region_Staple(dna, i, region, n_region)
    type(DNAType),    intent(inout) :: dna
    type(RegionType), intent(inout) :: region(:)
    integer,          intent(in)    :: i
    integer,          intent(inout) :: n_region

    integer :: j, k, base, across, sta_base, cen_base, end_base, cen_pos, len_cen
    logical :: b_region, b_vertex

    ! Find the starting point (down == -1)
    base = Mani_Go_Start_Base(dna, i)

    n_region = 0
    b_region = .false.
    b_vertex = .false.

    do j = 1, dna.strand(i).n_base

        ! Find across ID
        across = dna.top(base).across

        if(across == -1) then

            ! For Tn loop
            b_region = .true.
            b_vertex = .true.

            ! Update base to go up
            base = dna.top(base).up
            cycle
        else
            ! Check dn, up, xover and across's over
            ! Starting base always has -1 (down = -1)
            if( dna.top(base).dn      == -1 .or. &
                dna.top(base).up      == -1 .or. &
                dna.top(base).xover   /= -1 .or. &
                dna.top(across).xover /= -1 ) then

            b_region = .true.

            ! Update base to go up
            base = dna.top(base).up
            cycle
            end if
        end if

        if(b_region == .true.) then

            ! Make new region
            b_region = .false.
            n_region = n_region + 1

            ! Set region data
            region(n_region).sta_base = base
            region(n_region).length   = 1
            region(n_region).sta_pos  = j
            region(n_region).end_pos  = j

            ! Set vertex region
            if(b_vertex == .true.) then
                region(n_region).types   = 1
                region(n_region-1).types = 1
                b_vertex = .false.
            else
                region(n_region).types = 2
            end if
        else

            ! Increase length
            region(n_region).length  = region(n_region).length + 1
            region(n_region).end_pos = j
        end if

        ! Update base to go up
        base = dna.top(base).up
    end do

    ! Set center/end base and position
    do j = 1, n_region

        len_cen = (region(j).length+1)/2 - 1

        ! If the vetrex region
        if(region(j).types == 1) then
            if(dna.top(dna.top(region(j).sta_base).dn).across /= -1) then
                len_cen = len_cen - 1
            else
                len_cen = len_cen + 1
            end if
        end if

        ! Find center base
        cen_base = region(j).sta_base
        do k = 1, len_cen
            cen_base = dna.top(cen_base).up
        end do
        region(j).cen_base = cen_base

        ! Set center position
        region(j).cen_pos = region(j).sta_pos + len_cen

        ! Find end base
        end_base = region(j).sta_base
        do k = 1, region(j).length - 1
            end_base = dna.top(end_base).up
        end do
        region(j).end_base = end_base

        ! Set end position
        region(j).end_pos = region(j).sta_pos + (region(j).length - 1)
    end do
end subroutine SeqDesign_Build_Region_Staple

! -----------------------------------------------------------------------------

! Build staple region without crossovers and unpaired nucleotides
subroutine SeqDesign_Build_Region_Staple_1(dna, i, region, n_region)
    type(DNAType),    intent(inout) :: dna
    type(RegionType), intent(inout) :: region(:)
    integer,          intent(in)    :: i
    integer,          intent(inout) :: n_region

    integer :: j, k, base, across, sta_base, cen_base, end_base, cen_pos, len_cen
    logical :: b_region, b_vertex

    ! Find the starting point (down == -1)
    base = Mani_Go_Start_Base(dna, i)

    n_region = 0
    b_region = .false.
    b_vertex = .false.

    do j = 1, dna.strand(i).n_base

        ! Find across ID
        across = dna.top(base).across

        if(across == -1) then

            ! For Tn loop
            b_region = .true.
            b_vertex = .true.

            ! Update base to go up
            base = dna.top(base).up
            cycle
        else
            ! Check dn, up, xover and across's over
            ! Starting base always has -1 (down = -1)
            if( dna.top(base).dn      == -1 .or. &
                dna.top(base).up      == -1 .or. &
                dna.top(base).xover   /= -1 .or. &
                dna.top(across).xover /= -1 ) then

            b_region = .true.

            ! Update base to go up
            base = dna.top(base).up
            cycle
            end if
        end if

        if(b_region == .true.) then

            ! Make new region
            b_region = .false.
            n_region = n_region + 1

            ! Set region data
            region(n_region).sta_base = base
            region(n_region).length   = 1
            region(n_region).sta_pos  = j
            region(n_region).end_pos  = j

            ! Set vertex region
            if(b_vertex == .true.) then
                region(n_region).types   = 1
                region(n_region-1).types = 1
                b_vertex = .false.
            else
                region(n_region).types = 2
            end if
        else

            ! Increase length
            region(n_region).length  = region(n_region).length + 1
            region(n_region).end_pos = j
        end if

        ! Update base to go up
        base = dna.top(base).up
    end do

    ! Set center/end base and position
    do j = 1, n_region

        len_cen = (region(j).length+1)/2 - 1

        ! If the vetrex region
        if(region(j).types == 1) then
            if(dna.top(dna.top(region(j).sta_base).dn).across /= -1) then
                len_cen = len_cen - 1
            else
                len_cen = len_cen + 1
            end if
        end if

        ! Find center base
        cen_base = region(j).sta_base
        do k = 1, len_cen
            cen_base = dna.top(cen_base).up
        end do
        region(j).cen_base = cen_base

        ! Set center position
        region(j).cen_pos = region(j).sta_pos + len_cen

        ! Find end base
        end_base = region(j).sta_base
        do k = 1, region(j).length - 1
            end_base = dna.top(end_base).up
        end do
        region(j).end_base = end_base

        ! Set end position
        region(j).end_pos = region(j).sta_pos + (region(j).length - 1)

        ! --------------------------------------------------
        ! Set nucleotide of 14nt seeds
        ! --------------------------------------------------
        if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
            (region(j).types == 2 .and. region(j).length + 2 >= 14)) then

        sta_base = region(j).sta_base
        dna.top(dna.top(sta_base).dn).b_14nt = .true.
        dna.top(sta_base).b_14nt             = .true.
        dna.top(dna.top(sta_base).dn).status = "S"
        dna.top(sta_base).status             = "S"

        do k = 1, region(j).length - 1
            sta_base = dna.top(sta_base).up
            dna.top(sta_base).b_14nt = .true.
            dna.top(sta_base).status = "S"
        end do

        if(dna.top(dna.top(sta_base).up).across /= -1) then
            dna.top(dna.top(sta_base).up).b_14nt = .true.
            dna.top(dna.top(sta_base).up).status = "S"
        end if
            end if

        ! --------------------------------------------------
        ! Set nucleotide of 4nt dsDNA domain
        ! --------------------------------------------------
        if( (region(j).types == 1 .and. region(j).length + 1 <= 4) .or. &
            (region(j).types == 2 .and. region(j).length + 2 <= 4)) then

        sta_base = region(j).sta_base
        dna.top(dna.top(sta_base).dn).status = "F"
        dna.top(sta_base).status             = "F"

        do k = 1, region(j).length - 1
            sta_base = dna.top(sta_base).up
            dna.top(sta_base).status = "F"
        end do

        if(dna.top(dna.top(sta_base).up).across /= -1) then
            dna.top(dna.top(sta_base).up).status = "F"
        end if
        end if
    end do
end subroutine SeqDesign_Build_Region_Staple_1

! -----------------------------------------------------------------------------

! Build sequence design with non-circular staple strands
subroutine SeqDesign_Build_Sequence_Design(prob, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type :: RegionType
        integer :: types        ! 1 - vertex, 2 - edge
        integer :: length       ! total length at the region not including boundary

        integer :: strt_pos     ! start position
        integer :: end_pos      ! end position
        integer :: cntr_pos     ! center position

        integer :: strt_base    ! start base
        integer :: cntr_base    ! center base
    end type RegionType

    type(RegionType), allocatable :: region(:)

    integer :: i, j, jj, k, base, across, cntr_base, pre_base, up_base, pre_region
    integer :: n_region, cn_tn, length, final_length, cntr_pos, pre_pos, begin_pos
    logical :: b_region, b_cut, b_ext, b_vertex

    ! Loop to make short staple strand
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).type1 == "scaf") cycle

        ! For short staple (less than para_max_break_stap, 60)
        if(dna.strand(i).n_base < para_max_break_stap) cycle

        ! Count bases if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! --------------------------------------------------
        ! For vertex for DX tile
        ! --------------------------------------------------
        !if(cn_tn > 0 .and. prob.sel_edge_sec == 1 .and. dna.strand(i).n_base < 80) cycle

        ! Find the starting point (down == -1)
        base = Mani_Go_Start_Base(dna, i)

        ! --------------------------------------------------
        ! Build region array
        ! --------------------------------------------------
        ! Tn loop |<---------->| 12 bp
        !     ~~~~~============*===========*=====*  : base pair
        !                      |<--------->|
        !               region 11 not including crossover and nick region
        n_region = 0
        allocate(region(dna.strand(i).n_base))

        do j = 1, dna.strand(i).n_base

            ! Find across ID
            across = dna.top(base).across

            if(across == -1) then

                ! For Tn loop
                b_region = .true.
                b_vertex = .true.

                ! Update base to go up
                base = dna.top(base).up
                cycle
            else
                ! Check down, up, xover and across
                ! Starting base always has -1 (down = -1)
                if( dna.top(base).dn      == -1 .or. &
                    dna.top(base).up      == -1 .or. &
                    dna.top(base).xover   /= -1 .or. &
                    dna.top(across).xover /= -1 ) then

                    b_region = .true.

                    ! Update base to go up
                    base = dna.top(base).up
                    cycle
                end if
            end if

            if(b_region == .true.) then

                ! Make new region
                b_region = .false.
                n_region = n_region + 1

                ! Set region data
                region(n_region).strt_base = base
                region(n_region).length    = 1
                region(n_region).strt_pos  = j
                region(n_region).end_pos   = j
                region(n_region).types     = 2

                ! Set vertex region
                if(b_vertex == .true.) then
                    region(n_region).types   = 1
                    region(n_region-1).types = 1
                    b_vertex = .false.
                end if
            else

                ! Increase length
                region(n_region).length  = region(n_region).length + 1
                region(n_region).end_pos = j

                ! Set vertex region
                if(b_vertex == .true.) then
                    region(n_region).types   = 1
                    region(n_region-1).types = 1
                    b_vertex = .false.
                end if
            end if

            ! Update base to go up
            base = dna.top(base).up
        end do

        ! Set centered base and position
        do j = 1, n_region

            ! Find centered base
            cntr_base = region(j).strt_base
            do k = 1, (region(j).length + 1) / 2 - 1
                cntr_base = dna.top(cntr_base).up
            end do
            region(j).cntr_base = cntr_base

            ! Set centered position
            region(j).cntr_pos = region(j).strt_pos + ((region(j).length+1)/2-1)
        end do

        ! Print information on strand and region
        write(p_redir, "(i, 2(a, i5))"), i, &
            " - strd, # nts:", dna.strand(i).n_base, &
            ", # nts@Tn:", cn_tn
        do j = 1, n_region
            write(p_redir, "(i7, a$ )"), j, "-regn"
            write(p_redir, "(a,  i2$)"), ", T:",     region(j).types
            write(p_redir, "(a,  i4$)"), ", len:",   region(j).length
            write(p_redir, "(a,  i4$)"), ", str P:", region(j).strt_pos
            write(p_redir, "(a,  i4$)"), ", cen P:", region(j).cntr_pos
            write(p_redir, "(a,  i4 )"), ", end P:", region(j).end_pos
        end do
        write(p_redir, "(a)")
        write(p_redir, "(a)"), "      ----- Make nick ------"

        ! --------------------------------------------------
        ! Cut staple strand to make it short
        ! --------------------------------------------------
        begin_pos  = 0
        pre_region = 0
        b_cut      = .false.
        do j = 1, n_region

            if( para_break_stap_rule == "min" .and. &
                region(j).length < para_gap_xover_nick*2+2 ) cycle

            ! Find centered base and pos
            cntr_base = region(j).cntr_base
            cntr_pos  = region(j).cntr_pos
            length    = cntr_pos - begin_pos

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length >= para_min_break_stap) then
                b_cut = .true.
            end if

            ! Cutting method # 1 (minimum length staples)
            ! If minimum cutting length exceed, cut strand at the centered region
            if(para_break_stap_rule == "min" .and. b_cut == .true.) then

                pre_base = region(j).cntr_base
                pre_pos  = region(j).cntr_pos

                ! Skip the end iteration
                if(dna.strand(i).n_base - pre_pos < para_min_break_stap) then
                    write(p_redir, "(i7, a, i5)"), j, "-final rgn, len:", &
                        dna.strand(i).n_base - begin_pos
                    exit
                end if

                write(p_redir, "(i7, a, 2(a, i5))"), j, "->cut rgn", &
                    ", len:",     pre_pos-begin_pos, &
                    ", rmn len:", dna.strand(i).n_base - pre_pos

                ! Add # of staple
                dna.n_stap = dna.n_stap + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(pre_base).up
                dna.top(pre_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update new count and flag
                begin_pos = pre_pos
                b_cut     = .false.
            end if

            ! Cutting method # 2 and 3 (maximum or optimal length staples)
            if((para_break_stap_rule == "max" .or. para_break_stap_rule == "mid") .and. b_cut == .true.) then

                if( (para_break_stap_rule == "max" .and. length >= para_max_break_stap) .or. &
                    (para_break_stap_rule == "mid" .and. length >= para_mid_break_stap) ) then

                    ! To avoid small region in previous region
                    jj    = j - 1
                    b_ext = .false.

                    ! Find previous region
                    do
                        ! If the region length is longer than paramter value (default : 8)
                        if(region(jj).length >= para_gap_xover_nick*2+2) then
                            pre_base = region(jj).cntr_base
                            pre_pos  = region(jj).cntr_pos

                            ! Check remained length for the end of the strand
                            if(dna.strand(i).n_base - pre_pos >= para_min_break_stap) exit
                        end if

                        !Go back previous region
                        jj = jj - 1

                        ! Exception that the region cutted already
                        if(jj == pre_region) then
                            b_ext = .true.
                            exit
                        end if
                    end do

                    ! Make staple strand that exceeds para_max_break_stap/para_mid_break_stap
                    if(b_ext == .true.) then
                        if(j == n_region) then
                            write(p_redir, "(i7, a, i5)"), j, "-final rgn1, len:", &
                                dna.strand(i).n_base - begin_pos
                        end if
                        cycle
                    end if

                    ! If the cutted length is smaller than para_min_break_stap
                    if(pre_pos-begin_pos < para_min_break_stap) cycle

                    write(p_redir, "(i7, a, 3(a, i5))"), jj, "->cut rgn", &
                        ", len:",     pre_pos-begin_pos, &
                        ", cut pos:", pre_pos, &
                        ", rmn len:", dna.strand(i).n_base - pre_pos

                    ! Add # of staple
                    dna.n_stap = dna.n_stap + 1

                    ! Cut staple and make new connectivity
                    up_base = dna.top(pre_base).up
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1
                    pre_region           = jj

                    ! Update starting position and flag
                    begin_pos = pre_pos
                    b_cut     = .false.
                end if
            end if

            ! If remained length is smaller than para_max_break_stap, exit this loop
            if( (para_break_stap_rule == "max" .and. dna.strand(i).n_base-begin_pos < para_max_break_stap) .or. &
                (para_break_stap_rule == "mid" .and. dna.strand(i).n_base-begin_pos < para_mid_break_stap) ) then
                write(p_redir, "(i7, a, i5)"), j, "-final rgn2, len:", &
                    dna.strand(i).n_base - begin_pos
                exit
            end if

            ! Check final staple and print information (for only mid/max)
            if(j == n_region) then

                ! Final staple length
                final_length = dna.strand(i).n_base - begin_pos

                ! If the last staple exceeds maximum length
                if( (para_break_stap_rule == "max" .and. final_length >= para_max_break_stap) .or. &
                    (para_break_stap_rule == "mid" .and. final_length >= para_mid_break_stap) ) then

                    jj = j
                    do
                        ! Exception
                        if(jj == 0) then
                            deallocate(region)
                            return
                        end if

                        if(region(jj).length >= para_gap_xover_nick*2+2) then
                            pre_base = region(jj).cntr_base
                            pre_pos  = region(jj).cntr_pos
                            if(dna.strand(i).n_base - pre_pos >= para_min_break_stap) exit
                        end if
                        jj = jj - 1
                    end do

                    ! If the base has the same position with previous nick position
                    if(dna.top(pre_base).up == -1) then
                        write(p_redir, "(a)"), "   WARNING : This strand length exceeds para_max_break_stap"
                        cycle
                    end if

                    ! If the two cutted staples are smaller than para_min_break_stap
                    if( (final_length - (dna.strand(i).n_base-pre_pos) < para_min_break_stap) .or. &
                        (dna.strand(i).n_base - pre_pos < para_min_break_stap) ) then
                        write(p_redir, "(i7, a, i5)"), j, "-final rgn3, len:", &
                            dna.strand(i).n_base - begin_pos
                        cycle
                    end if

                    ! Add # of staple
                    up_base = dna.top(pre_base).up
                    dna.n_stap = dna.n_stap + 1
                    dna.top(pre_base).up = -1
                    dna.top(up_base).dn  = -1

                    write(p_redir, "(i7, a, i5)"), jj, "-final rgn4, len:", &
                        final_length - (dna.strand(i).n_base - pre_pos)
                    write(p_redir, "(i7, a, i5)"), jj, "-final rgn5, len:", &
                        dna.strand(i).n_base - pre_pos
                else
                    write(p_redir, "(i7, a, i5)"), j, "-final rgn6, len:", &
                        dna.strand(i).n_base - begin_pos
                end if
            end if
        end do
        write(p_redir, "(a)")

        deallocate(region)
    end do
end subroutine SeqDesign_Build_Sequence_Design

! -----------------------------------------------------------------------------

! Make the linear scaffold by the nick
subroutine SeqDesign_Make_Linear_Scaf_Nick_Outside(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Count # of crossovers based on edge section
    type :: CroLType
        integer :: n_xover
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, node, sec, across, edge, base, dn_base, strt_base
    integer :: max_strt_base, max_count, count, max_edge
    logical :: b_inside

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    do i = 1, geom.n_croL
        croL(i).n_xover = 0
    end do

    ! Find the number of crossovers and bps in terms of croL
    do i = 1, dna.strand(1).n_base
        base = dna.strand(1).base(i)
        node = dna.top(base).node

        ! For unpaired nucleotide
        if(node == -1) cycle

        edge = mesh.node(node).croL

        ! If the base is crossover
        if(dna.top(base).xover /= -1) then
            croL(edge).n_xover = croL(edge).n_xover + 1
        end if
    end do

    !do i = 1, geom.n_croL
    !    print *, i, "- : ", croL(i).n_xover
    !end do

    ! Loop for strands
    do i = 1, dna.n_strand

        ! only scaffold
        if(dna.strand(i).type1 == "stap") cycle

        b_inside = .false.
100     continue

        ! Set first base(going down) on the crossover
        base = dna.strand(i).base(1)
        do j = 1, dna.strand(i).n_base
            node = dna.top(base).node
            if(mesh.node(node).dn == -1) exit
            base = dna.top(base).dn
        end do

        ! Find maximum region in only non-crossover edges
        ! Scaffold nick will be located outside only for the 2D geometry
        max_strt_base = base
        max_count     = 0
        count         = 0
        do j = 1, dna.strand(i).n_base

            ! To avoid node ID is negative
            do
                if(dna.top(base).node /= -1) then

                    node = dna.top(base).node
                    sec  = mesh.node(node).sec

                    ! -------------------------------------
                    ! To replece nick at the bottom section
                    ! -------------------------------------
                    if(geom.sec.posR(sec+1) == 1) exit
                end if
                base = dna.top(base).up
            end do

            node = dna.top(base).node
            edge = mesh.node(node).croL

            ! Find max region only in non-crossover edges
            if(croL(edge).n_xover == 0) then
                node   = dna.top(base).node
                across = dna.top(base).across

                if( dna.top(base).xover     /= -1 .or. &    ! If crossover
                    dna.top(across).xover   /= -1 .or. &
                    dna.top(across).up      == -1 .or. &    ! If staple nick
                    dna.top(across).dn      == -1 .or. &
                    mesh.node(node).up      == -1 .or. &    ! If single crossover
                    mesh.node(node).dn      == -1 .or. &
                    mesh.node(node).mitered /= -1 .or. &    ! If mitered strand

                    ! To avoid placing inside
                    ( b_inside == .false. .and. geom.iniL(mesh.node(node).iniL).neiF(1) /= -1 .and. &
                      geom.iniL(mesh.node(node).iniL).neiF(2) /= -1 ) .or. &

                    ( b_inside == .false. .and. geom.iniL(mesh.node(node).iniL).neiF(2) == -1 .and. &     ! To place the unparied remaining scaffold outside
                      geom.croL(mesh.node(node).croL).sec == 0 ) .or. &

                    ( b_inside == .false. .and. geom.iniL(mesh.node(node).iniL).neiF(1) == -1 .and. &     ! To place the unparied remaining scaffold outside
                      geom.croL(mesh.node(node).croL).sec == 1 ) ) then

                    ! If the region exceeds max_count
                    if(max_count < count) then
                        max_count     = count
                        max_edge      = edge
                        max_strt_base = strt_base
                    end if

                    ! Set new starting base
                    count     = 0
                    strt_base = base
                else
                    count = count + 1
                end if
            end if

            ! Go to the upper base
            base = dna.top(base).up
        end do

        !max_count = max_count + 1

        ! If no nick position, set inside
        if(max_count == 0) then
            b_inside = .true.
            goto 100
        end if

        ! Set base in the middle at the maximum region
        base = dna.top(max_strt_base).id
        do j = 1, max_count / 2 + 1
            base = dna.top(base).up
        end do
        dn_base = dna.top(base).dn

        ! Disconnect between current and downward bases
        dna.top(base).dn    = -1
        dna.top(dn_base).up = -1
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.4. Make the linear scaffold by the nick"
    write(p_redir, "(a)"), "   * The edge number at the nick : "//trim(adjustl(Int2Str(max_edge)))
    write(p_redir, "(a)")

    ! Deallocate memory
    deallocate(croL)
end subroutine SeqDesign_Make_Linear_Scaf_Nick_Outside

! -----------------------------------------------------------------------------

! Make the linear scaffold by the nick
subroutine SeqDesign_Make_Linear_Scaf_Nick(geom, mesh, dna)
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    ! Count # of crossovers based on edge section
    type :: CroLType
        integer :: n_xover
    end type CroLType

    type(CroLType), allocatable :: croL(:)

    integer :: i, j, node, sec, across, edge, base, dn_base, strt_base
    integer :: max_strt_base, max_count, count, max_edge

    ! Allocate and initialize croL data
    allocate(croL(geom.n_croL))
    do i = 1, geom.n_croL
        croL(i).n_xover = 0
    end do

    ! Find the number of crossovers and bps in terms of croL
    do i = 1, dna.strand(1).n_base
        base = dna.strand(1).base(i)
        node = dna.top(base).node

        ! For unpaired nucleotide
        if(node == -1) cycle

        edge = mesh.node(node).croL

        ! If the base is crossover
        if(dna.top(base).xover /= -1) then
            croL(edge).n_xover = croL(edge).n_xover + 1
        end if
    end do

    ! Loop to make nick position in scaffold strand
    do i = 1, dna.n_strand

        ! Only for scaffold strand
        if(dna.strand(i).type1 == "stap") cycle

        ! Find first base(going down) that has single crossover
        base = dna.strand(i).base(1)
        do j = 1, dna.strand(i).n_base
            node = dna.top(base).node
            if(mesh.node(node).dn == -1) exit
            base = dna.top(base).dn
        end do

        ! Find maximum region in only non-crossover edges
        max_strt_base = base
        max_count     = 0
        count         = 0
        do j = 1, dna.strand(i).n_base

            ! To avoid node ID is negative
            do
                if(dna.top(base).node /= -1) then

                    node = dna.top(base).node
                    sec  = mesh.node(node).sec

                    ! To replece nick at the bottom section
                    if(geom.sec.posR(sec+1) == para_pos_nick_scaf) exit
                end if
                base = dna.top(base).up
            end do

            node = dna.top(base).node
            edge = mesh.node(node).croL

            ! Find max region only in non-crossover edges
            if(croL(edge).n_xover == 0) then
                node   = dna.top(base).node
                across = dna.top(base).across

                if( dna.top(base).xover   /= -1 .or. &  ! If there is crossover
                    dna.top(across).xover /= -1 .or. &
                    dna.top(across).up    == -1 .or. &  ! If there is staple nick
                    dna.top(across).dn    == -1 .or. &
                    mesh.node(node).up    == -1 .or. &  ! If there is single crossover
                    mesh.node(node).dn    == -1 ) then

                    ! If the region exceeds max_count
                    if(max_count < count) then
                        max_count     = count
                        max_edge      = edge
                        max_strt_base = strt_base
                    end if

                    ! Set new starting base
                    count     = 0
                    strt_base = base
                else
                    count = count + 1
                end if
            end if

            ! Go to the upper base
            base = dna.top(base).up
        end do

        ! Set base in the middle at the maximum region
        base = dna.top(max_strt_base).id
        do j = 1, max_count / 2 + 1
            base = dna.top(base).up
        end do
        dn_base = dna.top(base).dn

        ! Disconnect between current and downward bases
        dna.top(base).dn    = -1
        dna.top(dn_base).up = -1
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.4. Make the linear scaffold by the nick"
    write(p_redir, "(a)"), "   * The edge number at the nick : "//trim(adjustl(Int2Str(max_edge)))
    write(p_redir, "(a)")

    ! Deallocate memory
    deallocate(croL)
end subroutine SeqDesign_Make_Linear_Scaf_Nick

! -----------------------------------------------------------------------------

! Make short scaffold strand
subroutine SeqDesign_Make_Short_Scaf(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type :: RegionType
        integer :: length       ! total length at the region not including boundary

        integer :: strt_pos     ! start position
        integer :: end_pos      ! end position
        integer :: cntr_pos     ! center position

        integer :: strt_base    ! start base
        integer :: cntr_base    ! center base
    end type RegionType

    type(RegionType), allocatable :: region(:)

    integer :: i, j, k, base, across, cntr_base, pre_base, up_base
    integer :: n_region, length, cntr_pos, pre_pos, begin_pos
    logical :: b_region, b_cut

    ! Exception current subroutine
    if(para_max_break_scaf == 0) then
        return
    end if

    ! Loop to make short scaffold strand
    do i = 1, dna.n_strand

        ! Only for scaffold strand
        if(dna.strand(i).type1 == "stap") cycle

        ! Find the starting point (down == -1)
        base = Mani_Go_Start_Base(dna, i)

        ! --------------------------------------------------
        ! Build region array
        ! --------------------------------------------------
        ! Tn loop |<---------->| 12 bp
        !     ~~~~~============*===========*=====*  : base pair
        !                      |<--------->|
        !               region 11 not including crossover and nick region
        n_region = 0
        allocate(region(dna.strand(i).n_base))

        do j = 1, dna.strand(i).n_base

            ! Find across ID
            across = dna.top(base).across

            if(across == -1) then
                ! For Tn loop
                b_region = .true.

                ! Update base to go up
                base = dna.top(base).up
                cycle
            else
                ! Check down, up, xover and across
                ! Starting base always has -1 (down = -1)
                if( dna.top(base).dn      == -1 .or. &
                    dna.top(base).up      == -1 .or. &
                    dna.top(base).xover   /= -1 .or. &
                    dna.top(across).xover /= -1 ) then

                    b_region = .true.

                    ! Update base to go up
                    base = dna.top(base).up
                    cycle
                end if
            end if

            if(b_region == .true.) then
                ! Make new region
                b_region = .false.
                n_region = n_region + 1

                ! Set region data
                region(n_region).strt_base = base
                region(n_region).length    = 1
                region(n_region).strt_pos  = j
                region(n_region).end_pos   = j
            else
                ! Increase length
                region(n_region).length  = region(n_region).length + 1
                region(n_region).end_pos = j
            end if

            ! Update base to go up
            base = dna.top(base).up
        end do

        ! Set centered base and position
        do j = 1, n_region

            ! Find centered base
            cntr_base = region(j).strt_base
            do k = 1, (region(j).length+1)/2-1
                cntr_base = dna.top(cntr_base).up
            end do
            region(j).cntr_base = cntr_base

            ! Set centered position
            region(j).cntr_pos = region(j).strt_pos + ((region(j).length+1)/2-1)
        end do

        ! Print information on strand and region
        !write(0, "(i6, a, i5)"), i, &
        !    " - strd, # nts:", dna.strand(i).n_base
        !write(0, "(a)")

        !do j = 1, n_region
        !    write(0, "(i8, a, 4(a, i5))"), j, "-rgn", &
        !        ", len:",   region(j).length,   &
        !        ", str P:", region(j).strt_pos, &
        !        ", ctr P:", region(j).cntr_pos, &
        !        ", end P:", region(j).end_pos
        !end do
        !write(0, "(a)")
        !call Space(0, 19)
        !write(0, "(a)"), "----- Make nick ----"
        !write(0, "(a)")

        ! --------------------------------------------------
        ! Cut strand to make it short
        ! --------------------------------------------------
        begin_pos  = 0
        b_cut      = .false.
        do j = 1, n_region

            ! Find centered base and pos
            cntr_base = region(j).cntr_base
            cntr_pos  = region(j).cntr_pos
            length    = cntr_pos - begin_pos

            ! If the length exceeds the minimum cutting length
            if(b_cut == .false. .and. length > para_max_break_scaf) then
                b_cut = .true.
            end if

            ! Cutting method # 1 (minimum length staples)
            ! If minimum cutting length exceed, cut strand at the centered region
            if(b_cut == .true.) then

                pre_base = region(j).cntr_base
                pre_pos  = region(j).cntr_pos

                ! Skip the end iteration
                if(j == n_region) cycle

                write(p_redir, "(i16, a, 2(a, i5))"), j, "->cut rgn", &
                    ", len:",     pre_pos-begin_pos, &
                    ", rmn len:", dna.strand(i).n_base - pre_pos

                ! Add # of staple
                dna.n_scaf = dna.n_scaf + 1

                ! Cut staple and make new connectivity
                up_base = dna.top(pre_base).up
                dna.top(pre_base).up = -1
                dna.top(up_base).dn  = -1

                ! Update new count and flag
                begin_pos = pre_pos
                b_cut     = .false.
            end if

            ! Print final strand
            if(j == n_region) then
                write(p_redir, "(i16, a, i4)"), j, " -1 cut rgn, len:", &
                    dna.strand(i).n_base - pre_pos
            end if
        end do
        write(p_redir, "(a)")

        deallocate(region)
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.5. Make short scaffold"
    write(p_redir, "(a)"), "   * Maximum # of bases per scaffold : "//trim(adjustl(Int2Str(para_max_break_scaf)))
    write(p_redir, "(a)"), "   * # of modified scaffold : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a)"), "   * # of modified scaffold : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a)")
end subroutine SeqDesign_Make_Short_Scaf

! -----------------------------------------------------------------------------

! Move cur_base to avoid node without ID and crossovers
function SeqDesign_Avoid_Barrier(mesh, dna, base, gap) result(cur_base)
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer,        intent(in) :: base
    integer,        intent(in) :: gap

    integer :: i, cur_base, chk_base, node, xover
    logical :: nofront, noback

    cur_base = base

    ! --------------------------------------------------
    ! Move cur_base to avoid node without ID
    ! --------------------------------------------------
    do
        ! Check front region if there is node without ID
        nofront  = .true.
        chk_base = cur_base

        do i = 1, gap
            node = dna.top(chk_base).node

            ! For staple strand
            if(node == -1) then
                noback = .false.
                exit
            end if

            ! For scaffold strand
            if(mesh.node(node).up == -1 .or. mesh.node(node).dn == -1) then
                nofront = .false.
                exit
            end if
            chk_base = dna.top(chk_base).up
        end do

        ! Check backward region if there is node without ID
        noback   = .true.
        chk_base = cur_base

        do i = 1, gap
            node = dna.top(chk_base).node

            ! For staple strand
            if(node == -1) then
                noback = .false.
                exit
            end if

            ! For scaffold strand
            if(mesh.node(node).up == -1 .or. mesh.node(node).dn == -1) then
                noback = .false.
                exit
            end if
            chk_base = dna.top(chk_base).dn
        end do

        ! If two critera was satified
        if(nofront == .true. .and. noback == .true.) exit

        ! Update the cur_base
        cur_base = dna.top(cur_base).dn
    end do

    ! --------------------------------------------------
    ! Move cur_base to crossover if there is crossover nearby
    ! --------------------------------------------------
    do
        ! Check front region whether there is crossovers
        nofront  = .true.
        chk_base = cur_base

        do i = 1, gap
            xover = dna.top(chk_base).xover
            if(xover /= -1) then
                exit
            end if
            chk_base = dna.top(chk_base).up
        end do

        if(xover /= -1) then
            cur_base = dna.top(chk_base).up
            if(dna.top(cur_base).xover == -1) then
                cur_base = dna.top(cur_base).dn
            end if
            exit
        end if

        ! Check backward region whether there is crossovers
        noback   = .true.
        chk_base = cur_base

        do i = 1, gap
            xover = dna.top(chk_base).xover
            if(xover /= -1) then
                exit
            end if
            chk_base = dna.top(chk_base).dn
        end do

        if(xover /= -1) then
            cur_base = chk_base
            exit
        end if

        ! If two critera was satified
        if(nofront == .true. .and. noback == .true.) exit

        ! Update the cur_base
        cur_base = dna.top(cur_base).dn
    end do
end function SeqDesign_Avoid_Barrier

! -----------------------------------------------------------------------------

! Make short strand
subroutine SeqDesign_Make_Short_Strand(mesh, dna)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, max_strand, rem_count, cur_base, down_base
    logical :: b_last

    ! Loop for all strand
    do i = 1, dna.n_strand

        ! Set the maximum number of bases per strand
        if(dna.strand(i).type1 == "scaf") then
            max_strand = para_max_break_scaf
        else if(dna.strand(i).type1 == "stap") then
            max_strand = para_max_break_stap
        end if

        if(max_strand == 0) cycle

        ! If it exceeds certain length, cut strand
        if(dna.strand(i).n_base > max_strand) then

            ! Find starting point
            cur_base = dna.strand(i).base(1)
            do
                if(dna.top(cur_base).dn == -1) exit
                cur_base = dna.top(cur_base).dn
            end do

            ! Cut strand
            do
                ! Count the remaining number of bases
                rem_count = SeqDesign_Count_Remainder(dna, cur_base)

                ! Go cur_base up
                if(rem_count > 2 * max_strand) then
                    b_last = .false.
                    do j = 1, max_strand
                        cur_base = dna.top(cur_base).up
                    end do
                else
                    b_last = .true.
                    do j = 1, rem_count / 2
                        cur_base = dna.top(cur_base).up
                    end do
                end if

                ! Move cur_base to avoid node without ID and crossovers
                cur_base  = SeqDesign_Avoid_Barrier(mesh, dna, cur_base, 5)
                down_base = dna.top(cur_base).dn

                ! Disconnect between current and downward bases
                dna.top(cur_base).dn  = -1
                dna.top(down_base).up = -1

                ! If the cur_base is crossover, disconnect the crossover
                if(dna.top(cur_base).xover == dna.top(down_base).id) then
                    dna.top(cur_base).xover  = -1
                    dna.top(down_base).xover = -1
                end if

                if(dna.strand(i).type1 == "scaf") then
                    dna.n_scaf = dna.n_scaf + 1
                else if(dna.strand(i).type1 == "stap") then
                    dna.n_stap = dna.n_stap + 1
                end if

                if(b_last == .true.) exit
            end do
        end if
    end do

    ! Print progress
    write(p_redir, "(a )"), "6.4. Make short staples"
    write(p_redir, "(a$)"), "* Maximum # of bases per scaffold : "
    write(p_redir, "(a )"), trim(adjustl(Int2Str(para_max_break_scaf)))
    write(p_redir, "(a$)"), "* Maximum number of bases per staple : "
    write(p_redir, "(a )"), trim(adjustl(Int2Str(para_max_break_stap)))
    write(p_redir, "(a$)"), "* # of modified scaffold : "
    write(p_redir, "(a )"), trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a$)"), "* # of modified staples : "
    write(p_redir, "(a )"), trim(adjustl(Int2Str(dna.n_stap)))
end subroutine SeqDesign_Make_Short_Strand

! -----------------------------------------------------------------------------

! Count the number of remainder bases, strand should be non-circular
function SeqDesign_Count_Remainder(dna, cur_base) result(count)
    type(DNAType), intent(in) :: dna
    integer, intent(in) :: cur_base

    integer :: i, count, base

    count = 1
    base  = cur_base

    ! Find the number of remaining bases
    do
        base = dna.top(base).up
        if(base == -1) then
            exit
        else
            count = count + 1
        end if
    end do
end function SeqDesign_Count_Remainder

! -----------------------------------------------------------------------------

! Rebuild strand data from dnaTop data
subroutine SeqDesign_Rebuild_Strand(dna)
    type(DNAType), intent(inout) :: dna

    logical, allocatable :: b_visit(:)
    type(ListBase), pointer :: list_base
    type(ListBase), pointer :: ptr_base
    type(TopType) :: cur_base

    integer :: i, j, n_strand, n_base, int_base
    logical :: b_end

    ! Deallocate previous strand data
    deallocate(dna.strand)

    ! Allocate and initialize strand data
    dna.n_strand = dna.n_scaf + dna.n_stap
    allocate(dna.strand(dna.n_strand))
    call Mani_Init_StrandType(dna.strand, dna.n_strand)

    ! Set strand type
    do i = 1, dna.n_strand
        if(i <= dna.n_scaf) then
            dna.strand(i).type1 = "scaf"
        else
            dna.strand(i).type1 = "stap"
        end if
    end do

    ! Nullify the linked lists
    nullify(list_base, ptr_base)

    ! Allocate and initilize the b_visit data
    allocate(b_visit(dna.n_top))
    do i = 1, dna.n_top
        b_visit(i) = .false.
    end do
    
    n_strand = 0
    do
        ! Check b_visit if there is 0 (0 means not visiting base)
        b_end = .true.
        do i = 1, dna.n_top
            if(b_visit(dna.top(i).id) == .false.) then
                ! All visit of all bases
                b_end = .false.
                exit
            end if
        end do
        if(b_end == .true.) exit

        ! Increase the number of strands
        n_strand = n_strand + 1
        cur_base = dna.top(i)
        int_base = cur_base.id

        ! Find the first base in the current strand
        do
            if(cur_base.dn == -1) then
                ! cur_base is at the 3'-end of the strand
                dna.strand(n_strand).b_circular = .false.
                exit
            else if(cur_base.dn == int_base) then
                ! cur_base goes back to the starting point
                dna.strand(n_strand).b_circular = .true.
                cur_base = dna.top(int_base)
                exit
            end if

            cur_base = dna.top(cur_base.dn)

            if(b_visit(cur_base.id)) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Reached a visited base.                          |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if
        end do

        ! Walk through the current strand
        n_base = 1

        ! Insert data to linked list
        allocate(ptr_base)
        ptr_base%id = cur_base.id
        list_base => List_Insert_Base(list_base, ptr_base)

        dna.top(cur_base.id).strand  = n_strand
        dna.top(cur_base.id).address = n_base
        dna.top(cur_base.id).b_14nt  = .false.
        b_visit(cur_base.id)         = .true.

        ! Loop to add a new base into current strand
        do
            ! Check for going out loop
            if(dna.strand(n_strand).b_circular == .false.) then
                ! If non-circular and upper ID is equal to -1
                if(cur_base.up == -1) exit
            else
                ! If circular and base ID is equal to init base ID
                if(cur_base.up == int_base) exit
            end if

            cur_base = dna.top(cur_base.up)

            if(b_visit(cur_base.id)) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Reached a visited base.                          |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            n_base = n_base + 1

            ! Insert data to linked list
            allocate(ptr_base)
            ptr_base%id = cur_base.id
            list_base => List_Insert_Base(list_base, ptr_base)

            dna.top(cur_base.id).strand  = n_strand
            dna.top(cur_base.id).address = n_base
            dna.top(cur_base.id).b_14nt  = .false.
            b_visit(cur_base.id)         = .true.
        end do

        dna.strand(n_strand).n_base = n_base

        ! Allocate array using the linked list
        allocate(dna.strand(n_strand).base(n_base))

        ! Put strand data from linked list
        ptr_base => list_base
        do i = 1, n_base
            dna.strand(n_strand).base(n_base+1-i) = ptr_base%id
            ptr_base => ptr_base%next
        end do

    end do

    ! Deallocate b_visit data
    deallocate(b_visit)

    ! Delete linked list allocated
    call List_Delete_Base(list_base)
    !call List_Delete_Base(ptr_base)

    ! Write sequence data
    dna.len_min_stap =  10000
    dna.len_max_stap = -10000
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "stap") then
            if(dna.strand(i).n_base < dna.len_min_stap) dna.len_min_stap = dna.strand(i).n_base
            if(dna.strand(i).n_base > dna.len_max_stap) dna.len_max_stap = dna.strand(i).n_base
        end if
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.5. Rebuild strand data"
    write(p_redir, "(a)"), "   * Total # of strands        : "//trim(adjustl(Int2Str(dna.n_strand)))
    write(p_redir, "(a)"), "   * # of scaffold strands     : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a)"), "   * # of staple strands       : "//trim(adjustl(Int2Str(dna.n_stap)))
    write(p_redir, "(a)"), "   * The minimum staple length : "//trim(adjustl(Int2Str(dna.len_min_stap)))
    write(p_redir, "(a)"), "   * The maximum staple length : "//trim(adjustl(Int2Str(dna.len_max_stap)))
    write(p_redir, "(a)"), "   * Detailed info. on DNA strand data"

    ! Print progress in detail
    if(p_detail == .true.) then
        do i = 1, dna.n_strand
            write(p_redir, "(i11, a$)"), i, " strand -> "
            write(p_redir, "(a, i7$ )"), "# of bases : ", dna.strand(i).n_base
            write(p_redir, "(a, l, a)"), ", circular : ", dna.strand(i).b_circular, ", bases : "

            ! Print bases numbering
            write(p_redir, "(i20, a$)"), dna.strand(i).base(1), "->"
            do j = 2, dna.strand(i).n_base - 1
                if(mod(j, 20)== 0) then
                    write(p_redir, "(i7, a  )"), dna.strand(i).base(j), "->"
                else if(mod(j, 20)== 1) then
                    write(p_redir, "(i20, a$)"), dna.strand(i).base(j), "->"
                else
                    write(p_redir, "(i7, a$ )"), dna.strand(i).base(j), "->"
                end if
            end do
            write(p_redir, "(i7)"), dna.strand(i).base(dna.strand(i).n_base)
            write(p_redir, "(a )")
        end do
    end if
    write(p_redir, "(a)")
end subroutine SeqDesign_Rebuild_Strand

! -----------------------------------------------------------------------------

! List in long length order of the staple
subroutine SeqDesign_Order_Staple(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i

    ! dna.n_stap(i, 1) - original ID
    ! dna.n_stap(i, 2) - short order ID
    allocate(dna.order_stap(dna.n_stap, 2))

    do i = 1, dna.n_stap
        dna.order_stap(i, 1) = i + dna.n_scaf
        dna.order_stap(i, 2) = dna.strand(i+dna.n_scaf).n_base
    end do

    ! Sort
    call Sort2(dna.order_stap(:, 2), dna.order_stap(:, 1), dna.n_stap)

    ! Print progress
    write(p_redir, "(a)"), "     Short length order of staples"
    write(p_redir, "(a)"), "     #ID    # ori. ID     Staple length"
    write(p_redir, "(a)"), "     ---    ---------     -------------"
    do i = 1, dna.n_stap
        write(p_redir, "(i7$ )"), i
        write(p_redir, "(i12$)"), dna.order_stap(i, 1)
        write(p_redir, "(i14 )"), dna.order_stap(i, 2)
    end do
    write(p_redir, "(a)")
end subroutine SeqDesign_Order_Staple

! -----------------------------------------------------------------------------

! Print 14nt region with various representations
subroutine SeqDesign_Print_14nt_Region_Simple(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna

    type(RegionType), allocatable :: region(:)
    integer,          allocatable :: region_14nt(:)
    integer,          allocatable :: region_4nt(:)
    integer,          allocatable :: scaf_base2idx(:)
    type(GraphType)               :: graph

    integer :: n_nt_14nt, n_nt_4nt, cn_tn, n_win, n_14nt, n_4nt, n_only_4nt, n_sec_14nt, n_sec_4nt
    integer :: len_ave_stap, one_14nt, two_14nt, three_14nt, other_14nt, n_count(30)
    integer :: i, j, k, base, n_region, n_edge, tot_14nt, tot_4nt
    integer :: across1, across2, node, iniL
    integer :: cir_rgb(3), cir_size
    logical :: b_14nt, b_4nt
    character(200) :: path

    ! --------------------------------------------------
    ! Intialize and setup the data structures
    ! --------------------------------------------------
    n_14nt   = 0; n_4nt    = 0; n_only_4nt = 0
    one_14nt = 0; two_14nt = 0; three_14nt = 0; other_14nt = 0

    n_nt_14nt = 0; n_nt_4nt = 0
    len_ave_stap = 0

    ! Initialize the variable for couting regions
    allocate(region_14nt(dna.n_strand))
    allocate(region_4nt(dna.n_strand))
    region_14nt(:) = 0
    region_4nt(:)  = 0
    n_count(1:30)  = 0

    ! Count total edges for circular graph
    n_edge = SeqDesign_CirGraph_Count_Edge(dna)

    ! Allocate and initialize variables for the circular graph
    call SeqDesign_CirGraph_Init_Variable(dna, graph, n_edge)

    ! Build scaffold index for bar and circular graph
    allocate(scaf_base2idx(dna.n_base_scaf))
    base = Mani_Go_Start_Base(dna, 1)
    do i = 1, dna.n_base_scaf

        scaf_base2idx(base)   = i       ! From base to index
        graph.base2node(base) = i       ! From base to index
        graph.node2base(i)    = base    ! From index to base

        ! Node type depending on edge number
        if(dna.top(base).node == -1) then
            graph.node(base) = 0
        else
            graph.node(base) = mesh.node(dna.top(base).node).iniL
        end if

        base = dna.top(base).up
    end do

    ! Print progress
    write(p_redir, "(a)"), "  6.6. Print 14nt dsDNA domains"
    write(p_redir, "(a)"), "   * # of staple strands     : "//trim(adjustl(Int2Str(dna.n_stap)))
    write(p_redir, "(a)"), "   * # of total 14nt domains : "//trim(adjustl(Int2Str(ubound(graph.edge, 1))))
    write(p_redir, "(a)"), "   * Detailed info. on the dsDNA domain"
    write(p_redir, "(a)")

    ! --------------------------------------------------
    ! Calculate regions and # of staples
    ! --------------------------------------------------
    n_edge = 0
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).type1 == "scaf") then
            ! --------------------------------------------------
            ! For scaffold strand in bar graph - Graph #3
            ! --------------------------------------------------
            dna.strand(i).type2 = "vertex"
            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                node = dna.top(base).node
                if(node /= -1) then
                    iniL = mesh.node(node).iniL
                else
                    iniL = 0
                end if
            end do
            cycle
        end if

        ! Count staple length
        len_ave_stap = len_ave_stap + dna.strand(i).n_base

        ! Count unpaired staple nucleotides if it contains Tn loop
        cn_tn = 0
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            if(dna.top(base).across == -1) cn_tn = cn_tn + 1
        end do

        ! Build region of staple strands
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple_1(dna, i, region, n_region)

        ! --------------------------------------------------
        ! Count 14nt regions - Graph #2
        ! --------------------------------------------------
        n_sec_14nt = 0;   n_sec_4nt  = 0
        b_14nt = .false.; b_4nt= .false.

        do j = 1, n_region

            ! Check strand type 2
            if(region(j).types == 1) dna.strand(i).type2 = "vertex"

            if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
                (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then

                ! Just check one 14nt seeds
                if(b_14nt == .false.) then
                    n_14nt = n_14nt + 1
                    b_14nt = .true.
                end if

                ! All 14nt seeds in the region
                n_sec_14nt = n_sec_14nt + 1
                
                ! Count # of nucleotides in 14nt seed domains
                if(region(j).types == 1) n_nt_14nt = n_nt_14nt + region(j).length + 1
                if(region(j).types == 2) n_nt_14nt = n_nt_14nt + region(j).length + 2
            end if

            ! Global 14nt region couting
            if(region(j).length + 2 <= 30) then     ! There are no region longer than 30bp length
                if(region(j).types == 1) n_count(region(j).length + 1) = n_count(region(j).length + 1) + 1
                if(region(j).types == 2) n_count(region(j).length + 2) = n_count(region(j).length + 2) + 1

                ! Find two crossovers region, 2nt region
                if(j == 1) cycle 
                if(dna.top(dna.top(region(j  ).sta_base).dn).across == -1) cycle
                if(dna.top(dna.top(region(j-1).end_base).up).across == -1) cycle

                n_win = (region(j).sta_pos -1) - (region(j-1).end_pos + 1) - 1
                if(n_win >= 1) n_count(n_win) = n_count(n_win) + 1
            end if
        end do

        ! Count 4nt regions
        do j = 1, n_region
            if( (region(j).types == 1 .and. region(j).length + 1 <= 4) .or. &
                (region(j).types == 2 .and. region(j).length + 2 <= 4) ) then

                ! Just check one 4nt seeds
                if(b_4nt == .false.) then
                    n_4nt = n_4nt + 1
                    b_4nt = .true.
                end if

                ! # of 4nt seeds in the region
                n_sec_4nt = n_sec_4nt + 1

                ! Count # of nucleotides in 4nt seed domains
                if(region(j).types == 1) n_nt_4nt = n_nt_4nt + region(j).length + 1
                if(region(j).types == 2) n_nt_4nt = n_nt_4nt + region(j).length + 2
            end if
        end do

        if(n_sec_14nt == 1) region_14nt(i) = 1
        if(n_sec_14nt == 2) region_14nt(i) = 2
        if(n_sec_14nt == 3) region_14nt(i) = 3
        if(n_sec_14nt == 4) region_14nt(i) = 4

        if(n_sec_4nt  == 1) region_4nt(i)  = 1
        if(n_sec_4nt  == 2) region_4nt(i)  = 2
        if(n_sec_4nt  == 3) region_4nt(i)  = 3
        if(n_sec_4nt  == 4) region_4nt(i)  = 4

        ! --------------------------------------------------
        ! For staple strand in bar graph - Graph #3
        ! --------------------------------------------------
        ! Find the starting point (down == -1)
        base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base
            across1 = dna.top(base).across
            base = dna.top(base).up
        end do

        ! --------------------------------------------------
        ! Build node and edge for circular graph - Graph #4
        ! --------------------------------------------------
        n_sec_14nt = 0
        do j = 1, n_region

            ! Set node type with different color
            base = dna.top(region(j).cen_base).across

            ! Set node type depending on region size
            if( (region(j).types == 1 .and. region(j).length + 1 >= 14) .or. &
                (region(j).types == 2 .and. region(j).length + 2 >= 14) ) then

                n_sec_14nt = n_sec_14nt + 1
                !if(n_sec_14nt == 1) graph.node(base) = 1
                !if(n_sec_14nt == 2) graph.node(base) = 2
            else if( (region(j).types == 1 .and. region(j).length + 1 <= 4) .or. &
                     (region(j).types == 2 .and. region(j).length + 2 <= 4) ) then
                !graph.node(base) = 3
            else
                !graph.node(base) = 0
            end if

            ! Set Edge data
            if(j /= n_region) then

                n_edge = n_edge + 1

                ! Set edge in graph
                graph.edge(n_edge, 1) = dna.top(region(j  ).cen_base).across
                graph.edge(n_edge, 2) = dna.top(region(j+1).cen_base).across
            end if
        end do

        ! Update variables
        dna.strand(i).n_14nt = n_sec_14nt
        dna.strand(i).n_4nt  = n_sec_4nt
        dna.len_ave_stap     = dble(len_ave_stap)/dble(dna.n_stap)

        ! Check 14nt only region
        if(b_14nt == .false. .and. b_4nt == .true.) then
            n_only_4nt = n_only_4nt + 1
        end if

        ! Print information on region of staple strands
        write(p_redir, "(i6$   )"), i
        write(p_redir, "(a, i4$)"), "-stap, # nts:", dna.strand(i).n_base
        write(p_redir, "(a, i2$)"), ", # nts@Tn:", cn_tn
        write(p_redir, "(a, l$ )"), ", 14nt("//trim(adjustl(Int2Str(n_sec_14nt)))//")-", b_14nt
        write(p_redir, "(a, l  )"), ", 4nt("//trim(adjustl(Int2Str(n_sec_4nt)))//")-", b_4nt

        do j = 1, n_region
            write(p_redir, "(i8, a$)"), j, "-rgn"

            if(region(j).types == 1) write(p_redir, "(a$)"), " [vertex], "
            if(region(j).types == 2) write(p_redir, "(a$)"), " [ edge ], "

            if(region(j).types == 1) write(p_redir, "(a, 2(i2, a)$)"), "len:", region(j).length, "(", region(j).length+1, ")"
            if(region(j).types == 2) write(p_redir, "(a, 2(i2, a)$)"), "len:", region(j).length, "(", region(j).length+2, ")"

            write(p_redir, "(a, i3$)"), ", str P:", region(j).sta_pos
            write(p_redir, "(a, i3 )"), ", end P:",   region(j).end_pos
        end do
        write(p_redir, "(a)")

        ! Deallocate memory
        deallocate(region)
    end do

    ! Update global nucleotides in 14nt and 4nt dsDNA domains
    dna.n_nt_14nt = n_nt_14nt
    dna.n_nt_4nt  = n_nt_4nt

    ! --------------------------------------------------
    ! Print 14nt regions - Graph #2
    ! --------------------------------------------------
    tot_4nt  = 0
    tot_14nt = 0
    do i = 1, dna.n_strand
        if(region_14nt(i) == 1) one_14nt   = one_14nt   + 1
        if(region_14nt(i) == 2) two_14nt   = two_14nt   + 1
        if(region_14nt(i) == 3) three_14nt = three_14nt + 1
        if(region_14nt(i) == 4) other_14nt = other_14nt + 1

        tot_14nt = tot_14nt + region_14nt(i)
        tot_4nt  = tot_4nt  + region_4nt(i)
    end do

    ! Update n_14nt and n_4nt seeds in dna data
    dna.n_14nt       = n_14nt
    dna.n_s14nt      = two_14nt
    dna.n_4nt        = n_4nt
    dna.n_only_4nt   = n_only_4nt
    dna.n_tot_region = n_edge
    dna.n_tot_14nt   = tot_14nt
    dna.n_tot_4nt    = tot_4nt

    ! Print progress
    write(p_redir, "(a)"), "   * # of staples with 14nt dsDNA domains : "//trim(adjustl(Int2Str(n_14nt)))
    write(p_redir, "(a)"), "     * One   14nt dsDNA domain  : "//trim(adjustl(Int2Str(one_14nt)))
    write(p_redir, "(a)"), "     * Two   14nt dsDNA domains : "//trim(adjustl(Int2Str(two_14nt)))
    write(p_redir, "(a)"), "     * Three 14nt dsDNA domains : "//trim(adjustl(Int2Str(three_14nt)))
    write(p_redir, "(a)"), "     * More  14nt dsDNA domains : "//trim(adjustl(Int2Str(other_14nt)))
    write(p_redir, "(a)")
    write(p_redir, "(a)"), "   * Length of the dsDNA domain (bp)"
    do k = 1, 25
        write(p_redir, "(a, i3, a)"), "      * ", k, &
            " bp dsDNA domain : "//trim(adjustl(Int2Str(n_count(k))))
    end do
    write(p_redir, "(a)")

    ! Deallocate memory
    deallocate(region_14nt)
    deallocate(region_4nt)
    deallocate(scaf_base2idx)
    deallocate(graph.node2base)
    deallocate(graph.base2node)
    deallocate(graph.node)
    deallocate(graph.edge)
end subroutine SeqDesign_Print_14nt_Region_Simple

! -----------------------------------------------------------------------------

! Count the total number of regions for all staples
function SeqDesign_CirGraph_Count_Edge(dna) result(n_edge)
    type(DNAType), intent(inout) :: dna

    type(RegionType), allocatable :: region(:)
    integer :: i, n_region, n_edge

    ! Count the number of edges
    n_edge = 0
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).type1 == "scaf") cycle

        ! Find region
        allocate(region(dna.strand(i).n_base))
        call SeqDesign_Build_Region_Staple(dna, i, region, n_region)

        ! Increse edge number
        n_edge = n_edge + (n_region - 1)
        deallocate(region)
    end do
end function SeqDesign_CirGraph_Count_Edge

! -----------------------------------------------------------------------------

! Allocate and initialize variables for the circular graph
subroutine SeqDesign_CirGraph_Init_Variable(dna, graph, n_edge)
    type(DNAType),   intent(in)    :: dna
    type(GraphType), intent(inout) :: graph
    integer,         intent(in)    :: n_edge

    allocate(graph.node2base(dna.n_base_scaf))
    allocate(graph.base2node(dna.n_base_scaf))
    allocate(graph.node(dna.n_base_scaf))
    allocate(graph.edge(n_edge, 3))

    graph.node2base(1:dna.n_base_scaf) = 0
    graph.base2node(1:dna.n_base_scaf) = 0
    graph.node(1:dna.n_base_scaf)      = 0
    graph.edge(1:n_edge, 1) = 0
    graph.edge(1:n_edge, 2) = 0
    graph.edge(1:n_edge, 3) = 0
end subroutine SeqDesign_CirGraph_Init_Variable

! -----------------------------------------------------------------------------

! Assign DNA sequence according to para_set_seq_scaf
subroutine SeqDesign_Assign_Sequence(prob, dna)
    type(ProbType), intent(in)    :: prob
    type(DNAType),  intent(inout) :: dna

    integer :: i

    ! Print progress
    write(p_redir, "(a)"), "  6.7. Set staple sequences"
    write(p_redir, "(a)"), "   * Scaffold sequence        : "//trim(para_scaf_seq)
    write(p_redir, "(a)"), "   * # of scaffolds           : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a)"), "   * # of nts in the scaffold : "//trim(adjustl(Int2Str(dna.n_base_scaf)))

    if(trim(para_scaf_seq) == "m13") then

        ! Set M13mp18 DNA sequence
        call SeqDesign_Set_M13mp18(dna)
    else if(trim(para_scaf_seq) == "user") then

        ! Import sequence from txt file
        call SeqDesign_Import_Sequence(prob, dna)
    else if(trim(para_scaf_seq) == "rand") then

        ! Set random sequence
        call SeqDesign_Set_Rand_Sequence(dna)
    else
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Scaffold sequences are not assigned.             |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Print progress
    write(p_redir, "(a)"), "   * Detailed info. on DNAtop data"

    ! Print detail information
    if(p_detail == .true.) then
        do i = 1, dna.n_top
            write(p_redir, "(i10, a$)"), dna.top(i).id, " base ==>"
            write(p_redir, "(a,  i6$)"), "  up: ",   dna.top(i).up
            write(p_redir, "(a,  i6$)"), ", dn: ",   dna.top(i).dn
            write(p_redir, "(a,  i6$)"), ", ac: ",   dna.top(i).across
            write(p_redir, "(a,  i6$)"), ", xo: ",   dna.top(i).xover
            write(p_redir, "(a,  i6$)"), ", nde: ",  dna.top(i).node
            write(p_redir, "(a,  a3$)"), ", seq: ",  dna.top(i).seq
            write(p_redir, "(a,  i3$)"), ", strd: ", dna.top(i).strand
            write(p_redir, "(a,  i6 )"), ", addr: ", dna.top(i).address
        end do
    end if
    write(p_redir, "(a)")
end subroutine SeqDesign_Assign_Sequence

! -----------------------------------------------------------------------------

! Set M13mp18 DNA sequence
subroutine SeqDesign_Set_M13mp18(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i, j, count, n_base_scaf, base, across
    integer, parameter :: len_M13   = 7249
    integer, parameter :: len_lamda = 48502
    character(len_M13) :: M13_seq
    character(len_lamda) :: lamda_seq

    ! Check the number of bases in scaffold strands
    n_base_scaf = 0
    do i = 1, dna.n_scaf
        n_base_scaf = n_base_scaf + dna.strand(i).n_base
    end do

    if(n_base_scaf /= dna.n_base_scaf) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | # of nts in scaffold is not consistent.          |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Set M13mp18 or Lamda sequence
    if(n_base_scaf <= len_M13) then
        M13_seq   = SeqDesign_Get_M13mp18(len_M13)
    else if(n_base_scaf > len_M13 .and. n_base_scaf <= len_Lamda) then
        lamda_seq = SeqDesign_Get_Lamda(len_Lamda)
    else
    end if

    ! Set scaffold sequence as full M13 scaffold
    do i = 1, dna.n_strand

        ! Set sequence for scaffold
        if(dna.strand(i).type1 /= "scaf") cycle

        ! Find the starting point (down == -1)
        base   = Mani_Go_Start_Base(dna, i)
        across = dna.top(base).across

        do j = 1, dna.strand(i).n_base

            count = j+para_set_start_scaf-1
            if(count >= len_M13) then
                count = mod(count, len_M13)
                if(count == 0) count = len_M13
            end if

            ! Assign M13mp18 sequence
            if(n_base_scaf <= len_M13) then
                dna.top(base).seq = M13_seq(count:count)
            else if(n_base_scaf > len_M13 .and. n_base_scaf <= len_Lamda) then
                dna.top(base).seq = Lamda_seq(count:count)
            else
                dna.top(base).seq = SeqDesign_Get_Rand_Sequence()
            end if

            ! Set complementary sequence
            if(across /= -1) then
                dna.top(across).seq = SeqDesign_Get_Comp_Sequence(dna.top(base).seq)
            end if

            ! Update base
            if(j /= dna.strand(i).n_base) then
                base   = dna.top(base).up
                across = dna.top(base).across
            end if
        end do
    end do
end subroutine SeqDesign_Set_M13mp18

! -----------------------------------------------------------------------------

! Get M13mp18 DNA sequence
function SeqDesign_Get_M13mp18(len_M13) result(M13_seq)
    integer, intent(in) :: len_M13

    character(len_M13) :: M13_seq, M13_seq_old
    integer :: i, count

    ! Initialize M13 sequence data as "N"
    do i = 1, len_M13
        M13_seq(i:i)     = "N"
        M13_seq_old(i:i) = "N"
    end do

    ! Set scaffold sequence as full M13 scaffold, 7249 nucleotides
    ! M13mp18 [length=7249] [version=09-MAY-2008] [topology=circular] Cloning vector M13mp18
    M13_seq = &
        "AATGCTACTACTATTAGTAGAATTGATGCCACCTTTTCAGCTCGCGCCCC"//&
        "AAATGAAAATATAGCTAAACAGGTTATTGACCATTTGCGAAATGTATCTA"//&   !    1 -  100
        "ATGGTCAAACTAAATCTACTCGTTCGCAGAATTGGGAATCAACTGTTATA"//&
        "TGGAATGAAACTTCCAGACACCGTACTTTAGTTGCATATTTAAAACATGT"//&   !  101 -  200
        "TGAGCTACAGCATTATATTCAGCAATTAAGCTCTAAGCCATCCGCAAAAA"//&
        "TGACCTCTTATCAAAAGGAGCAATTAAAGGTACTCTCTAATCCTGACCTG"//&   !  201 -  300
        "TTGGAGTTTGCTTCCGGTCTGGTTCGCTTTGAAGCTCGAATTAAAACGCG"//&
        "ATATTTGAAGTCTTTCGGGCTTCCTCTTAATCTTTTTGATGCAATCCGCT"//&   !  301 -  400

        "TTGCTTCTGACTATAATAGTCAGGGTAAAGACCTGATTTTTGATTTATGG"//&
        "TCATTCTCGTTTTCTGAACTGTTTAAAGCATTTGAGGGGGATTCAATGAA"//&   !  401 -  500
        "TATTTATGACGATTCCGCAGTATTGGACGCTATCCAGTCTAAACATTTTA"//&
        "CTATTACCCCCTCTGGCAAAACTTCTTTTGCAAAAGCCTCTCGCTATTTT"//&   !  501 -  600
        "GGTTTTTATCGTCGTCTGGTAAACGAGGGTTATGATAGTGTTGCTCTTAC"//&
        "TATGCCTCGTAATTCCTTTTGGCGTTATGTATCTGCATTAGTTGAATGTG"//&   !  601 -  700
        "GTATTCCTAAATCTCAACTGATGAATCTTTCTACCTGTAATAATGTTGTT"//&
        "CCGTTAGTTCGTTTTATTAACGTAGATTTTTCTTCCCAACGTCCTGACTG"//&   !  701 -  800

        "GTATAATGAGCCAGTTCTTAAAATCGCATAAGGTAATTCACAATGATTAA"//&
        "AGTTGAAATTAAACCATCTCAAGCCCAATTTACTACTCGTTCTGGTGTTT"//&   !  801 -  900
        "CTCGTCAGGGCAAGCCTTATTCACTGAATGAGCAGCTTTGTTACGTTGAT"//&
        "TTGGGTAATGAATATCCGGTTCTTGTCAAGATTACTCTTGATGAAGGTCA"//&   !  901 - 1000
        "GCCAGCCTATGCGCCTGGTCTGTACACCGTTCATCTGTCCTCTTTCAAAG"//&
        "TTGGTCAGTTCGGTTCCCTTATGATTGACCGTCTGCGCCTCGTTCCGGCT"//&   ! 1001 - 1100
        "AAGTAACATGGAGCAGGTCGCGGATTTCGACACAATTTATCAGGCGATGA"//&
        "TACAAATCTCCGTTGTACTTTGTTTCGCGCTTGGTATAATCGCTGGGGGT"//&   ! 1101 - 1200

        "CAAAGATGAGTGTTTTAGTGTATTCTTTTGCCTCTTTCGTTTTAGGTTGG"//&
        "TGCCTTCGTAGTGGCATTACGTATTTTACCCGTTTAATGGAAACTTCCTC"//&   ! 1201 - 1300
        "ATGAAAAAGTCTTTAGTCCTCAAAGCCTCTGTAGCCGTTGCTACCCTCGT"//&
        "TCCGATGCTGTCTTTCGCTGCTGAGGGTGACGATCCCGCAAAAGCGGCCT"//&   ! 1301 - 1400
        "TTAACTCCCTGCAAGCCTCAGCGACCGAATATATCGGTTATGCGTGGGCG"//&
        "ATGGTTGTTGTCATTGTCGGCGCAACTATCGGTATCAAGCTGTTTAAGAA"//&   ! 1401 - 1500
        "ATTCACCTCGAAAGCAAGCTGATAAACCGATACAATTAAAGGCTCCTTTT"//&
        "GGAGCCTTTTTTTTGGAGATTTTCAACGTGAAAAAATTATTATTCGCAAT"//&   ! 1501 - 1600

        "TCCTTTAGTTGTTCCTTTCTATTCTCACTCCGCTGAAACTGTTGAAAGTT"//&
        "GTTTAGCAAAATCCCATACAGAAAATTCATTTACTAACGTCTGGAAAGAC"//&   ! 1601 - 1700
        "GACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTGTCTGTGGAATGC"//&
        "TACAGGCGTTGTAGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACAT"//&   ! 1701 - 1800
        "GGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAGGGTGGTGGCTCTGAG"//&
        "GGTGGCGGTTCTGAGGGTGGCGGTTCTGAGGGTGGCGGTACTAAACCTCC"//&   ! 1801 - 1900
        "TGAGTACGGTGATACACCTATTCCGGGCTATACTTATATCAACCCTCTCG"//&
        "ACGGCACTTATCCGCCTGGTACTGAGCAAAACCCCGCTAATCCTAATCCT"//&   ! 1901 - 2000

        "TCTCTTGAGGAGTCTCAGCCTCTTAATACTTTCATGTTTCAGAATAATAG"//&
        "GTTCCGAAATAGGCAGGGGGCATTAACTGTTTATACGGGCACTGTTACTC"//&   ! 2001 - 2100
        "AAGGCACTGACCCCGTTAAAACTTATTACCAGTACACTCCTGTATCATCA"//&
        "AAAGCCATGTATGACGCTTACTGGAACGGTAAATTCAGAGACTGCGCTTT"//&   ! 2101 - 2200
        "CCATTCTGGCTTTAATGAGGATTTATTTGTTTGTGAATATCAAGGCCAAT"//&
        "CGTCTGACCTGCCTCAACCTCCTGTCAATGCTGGCGGCGGCTCTGGTGGT"//&   ! 2201 - 2300
        "GGTTCTGGTGGCGGCTCTGAGGGTGGTGGCTCTGAGGGTGGCGGTTCTGA"//&
        "GGGTGGCGGCTCTGAGGGAGGCGGTTCCGGTGGTGGCTCTGGTTCCGGTG"//&   ! 2301 - 2400

        "ATTTTGATTATGAAAAGATGGCAAACGCTAATAAGGGGGCTATGACCGAA"//&
        "AATGCCGATGAAAACGCGCTACAGTCTGACGCTAAAGGCAAACTTGATTC"//&   ! 2401 - 2500
        "TGTCGCTACTGATTACGGTGCTGCTATCGATGGTTTCATTGGTGACGTTT"//&
        "CCGGCCTTGCTAATGGTAATGGTGCTACTGGTGATTTTGCTGGCTCTAAT"//&   ! 2501 - 2600
        "TCCCAAATGGCTCAAGTCGGTGACGGTGATAATTCACCTTTAATGAATAA"//&
        "TTTCCGTCAATATTTACCTTCCCTCCCTCAATCGGTTGAATGTCGCCCTT"//&   ! 2601 - 2700
        "TTGTCTTTGGCGCTGGTAAACCATATGAATTTTCTATTGATTGTGACAAA"//&
        "ATAAACTTATTCCGTGGTGTCTTTGCGTTTCTTTTATATGTTGCCACCTT"//&   ! 2701 - 2800

        "TATGTATGTATTTTCTACGTTTGCTAACATACTGCGTAATAAGGAGTCTT"//&
        "AATCATGCCAGTTCTTTTGGGTATTCCGTTATTATTGCGTTTCCTCGGTT"//&   ! 2801 - 2900
        "TCCTTCTGGTAACTTTGTTCGGCTATCTGCTTACTTTTCTTAAAAAGGGC"//&
        "TTCGGTAAGATAGCTATTGCTATTTCATTGTTTCTTGCTCTTATTATTGG"//&   ! 2901 - 3000
        "GCTTAACTCAATTCTTGTGGGTTATCTCTCTGATATTAGCGCTCAATTAC"//&
        "CCTCTGACTTTGTTCAGGGTGTTCAGTTAATTCTCCCGTCTAATGCGCTT"//&   ! 3001 - 3100
        "CCCTGTTTTTATGTTATTCTCTCTGTAAAGGCTGCTATTTTCATTTTTGA"//&
        "CGTTAAACAAAAAATCGTTTCTTATTTGGATTGGGATAAATAATATGGCT"//&   ! 3101 - 3200

        "GTTTATTTTGTAACTGGCAAATTAGGCTCTGGAAAGACGCTCGTTAGCGT"//&
        "TGGTAAGATTCAGGATAAAATTGTAGCTGGGTGCAAAATAGCAACTAATC"//&   ! 3201 - 3300
        "TTGATTTAAGGCTTCAAAACCTCCCGCAAGTCGGGAGGTTCGCTAAAACG"//&
        "CCTCGCGTTCTTAGAATACCGGATAAGCCTTCTATATCTGATTTGCTTGC"//&   ! 3301 - 3400
        "TATTGGGCGCGGTAATGATTCCTACGATGAAAATAAAAACGGCTTGCTTG"//&
        "TTCTCGATGAGTGCGGTACTTGGTTTAATACCCGTTCTTGGAATGATAAG"//&   ! 3401 - 3500
        "GAAAGACAGCCGATTATTGATTGGTTTCTACATGCTCGTAAATTAGGATG"//&
        "GGATATTATTTTTCTTGTTCAGGACTTATCTATTGTTGATAAACAGGCGC"//&   ! 3501 - 3600

        "GTTCTGCATTAGCTGAACATGTTGTTTATTGTCGTCGTCTGGACAGAATT"//&
        "ACTTTACCTTTTGTCGGTACTTTATATTCTCTTATTACTGGCTCGAAAAT"//&   ! 3601 - 3700
        "GCCTCTGCCTAAATTACATGTTGGCGTTGTTAAATATGGCGATTCTCAAT"//&
        "TAAGCCCTACTGTTGAGCGTTGGCTTTATACTGGTAAGAATTTGTATAAC"//&   ! 3701 - 3800
        "GCATATGATACTAAACAGGCTTTTTCTAGTAATTATGATTCCGGTGTTTA"//&
        "TTCTTATTTAACGCCTTATTTATCACACGGTCGGTATTTCAAACCATTAA"//&   ! 3801 - 3900
        "ATTTAGGTCAGAAGATGAAATTAACTAAAATATATTTGAAAAAGTTTTCT"//&
        "CGCGTTCTTTGTCTTGCGATTGGATTTGCATCAGCATTTACATATAGTTA"//&   ! 3901 - 4000

        "TATAACCCAACCTAAGCCGGAGGTTAAAAAGGTAGTCTCTCAGACCTATG"//&
        "ATTTTGATAAATTCACTATTGACTCTTCTCAGCGTCTTAATCTAAGCTAT"//&   ! 4001 - 4100
        "CGCTATGTTTTCAAGGATTCTAAGGGAAAATTAATTAATAGCGACGATTT"//&
        "ACAGAAGCAAGGTTATTCACTCACATATATTGATTTATGTACTGTTTCCA"//&   ! 4101 - 4200
        "TTAAAAAAGGTAATTCAAATGAAATTGTTAAATGTAATTAATTTTGTTTT"//&
        "CTTGATGTTTGTTTCATCATCTTCTTTTGCTCAGGTAATTGAAATGAATA"//&   ! 4201 - 4300
        "ATTCGCCTCTGCGCGATTTTGTAACTTGGTATTCAAAGCAATCAGGCGAA"//&
        "TCCGTTATTGTTTCTCCCGATGTAAAAGGTACTGTTACTGTATATTCATC"//&   ! 4301 - 4400

        "TGACGTTAAACCTGAAAATCTACGCAATTTCTTTATTTCTGTTTTACGTG"//&
        "CAAATAATTTTGATATGGTAGGTTCTAACCCTTCCATTATTCAGAAGTAT"//&   ! 4401 - 4500
        "AATCCAAACAATCAGGATTATATTGATGAATTGCCATCATCTGATAATCA"//&
        "GGAATATGATGATAATTCCGCTCCTTCTGGTGGTTTCTTTGTTCCGCAAA"//&   ! 4501 - 4600
        "ATGATAATGTTACTCAAACTTTTAAAATTAATAACGTTCGGGCAAAGGAT"//&
        "TTAATACGAGTTGTCGAATTGTTTGTAAAGTCTAATACTTCTAAATCCTC"//&   ! 4601 - 4700
        "AAATGTATTATCTATTGACGGCTCTAATCTATTAGTTGTTAGTGCTCCTA"//&
        "AAGATATTTTAGATAACCTTCCTCAATTCCTTTCAACTGTTGATTTGCCA"//&   ! 4701 - 4800

        "ACTGACCAGATATTGATTGAGGGTTTGATATTTGAGGTTCAGCAAGGTGA"//&
        "TGCTTTAGATTTTTCATTTGCTGCTGGCTCTCAGCGTGGCACTGTTGCAG"//&   ! 4801 - 4900
        "GCGGTGTTAATACTGACCGCCTCACCTCTGTTTTATCTTCTGCTGGTGGT"//&
        "TCGTTCGGTATTTTTAATGGCGATGTTTTAGGGCTATCAGTTCGCGCATT"//&   ! 4901 - 5000
        "AAAGACTAATAGCCATTCAAAAATATTGTCTGTGCCACGTATTCTTACGC"//&
        "TTTCAGGTCAGAAGGGTTCTATCTCTGTTGGCCAGAATGTCCCTTTTATT"//&   ! 5001 - 5100
        "ACTGGTCGTGTGACTGGTGAATCTGCCAATGTAAATAATCCATTTCAGAC"//&
        "GATTGAGCGTCAAAATGTAGGTATTTCCATGAGCGTTTTTCCTGTTGCAA"//&   ! 5101 - 5200

        "TGGCTGGCGGTAATATTGTTCTGGATATTACCAGCAAGGCCGATAGTTTG"//&
        "AGTTCTTCTACTCAGGCAAGTGATGTTATTACTAATCAAAGAAGTATTGC"//&   ! 5201 - 5300
        "TACAACGGTTAATTTGCGTGATGGACAGACTCTTTTACTCGGTGGCCTCA"//&
        "CTGATTATAAAAACACTTCTCAGGATTCTGGCGTACCGTTCCTGTCTAAA"//&   ! 5301 - 5400
        "ATCCCTTTAATCGGCCTCCTGTTTAGCTCCCGCTCTGATTCTAACGAGGA"//&
        "AAGCACGTTATACGTGCTCGTCAAAGCAACCATAGTACGCGCCCTGTAGC"//&   ! 5401 - 5500
        "GGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTAC"//&
        "ACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTC"//&   ! 5501 - 5600

        "TCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCT"//&
        "TTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGA"//&   ! 5601 - 5700
        "TTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTC"//&
        "GCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAA"//&   ! 5701 - 5800
        "ACTGGAACAACACTCAACCCTATCTCGGGCTATTCTTTTGATTTATAAGG"//&
        "GATTTTGCCGATTTCGGAACCACCATCAAACAGGATTTTCGCCTGCTGGG"//&   ! 5801 - 5900
        "GCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGA"//&
        "AGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTG"//&   ! 5901 - 6000

        "GCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAAT"//&
        "GCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAAC"//&   ! 6001 - 6100
        "GCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTT"//&
        "TATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTC"//&   ! 6101 - 6200
        "ACACAGGAAACAGCTATGACCATGATTACGAATTCGAGCTCGGTACCCGG"//&
        "GGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCG"//&   ! 6201 - 6300
        "TTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGC"//&
        "CTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCG"//&   ! 6301 - 6400

        "CACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCT"//&
        "TTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGC"//&   ! 6401 - 6500
        "GATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCA"//&
        "CGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCA"//&   ! 6501 - 6600
        "ATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACA"//&
        "TTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTT"//&   ! 6601 - 6700
        "TGATGGCGTTCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAA"//&
        "TGCGAATTTTAACAAAATATTAACGTTTACAATTTAAATATTTGCTTATA"//&   ! 6701 - 6800

        "CAATCTTCCTGTTTTTGGGGCTTTTCTGATTATCAACCGGGGTACATATG"//&
        "ATTGACATGCTAGTTTTACGATTACCGTTCATCGATTCTCTTGTTTGCTC"//&   ! 6801 - 6900
        "CAGACTCTCAGGCAATGACCTGATAGCCTTTGTAGATCTCTCAAAAATAG"//&
        "CTACCCTCTCCGGCATTAATTTATCAGCTAGAACGGTTGAATATCATATT"//&   ! 6901 - 7000
        "GATGGTGATTTGACTGTCTCCGGCCTTTCTCACCCTTTTGAATCTTTACC"//&
        "TACACATTACTCAGGCATTGCATTTAAAATATATGAGGGTTCTAAAAATT"//&   ! 7001 - 7100
        "TTTATCCTTGCGTTGAAATAAAGGCTTCTCCCGCAAAAGTATTACAGGGT"//&
        "CATAATGTTTTTGGTACAACCGATTTAGCTTTATGCTCTGAGGCTTTATT"//&   ! 7101 - 7200
        "GCTTAATTTTGCTAATTCTTTGCCTTGCCTGTATGATTTATTGGATGTT"       ! 7200 - 7249

    M13_seq_old(1:3652)    = "AATGC"
    M13_seq_old(3653:7249) = "TTTAC"

    ! Compare M13 with old version
    !count = 0
    !do i = 1, len_M13
    !    if(M13_seq(i:i) /= M13_seq_old(i:i)) then
    !        count = count + 1
    !        write(0, "(a)"), trim(adjustl(Int2Str(count)))//" th, pos : "&
    !            //trim(adjustl(Int2Str(i)))//", current : "&
    !            //M13_seq(i:i)//", previous : "//M13_seq_old(i:i)
    !    end if
    !end do
    !stop

    ! Check M13 sequence data
    do i = 1, len_M13
        if(M13_seq(i:i) == "N") then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | Wrong sequences in M13mp18.                      |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if
    end do
end function SeqDesign_Get_M13mp18

! -----------------------------------------------------------------------------

! Get M13mp18 DNA sequence
function SeqDesign_Get_Lamda(len_Lamda) result(Lamda_seq)
    integer, intent(in) :: len_Lamda

    character(len_Lamda) :: Lamda_seq
    integer :: i, count

    ! Initialize Lamda sequence data as "N"
    do i = 1, len_Lamda
        Lamda_seq(i:i) = "N"
    end do

    ! Set scaffold sequence as full Lamda scaffold, 48502 nucleotides
    ! https://www.neb.com/tools-and-resources/interactive-tools/dna-sequences-and-maps-tool
    Lamda_seq = &
        "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG"//&
        "TCATAACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGGC"//&
        "TTTTTGGCCTCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCA"//&
        "GCTGGCTGACATTTTCGGTGCGAGTATCCGTACCATTCAGAACTGGCAGGAACAGGGAATGCCCGTTCTG"//&
        "CGAGGCGGTGGCAAGGGTAATGAGGTGCTTTATGACTCTGCCGCCGTCATAAAATGGTATGCCGAAAGGG"//&
        "ATGCTGAAATTGAGAACGAAAAGCTGCGCCGGGAGGTTGAAGAACTGCGGCAGGCCAGCGAGGCAGATCT"//&
        "CCAGCCAGGAACTATTGAGTACGAACGCCATCGACTTACGCGTGCGCAGGCCGACGCACAGGAACTGAAG"//&
        "AATGCCAGAGACTCCGCTGAAGTGGTGGAAACCGCATTCTGTACTTTCGTGCTGTCGCGGATCGCAGGTG"//&
        "AAATTGCCAGTATTCTCGACGGGCTCCCCCTGTCGGTGCAGCGGCGTTTTCCGGAACTGGAAAACCGACA"//&
        "TGTTGATTTCCTGAAACGGGATATCATCAAAGCCATGAACAAAGCAGCCGCGCTGGATGAACTGATACCG"//&
        "GGGTTGCTGAGTGAATATATCGAACAGTCAGGTTAACAGGCTGCGGCATTTTGTCCGCGCCGGGCTTCGC"//&
        "TCACTGTTCAGGCCGGAGCCACAGACCGCCGTTGAATGGGCGGATGCTAATTACTATCTCCCGAAAGAAT"//&
        "CCGCATACCAGGAAGGGCGCTGGGAAACACTGCCCTTTCAGCGGGCCATCATGAATGCGATGGGCAGCGA"//&
        "CTACATCCGTGAGGTGAATGTGGTGAAGTCTGCCCGTGTCGGTTATTCCAAAATGCTGCTGGGTGTTTAT"//&
        "GCCTACTTTATAGAGCATAAGCAGCGCAACACCCTTATCTGGTTGCCGACGGATGGTGATGCCGAGAACT"//&
        "TTATGAAAACCCACGTTGAGCCGACTATTCGTGATATTCCGTCGCTGCTGGCGCTGGCCCCGTGGTATGG"//&
        "CAAAAAGCACCGGGATAACACGCTCACCATGAAGCGTTTCACTAATGGGCGTGGCTTCTGGTGCCTGGGC"//&
        "GGTAAAGCGGCAAAAAACTACCGTGAAAAGTCGGTGGATGTGGCGGGTTATGATGAACTTGCTGCTTTTG"//&
        "ATGATGATATTGAACAGGAAGGCTCTCCGACGTTCCTGGGTGACAAGCGTATTGAAGGCTCGGTCTGGCC"//&
        "AAAGTCCATCCGTGGCTCCACGCCAAAAGTGAGAGGCACCTGTCAGATTGAGCGTGCAGCCAGTGAATCC"//&
        "CCGCATTTTATGCGTTTTCATGTTGCCTGCCCGCATTGCGGGGAGGAGCAGTATCTTAAATTTGGCGACA"//&
        "AAGAGACGCCGTTTGGCCTCAAATGGACGCCGGATGACCCCTCCAGCGTGTTTTATCTCTGCGAGCATAA"//&
        "TGCCTGCGTCATCCGCCAGCAGGAGCTGGACTTTACTGATGCCCGTTATATCTGCGAAAAGACCGGGATC"//&
        "TGGACCCGTGATGGCATTCTCTGGTTTTCGTCATCCGGTGAAGAGATTGAGCCACCTGACAGTGTGACCT"//&
        "TTCACATCTGGACAGCGTACAGCCCGTTCACCACCTGGGTGCAGATTGTCAAAGACTGGATGAAAACGAA"//&
        "AGGGGATACGGGAAAACGTAAAACCTTCGTAAACACCACGCTCGGTGAGACGTGGGAGGCGAAAATTGGC"//&
        "GAACGTCCGGATGCTGAAGTGATGGCAGAGCGGAAAGAGCATTATTCAGCGCCCGTTCCTGACCGTGTGG"//&
        "CTTACCTGACCGCCGGTATCGACTCCCAGCTGGACCGCTACGAAATGCGCGTATGGGGATGGGGGCCGGG"//&
        "TGAGGAAAGCTGGCTGATTGACCGGCAGATTATTATGGGCCGCCACGACGATGAACAGACGCTGCTGCGT"//&
        "GTGGATGAGGCCATCAATAAAACCTATACCCGCCGGAATGGTGCAGAAATGTCGATATCCCGTATCTGCT"//&
        "GGGATACTGGCGGGATTGACCCGACCATTGTGTATGAACGCTCGAAAAAACATGGGCTGTTCCGGGTGAT"//&
        "CCCCATTAAAGGGGCATCCGTCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGG"//&
        "GTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAACCGCTTCACACTGACGCCGG"//&
        "AAGGGGATGAACCGCTTCCCGGTGCCGTTCACTTCCCGAATAACCCGGATATTTTTGATCTGACCGAAGC"//&
        "GCAGCAGCTGACTGCTGAAGAGCAGGTCGAAAAATGGGTGGATGGCAGGAAAAAAATACTGTGGGACAGC"//&
        "AAAAAGCGACGCAATGAGGCACTCGACTGCTTCGTTTATGCGCTGGCGGCGCTGCGCATCAGTATTTCCC"//&
        "GCTGGCAGCTGGATCTCAGTGCGCTGCTGGCGAGCCTGCAGGAAGAGGATGGTGCAGCAACCAACAAGAA"//&
        "AACACTGGCAGATTACGCCCGTGCCTTATCCGGAGAGGATGAATGACGCGACAGGAAGAACTTGCCGCTG"//&
        "CCCGTGCGGCACTGCATGACCTGATGACAGGTAAACGGGTGGCAACAGTACAGAAAGACGGACGAAGGGT"//&
        "GGAGTTTACGGCCACTTCCGTGTCTGACCTGAAAAAATATATTGCAGAGCTGGAAGTGCAGACCGGCATG"//&
        "ACACAGCGACGCAGGGGACCTGCAGGATTTTATGTATGAAAACGCCCACCATTCCCACCCTTCTGGGGCC"//&
        "GGACGGCATGACATCGCTGCGCGAATATGCCGGTTATCACGGCGGTGGCAGCGGATTTGGAGGGCAGTTG"//&
        "CGGTCGTGGAACCCACCGAGTGAAAGTGTGGATGCAGCCCTGTTGCCCAACTTTACCCGTGGCAATGCCC"//&
        "GCGCAGACGATCTGGTACGCAATAACGGCTATGCCGCCAACGCCATCCAGCTGCATCAGGATCATATCGT"//&
        "CGGGTCTTTTTTCCGGCTCAGTCATCGCCCAAGCTGGCGCTATCTGGGCATCGGGGAGGAAGAAGCCCGT"//&
        "GCCTTTTCCCGCGAGGTTGAAGCGGCATGGAAAGAGTTTGCCGAGGATGACTGCTGCTGCATTGACGTTG"//&
        "AGCGAAAACGCACGTTTACCATGATGATTCGGGAAGGTGTGGCCATGCACGCCTTTAACGGTGAACTGTT"//&
        "CGTTCAGGCCACCTGGGATACCAGTTCGTCGCGGCTTTTCCGGACACAGTTCCGGATGGTCAGCCCGAAG"//&
        "CGCATCAGCAACCCGAACAATACCGGCGACAGCCGGAACTGCCGTGCCGGTGTGCAGATTAATGACAGCG"//&
        "GTGCGGCGCTGGGATATTACGTCAGCGAGGACGGGTATCCTGGCTGGATGCCGCAGAAATGGACATGGAT"//&
        "ACCCCGTGAGTTACCCGGCGGGCGCGCCTCGTTCATTCACGTTTTTGAACCCGTGGAGGACGGGCAGACT"//&
        "CGCGGTGCAAATGTGTTTTACAGCGTGATGGAGCAGATGAAGATGCTCGACACGCTGCAGAACACGCAGC"//&
        "TGCAGAGCGCCATTGTGAAGGCGATGTATGCCGCCACCATTGAGAGTGAGCTGGATACGCAGTCAGCGAT"//&
        "GGATTTTATTCTGGGCGCGAACAGTCAGGAGCAGCGGGAAAGGCTGACCGGCTGGATTGGTGAAATTGCC"//&
        "GCGTATTACGCCGCAGCGCCGGTCCGGCTGGGAGGCGCAAAAGTACCGCACCTGATGCCGGGTGACTCAC"//&
        "TGAACCTGCAGACGGCTCAGGATACGGATAACGGCTACTCCGTGTTTGAGCAGTCACTGCTGCGGTATAT"//&
        "CGCTGCCGGGCTGGGTGTCTCGTATGAGCAGCTTTCCCGGAATTACGCCCAGATGAGCTACTCCACGGCA"//&
        "CGGGCCAGTGCGAACGAGTCGTGGGCGTACTTTATGGGGCGGCGAAAATTCGTCGCATCCCGTCAGGCGA"//&
        "GCCAGATGTTTCTGTGCTGGCTGGAAGAGGCCATCGTTCGCCGCGTGGTGACGTTACCTTCAAAAGCGCG"//&
        "CTTCAGTTTTCAGGAAGCCCGCAGTGCCTGGGGGAACTGCGACTGGATAGGCTCCGGTCGTATGGCCATC"//&
        "GATGGTCTGAAAGAAGTTCAGGAAGCGGTGATGCTGATAGAAGCCGGACTGAGTACCTACGAGAAAGAGT"//&
        "GCGCAAAACGCGGTGACGACTATCAGGAAATTTTTGCCCAGCAGGTCCGTGAAACGATGGAGCGCCGTGC"//&
        "AGCCGGTCTTAAACCGCCCGCCTGGGCGGCTGCAGCATTTGAATCCGGGCTGCGACAATCAACAGAGGAG"//&
        "GAGAAGAGTGACAGCAGAGCTGCGTAATCTCCCGCATATTGCCAGCATGGCCTTTAATGAGCCGCTGATG"//&
        "CTTGAACCCGCCTATGCGCGGGTTTTCTTTTGTGCGCTTGCAGGCCAGCTTGGGATCAGCAGCCTGACGG"//&
        "ATGCGGTGTCCGGCGACAGCCTGACTGCCCAGGAGGCACTCGCGACGCTGGCATTATCCGGTGATGATGA"//&
        "CGGACCACGACAGGCCCGCAGTTATCAGGTCATGAACGGCATCGCCGTGCTGCCGGTGTCCGGCACGCTG"//&
        "GTCAGCCGGACGCGGGCGCTGCAGCCGTACTCGGGGATGACCGGTTACAACGGCATTATCGCCCGTCTGC"//&
        "AACAGGCTGCCAGCGATCCGATGGTGGACGGCATTCTGCTCGATATGGACACGCCCGGCGGGATGGTGGC"//&
        "GGGGGCATTTGACTGCGCTGACATCATCGCCCGTGTGCGTGACATAAAACCGGTATGGGCGCTTGCCAAC"//&
        "GACATGAACTGCAGTGCAGGTCAGTTGCTTGCCAGTGCCGCCTCCCGGCGTCTGGTCACGCAGACCGCCC"//&
        "GGACAGGCTCCATCGGCGTCATGATGGCTCACAGTAATTACGGTGCTGCGCTGGAGAAACAGGGTGTGGA"//&
        "AATCACGCTGATTTACAGCGGCAGCCATAAGGTGGATGGCAACCCCTACAGCCATCTTCCGGATGACGTC"//&
        "CGGGAGACACTGCAGTCCCGGATGGACGCAACCCGCCAGATGTTTGCGCAGAAGGTGTCGGCATATACCG"//&
        "GCCTGTCCGTGCAGGTTGTGCTGGATACCGAGGCTGCAGTGTACAGCGGTCAGGAGGCCATTGATGCCGG"//&
        "ACTGGCTGATGAACTTGTTAACAGCACCGATGCGATCACCGTCATGCGTGATGCACTGGATGCACGTAAA"//&
        "TCCCGTCTCTCAGGAGGGCGAATGACCAAAGAGACTCAATCAACAACTGTTTCAGCCACTGCTTCGCAGG"//&
        "CTGACGTTACTGACGTGGTGCCAGCGACGGAGGGCGAGAACGCCAGCGCGGCGCAGCCGGACGTGAACGC"//&
        "GCAGATCACCGCAGCGGTTGCGGCAGAAAACAGCCGCATTATGGGGATCCTCAACTGTGAGGAGGCTCAC"//&
        "GGACGCGAAGAACAGGCACGCGTGCTGGCAGAAACCCCCGGTATGACCGTGAAAACGGCCCGCCGCATTC"//&
        "TGGCCGCAGCACCACAGAGTGCACAGGCGCGCAGTGACACTGCGCTGGATCGTCTGATGCAGGGGGCACC"//&
        "GGCACCGCTGGCTGCAGGTAACCCGGCATCTGATGCCGTTAACGATTTGCTGAACACACCAGTGTAAGGG"//&
        "ATGTTTATGACGAGCAAAGAAACCTTTACCCATTACCAGCCGCAGGGCAACAGTGACCCGGCTCATACCG"//&
        "CAACCGCGCCCGGCGGATTGAGTGCGAAAGCGCCTGCAATGACCCCGCTGATGCTGGACACCTCCAGCCG"//&
        "TAAGCTGGTTGCGTGGGATGGCACCACCGACGGTGCTGCCGTTGGCATTCTTGCGGTTGCTGCTGACCAG"//&
        "ACCAGCACCACGCTGACGTTCTACAAGTCCGGCACGTTCCGTTATGAGGATGTGCTCTGGCCGGAGGCTG"//&
        "CCAGCGACGAGACGAAAAAACGGACCGCGTTTGCCGGAACGGCAATCAGCATCGTTTAACTTTACCCTTC"//&
        "ATCACTAAAGGCCGCCTGTGCGGCTTTTTTTACGGGATTTTTTTATGTCGATGTACACAACCGCCCAACT"//&
        "GCTGGCGGCAAATGAGCAGAAATTTAAGTTTGATCCGCTGTTTCTGCGTCTCTTTTTCCGTGAGAGCTAT"//&
        "CCCTTCACCACGGAGAAAGTCTATCTCTCACAAATTCCGGGACTGGTAAACATGGCGCTGTACGTTTCGC"//&
        "CGATTGTTTCCGGTGAGGTTATCCGTTCCCGTGGCGGCTCCACCTCTGAATTTACGCCGGGATATGTCAA"//&
        "GCCGAAGCATGAAGTGAATCCGCAGATGACCCTGCGTCGCCTGCCGGATGAAGATCCGCAGAATCTGGCG"//&
        "GACCCGGCTTACCGCCGCCGTCGCATCATCATGCAGAACATGCGTGACGAAGAGCTGGCCATTGCTCAGG"//&
        "TCGAAGAGATGCAGGCAGTTTCTGCCGTGCTTAAGGGCAAATACACCATGACCGGTGAAGCCTTCGATCC"//&
        "GGTTGAGGTGGATATGGGCCGCAGTGAGGAGAATAACATCACGCAGTCCGGCGGCACGGAGTGGAGCAAG"//&
        "CGTGACAAGTCCACGTATGACCCGACCGACGATATCGAAGCCTACGCGCTGAACGCCAGCGGTGTGGTGA"//&
        "ATATCATCGTGTTCGATCCGAAAGGCTGGGCGCTGTTCCGTTCCTTCAAAGCCGTCAAGGAGAAGCTGGA"//&
        "TACCCGTCGTGGCTCTAATTCCGAGCTGGAGACAGCGGTGAAAGACCTGGGCAAAGCGGTGTCCTATAAG"//&
        "GGGATGTATGGCGATGTGGCCATCGTCGTGTATTCCGGACAGTACGTGGAAAACGGCGTCAAAAAGAACT"//&
        "TCCTGCCGGACAACACGATGGTGCTGGGGAACACTCAGGCACGCGGTCTGCGCACCTATGGCTGCATTCA"//&
        "GGATGCGGACGCACAGCGCGAAGGCATTAACGCCTCTGCCCGTTACCCGAAAAACTGGGTGACCACCGGC"//&
        "GATCCGGCGCGTGAGTTCACCATGATTCAGTCAGCACCGCTGATGCTGCTGGCTGACCCTGATGAGTTCG"//&
        "TGTCCGTACAACTGGCGTAATCATGGCCCTTCGGGGCCATTGTTTCTCTGTGGAGGAGTCCATGACGAAA"//&
        "GATGAACTGATTGCCCGTCTCCGCTCGCTGGGTGAACAACTGAACCGTGATGTCAGCCTGACGGGGACGA"//&
        "AAGAAGAACTGGCGCTCCGTGTGGCAGAGCTGAAAGAGGAGCTTGATGACACGGATGAAACTGCCGGTCA"//&
        "GGACACCCCTCTCAGCCGGGAAAATGTGCTGACCGGACATGAAAATGAGGTGGGATCAGCGCAGCCGGAT"//&
        "ACCGTGATTCTGGATACGTCTGAACTGGTCACGGTCGTGGCACTGGTGAAGCTGCATACTGATGCACTTC"//&
        "ACGCCACGCGGGATGAACCTGTGGCATTTGTGCTGCCGGGAACGGCGTTTCGTGTCTCTGCCGGTGTGGC"//&
        "AGCCGAAATGACAGAGCGCGGCCTGGCCAGAATGCAATAACGGGAGGCGCTGTGGCTGATTTCGATAACC"//&
        "TGTTCGATGCTGCCATTGCCCGCGCCGATGAAACGATACGCGGGTACATGGGAACGTCAGCCACCATTAC"//&
        "ATCCGGTGAGCAGTCAGGTGCGGTGATACGTGGTGTTTTTGATGACCCTGAAAATATCAGCTATGCCGGA"//&
        "CAGGGCGTGCGCGTTGAAGGCTCCAGCCCGTCCCTGTTTGTCCGGACTGATGAGGTGCGGCAGCTGCGGC"//&
        "GTGGAGACACGCTGACCATCGGTGAGGAAAATTTCTGGGTAGATCGGGTTTCGCCGGATGATGGCGGAAG"//&
        "TTGTCATCTCTGGCTTGGACGGGGCGTACCGCCTGCCGTTAACCGTCGCCGCTGAAAGGGGGATGTATGG"//&
        "CCATAAAAGGTCTTGAGCAGGCCGTTGAAAACCTCAGCCGTATCAGCAAAACGGCGGTGCCTGGTGCCGC"//&
        "CGCAATGGCCATTAACCGCGTTGCTTCATCCGCGATATCGCAGTCGGCGTCACAGGTTGCCCGTGAGACA"//&
        "AAGGTACGCCGGAAACTGGTAAAGGAAAGGGCCAGGCTGAAAAGGGCCACGGTCAAAAATCCGCAGGCCA"//&
        "GAATCAAAGTTAACCGGGGGGATTTGCCCGTAATCAAGCTGGGTAATGCGCGGGTTGTCCTTTCGCGCCG"//&
        "CAGGCGTCGTAAAAAGGGGCAGCGTTCATCCCTGAAAGGTGGCGGCAGCGTGCTTGTGGTGGGTAACCGT"//&
        "CGTATTCCCGGCGCGTTTATTCAGCAACTGAAAAATGGCCGGTGGCATGTCATGCAGCGTGTGGCTGGGA"//&
        "AAAACCGTTACCCCATTGATGTGGTGAAAATCCCGATGGCGGTGCCGCTGACCACGGCGTTTAAACAAAA"//&
        "TATTGAGCGGATACGGCGTGAACGTCTTCCGAAAGAGCTGGGCTATGCGCTGCAGCATCAACTGAGGATG"//&
        "GTAATAAAGCGATGAAACATACTGAACTCCGTGCAGCCGTACTGGATGCACTGGAGAAGCATGACACCGG"//&
        "GGCGACGTTTTTTGATGGTCGCCCCGCTGTTTTTGATGAGGCGGATTTTCCGGCAGTTGCCGTTTATCTC"//&
        "ACCGGCGCTGAATACACGGGCGAAGAGCTGGACAGCGATACCTGGCAGGCGGAGCTGCATATCGAAGTTT"//&
        "TCCTGCCTGCTCAGGTGCCGGATTCAGAGCTGGATGCGTGGATGGAGTCCCGGATTTATCCGGTGATGAG"//&
        "CGATATCCCGGCACTGTCAGATTTGATCACCAGTATGGTGGCCAGCGGCTATGACTACCGGCGCGACGAT"//&
        "GATGCGGGCTTGTGGAGTTCAGCCGATCTGACTTATGTCATTACCTATGAAATGTGAGGACGCTATGCCT"//&
        "GTACCAAATCCTACAATGCCGGTGAAAGGTGCCGGGACCACCCTGTGGGTTTATAAGGGGAGCGGTGACC"//&
        "CTTACGCGAATCCGCTTTCAGACGTTGACTGGTCGCGTCTGGCAAAAGTTAAAGACCTGACGCCCGGCGA"//&
        "ACTGACCGCTGAGTCCTATGACGACAGCTATCTCGATGATGAAGATGCAGACTGGACTGCGACCGGGCAG"//&
        "GGGCAGAAATCTGCCGGAGATACCAGCTTCACGCTGGCGTGGATGCCCGGAGAGCAGGGGCAGCAGGCGC"//&
        "TGCTGGCGTGGTTTAATGAAGGCGATACCCGTGCCTATAAAATCCGCTTCCCGAACGGCACGGTCGATGT"//&
        "GTTCCGTGGCTGGGTCAGCAGTATCGGTAAGGCGGTGACGGCGAAGGAAGTGATCACCCGCACGGTGAAA"//&
        "GTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCGCAGCACGGTAACAGCGGCAACCGGCATGACCG"//&
        "TGACGCCTGCCAGCACCTCGGTGGTGAAAGGGCAGAGCACCACGCTGACCGTGGCCTTCCAGCCGGAGGG"//&
        "CGTAACCGACAAGAGCTTTCGTGCGGTGTCTGCGGATAAAACAAAAGCCACCGTGTCGGTCAGTGGTATG"//&
        "ACCATCACCGTGAACGGCGTTGCTGCAGGCAAGGTCAACATTCCGGTTGTATCCGGTAATGGTGAGTTTG"//&
        "CTGCGGTTGCAGAAATTACCGTCACCGCCAGTTAATCCGGAGAGTCAGCGATGTTCCTGAAAACCGAATC"//&
        "ATTTGAACATAACGGTGTGACCGTCACGCTTTCTGAACTGTCAGCCCTGCAGCGCATTGAGCATCTCGCC"//&
        "CTGATGAAACGGCAGGCAGAACAGGCGGAGTCAGACAGCAACCGGAAGTTTACTGTGGAAGACGCCATCA"//&
        "GAACCGGCGCGTTTCTGGTGGCGATGTCCCTGTGGCATAACCATCCGCAGAAGACGCAGATGCCGTCCAT"//&
        "GAATGAAGCCGTTAAACAGATTGAGCAGGAAGTGCTTACCACCTGGCCCACGGAGGCAATTTCTCATGCT"//&
        "GAAAACGTGGTGTACCGGCTGTCTGGTATGTATGAGTTTGTGGTGAATAATGCCCCTGAACAGACAGAGG"//&
        "ACGCCGGGCCCGCAGAGCCTGTTTCTGCGGGAAAGTGTTCGACGGTGAGCTGAGTTTTGCCCTGAAACTG"//&
        "GCGCGTGAGATGGGGCGACCCGACTGGCGTGCCATGCTTGCCGGGATGTCATCCACGGAGTATGCCGACT"//&
        "GGCACCGCTTTTACAGTACCCATTATTTTCATGATGTTCTGCTGGATATGCACTTTTCCGGGCTGACGTA"//&
        "CACCGTGCTCAGCCTGTTTTTCAGCGATCCGGATATGCATCCGCTGGATTTCAGTCTGCTGAACCGGCGC"//&
        "GAGGCTGACGAAGAGCCTGAAGATGATGTGCTGATGCAGAAAGCGGCAGGGCTTGCCGGAGGTGTCCGCT"//&
        "TTGGCCCGGACGGGAATGAAGTTATCCCCGCTTCCCCGGATGTGGCGGACATGACGGAGGATGACGTAAT"//&
        "GCTGATGACAGTATCAGAAGGGATCGCAGGAGGAGTCCGGTATGGCTGAACCGGTAGGCGATCTGGTCGT"//&
        "TGATTTGAGTCTGGATGCGGCCAGATTTGACGAGCAGATGGCCAGAGTCAGGCGTCATTTTTCTGGTACG"//&
        "GAAAGTGATGCGAAAAAAACAGCGGCAGTCGTTGAACAGTCGCTGAGCCGACAGGCGCTGGCTGCACAGA"//&
        "AAGCGGGGATTTCCGTCGGGCAGTATAAAGCCGCCATGCGTATGCTGCCTGCACAGTTCACCGACGTGGC"//&
        "CACGCAGCTTGCAGGCGGGCAAAGTCCGTGGCTGATCCTGCTGCAACAGGGGGGGCAGGTGAAGGACTCC"//&
        "TTCGGCGGGATGATCCCCATGTTCAGGGGGCTTGCCGGTGCGATCACCCTGCCGATGGTGGGGGCCACCT"//&
        "CGCTGGCGGTGGCGACCGGTGCGCTGGCGTATGCCTGGTATCAGGGCAACTCAACCCTGTCCGATTTCAA"//&
        "CAAAACGCTGGTCCTTTCCGGCAATCAGGCGGGACTGACGGCAGATCGTATGCTGGTCCTGTCCAGAGCC"//&
        "GGGCAGGCGGCAGGGCTGACGTTTAACCAGACCAGCGAGTCACTCAGCGCACTGGTTAAGGCGGGGGTAA"//&
        "GCGGTGAGGCTCAGATTGCGTCCATCAGCCAGAGTGTGGCGCGTTTCTCCTCTGCATCCGGCGTGGAGGT"//&
        "GGACAAGGTCGCTGAAGCCTTCGGGAAGCTGACCACAGACCCGACGTCGGGGCTGACGGCGATGGCTCGC"//&
        "CAGTTCCATAACGTGTCGGCGGAGCAGATTGCGTATGTTGCTCAGTTGCAGCGTTCCGGCGATGAAGCCG"//&
        "GGGCATTGCAGGCGGCGAACGAGGCCGCAACGAAAGGGTTTGATGACCAGACCCGCCGCCTGAAAGAGAA"//&
        "CATGGGCACGCTGGAGACCTGGGCAGACAGGACTGCGCGGGCATTCAAATCCATGTGGGATGCGGTGCTG"//&
        "GATATTGGTCGTCCTGATACCGCGCAGGAGATGCTGATTAAGGCAGAGGCTGCGTATAAGAAAGCAGACG"//&
        "ACATCTGGAATCTGCGCAAGGATGATTATTTTGTTAACGATGAAGCGCGGGCGCGTTACTGGGATGATCG"//&
        "TGAAAAGGCCCGTCTTGCGCTTGAAGCCGCCCGAAAGAAGGCTGAGCAGCAGACTCAACAGGACAAAAAT"//&
        "GCGCAGCAGCAGAGCGATACCGAAGCGTCACGGCTGAAATATACCGAAGAGGCGCAGAAGGCTTACGAAC"//&
        "GGCTGCAGACGCCGCTGGAGAAATATACCGCCCGTCAGGAAGAACTGAACAAGGCACTGAAAGACGGGAA"//&
        "AATCCTGCAGGCGGATTACAACACGCTGATGGCGGCGGCGAAAAAGGATTATGAAGCGACGCTGAAAAAG"//&
        "CCGAAACAGTCCAGCGTGAAGGTGTCTGCGGGCGATCGTCAGGAAGACAGTGCTCATGCTGCCCTGCTGA"//&
        "CGCTTCAGGCAGAACTCCGGACGCTGGAGAAGCATGCCGGAGCAAATGAGAAAATCAGCCAGCAGCGCCG"//&
        "GGATTTGTGGAAGGCGGAGAGTCAGTTCGCGGTACTGGAGGAGGCGGCGCAACGTCGCCAGCTGTCTGCA"//&
        "CAGGAGAAATCCCTGCTGGCGCATAAAGATGAGACGCTGGAGTACAAACGCCAGCTGGCTGCACTTGGCG"//&
        "ACAAGGTTACGTATCAGGAGCGCCTGAACGCGCTGGCGCAGCAGGCGGATAAATTCGCACAGCAGCAACG"//&
        "GGCAAAACGGGCCGCCATTGATGCGAAAAGCCGGGGGCTGACTGACCGGCAGGCAGAACGGGAAGCCACG"//&
        "GAACAGCGCCTGAAGGAACAGTATGGCGATAATCCGCTGGCGCTGAATAACGTCATGTCAGAGCAGAAAA"//&
        "AGACCTGGGCGGCTGAAGACCAGCTTCGCGGGAACTGGATGGCAGGCCTGAAGTCCGGCTGGAGTGAGTG"//&
        "GGAAGAGAGCGCCACGGACAGTATGTCGCAGGTAAAAAGTGCAGCCACGCAGACCTTTGATGGTATTGCA"//&
        "CAGAATATGGCGGCGATGCTGACCGGCAGTGAGCAGAACTGGCGCAGCTTCACCCGTTCCGTGCTGTCCA"//&
        "TGATGACAGAAATTCTGCTTAAGCAGGCAATGGTGGGGATTGTCGGGAGTATCGGCAGCGCCATTGGCGG"//&
        "GGCTGTTGGTGGCGGCGCATCCGCGTCAGGCGGTACAGCCATTCAGGCCGCTGCGGCGAAATTCCATTTT"//&
        "GCAACCGGAGGATTTACGGGAACCGGCGGCAAATATGAGCCAGCGGGGATTGTTCACCGTGGTGAGTTTG"//&
        "TCTTCACGAAGGAGGCAACCAGCCGGATTGGCGTGGGGAATCTTTACCGGCTGATGCGCGGCTATGCCAC"//&
        "CGGCGGTTATGTCGGTACACCGGGCAGCATGGCAGACAGCCGGTCGCAGGCGTCCGGGACGTTTGAGCAG"//&
        "AATAACCATGTGGTGATTAACAACGACGGCACGAACGGGCAGATAGGTCCGGCTGCTCTGAAGGCGGTGT"//&
        "ATGACATGGCCCGCAAGGGTGCCCGTGATGAAATTCAGACACAGATGCGTGATGGTGGCCTGTTCTCCGG"//&
        "AGGTGGACGATGAAGACCTTCCGCTGGAAAGTGAAACCCGGTATGGATGTGGCTTCGGTCCCTTCTGTAA"//&
        "GAAAGGTGCGCTTTGGTGATGGCTATTCTCAGCGAGCGCCTGCCGGGCTGAATGCCAACCTGAAAACGTA"//&
        "CAGCGTGACGCTTTCTGTCCCCCGTGAGGAGGCCACGGTACTGGAGTCGTTTCTGGAAGAGCACGGGGGC"//&
        "TGGAAATCCTTTCTGTGGACGCCGCCTTATGAGTGGCGGCAGATAAAGGTGACCTGCGCAAAATGGTCGT"//&
        "CGCGGGTCAGTATGCTGCGTGTTGAGTTCAGCGCAGAGTTTGAACAGGTGGTGAACTGATGCAGGATATC"//&
        "CGGCAGGAAACACTGAATGAATGCACCCGTGCGGAGCAGTCGGCCAGCGTGGTGCTCTGGGAAATCGACC"//&
        "TGACAGAGGTCGGTGGAGAACGTTATTTTTTCTGTAATGAGCAGAACGAAAAAGGTGAGCCGGTCACCTG"//&
        "GCAGGGGCGACAGTATCAGCCGTATCCCATTCAGGGGAGCGGTTTTGAACTGAATGGCAAAGGCACCAGT"//&
        "ACGCGCCCCACGCTGACGGTTTCTAACCTGTACGGTATGGTCACCGGGATGGCGGAAGATATGCAGAGTC"//&
        "TGGTCGGCGGAACGGTGGTCCGGCGTAAGGTTTACGCCCGTTTTCTGGATGCGGTGAACTTCGTCAACGG"//&
        "AAACAGTTACGCCGATCCGGAGCAGGAGGTGATCAGCCGCTGGCGCATTGAGCAGTGCAGCGAACTGAGC"//&
        "GCGGTGAGTGCCTCCTTTGTACTGTCCACGCCGACGGAAACGGATGGCGCTGTTTTTCCGGGACGTATCA"//&
        "TGCTGGCCAACACCTGCACCTGGACCTATCGCGGTGACGAGTGCGGTTATAGCGGTCCGGCTGTCGCGGA"//&
        "TGAATATGACCAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTC"//&
        "CGCAATAACGTCGGCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCGCAGTAAATCCCATGACACA"//&
        "GACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGC"//&
        "ACGCCGGAGGGGGAAAGATATTTCCCCTGCGTGAATATCTCCGGTGAGCCGGAGGCTATTTCCGTATGTC"//&
        "GCCGGAAGACTGGCTGCAGGCAGAAATGCAGGGTGAGATTGTGGCGCTGGTCCACAGCCACCCCGGTGGT"//&
        "CTGCCCTGGCTGAGTGAGGCCGACCGGCGGCTGCAGGTGCAGAGTGATTTGCCGTGGTGGCTGGTCTGCC"//&
        "GGGGGACGATTCATAAGTTCCGCTGTGTGCCGCATCTCACCGGGCGGCGCTTTGAGCACGGTGTGACGGA"//&
        "CTGTTACACACTGTTCCGGGATGCTTATCATCTGGCGGGGATTGAGATGCCGGACTTTCATCGTGAGGAT"//&
        "GACTGGTGGCGTAACGGCCAGAATCTCTATCTGGATAATCTGGAGGCGACGGGGCTGTATCAGGTGCCGT"//&
        "TGTCAGCGGCACAGCCGGGCGATGTGCTGCTGTGCTGTTTTGGTTCATCAGTGCCGAATCACGCCGCAAT"//&
        "TTACTGCGGCGACGGCGAGCTGCTGCACCATATTCCTGAACAACTGAGCAAACGAGAGAGGTACACCGAC"//&
        "AAATGGCAGCGACGCACACACTCCCTCTGGCGTCACCGGGCATGGCGCGCATCTGCCTTTACGGGGATTT"//&
        "ACAACGATTTGGTCGCCGCATCGACCTTCGTGTGAAAACGGGGGCTGAAGCCATCCGGGCACTGGCCACA"//&
        "CAGCTCCCGGCGTTTCGTCAGAAACTGAGCGACGGCTGGTATCAGGTACGGATTGCCGGGCGGGACGTCA"//&
        "GCACGTCCGGGTTAACGGCGCAGTTACATGAGACTCTGCCTGATGGCGCTGTAATTCATATTGTTCCCAG"//&
        "AGTCGCCGGGGCCAAGTCAGGTGGCGTATTCCAGATTGTCCTGGGGGCTGCCGCCATTGCCGGATCATTC"//&
        "TTTACCGCCGGAGCCACCCTTGCAGCATGGGGGGCAGCCATTGGGGCCGGTGGTATGACCGGCATCCTGT"//&
        "TTTCTCTCGGTGCCAGTATGGTGCTCGGTGGTGTGGCGCAGATGCTGGCACCGAAAGCCAGAACTCCCCG"//&
        "TATACAGACAACGGATAACGGTAAGCAGAACACCTATTTCTCCTCACTGGATAACATGGTTGCCCAGGGC"//&
        "AATGTTCTGCCTGTTCTGTACGGGGAAATGCGCGTGGGGTCACGCGTGGTTTCTCAGGAGATCAGCACGG"//&
        "CAGACGAAGGGGACGGTGGTCAGGTTGTGGTGATTGGTCGCTGATGCAAAATGTTTTATGTGAAACCGCC"//&
        "TGCGGGCGGTTTTGTCATTTATGGAGCGTGAGGAATGGGTAAAGGAAGCAGTAAGGGGCATACCCCGCGC"//&
        "GAAGCGAAGGACAACCTGAAGTCCACGCAGTTGCTGAGTGTGATCGATGCCATCAGCGAAGGGCCGATTG"//&
        "AAGGTCCGGTGGATGGCTTAAAAAGCGTGCTGCTGAACAGTACGCCGGTGCTGGACACTGAGGGGAATAC"//&
        "CAACATATCCGGTGTCACGGTGGTGTTCCGGGCTGGTGAGCAGGAGCAGACTCCGCCGGAGGGATTTGAA"//&
        "TCCTCCGGCTCCGAGACGGTGCTGGGTACGGAAGTGAAATATGACACGCCGATCACCCGCACCATTACGT"//&
        "CTGCAAACATCGACCGTCTGCGCTTTACCTTCGGTGTACAGGCACTGGTGGAAACCACCTCAAAGGGTGA"//&
        "CAGGAATCCGTCGGAAGTCCGCCTGCTGGTTCAGATACAACGTAACGGTGGCTGGGTGACGGAAAAAGAC"//&
        "ATCACCATTAAGGGCAAAACCACCTCGCAGTATCTGGCCTCGGTGGTGATGGGTAACCTGCCGCCGCGCC"//&
        "CGTTTAATATCCGGATGCGCAGGATGACGCCGGACAGCACCACAGACCAGCTGCAGAACAAAACGCTCTG"//&
        "GTCGTCATACACTGAAATCATCGATGTGAAACAGTGCTACCCGAACACGGCACTGGTCGGCGTGCAGGTG"//&
        "GACTCGGAGCAGTTCGGCAGCCAGCAGGTGAGCCGTAATTATCATCTGCGCGGGCGTATTCTGCAGGTGC"//&
        "CGTCGAACTATAACCCGCAGACGCGGCAATACAGCGGTATCTGGGACGGAACGTTTAAACCGGCATACAG"//&
        "CAACAACATGGCCTGGTGTCTGTGGGATATGCTGACCCATCCGCGCTACGGCATGGGGAAACGTCTTGGT"//&
        "GCGGCGGATGTGGATAAATGGGCGCTGTATGTCATCGGCCAGTACTGCGACCAGTCAGTGCCGGACGGCT"//&
        "TTGGCGGCACGGAGCCGCGCATCACCTGTAATGCGTACCTGACCACACAGCGTAAGGCGTGGGATGTGCT"//&
        "CAGCGATTTCTGCTCGGCGATGCGCTGTATGCCGGTATGGAACGGGCAGACGCTGACGTTCGTGCAGGAC"//&
        "CGACCGTCGGATAAGACGTGGACCTATAACCGCAGTAATGTGGTGATGCCGGATGATGGCGCGCCGTTCC"//&
        "GCTACAGCTTCAGCGCCCTGAAGGACCGCCATAATGCCGTTGAGGTGAACTGGATTGACCCGAACAACGG"//&
        "CTGGGAGACGGCGACAGAGCTTGTTGAAGATACGCAGGCCATTGCCCGTTACGGTCGTAATGTTACGAAG"//&
        "ATGGATGCCTTTGGCTGTACCAGCCGGGGGCAGGCACACCGCGCCGGGCTGTGGCTGATTAAAACAGAAC"//&
        "TGCTGGAAACGCAGACCGTGGATTTCAGCGTCGGCGCAGAAGGGCTTCGCCATGTACCGGGCGATGTTAT"//&
        "TGAAATCTGCGATGATGACTATGCCGGTATCAGCACCGGTGGTCGTGTGCTGGCGGTGAACAGCCAGACC"//&
        "CGGACGCTGACGCTCGACCGTGAAATCACGCTGCCATCCTCCGGTACCGCGCTGATAAGCCTGGTTGACG"//&
        "GAAGTGGCAATCCGGTCAGCGTGGAGGTTCAGTCCGTCACCGACGGCGTGAAGGTAAAAGTGAGCCGTGT"//&
        "TCCTGACGGTGTTGCTGAATACAGCGTATGGGAGCTGAAGCTGCCGACGCTGCGCCAGCGACTGTTCCGC"//&
        "TGCGTGAGTATCCGTGAGAACGACGACGGCACGTATGCCATCACCGCCGTGCAGCATGTGCCGGAAAAAG"//&
        "AGGCCATCGTGGATAACGGGGCGCACTTTGACGGCGAACAGAGTGGCACGGTGAATGGTGTCACGCCGCC"//&
        "AGCGGTGCAGCACCTGACCGCAGAAGTCACTGCAGACAGCGGGGAATATCAGGTGCTGGCGCGATGGGAC"//&
        "ACACCGAAGGTGGTGAAGGGCGTGAGTTTCCTGCTCCGTCTGACCGTAACAGCGGACGACGGCAGTGAGC"//&
        "GGCTGGTCAGCACGGCCCGGACGACGGAAACCACATACCGCTTCACGCAACTGGCGCTGGGGAACTACAG"//&
        "GCTGACAGTCCGGGCGGTAAATGCGTGGGGGCAGCAGGGCGATCCGGCGTCGGTATCGTTCCGGATTGCC"//&
        "GCACCGGCAGCACCGTCGAGGATTGAGCTGACGCCGGGCTATTTTCAGATAACCGCCACGCCGCATCTTG"//&
        "CCGTTTATGACCCGACGGTACAGTTTGAGTTCTGGTTCTCGGAAAAGCAGATTGCGGATATCAGACAGGT"//&
        "TGAAACCAGCACGCGTTATCTTGGTACGGCGCTGTACTGGATAGCCGCCAGTATCAATATCAAACCGGGC"//&
        "CATGATTATTACTTTTATATCCGCAGTGTGAACACCGTTGGCAAATCGGCATTCGTGGAGGCCGTCGGTC"//&
        "GGGCGAGCGATGATGCGGAAGGTTACCTGGATTTTTTCAAAGGCAAGATAACCGAATCCCATCTCGGCAA"//&
        "GGAGCTGCTGGAAAAAGTCGAGCTGACGGAGGATAACGCCAGCAGACTGGAGGAGTTTTCGAAAGAGTGG"//&
        "AAGGATGCCAGTGATAAGTGGAATGCCATGTGGGCTGTCAAAATTGAGCAGACCAAAGACGGCAAACATT"//&
        "ATGTCGCGGGTATTGGCCTCAGCATGGAGGACACGGAGGAAGGCAAACTGAGCCAGTTTCTGGTTGCCGC"//&
        "CAATCGTATCGCATTTATTGACCCGGCAAACGGGAATGAAACGCCGATGTTTGTGGCGCAGGGCAACCAG"//&
        "ATATTCATGAACGACGTGTTCCTGAAGCGCCTGACGGCCCCCACCATTACCAGCGGCGGCAATCCTCCGG"//&
        "CCTTTTCCCTGACACCGGACGGAAAGCTGACCGCTAAAAATGCGGATATCAGTGGCAGTGTGAATGCGAA"//&
        "CTCCGGGACGCTCAGTAATGTGACGATAGCTGAAAACTGTACGATAAACGGTACGCTGAGGGCGGAAAAA"//&
        "ATCGTCGGGGACATTGTAAAGGCGGCGAGCGCGGCTTTTCCGCGCCAGCGTGAAAGCAGTGTGGACTGGC"//&
        "CGTCAGGTACCCGTACTGTCACCGTGACCGATGACCATCCTTTTGATCGCCAGATAGTGGTGCTTCCGCT"//&
        "GACGTTTCGCGGAAGTAAGCGTACTGTCAGCGGCAGGACAACGTATTCGATGTGTTATCTGAAAGTACTG"//&
        "ATGAACGGTGCGGTGATTTATGATGGCGCGGCGAACGAGGCGGTACAGGTGTTCTCCCGTATTGTTGACA"//&
        "TGCCAGCGGGTCGGGGAAACGTGATCCTGACGTTCACGCTTACGTCCACACGGCATTCGGCAGATATTCC"//&
        "GCCGTATACGTTTGCCAGCGATGTGCAGGTTATGGTGATTAAGAAACAGGCGCTGGGCATCAGCGTGGTC"//&
        "TGAGTGTGTTACAGAGGTTCGTCCGGGAACGGGCGTTTTATTATAAAACAGTGAGAGGTGAACGATGCGT"//&
        "AATGTGTGTATTGCCGTTGCTGTCTTTGCCGCACTTGCGGTGACAGTCACTCCGGCCCGTGCGGAAGGTG"//&
        "GACATGGTACGTTTACGGTGGGCTATTTTCAAGTGAAACCGGGTACATTGCCGTCGTTGTCGGGCGGGGA"//&
        "TACCGGTGTGAGTCATCTGAAAGGGATTAACGTGAAGTACCGTTATGAGCTGACGGACAGTGTGGGGGTG"//&
        "ATGGCTTCCCTGGGGTTCGCCGCGTCGAAAAAGAGCAGCACAGTGATGACCGGGGAGGATACGTTTCACT"//&
        "ATGAGAGCCTGCGTGGACGTTATGTGAGCGTGATGGCCGGACCGGTTTTACAAATCAGTAAGCAGGTCAG"//&
        "TGCGTACGCCATGGCCGGAGTGGCTCACAGTCGGTGGTCCGGCAGTACAATGGATTACCGTAAGACGGAA"//&
        "ATCACTCCCGGGTATATGAAAGAGACGACCACTGCCAGGGACGAAAGTGCAATGCGGCATACCTCAGTGG"//&
        "CGTGGAGTGCAGGTATACAGATTAATCCGGCAGCGTCCGTCGTTGTTGATATTGCTTATGAAGGCTCCGG"//&
        "CAGTGGCGACTGGCGTACTGACGGATTCATCGTTGGGGTCGGTTATAAATTCTGATTAGCCAGGTAACAC"//&
        "AGTGTTATGACAGCCCGCCGGAACCGGTGGGCTTTTTTGTGGGGTGAATATGGCAGTAAAGATTTCAGGA"//&
        "GTCCTGAAAGACGGCACAGGAAAACCGGTACAGAACTGCACCATTCAGCTGAAAGCCAGACGTAACAGCA"//&
        "CCACGGTGGTGGTGAACACGGTGGGCTCAGAGAATCCGGATGAAGCCGGGCGTTACAGCATGGATGTGGA"//&
        "GTACGGTCAGTACAGTGTCATCCTGCAGGTTGACGGTTTTCCACCATCGCACGCCGGGACCATCACCGTG"//&
        "TATGAAGATTCACAACCGGGGACGCTGAATGATTTTCTCTGTGCCATGACGGAGGATGATGCCCGGCCGG"//&
        "AGGTGCTGCGTCGTCTTGAACTGATGGTGGAAGAGGTGGCGCGTAACGCGTCCGTGGTGGCACAGAGTAC"//&
        "GGCAGACGCGAAGAAATCAGCCGGCGATGCCAGTGCATCAGCTGCTCAGGTCGCGGCCCTTGTGACTGAT"//&
        "GCAACTGACTCAGCACGCGCCGCCAGCACGTCCGCCGGACAGGCTGCATCGTCAGCTCAGGAAGCGTCCT"//&
        "CCGGCGCAGAAGCGGCATCAGCAAAGGCCACTGAAGCGGAAAAAAGTGCCGCAGCCGCAGAGTCCTCAAA"//&
        "AAACGCGGCGGCCACCAGTGCCGGTGCGGCGAAAACGTCAGAAACGAATGCTGCAGCGTCACAACAATCA"//&
        "GCCGCCACGTCTGCCTCCACCGCGGCCACGAAAGCGTCAGAGGCCGCCACTTCAGCACGAGATGCGGTGG"//&
        "CCTCAAAAGAGGCAGCAAAATCATCAGAAACGAACGCATCATCAAGTGCCGGTCGTGCAGCTTCCTCGGC"//&
        "AACGGCGGCAGAAAATTCTGCCAGGGCGGCAAAAACGTCCGAGACGAATGCCAGGTCATCTGAAACAGCA"//&
        "GCGGAACGGAGCGCCTCTGCCGCGGCAGACGCAAAAACAGCGGCGGCGGGGAGTGCGTCAACGGCATCCA"//&
        "CGAAGGCGACAGAGGCTGCGGGAAGTGCGGTATCAGCATCGCAGAGCAAAAGTGCGGCAGAAGCGGCGGC"//&
        "AATACGTGCAAAAAATTCGGCAAAACGTGCAGAAGATATAGCTTCAGCTGTCGCGCTTGAGGATGCGGAC"//&
        "ACAACGAGAAAGGGGATAGTGCAGCTCAGCAGTGCAACCAACAGCACGTCTGAAACGCTTGCTGCAACGC"//&
        "CAAAGGCGGTTAAGGTGGTAATGGATGAAACGAACAGAAAAGCCCACTGGACAGTCCGGCACTGACCGGA"//&
        "ACGCCAACAGCACCAACCGCGCTCAGGGGAACAAACAATACCCAGATTGCGAACACCGCTTTTGTACTGG"//&
        "CCGCGATTGCAGATGTTATCGACGCGTCACCTGACGCACTGAATACGCTGAATGAACTGGCCGCAGCGCT"//&
        "CGGGAATGATCCAGATTTTGCTACCACCATGACTAACGCGCTTGCGGGTAAACAACCGAAGAATGCGACA"//&
        "CTGACGGCGCTGGCAGGGCTTTCCACGGCGAAAAATAAATTACCGTATTTTGCGGAAAATGATGCCGCCA"//&
        "GCCTGACTGAACTGACTCAGGTTGGCAGGGATATTCTGGCAAAAAATTCCGTTGCAGATGTTCTTGAATA"//&
        "CCTTGGGGCCGGTGAGAATTCGGCCTTTCCGGCAGGTGCGCCGATCCCGTGGCCATCAGATATCGTTCCG"//&
        "TCTGGCTACGTCCTGATGCAGGGGCAGGCGTTTGACAAATCAGCCTACCCAAAACTTGCTGTCGCGTATC"//&
        "CATCGGGTGTGCTTCCTGATATGCGAGGCTGGACAATCAAGGGGAAACCCGCCAGCGGTCGTGCTGTATT"//&
        "GTCTCAGGAACAGGATGGAATTAAGTCGCACACCCACAGTGCCAGTGCATCCGGTACGGATTTGGGGACG"//&
        "AAAACCACATCGTCGTTTGATTACGGGACGAAAACAACAGGCAGTTTCGATTACGGCACCAAATCGACGA"//&
        "ATAACACGGGGGCTCATGCTCACAGTCTGAGCGGTTCAACAGGGGCCGCGGGTGCTCATGCCCACACAAG"//&
        "TGGTTTAAGGATGAACAGTTCTGGCTGGAGTCAGTATGGAACAGCAACCATTACAGGAAGTTTATCCACA"//&
        "GTTAAAGGAACCAGCACACAGGGTATTGCTTATTTATCGAAAACGGACAGTCAGGGCAGCCACAGTCACT"//&
        "CATTGTCCGGTACAGCCGTGAGTGCCGGTGCACATGCGCATACAGTTGGTATTGGTGCGCACCAGCATCC"//&
        "GGTTGTTATCGGTGCTCATGCCCATTCTTTCAGTATTGGTTCACACGGACACACCATCACCGTTAACGCT"//&
        "GCGGGTAACGCGGAAAACACCGTCAAAAACATTGCATTTAACTATATTGTGAGGCTTGCATAATGGCATT"//&
        "CAGAATGAGTGAACAACCACGGACCATAAAAATTTATAATCTGCTGGCCGGAACTAATGAATTTATTGGT"//&
        "GAAGGTGACGCATATATTCCGCCTCATACCGGTCTGCCTGCAAACAGTACCGATATTGCACCGCCAGATA"//&
        "TTCCGGCTGGCTTTGTGGCTGTTTTCAACAGTGATGAGGCATCGTGGCATCTCGTTGAAGACCATCGGGG"//&
        "TAAAACCGTCTATGACGTGGCTTCCGGCGACGCGTTATTTATTTCTGAACTCGGTCCGTTACCGGAAAAT"//&
        "TTTACCTGGTTATCGCCGGGAGGGGAATATCAGAAGTGGAACGGCACAGCCTGGGTGAAGGATACGGAAG"//&
        "CAGAAAAACTGTTCCGGATCCGGGAGGCGGAAGAAACAAAAAAAAGCCTGATGCAGGTAGCCAGTGAGCA"//&
        "TATTGCGCCGCTTCAGGATGCTGCAGATCTGGAAATTGCAACGAAGGAAGAAACCTCGTTGCTGGAAGCC"//&
        "TGGAAGAAGTATCGGGTGTTGCTGAACCGTGTTGATACATCAACTGCACCTGATATTGAGTGGCCTGCTG"//&
        "TCCCTGTTATGGAGTAATCGTTTTGTGATATGCCGCAGAAACGTTGTATGAAATAACGTTCTGCGGTTAG"//&
        "TTAGTATATTGTAAAGCTGAGTATTGGTTTATTTGGCGATTATTATCTTCAGGAGAATAATGGAAGTTCT"//&
        "ATGACTCAATTGTTCATAGTGTTTACATCACCGCCAATTGCTTTTAAGACTGAACGCATGAAATATGGTT"//&
        "TTTCGTCATGTTTTGAGTCTGCTGTTGATATTTCTAAAGTCGGTTTTTTTTCTTCGTTTTCTCTAACTAT"//&
        "TTTCCATGAAATACATTTTTGATTATTATTTGAATCAATTCCAATTACCTGAAGTCTTTCATCTATAATT"//&
        "GGCATTGTATGTATTGGTTTATTGGAGTAGATGCTTGCTTTTCTGAGCCATAGCTCTGATATCCAAATGA"//&
        "AGCCATAGGCATTTGTTATTTTGGCTCTGTCAGCTGCATAACGCCAAAAAATATATTTATCTGCTTGATC"//&
        "TTCAAATGTTGTATTGATTAAATCAATTGGATGGAATTGTTTATCATAAAAAATTAATGTTTGAATGTGA"//&
        "TAACCGTCCTTTAAAAAAGTCGTTTCTGCAAGCTTGGCTGTATAGTCAACTAACTCTTCTGTCGAAGTGA"//&
        "TATTTTTAGGCTTATCTACCAGTTTTAGACGCTCTTTAATATCTTCAGGAATTATTTTATTGTCATATTG"//&
        "TATCATGCTAAATGACAATTTGCTTATGGAGTAATCTTTTAATTTTAAATAAGTTATTCTCCTGGCTTCA"//&
        "TCAAATAAAGAGTCGAATGATGTTGGCGAAATCACATCGTCACCCATTGGATTGTTTATTTGTATGCCAA"//&
        "GAGAGTTACAGCAGTTATACATTCTGCCATAGATTATAGCTAAGGCATGTAATAATTCGTAATCTTTTAG"//&
        "CGTATTAGCGACCCATCGTCTTTCTGATTTAATAATAGATGATTCAGTTAAATATGAAGGTAATTTCTTT"//&
        "TGTGCAAGTCTGACTAACTTTTTTATACCAATGTTTAACATACTTTCATTTGTAATAAACTCAATGTCAT"//&
        "TTTCTTCAATGTAAGATGAAATAAGAGTAGCCTTTGCCTCGCTATACATTTCTAAATCGCCTTGTTTTTC"//&
        "TATCGTATTGCGAGAATTTTTAGCCCAAGCCATTAATGGATCATTTTTCCATTTTTCAATAACATTATTG"//&
        "TTATACCAAATGTCATATCCTATAATCTGGTTTTTGTTTTTTTGAATAATAAATGTTACTGTTCTTGCGG"//&
        "TTTGGAGGAATTGATTCAAATTCAAGCGAAATAATTCAGGGTCAAAATATGTATCAATGCAGCATTTGAG"//&
        "CAAGTGCGATAAATCTTTAAGTCTTCTTTCCCATGGTTTTTTAGTCATAAAACTCTCCATTTTGATAGGT"//&
        "TGCATGCTAGATGCTGATATATTTTAGAGGTGATAAAATTAACTGCTTAACTGTCAATGTAATACAAGTT"//&
        "GTTTGATCTTTGCAATGATTCTTATCAGAAACCATATAGTAAATTAGTTACACAGGAAATTTTTAATATT"//&
        "ATTATTATCATTCATTATGTATTAAAATTAGAGTTGTGGCTTGGCTCTGCTAACACGTTGCTCATAGGAG"//&
        "ATATGGTAGAGCCGCAGACACGTCGTATGCAGGAACGTGCTGCGGCTGGCTGGTGAACTTCCGATAGTGC"//&
        "GGGTGTTGAATGATTTCCAGTTGCTACCGATTTTACATATTTTTTGCATGAGAGAATTTGTACCACCTCC"//&
        "CACCGACCATCTATGACTGTACGCCACTGTCCCTAGGACTGCTATGTGCCGGAGCGGACATTACAAACGT"//&
        "CCTTCTCGGTGCATGCCACTGTTGCCAATGACCTGCCTAGGAATTGGTTAGCAAGTTACTACCGGATTTT"//&
        "GTAAAAACAGCCCTCCTCATATAAAAAGTATTCGTTCACTTCCGATAAGCGTCGTAATTTTCTATCTTTC"//&
        "ATCATATTCTAGATCCCTCTGAAAAAATCTTCCGAGTTTGCTAGGCACTGATACATAACTCTTTTCCAAT"//&
        "AATTGGGGAAGTCATTCAAATCTATAATAGGTTTCAGATTTGCTTCAATAAATTCTGACTGTAGCTGCTG"//&
        "AAACGTTGCGGTTGAACTATATTTCCTTATAACTTTTACGAAAGAGTTTCTTTGAGTAATCACTTCACTC"//&
        "AAGTGCTTCCCTGCCTCCAAACGATACCTGTTAGCAATATTTAATAGCTTGAAATGATGAAGAGCTCTGT"//&
        "GTTTGTCTTCCTGCCTCCAGTTCGCCGGGCATTCAACATAAAAACTGATAGCACCCGGAGTTCCGGAAAC"//&
        "GAAATTTGCATATACCCATTGCTCACGAAAAAAAATGTCCTTGTCGATATAGGGATGAATCGCTTGGTGT"//&
        "ACCTCATCTACTGCGAAAACTTGACCTTTCTCTCCCATATTGCAGTCGCGGCACGATGGAACTAAATTAA"//&
        "TAGGCATCACCGAAAATTCAGGATAATGTGCAATAGGAAGAAAATGATCTATATTTTTTGTCTGTCCTAT"//&
        "ATCACCACAAAATGGACATTTTTCACCTGATGAAACAAGCATGTCATCGTAATATGTTCTAGCGGGTTTG"//&
        "TTTTTATCTCGGAGATTATTTTCATAAAGCTTTTCTAATTTAACCTTTGTCAGGTTACCAACTACTAAGG"//&
        "TTGTAGGCTCAAGAGGGTGTGTCCTGTCGTAGGTAAATAACTGACCTGTCGAGCTTAATATTCTATATTG"//&
        "TTGTTCTTTCTGCAAAAAAGTGGGGAAGTGAGTAATGAAATTATTTCTAACATTTATCTGCATCATACCT"//&
        "TCCGAGCATTTATTAAGCATTTCGCTATAAGTTCTCGCTGGAAGAGGTAGTTTTTTCATTGTACTTTACC"//&
        "TTCATCTCTGTTCATTATCATCGCTTTTAAAACGGTTCGACCTTCTAATCCTATCTGACCATTATAATTT"//&
        "TTTAGAATGGTTTCATAAGAAAGCTCTGAATCAACGGACTGCGATAATAAGTGGTGGTATCCAGAATTTG"//&
        "TCACTTCAAGTAAAAACACCTCACGAGTTAAAACACCTAAGTTCTCACCGAATGTCTCAATATCCGGACG"//&
        "GATAATATTTATTGCTTCTCTTGACCGTAGGACTTTCCACATGCAGGATTTTGGAACCTCTTGCAGTACT"//&
        "ACTGGGGAATGAGTTGCAATTATTGCTACACCATTGCGTGCATCGAGTAAGTCGCTTAATGTTCGTAAAA"//&
        "AAGCAGAGAGCAAAGGTGGATGCAGATGAACCTCTGGTTCATCGAATAAAACTAATGACTTTTCGCCAAC"//&
        "GACATCTACTAATCTTGTGATAGTAAATAAAACAATTGCATGTCCAGAGCTCATTCGAAGCAGATATTTC"//&
        "TGGATATTGTCATAAAACAATTTAGTGAATTTATCATCGTCCACTTGAATCTGTGGTTCATTACGTCTTA"//&
        "ACTCTTCATATTTAGAAATGAGGCTGATGAGTTCCATATTTGAAAAGTTTTCATCACTACTTAGTTTTTT"//&
        "GATAGCTTCAAGCCAGAGTTGTCTTTTTCTATCTACTCTCATACAACCAATAAATGCTGAAATGAATTCT"//&
        "AAGCGGAGATCGCCTAGTGATTTTAAACTATTGCTGGCAGCATTCTTGAGTCCAATATAAAAGTATTGTG"//&
        "TACCTTTTGCTGGGTCAGGTTGTTCTTTAGGAGGAGTAAAAGGATCAAATGCACTAAACGAAACTGAAAC"//&
        "AAGCGATCGAAAATATCCCTTTGGGATTCTTGACTCGATAAGTCTATTATTTTCAGAGAAAAAATATTCA"//&
        "TTGTTTTCTGGGTTGGTGATTGCACCAATCATTCCATTCAAAATTGTTGTTTTACCACACCCATTCCGCC"//&
        "CGATAAAAGCATGAATGTTCGTGCTGGGCATAGAATTAACCGTCACCTCAAAAGGTATAGTTAAATCACT"//&
        "GAATCCGGGAGCACTTTTTCTATTAAATGAAAAGTGGAAATCTGACAATTCTGGCAAACCATTTAACACA"//&
        "CGTGCGAACTGTCCATGAATTTCTGAAAGAGTTACCCCTCTAAGTAATGAGGTGTTAAGGACGCTTTCAT"//&
        "TTTCAATGTCGGCTAATCGATTTGGCCATACTACTAAATCCTGAATAGCTTTAAGAAGGTTATGTTTAAA"//&
        "ACCATCGCTTAATTTGCTGAGATTAACATAGTAGTCAATGCTTTCACCTAAGGAAAAAAACATTTCAGGG"//&
        "AGTTGACTGAATTTTTTATCTATTAATGAATAAGTGCTTACTTCTTCTTTTTGACCTACAAAACCAATTT"//&
        "TAACATTTCCGATATCGCATTTTTCACCATGCTCATCAAAGACAGTAAGATAAAACATTGTAACAAAGGA"//&
        "ATAGTCATTCCAACCATCTGCTCGTAGGAATGCCTTATTTTTTTCTACTGCAGGAATATACCCGCCTCTT"//&
        "TCAATAACACTAAACTCCAACATATAGTAACCCTTAATTTTATTAAAATAACCGCAATTTATTTGGCGGC"//&
        "AACACAGGATCTCTCTTTTAAGTTACTCTCTATTACATACGTTTTCCATCTAAAAATTAGTAGTATTGAA"//&
        "CTTAACGGGGCATCGTATTGTAGTTTTCCATATTTAGCTTTCTGCTTCCTTTTGGATAACCCACTGTTAT"//&
        "TCATGTTGCATGGTGCACTGTTTATACCAACGATATAGTCTATTAATGCATATATAGTATCGCCGAACGA"//&
        "TTAGCTCTTCAGGCTTCTGAAGAAGCGTTTCAAGTACTAATAAGCCGATAGATAGCCACGGACTTCGTAG"//&
        "CCATTTTTCATAAGTGTTAACTTCCGCTCCTCGCTCATAACAGACATTCACTACAGTTATGGCGGAAAGG"//&
        "TATGCATGCTGGGTGTGGGGAAGTCGTGAAAGAAAAGAAGTCAGCTGCGTCGTTTGACATCACTGCTATC"//&
        "TTCTTACTGGTTATGCAGGTCGTAGTGGGTGGCACACAAAGCTTTGCACTGGATTGCGAGGCTTTGTGCT"//&
        "TCTCTGGAGTGCGACAGGTTTGATGACAAAAAATTAGCGCAAGAAGACAAAAATCACCTTGCGCTAATGC"//&
        "TCTGTTACAGGTCACTAATACCATCTAAGTAGTTGATTCATAGTGACTGCATATGTTGTGTTTTACAGTA"//&
        "TTATGTAGTCTGTTTTTTATGCAAAATCTAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGT"//&
        "TCAGCTTTTTTATACTAAGTTGGCATTATAAAAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCA"//&
        "CTATCAGTCAAAATAAAATCATTATTTGATTTCAATTTTGTCCCACTCCCTGCCTCTGTCATCACGATAC"//&
        "TGTGATGCCATGGTGTCCGACTTATGCCCGAGAAGATGTTGAGCAAACTTATCGCTTATCTGCTTCTCAT"//&
        "AGAGTCTTGCAGACAAACTGCGCAACTCGTGAAAGGTAGGCGGATCCCCTTCGAAGGAAAGACCTGATGC"//&
        "TTTTCGTGCGCGCATAAAATACCTTGATACTGTGCCGGATGAAAGCGGTTCGCGACGAGTAGATGCAATT"//&
        "ATGGTTTCTCCGCCAAGAATCTCTTTGCATTTATCAAGTGTTTCCTTCATTGATATTCCGAGAGCATCAA"//&
        "TATGCAATGCTGTTGGGATGGCAATTTTTACGCCTGTTTTGCTTTGCTCGACATAAAGATATCCATCTAC"//&
        "GATATCAGACCACTTCATTTCGCATAAATCACCAACTCGTTGCCCGGTAACAACAGCCAGTTCCATTGCA"//&
        "AGTCTGAGCCAACATGGTGATGATTCTGCTGCTTGATAAATTTTCAGGTATTCGTCAGCCGTAAGTCTTG"//&
        "ATCTCCTTACCTCTGATTTTGCTGCGCGAGTGGCAGCGACATGGTTTGTTGTTATATGGCCTTCAGCTAT"//&
        "TGCCTCTCGGAATGCATCGCTCAGTGTTGATCTGATTAACTTGGCTGACGCCGCCTTGCCCTCGTCTATG"//&
        "TATCCATTGAGCATTGCCGCAATTTCTTTTGTGGTGATGTCTTCAAGTGGAGCATCAGGCAGACCCCTCC"//&
        "TTATTGCTTTAATTTTGCTCATGTAATTTATGAGTGTCTTCTGCTTGATTCCTCTGCTGGCCAGGATTTT"//&
        "TTCGTAGCGATCAAGCCATGAATGTAACGTAACGGAATTATCACTGTTGATTCTCGCTGTCAGAGGCTTG"//&
        "TGTTTGTGTCCTGAAAATAACTCAATGTTGGCCTGTATAGCTTCAGTGATTGCGATTCGCCTGTCTCTGC"//&
        "CTAATCCAAACTCTTTACCCGTCCTTGGGTCCCTGTAGCAGTAATATCCATTGTTTCTTATATAAAGGTT"//&
        "AGGGGGTAAATCCCGGCGCTCATGACTTCGCCTTCTTCCCATTTCTGATCCTCTTCAAAAGGCCACCTGT"//&
        "TACTGGTCGATTTAAGTCAACCTTTACCGCTGATTCGTGGAACAGATACTCTCTTCCATCCTTAACCGGA"//&
        "GGTGGGAATATCCTGCATTCCCGAACCCATCGACGAACTGTTTCAAGGCTTCTTGGACGTCGCTGGCGTG"//&
        "CGTTCCACTCCTGAAGTGTCAAGTACATCGCAAAGTCTCCGCAATTACACGCAAGAAAAAACCGCCATCA"//&
        "GGCGGCTTGGTGTTCTTTCAGTTCTTCAATTCGAATATTGGTTACGTCTGCATGTGCTATCTGCGCCCAT"//&
        "ATCATCCAGTGGTCGTAGCAGTCGTTGATGTTCTCCGCTTCGATAACTCTGTTGAATGGCTCTCCATTCC"//&
        "ATTCTCCTGTGACTCGGAAGTGCATTTATCATCTCCATAAAACAAAACCCGCCGTAGCGAGTTCAGATAA"//&
        "AATAAATCCCCGCGAGTGCGAGGATTGTTATGTAATATTGGGTTTAATCATCTATATGTTTTGTACAGAG"//&
        "AGGGCAAGTATCGTTTCCACCGTACTCGTGATAATAATTTTGCACGGTATCAGTCATTTCTCGCACATTG"//&
        "CAGAATGGGGATTTGTCTTCATTAGACTTATAAACCTTCATGGAATATTTGTATGCCGACTCTATATCTA"//&
        "TACCTTCATCTACATAAACACCTTCGTGATGTCTGCATGGAGACAAGACACCGGATCTGCACAACATTGA"//&
        "TAACGCCCAATCTTTTTGCTCAGACTCTAACTCATTGATACTCATTTATAAACTCCTTGCAATGTATGTC"//&
        "GTTTCAGCTAAACGGTATCAGCAATGTTTATGTAAAGAAACAGTAAGATAATACTCAACCCGATGTTTGA"//&
        "GTACGGTCATCATCTGACACTACAGACTCTGGCATCGCTGTGAAGACGACGCGAAATTCAGCATTTTCAC"//&
        "AAGCGTTATCTTTTACAAAACCGATCTCACTCTCCTTTGATGCGAATGCCAGCGTCAGACATCATATGCA"//&
        "GATACTCACCTGCATCCTGAACCCATTGACCTCCAACCCCGTAATAGCGATGCGTAATGATGTCGATAGT"//&
        "TACTAACGGGTCTTGTTCGATTAACTGCCGCAGAAACTCTTCCAGGTCACCAGTGCAGTGCTTGATAACA"//&
        "GGAGTCTTCCCAGGATGGCGAACAACAAGAAACTGGTTTCCGTCTTCACGGACTTCGTTGCTTTCCAGTT"//&
        "TAGCAATACGCTTACTCCCATCCGAGATAACACCTTCGTAATACTCACGCTGCTCGTTGAGTTTTGATTT"//&
        "TGCTGTTTCAAGCTCAACACGCAGTTTCCCTACTGTTAGCGCAATATCCTCGTTCTCCTGGTCGCGGCGT"//&
        "TTGATGTATTGCTGGTTTCTTTCCCGTTCATCCAGCAGTTCCAGCACAATCGATGGTGTTACCAATTCAT"//&
        "GGAAAAGGTCTGCGTCAAATCCCCAGTCGTCATGCATTGCCTGCTCTGCCGCTTCACGCAGTGCCTGAGA"//&
        "GTTAATTTCGCTCACTTCGAACCTCTCTGTTTACTGATAAGTTCCAGATCCTCCTGGCAACTTGCACAAG"//&
        "TCCGACAACCCTGAACGACCAGGCGTCTTCGTTCATCTATCGGATCGCCACACTCACAACAATGAGTGGC"//&
        "AGATATAGCCTGGTGGTTCAGGCGGCGCATTTTTATTGCTGTGTTGCGCTGTAATTCTTCTATTTCTGAT"//&
        "GCTGAATCAATGATGTCTGCCATCTTTCATTAATCCCTGAACTGTTGGTTAATACGCTTGAGGGTGAATG"//&
        "CGAATAATAAAAAAGGAGCCTGTAGCTCCCTGATGATTTTGCTTTTCATGTTCATCGTTCCTTAAAGACG"//&
        "CCGTTTAACATGCCGATTGCCAGGCTTAAATGAGTCGGTGTGAATCCCATCAGCGTTACCGTTTCGCGGT"//&
        "GCTTCTTCAGTACGCTACGGCAAATGTCATCGACGTTTTTATCCGGAAACTGCTGTCTGGCTTTTTTTGA"//&
        "TTTCAGAATTAGCCTGACGGGCAATGCTGCGAAGGGCGTTTTCCTGCTGAGGTGTCATTGAACAAGTCCC"//&
        "ATGTCGGCAAGCATAAGCACACAGAATATGAAGCCCGCTGCCAGAAAAATGCATTCCGTGGTTGTCATAC"//&
        "CTGGTTTCTCTCATCTGCTTCTGCTTTCGCCACCATCATTTCCAGCTTTTGTGAAAGGGATGCGGCTAAC"//&
        "GTATGAAATTCTTCGTCTGTTTCTACTGGTATTGGCACAAACCTGATTCCAATTTGAGCAAGGCTATGTG"//&
        "CCATCTCGATACTCGTTCTTAACTCAACAGAAGATGCTTTGTGCATACAGCCCCTCGTTTATTATTTATC"//&
        "TCCTCAGCCAGCCGCTGTGCTTTCAGTGGATTTCGGATAACAGAAAGGCCGGGAAATACCCAGCCTCGCT"//&
        "TTGTAACGGAGTAGACGAAAGTGATTGCGCCTACCCGGATATTATCGTGAGGATGCGTCATCGCCATTGC"//&
        "TCCCCAAATACAAAACCAATTTCAGCCAGTGCCTCGTCCATTTTTTCGATGAACTCCGGCACGATCTCGT"//&
        "CAAAACTCGCCATGTACTTTTCATCCCGCTCAATCACGACATAATGCAGGCCTTCACGCTTCATACGCGG"//&
        "GTCATAGTTGGCAAAGTACCAGGCATTTTTTCGCGTCACCCACATGCTGTACTGCACCTGGGCCATGTAA"//&
        "GCTGACTTTATGGCCTCGAAACCACCGAGCCGGAACTTCATGAAATCCCGGGAGGTAAACGGGCATTTCA"//&
        "GTTCAAGGCCGTTGCCGTCACTGCATAAACCATCGGGAGAGCAGGCGGTACGCATACTTTCGTCGCGATA"//&
        "GATGATCGGGGATTCAGTAACATTCACGCCGGAAGTGAATTCAAACAGGGTTCTGGCGTCGTTCTCGTAC"//&
        "TGTTTTCCCCAGGCCAGTGCTTTAGCGTTAACTTCCGGAGCCACACCGGTGCAAACCTCAGCAAGCAGGG"//&
        "TGTGGAAGTAGGACATTTTCATGTCAGGCCACTTCTTTCCGGAGCGGGGTTTTGCTATCACGTTGTGAAC"//&
        "TTCTGAAGCGGTGATGACGCCGAGCCGTAATTTGTGCCACGCATCATCCCCCTGTTCGACAGCTCTCACA"//&
        "TCGATCCCGGTACGCTGCAGGATAATGTCCGGTGTCATGCTGCCACCTTCTGCTCTGCGGCTTTCTGTTT"//&
        "CAGGAATCCAAGAGCTTTTACTGCTTCGGCCTGTGTCAGTTCTGACGATGCACGAATGTCGCGGCGAAAT"//&
        "ATCTGGGAACAGAGCGGCAATAAGTCGTCATCCCATGTTTTATCCAGGGCGATCAGCAGAGTGTTAATCT"//&
        "CCTGCATGGTTTCATCGTTAACCGGAGTGATGTCGCGTTCCGGCTGACGTTCTGCAGTGTATGCAGTATT"//&
        "TTCGACAATGCGCTCGGCTTCATCCTTGTCATAGATACCAGCAAATCCGAAGGCCAGACGGGCACACTGA"//&
        "ATCATGGCTTTATGACGTAACATCCGTTTGGGATGCGACTGCCACGGCCCCGTGATTTCTCTGCCTTCGC"//&
        "GAGTTTTGAATGGTTCGCGGCGGCATTCATCCATCCATTCGGTAACGCAGATCGGATGATTACGGTCCTT"//&
        "GCGGTAAATCCGGCATGTACAGGATTCATTGTCCTGCTCAAAGTCCATGCCATCAAACTGCTGGTTTTCA"//&
        "TTGATGATGCGGGACCAGCCATCAACGCCCACCACCGGAACGATGCCATTCTGCTTATCAGGAAAGGCGT"//&
        "AAATTTCTTTCGTCCACGGATTAAGGCCGTACTGGTTGGCAACGATCAGTAATGCGATGAACTGCGCATC"//&
        "GCTGGCATCACCTTTAAATGCCGTCTGGCGAAGAGTGGTGATCAGTTCCTGTGGGTCGACAGAATCCATG"//&
        "CCGACACGTTCAGCCAGCTTCCCAGCCAGCGTTGCGAGTGCAGTACTCATTCGTTTTATACCTCTGAATC"//&
        "AATATCAACCTGGTGGTGAGCAATGGTTTCAACCATGTACCGGATGTGTTCTGCCATGCGCTCCTGAAAC"//&
        "TCAACATCGTCATCAAACGCACGGGTAATGGATTTTTTGCTGGCCCCGTGGCGTTGCAAATGATCGATGC"//&
        "ATAGCGATTCAAACAGGTGCTGGGGCAGGCCTTTTTCCATGTCGTCTGCCAGTTCTGCCTCTTTCTCTTC"//&
        "ACGGGCGAGCTGCTGGTAGTGACGCGCCCAGCTCTGAGCCTCAAGACGATCCTGAATGTAATAAGCGTTC"//&
        "ATGGCTGAACTCCTGAAATAGCTGTGAAAATATCGCCCGCGAAATGCCGGGCTGATTAGGAAAACAGGAA"//&
        "AGGGGGTTAGTGAATGCTTTTGCTTGATCTCAGTTTCAGTATTAATATCCATTTTTTATAAGCGTCGACG"//&
        "GCTTCACGAAACATCTTTTCATCGCCAATAAAAGTGGCGATAGTGAATTTAGTCTGGATAGCCATAAGTG"//&
        "TTTGATCCATTCTTTGGGACTCCTGGCTGATTAAGTATGTCGATAAGGCGTTTCCATCCGTCACGTAATT"//&
        "TACGGGTGATTCGTTCAAGTAAAGATTCGGAAGGGCAGCCAGCAACAGGCCACCCTGCAATGGCATATTG"//&
        "CATGGTGTGCTCCTTATTTATACATAACGAAAAACGCCTCGAGTGAAGCGTTATTGGTATGCGGTAAAAC"//&
        "CGCACTCAGGCGGCCTTGATAGTCATATCATCTGAATCAAATATTCCTGATGTATCGATATCGGTAATTC"//&
        "TTATTCCTTCGCTACCATCCATTGGAGGCCATCCTTCCTGACCATTTCCATCATTCCAGTCGAACTCACA"//&
        "CACAACACCATATGCATTTAAGTCGCTTGAAATTGCTATAAGCAGAGCATGTTGCGCCAGCATGATTAAT"//&
        "ACAGCATTTAATACAGAGCCGTGTTTATTGAGTCGGTATTCAGAGTCTGACCAGAAATTATTAATCTGGT"//&
        "GAAGTTTTTCCTCTGTCATTACGTCATGGTCGATTTCAATTTCTATTGATGCTTTCCAGTCGTAATCAAT"//&
        "GATGTATTTTTTGATGTTTGACATCTGTTCATATCCTCACAGATAAAAAATCGCCCTCACACTGGAGGGC"//&
        "AAAGAAGATTTCCAATAATCAGAACAAGTCGGCTCCTGTTTAGTTACGAGCGACATTGCTCCGTGTATTC"//&
        "ACTCGTTGGAATGAATACACAGTGCAGTGTTTATTCTGTTATTTATGCCAAAAATAAAGGCCACTATCAG"//&
        "GCAGCTTTGTTGTTCTGTTTACCAAGTTCTCTGGCAATCATTGCCGTCGTTCGTATTGCCCATTTATCGA"//&
        "CATATTTCCCATCTTCCATTACAGGAAACATTTCTTCAGGCTTAACCATGCATTCCGATTGCAGCTTGCA"//&
        "TCCATTGCATCGCTTGAATTGTCCACACCATTGATTTTTATCAATAGTCGTAGTCATACGGATAGTCCTG"//&
        "GTATTGTTCCATCACATCCTGAGGATGCTCTTCGAACTCTTCAAATTCTTCTTCCATATATCACCTTAAA"//&
        "TAGTGGATTGCGGTAGTAAAGATTGTGCCTGTCTTTTAACCACATCAGGCTCGGTGGTTCTCGTGTACCC"//&
        "CTACAGCGAGAAATCGGATAAACTATTACAACCCCTACAGTTTGATGAGTATAGAAATGGATCCACTCGT"//&
        "TATTCTCGGACGAGTGTTCAGTAATGAACCTCTGGAGAGAACCATGTATATGATCGTTATCTGGGTTGGA"//&
        "CTTCTGCTTTTAAGCCCAGATAACTGGCCTGAATATGTTAATGAGAGAATCGGTATTCCTCATGTGTGGC"//&
        "ATGTTTTCGTCTTTGCTCTTGCATTTTCGCTAGCAATTAATGTGCATCGATTATCAGCTATTGCCAGCGC"//&
        "CAGATATAAGCGATTTAAGCTAAGAAAACGCATTAAGATGCAAAACGATAAAGTGCGATCAGTAATTCAA"//&
        "AACCTTACAGAAGAGCAATCTATGGTTTTGTGCGCAGCCCTTAATGAAGGCAGGAAGTATGTGGTTACAT"//&
        "CAAAACAATTCCCATACATTAGTGAGTTGATTGAGCTTGGTGTGTTGAACAAAACTTTTTCCCGATGGAA"//&
        "TGGAAAGCATATATTATTCCCTATTGAGGATATTTACTGGACTGAATTAGTTGCCAGCTATGATCCATAT"//&
        "AATATTGAGATAAAGCCAAGGCCAATATCTAAGTAACTAGATAAGAGGAATCGATTTTCCCTTAATTTTC"//&
        "TGGCGTCCACTGCATGTTATGCCGCGTTCGCCAGGCTTGCTGTACCATGTGCGCTGATTCTTGCGCTCAA"//&
        "TACGTTGCAGGTTGCTTTCAATCTGTTTGTGGTATTCAGCCAGCACTGTAAGGTCTATCGGATTTAGTGC"//&
        "GCTTTCTACTCGTGATTTCGGTTTGCGATTCAGCGAGAGAATAGGGCGGTTAACTGGTTTTGCGCTTACC"//&
        "CCAACCAACAGGGGATTTGCTGCTTTCCATTGAGCCTGTTTCTCTGCGCGACGTTCGCGGCGGCGTGTTT"//&
        "GTGCATCCATCTGGATTCTCCTGTCAGTTAGCTTTGGTGGTGTGTGGCAGTTGTAGTCCTGAACGAAAAC"//&
        "CCCCCGCGATTGGCACATTGGCAGCTAATCCGGAATCGCACTTACGGCCAATGCTTCGTTTCGTATCACA"//&
        "CACCCCAAAGCCTTCTGCTTTGAATGCTGCCCTTCTTCAGGGCTTAATTTTTAAGAGCGTCACCTTCATG"//&
        "GTGGTCAGTGCGTCCTGCTGATGTGCTCAGTATCACCGCCAGTGGTATTTATGTCAACACCGCCAGAGAT"//&
        "AATTTATCACCGCAGATGGTTATCTGTATGTTTTTTATATGAATTTATTTTTTGCAGGGGGGCATTGTTT"//&
        "GGTAGGTGAGAGATCTGAATTGCTATGTTTAGTGAGTTGTATCTATTTATTTTTCAATAAATACAATTGG"//&
        "TTATGTGTTTTGGGGGCGATCGTGAGGCAAAGAAAACCCGGCGCTGAGGCCGGGTTATTCTTGTTCTCTG"//&
        "GTCAAATTATATAGTTGGAAAACAAGGATGCATATATGAATGAACGATGCAGAGGCAATGCCGATGGCGA"//&
        "TAGTGGGTATCATGTAGCCGCTTATGCTGGAAAGAAGCAATAACCCGCAGAAAAACAAAGCTCCAAGCTC"//&
        "AACAAAACTAAGGGCATAGACAATAACTACCGATGTCATATACCCATACTCTCTAATCTTGGCCAGTCGG"//&
        "CGCGTTCTGCTTCCGATTAGAAACGTCAAGGCAGCAATCAGGATTGCAATCATGGTTCCTGCATATGATG"//&
        "ACAATGTCGCCCCAAGACCATCTCTATGAGCTGAAAAAGAAACACCAGGAATGTAGTGGCGGAAAAGGAG"//&
        "ATAGCAAATGCTTACGATAACGTAAGGAATTATTACTATGTAAACACCAGGCATGATTCTGTTCCGCATA"//&
        "ATTACTCCTGATAATTAATCCTTAACTTTGCCCACCTGCCTTTTAAAACATTCCAGTATATCACTTTTCA"//&
        "TTCTTGCGTAGCAATATGCCATCTCTTCAGCTATCTCAGCATTGGTGACCTTGTTCAGAGGCGCTGAGAG"//&
        "ATGGCCTTTTTCTGATAGATAATGTTCTGTTAAAATATCTCCGGCCTCATCTTTTGCCCGCAGGCTAATG"//&
        "TCTGAAAATTGAGGTGACGGGTTAAAAATAATATCCTTGGCAACCTTTTTTATATCCCTTTTAAATTTTG"//&
        "GCTTAATGACTATATCCAATGAGTCAAAAAGCTCCCCTTCAATATCTGTTGCCCCTAAGACCTTTAATAT"//&
        "ATCGCCAAATACAGGTAGCTTGGCTTCTACCTTCACCGTTGTTCGGCCGATGAAATGCATATGCATAACA"//&
        "TCGTCTTTGGTGGTTCCCCTCATCAGTGGCTCTATCTGAACGCGCTCTCCACTGCTTAATGACATTCCTT"//&
        "TCCCGATTAAAAAATCTGTCAGATCGGATGTGGTCGGCCCGAAAACAGTTCTGGCAAAACCAATGGTGTC"//&
        "GCCTTCAACAAACAAAAAAGATGGGAATCCCAATGATTCGTCATCTGCGAGGCTGTTCTTAATATCTTCA"//&
        "ACTGAAGCTTTAGAGCGATTTATCTTCTGAACCAGACTCTTGTCATTTGTTTTGGTAAAGAGAAAAGTTT"//&
        "TTCCATCGATTTTATGAATATACAAATAATTGGAGCCAACCTGCAGGTGATGATTATCAGCCAGCAGAGA"//&
        "ATTAAGGAAAACAGACAGGTTTATTGAGCGCTTATCTTTCCCTTTATTTTTGCTGCGGTAAGTCGCATAA"//&
        "AAACCATTCTTCATAATTCAATCCATTTACTATGTTATGTTCTGAGGGGAGTGAAAATTCCCCTAATTCG"//&
        "ATGAAGATTCTTGCTCAATTGTTATCAGCTATGCGCCGACCAGAACACCTTGCCGATCAGCCAAACGTCT"//&
        "CTTCAGGCCACTGACTAGCGATAACTTTCCCCACAACGGAACAACTCTCATTGCATGGGATCATTGGGTA"//&
        "CTGTGGGTTTAGTGGTTGTAAAAACACCTGACCGCTATCCCTGATCAGTTTCTTGAAGGTAAACTCATCA"//&
        "CCCCCAAGTCTGGCTATGCAGAAATCACCTGGCTCAACAGCCTGCTCAGGGTCAACGAGAATTAACATTC"//&
        "CGTCAGGAAAGCTTGGCTTGGAGCCTGTTGGTGCGGTCATGGAATTACCTTCAACCTCAAGCCAGAATGC"//&
        "AGAATCACTGGCTTTTTTGGTTGTGCTTACCCATCTCTCCGCATCACCTTTGGTAAAGGTTCTAAGCTTA"//&
        "GGTGAGAACATCCCTGCCTGAACATGAGAAAAAACAGGGTACTCATACTCACTTCTAAGTGACGGCTGCA"//&
        "TACTAACCGCTTCATACATCTCGTAGATTTCTCTGGCGATTGAAGGGCTAAATTCTTCAACGCTAACTTT"//&
        "GAGAATTTTTGTAAGCAATGCGGCGTTATAAGCATTTAATGCATTGATGCCATTAAATAAAGCACCAACG"//&
        "CCTGACTGCCCCATCCCCATCTTGTCTGCGACAGATTCCTGGGATAAGCCAAGTTCATTTTTCTTTTTTT"//&
        "CATAAATTGCTTTAAGGCGACGTGCGTCCTCAAGCTGCTCTTGTGTTAATGGTTTCTTTTTTGTGCTCAT"//&
        "ACGTTAAATCTATCACCGCAAGGGATAAATATCTAACACCGTGCGTGTTGACTATTTTACCTCTGGCGGT"//&
        "GATAATGGTTGCATGTACTAAGGAGGTTGTATGGAACAACGCATAACCCTGAAAGATTATGCAATGCGCT"//&
        "TTGGGCAAACCAAGACAGCTAAAGATCTCGGCGTATATCAAAGCGCGATCAACAAGGCCATTCATGCAGG"//&
        "CCGAAAGATTTTTTTAACTATAAACGCTGATGGAAGCGTTTATGCGGAAGAGGTAAAGCCCTTCCCGAGT"//&
        "AACAAAAAAACAACAGCATAAATAACCCCGCTCTTACACATTCCAGCCCTGAAAAAGGGCATCAAATTAA"//&
        "ACCACACCTATGGTGTATGCATTTATTTGCATACATTCAATCAATTGTTATCTAAGGAAATACTTACATA"//&
        "TGGTTCGTGCAAACAAACGCAACGAGGCTCTACGAATCGAGAGTGCGTTGCTTAACAAAATCGCAATGCT"//&
        "TGGAACTGAGAAGACAGCGGAAGCTGTGGGCGTTGATAAGTCGCAGATCAGCAGGTGGAAGAGGGACTGG"//&
        "ATTCCAAAGTTCTCAATGCTGCTTGCTGTTCTTGAATGGGGGGTCGTTGACGACGACATGGCTCGATTGG"//&
        "CGCGACAAGTTGCTGCGATTCTCACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCA"//&
        "GATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCATTATGACAAATACAGCAAAAATACTCA"//&
        "ACTTCGGCAGAGGTAACTTTGCCGGACAGGAGCGTAATGTGGCAGATCTCGATGATGGTTACGCCAGACT"//&
        "ATCAAATATGCTGCTTGAGGCTTATTCGGGCGCAGATCTGACCAAGCGACAGTTTAAAGTGCTGCTTGCC"//&
        "ATTCTGCGTAAAACCTATGGGTGGAATAAACCAATGGACAGAATCACCGATTCTCAACTTAGCGAGATTA"//&
        "CAAAGTTACCTGTCAAACGGTGCAATGAAGCCAAGTTAGAACTCGTCAGAATGAATATTATCAAGCAGCA"//&
        "AGGCGGCATGTTTGGACCAAATAAAAACATCTCAGAATGGTGCATCCCTCAAAACGAGGGAAAATCCCCT"//&
        "AAAACGAGGGATAAAACATCCCTCAAATTGGGGGATTGCTATCCCTCAAAACAGGGGGACACAAAAGACA"//&
        "CTATTACAAAAGAAAAAAGAAAAGATTATTCGTCAGAGAATTCTGGCGAATCCTCTGACCAGCCAGAAAA"//&
        "CGACCTTTCTGTGGTGAAACCGGATGCTGCAATTCAGAGCGGCAGCAAGTGGGGGACAGCAGAAGACCTG"//&
        "ACCGCCGCAGAGTGGATGTTTGACATGGTGAAGACTATCGCACCATCAGCCAGAAAACCGAATTTTGCTG"//&
        "GGTGGGCTAACGATATCCGCCTGATGCGTGAACGTGACGGACGTAACCACCGCGACATGTGTGTGCTGTT"//&
        "CCGCTGGGCATGCCAGGACAACTTCTGGTCCGGTAACGTGCTGAGCCCGGCCAAACTCCGCGATAAGTGG"//&
        "ACCCAACTCGAAATCAACCGTAACAAGCAACAGGCAGGCGTGACAGCCAGCAAACCAAAACTCGACCTGA"//&
        "CAAACACAGACTGGATTTACGGGGTGGATCTATGAAAAACATCGCCGCACAGATGGTTAACTTTGACCGT"//&
        "GAGCAGATGCGTCGGATCGCCAACAACATGCCGGAACAGTACGACGAAAAGCCGCAGGTACAGCAGGTAG"//&
        "CGCAGATCATCAACGGTGTGTTCAGCCAGTTACTGGCAACTTTCCCGGCGAGCCTGGCTAACCGTGACCA"//&
        "GAACGAAGTGAACGAAATCCGTCGCCAGTGGGTTCTGGCTTTTCGGGAAAACGGGATCACCACGATGGAA"//&
        "CAGGTTAACGCAGGAATGCGCGTAGCCCGTCGGCAGAATCGACCATTTCTGCCATCACCCGGGCAGTTTG"//&
        "TTGCATGGTGCCGGGAAGAAGCATCCGTTACCGCCGGACTGCCAAACGTCAGCGAGCTGGTTGATATGGT"//&
        "TTACGAGTATTGCCGGAAGCGAGGCCTGTATCCGGATGCGGAGTCTTATCCGTGGAAATCAAACGCGCAC"//&
        "TACTGGCTGGTTACCAACCTGTATCAGAACATGCGGGCCAATGCGCTTACTGATGCGGAATTACGCCGTA"//&
        "AGGCCGCAGATGAGCTTGTCCATATGACTGCGAGAATTAACCGTGGTGAGGCGATCCCTGAACCAGTAAA"//&
        "ACAACTTCCTGTCATGGGCGGTAGACCTCTAAATCGTGCACAGGCTCTGGCGAAGATCGCAGAAATCAAA"//&
        "GCTAAGTTCGGACTGAAAGGAGCAAGTGTATGACGGGCAAAGAGGCAATTATTCATTACCTGGGGACGCA"//&
        "TAATAGCTTCTGTGCGCCGGACGTTGCCGCGCTAACAGGCGCAACAGTAACCAGCATAAATCAGGCCGCG"//&
        "GCTAAAATGGCACGGGCAGGTCTTCTGGTTATCGAAGGTAAGGTCTGGCGAACGGTGTATTACCGGTTTG"//&
        "CTACCAGGGAAGAACGGGAAGGAAAGATGAGCACGAACCTGGTTTTTAAGGAGTGTCGCCAGAGTGCCGC"//&
        "GATGAAACGGGTATTGGCGGTATATGGAGTTAAAAGATGACCATCTACATTACTGAGCTAATAACAGGCC"//&
        "TGCTGGTAATCGCAGGCCTTTTTATTTGGGGGAGAGGGAAGTCATGAAAAAACTAACCTTTGAAATTCGA"//&
        "TCTCCAGCACATCAGCAAAACGCTATTCACGCAGTACAGCAAATCCTTCCAGACCCAACCAAACCAATCG"//&
        "TAGTAACCATTCAGGAACGCAACCGCAGCTTAGACCAAAACAGGAAGCTATGGGCCTGCTTAGGTGACGT"//&
        "CTCTCGTCAGGTTGAATGGCATGGTCGCTGGCTGGATGCAGAAAGCTGGAAGTGTGTGTTTACCGCAGCA"//&
        "TTAAAGCAGCAGGATGTTGTTCCTAACCTTGCCGGGAATGGCTTTGTGGTAATAGGCCAGTCAACCAGCA"//&
        "GGATGCGTGTAGGCGAATTTGCGGAGCTATTAGAGCTTATACAGGCATTCGGTACAGAGCGTGGCGTTAA"//&
        "GTGGTCAGACGAAGCGAGACTGGCTCTGGAGTGGAAAGCGAGATGGGGAGACAGGGCTGCATGATAAATG"//&
        "TCGTTAGTTTCTCCGGTGGCAGGACGTCAGCATATTTGCTCTGGCTAATGGAGCAAAAGCGACGGGCAGG"//&
        "TAAAGACGTGCATTACGTTTTCATGGATACAGGTTGTGAACATCCAATGACATATCGGTTTGTCAGGGAA"//&
        "GTTGTGAAGTTCTGGGATATACCGCTCACCGTATTGCAGGTTGATATCAACCCGGAGCTTGGACAGCCAA"//&
        "ATGGTTATACGGTATGGGAACCAAAGGATATTCAGACGCGAATGCCTGTTCTGAAGCCATTTATCGATAT"//&
        "GGTAAAGAAATATGGCACTCCATACGTCGGCGGCGCGTTCTGCACTGACAGATTAAAACTCGTTCCCTTC"//&
        "ACCAAATACTGTGATGACCATTTCGGGCGAGGGAATTACACCACGTGGATTGGCATCAGAGCTGATGAAC"//&
        "CGAAGCGGCTAAAGCCAAAGCCTGGAATCAGATATCTTGCTGAACTGTCAGACTTTGAGAAGGAAGATAT"//&
        "CCTCGCATGGTGGAAGCAACAACCATTCGATTTGCAAATACCGGAACATCTCGGTAACTGCATATTCTGC"//&
        "ATTAAAAAATCAACGCAAAAAATCGGACTTGCCTGCAAAGATGAGGAGGGATTGCAGCGTGTTTTTAATG"//&
        "AGGTCATCACGGGATCCCATGTGCGTGACGGACATCGGGAAACGCCAAAGGAGATTATGTACCGAGGAAG"//&
        "AATGTCGCTGGACGGTATCGCGAAAATGTATTCAGAAAATGATTATCAAGCCCTGTATCAGGACATGGTA"//&
        "CGAGCTAAAAGATTCGATACCGGCTCTTGTTCTGAGTCATGCGAAATATTTGGAGGGCAGCTTGATTTCG"//&
        "ACTTCGGGAGGGAAGCTGCATGATGCGATGTTATCGGTGCGGTGAATGCAAAGAAGATAACCGCTTCCGA"//&
        "CCAAATCAACCTTACTGGAATCGATGGTGTCTCCGGTGTGAAAGAACACCAACAGGGGTGTTACCACTAC"//&
        "CGCAGGAAAAGGAGGACGTGTGGCGAGACAGCGACGAAGTATCACCGACATAATCTGCGAAAACTGCAAA"//&
        "TACCTTCCAACGAAACGCACCAGAAATAAACCCAAGCCAATCCCAAAAGAATCTGACGTAAAAACCTTCA"//&
        "ACTACACGGCTCACCTGTGGGATATCCGGTGGCTAAGACGTCGTGCGAGGAAAACAAGGTGATTGACCAA"//&
        "AATCGAAGTTACGAACAAGAAAGCGTCGAGCGAGCTTTAACGTGCGCTAACTGCGGTCAGAAGCTGCATG"//&
        "TGCTGGAAGTTCACGTGTGTGAGCACTGCTGCGCAGAACTGATGAGCGATCCGAATAGCTCGATGCACGA"//&
        "GGAAGAAGATGATGGCTAAACCAGCGCGAAGACGATGTAAAAACGATGAATGCCGGGAATGGTTTCACCC"//&
        "TGCATTCGCTAATCAGTGGTGGTGCTCTCCAGAGTGTGGAACCAAGATAGCACTCGAACGACGAAGTAAA"//&
        "GAACGCGAAAAAGCGGAAAAAGCAGCAGAGAAGAAACGACGACGAGAGGAGCAGAAACAGAAAGATAAAC"//&
        "TTAAGATTCGAAAACTCGCCTTAAAGCCCCGCAGTTACTGGATTAAACAAGCCCAACAAGCCGTAAACGC"//&
        "CTTCATCAGAGAAAGAGACCGCGACTTACCATGTATCTCGTGCGGAACGCTCACGTCTGCTCAGTGGGAT"//&
        "GCCGGACATTACCGGACAACTGCTGCGGCACCTCAACTCCGATTTAATGAACGCAATATTCACAAGCAAT"//&
        "GCGTGGTGTGCAACCAGCACAAAAGCGGAAATCTCGTTCCGTATCGCGTCGAACTGATTAGCCGCATCGG"//&
        "GCAGGAAGCAGTAGACGAAATCGAATCAAACCATAACCGCCATCGCTGGACTATCGAAGAGTGCAAGGCG"//&
        "ATCAAGGCAGAGTACCAACAGAAACTCAAAGACCTGCGAAATAGCAGAAGTGAGGCCGCATGACGTTCTC"//&
        "AGTAAAAACCATTCCAGACATGCTCGTTGAAACATACGGAAATCAGACAGAAGTAGCACGCAGACTGAAA"//&
        "TGTAGTCGCGGTACGGTCAGAAAATACGTTGATGATAAAGACGGGAAAATGCACGCCATCGTCAACGACG"//&
        "TTCTCATGGTTCATCGCGGATGGAGTGAAAGAGATGCGCTATTACGAAAAAATTGATGGCAGCAAATACC"//&
        "GAAATATTTGGGTAGTTGGCGATCTGCACGGATGCTACACGAACCTGATGAACAAACTGGATACGATTGG"//&
        "ATTCGACAACAAAAAAGACCTGCTTATCTCGGTGGGCGATTTGGTTGATCGTGGTGCAGAGAACGTTGAA"//&
        "TGCCTGGAATTAATCACATTCCCCTGGTTCAGAGCTGTACGTGGAAACCATGAGCAAATGATGATTGATG"//&
        "GCTTATCAGAGCGTGGAAACGTTAATCACTGGCTGCTTAATGGCGGTGGCTGGTTCTTTAATCTCGATTA"//&
        "CGACAAAGAAATTCTGGCTAAAGCTCTTGCCCATAAAGCAGATGAACTTCCGTTAATCATCGAACTGGTG"//&
        "AGCAAAGATAAAAAATATGTTATCTGCCACGCCGATTATCCCTTTGACGAATACGAGTTTGGAAAGCCAG"//&
        "TTGATCATCAGCAGGTAATCTGGAACCGCGAACGAATCAGCAACTCACAAAACGGGATCGTGAAAGAAAT"//&
        "CAAAGGCGCGGACACGTTCATCTTTGGTCATACGCCAGCAGTGAAACCACTCAAGTTTGCCAACCAAATG"//&
        "TATATCGATACCGGCGCAGTGTTCTGCGGAAACCTAACATTGATTCAGGTACAGGGAGAAGGCGCATGAG"//&
        "ACTCGAAAGCGTAGCTAAATTTCATTCGCCAAAAAGCCCGATGATGAGCGACTCACCACGGGCCACGGCT"//&
        "TCTGACTCTCTTTCCGGTACTGATGTGATGGCTGCTATGGGGATGGCGCAATCACAAGCCGGATTCGGTA"//&
        "TGGCTGCATTCTGCGGTAAGCACGAACTCAGCCAGAACGACAAACAAAAGGCTATCAACTATCTGATGCA"//&
        "ATTTGCACACAAGGTATCGGGGAAATACCGTGGTGTGGCAAAGCTTGAAGGAAATACTAAGGCAAAGGTA"//&
        "CTGCAAGTGCTCGCAACATTCGCTTATGCGGATTATTGCCGTAGTGCCGCGACGCCGGGGGCAAGATGCA"//&
        "GAGATTGCCATGGTACAGGCCGTGCGGTTGATATTGCCAAAACAGAGCTGTGGGGGAGAGTTGTCGAGAA"//&
        "AGAGTGCGGAAGATGCAAAGGCGTCGGCTATTCAAGGATGCCAGCAAGCGCAGCATATCGCGCTGTGACG"//&
        "ATGCTAATCCCAAACCTTACCCAACCCACCTGGTCACGCACTGTTAAGCCGCTGTATGACGCTCTGGTGG"//&
        "TGCAATGCCACAAAGAAGAGTCAATCGCAGACAACATTTTGAATGCGGTCACACGTTAGCAGCATGATTG"//&
        "CCACGGATGGCAACATATTAACGGCATGATATTGACTTATTGAATAAAATTGGGTAAATTTGACTCAACG"//&
        "ATGGGTTAATTCGCTCGTTGTGGTAGTGAGATGAAAAGAGGCGGCGCTTACTACCGATTCCGCCTAGTTG"//&
        "GTCACTTCGACGTATCGTCTGGAACTCCAACCATCGCAGGCAGAGAGGTCTGCAAAATGCAATCCCGAAA"//&
        "CAGTTCGCAGGTAATAGTTAGAGCCTGCATAACGGTTTCGGGATTTTTTATATCTGCACAACAGGTAAGA"//&
        "GCATTGAGTCGATAATCGTGAAGAGTCGGCGAGCCTGGTTAGCCAGTGCTCTTTCCGTTGTGCTGAATTA"//&
        "AGCGAATACCGGAAGCAGAACCGGATCACCAAATGCGTACAGGCGTCATCGCCGCCCAGCAACAGCACAA"//&
        "CCCAAACTGAGCCGTAGCCACTGTCTGTCCTGAATTCATTAGTAATAGTTACGCTGCGGCCTTTTACACA"//&
        "TGACCTTCGTGAAAGCGGGTGGCAGGAGGTCGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGG"//&
        "TCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAA"//&
        "TTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGC"//&
        "CATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGA"//&
        "TATAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTC"//&
        "GTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACAT"//&
        "CGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGA"//&
        "AATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAG"//&
        "AAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTC"//&
        "GCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTG"//&
        "GTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTG"//&
        "CAGCAGATTAAGGAGCGTGGCGCTTTACCTATGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTT"//&
        "GCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGAT"//&
        "TGCAAAATTCAAAGAAGCGGGCGGAACGGTCAGAGAGATTGATGTATGAGCAGAGTCACCGCGATTATCT"//&
        "CCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTAC"//&
        "CTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAG"//&
        "ATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAA"//&
        "ATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGT"//&
        "GCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGG"//&
        "GATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTA"//&
        "ATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGC"//&
        "GCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTT"//&
        "CTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTCTGCCCAATTCC"//&
        "AGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAAC"//&
        "GTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTAT"//&
        "TCCCGATGCTTTTTGAAGTTCGCAGAATCGTATGTGTAGAAAATTAAACAAACCCTAAACAATGAGTTGA"//&
        "AATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAAC"//&
        "TATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGG"//&
        "GTTTAGCGCGTACACGTATTGCATTATGCCAACGCCCCGGTGCTGACACGGAAGAAACCGGACGTTATGA"//&
        "TTTAGCGTGGAAAGATTTGTGTAGTGTTCTGAATGCTCTCAGTAAATAGTAATGAATTATCAAAGGTATA"//&
        "GTAATATCTTTTATGTTCATGGATATTTGTAACCCATCGGAAAACTCCTGCTTTAGCAAGATTTTCCCTG"//&
        "TATTGCTGAAATGTGATTTCTCTTGATTTCAACCTATCATAGGACGTTTCTATAAGATGCGTGTTTCTTG"//&
        "AGAATTTAACATTTACAACCTTTTTAAGTCCTTTTATTAACACGGTGTTATCGTTTTCTAACACGATGTG"//&
        "AATATTATCTGTGGCTAGATAGTAAATATAATGTGAGACGTTGTGACGTTTTAGTTCAGAATAAAACAAT"//&
        "TCACAGTCTAAATCTTTTCGCACTTGATCGAATATTTCTTTAAAAATGGCAACCTGAGCCATTGGTAAAA"//&
        "CCTTCCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAATCTGGTCTGACCTCCT"//&
        "TGTGTTTTGTTGATGATTTATGTCAAATATTAGGAATGTTTTCACTTAATAGTATTGGTTGCGTAACAAA"//&
        "GTGCGGTCCTGCTGGCATTCTGGAGGGAAATACAACCGACAGATGTATGTAAGGCCAACGTGCTCAAATC"//&
        "TTCATACAGAAAGATTTGAAGTAATATTTTAACCGCTAGATGAAGAGCAAGCGCATGGAGCGACAAAATG"//&
        "AATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAAAAAATATGCTTAATAGCACCA"//&
        "TTTCTATGAGTTACCCTGATGTTGTAATTGCATGTATAGAACATAAGGTGTCTCTGGAAGCATTCAGAGC"//&
        "AATTGAGGCAGCGTTGGTGAAGCACGATAATAATATGAAGGATTATTCCCTGGTGGTTGACTGATCACCA"//&
        "TAACTGCTAATCATTCAAACTATTTAGTCTGTGACAGAGCCAACACGCAGTCTGTCACTGTCAGGAAAGT"//&
        "GGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGA"//&
        "GGGAAGAACGCGGGATGTTCATTCTTCATCACTTTTAATTGATGTATATGCTCTCTTTTCTGACGTTAGT"//&
        "CTCCGACGGCAGGCTTCAATGACCCAGGCTGAGAAATTCCCGGACCCTTTTTGCTCAAGAGCGATGTTAA"//&
        "TTTGTTCAATCATTTGGTTAGGAAAGCGGATGTTGCGGGTTGTTGTTCTGCGGGTTCTGTTCTTCGTTGA"//&
        "CATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGTTAAGTTGATGC"//&
        "AGATCAATTAATACGATACCTGCGTCATAATTGATTATTTGACGTGGTTTGATGGCCTCCACGCACGTTG"//&
        "TGATATGTAGATGATAATCATTATCACTTTACGGGTCCTTTCCGGTGATCCGACAGGTTACG"

    ! Check Lamda sequence data
    do i = 1, len_Lamda
        if(Lamda_seq(i:i) == "N") then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error ========================================+"
            write(p_redir, "(a)"), " | Wrong assigned sequences in Lamda.               |"
            write(p_redir, "(a)"), " +==================================================+"
            stop
        end if
    end do
end function SeqDesign_Get_Lamda

! -----------------------------------------------------------------------------

! Import sequence from txt file
subroutine SeqDesign_Import_Sequence(prob, dna)
    type(ProbType), intent(in)    :: prob
    type(DNAType),  intent(inout) :: dna

    integer :: i, j, count, len_seq, n_base_scaf, base, across

    ! Get the sequence length
    len_seq = len(trim(prob.scaf_seq))

    if(len_seq < dna.n_base_scaf) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | User-defined sequences are short                 |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Print progress
    write(p_redir, "(a)"), "  6.6. Set sequence from file"
    write(p_redir, "(a)"), "   * # of sequence from file : "//trim(adjustl(Int2Str(len_seq)))
    write(p_redir, "(a)"), "   * # of scaffold           : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a)"), "   * # of nts in scaffold    : "//trim(adjustl(Int2Str(dna.n_base_scaf)))

    ! Check the number of bases in scaffold strands
    n_base_scaf = 0
    do i = 1, dna.n_scaf
        n_base_scaf = n_base_scaf + dna.strand(i).n_base
    end do

    if(n_base_scaf /= dna.n_base_scaf) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | # of nts in scaffold is not consistent.          |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Set scaffold sequence from file
    do i = 1, dna.n_strand

        ! Set sequence for scaffold
        if(dna.strand(i).type1 /= "scaf") cycle

        ! Find the starting point (down == -1)
        base   = Mani_Go_Start_Base(dna, i)
        across = dna.top(base).across

        do j = 1, dna.strand(i).n_base

            count = j+para_set_start_scaf-1
            if(count >= len_seq) then
                count = mod(count, len_seq)
                if(count == 0) count = len_seq
            end if

            ! Assign sequence from user-defined scaffold sequence
            dna.top(base).seq = prob.scaf_seq(count:count)

            ! Set complementary sequence
            if(across /= -1) then
                dna.top(across).seq = SeqDesign_Get_Comp_Sequence(dna.top(base).seq)
            end if

            ! Update base
            if(j /= dna.strand(i).n_base) then
                base   = dna.top(base).up
                across = dna.top(base).across
            end if
        end do
    end do
end subroutine SeqDesign_Import_Sequence

! -----------------------------------------------------------------------------

! Set random sequence
subroutine SeqDesign_Set_Rand_Sequence(dna)
    type(DNAType), intent(inout) :: dna

    integer :: i, j, len_seq, n_base_scaf, base, across

    ! Print progress
    write(p_redir, "(a )"), "6.8. Set random sequence"
    write(p_redir, "(a$)"), "* # of scaffolds       : "//trim(adjustl(Int2Str(dna.n_scaf)))
    write(p_redir, "(a$)"), "* # of nts in scaffold : "//trim(adjustl(Int2Str(dna.n_base_scaf)))
    write(p_redir, "(a )")

    ! Check the number of bases in scaffold strands
    n_base_scaf = 0
    do i = 1, dna.n_scaf
        n_base_scaf = n_base_scaf + dna.strand(i).n_base
    end do

    if(n_base_scaf /= dna.n_base_scaf) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | # of nts in scaffold is not consistent.          |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Set scaffold sequence from file
    do i = 1, dna.n_strand

        ! For scaffold strand
        if(dna.strand(i).type1 /= "scaf") cycle

        do j = 1, dna.strand(i).n_base

            base   = dna.strand(i).base(j)
            across = dna.top(base).across

            ! Assign random sequence
            dna.top(base).seq = SeqDesign_Get_Rand_Sequence()

            ! Set complementary sequence
            if(across /= -1) then
                dna.top(across).seq = SeqDesign_Get_Comp_Sequence(dna.top(base).seq)
            end if
        end do
    end do
end subroutine SeqDesign_Set_Rand_Sequence

! -----------------------------------------------------------------------------

! Write atom model by dnaTop and strand data
subroutine SeqDesign_Chimera_Atom(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    double precision :: pos_1(3), pos_2(3)
    integer :: i, j, base, up, xover, across
    logical :: f_axis
    character(200) :: path

    if(para_write_702 == .false.) return

    f_axis = para_chimera_axis

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 702, file = trim(path)//"_09_atomic_model.bild", form = "formatted")

    ! For all bases
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base

            if(dna.strand(i).type1 == "scaf") write(702, "(a)"), ".color steel blue"
            if(dna.strand(i).type1 == "stap") write(702, "(a)"), ".color orange"

            ! Draw bases
            base = dna.strand(i).base(j)
            write(702, "(a$    )"), ".sphere "
            write(702, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(702, "(1f9.3 )"), 0.15d0

            ! Draw backbones
            up = dna.top(base).up
            if(up /= -1) then
                write(702, "(a$    )"), ".cylinder "
                write(702, "(3f9.3$)"), dna.top(base).pos(1:3)
                write(702, "(3f9.3$)"), dna.top(up).pos(1:3)
                write(702, "(1f9.3 )"), 0.05d0
            end if

            ! Draw crossovers
            xover = dna.top(base).xover
            if(xover /= -1 .and. base < xover) then

                if(dna.strand(i).type1 == "scaf") write(702, "(a)"), ".color blue"
                if(dna.strand(i).type1 == "stap") write(702, "(a)"), ".color red"

                pos_1(:) = dna.top(base).pos(1:3)
                pos_2(:) = dna.top(xover).pos(1:3)

                write(702, "(a$    )"), ".cylinder "
                write(702, "(3f9.3$)"), pos_1(1:3)
                write(702, "(3f9.3$)"), pos_2(1:3)
                write(702, "(1f9.3 )"), 0.08d0
            end if

            ! Draw the Watson-Crick connections
            across = dna.top(base).across
            if(across > 0) then
                write(702, "(a     )"), ".color light gray"
                write(702, "(a$    )"), ".cylinder "
                write(702, "(3f9.3$)"), dna.top(base).pos(1:3)
                write(702, "(3f9.3$)"), dna.top(across).pos(1:3)
                write(702, "(1f9.3 )"), 0.025d0
            end if
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(702)
    close(unit = 702)

    ! ---------------------------------------------
    ! Write the file for Tecplot
    ! ---------------------------------------------
    if(para_tecplot == .false.) return

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit = 702, file = trim(path)//"_09_atomic_model.dat", form = "formatted")

    write(702, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'

    ! For bases in scaffold strands
    do i = 1, dna.n_strand

        write(702, "(a )"), 'VARIABLES = "X", "Y", "Z", "Weight"'
        write(702, "(a$)"), 'ZONE F = FEPOINT'
        write(702, "(a$)"), ', N='//trim(adjustl(Int2Str(dna.strand(i).n_base)))
        write(702, "(a$)"), ', E='//trim(adjustl(Int2Str(dna.strand(i).n_base - 1)))
        write(702, "(a )"), ', ET=LINESEG'

        ! Draw bases
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            write(702, "(4f9.3)"), dna.top(base).pos(1:3), 1.0d0
        end do

        ! Write elements
        do j = 1, dna.strand(i).n_base - 1
            write(702, "(2i7)"), j, j + 1
        end do
    end do

    write(702, "(a )"), 'VARIABLES = "X", "Y", "Z", "Weight"'
    write(702, "(a$)"), 'ZONE F = FEPOINT'
    write(702, "(a$)"), ', N='//trim(adjustl(Int2Str(dna.n_base_scaf*2)))
    write(702, "(a$)"), ', E='//trim(adjustl(Int2Str(dna.n_base_scaf)))
    write(702, "(a )"), ', ET=LINESEG'

    ! For bases in scaffold strands
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base
            if(dna.strand(i).type1 == "scaf") then
                base   = dna.strand(i).base(j)
                across = dna.top(base).across

                if(across > 0) then
                    write(702, "(4f9.3)"), dna.top(base).pos(1:3),   1.0d0
                    write(702, "(4f9.3)"), dna.top(across).pos(1:3), 1.0d0
                end if
            end if
        end do
    end do

    ! Write elements
    do i = 1, dna.n_base_scaf
        write(702, "(2i7)"), 2*i-1, 2*i
    end do

    close(unit = 702)
end subroutine SeqDesign_Chimera_Atom

! -----------------------------------------------------------------------------

! Write scaffold or staple route for Chimera
subroutine SeqDesign_Chimera_Route(prob, geom, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    double precision, allocatable :: base_scaf(:,:), base_stap(:,:)
    double precision :: pos_1(3), pos_2(3), vec(3)
    integer :: i, j, k, cur_base, up_base, down_base, cur_node, up_node
    integer :: n_base_scaf, n_base_stap, gap_layer = 30
    logical :: f_axis, f_opt
    character(200) :: path

    ! For option
    f_opt = .true.

    if(para_write_703 == .false.) return

    ! Exception for Tecplot drawing
    if(para_tecplot == .true.) then
        allocate(base_scaf (dna.n_base_scaf*2, 3 + 1 + 1))
        allocate(base_stap (dna.n_base_stap*2, 3 + 1 + 1))
        n_base_scaf = 0
        n_base_stap = 0
    end if

    ! Set option
    f_axis = para_chimera_axis

    ! File open for route step
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 703, file = trim(path)//"_10_routing_scaf.bild", form = "formatted")
    open(unit = 704, file = trim(path)//"_11_routing_stap.bild", form = "formatted")

    ! --------------------------------------------------
    ! For scaffold strand
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! Only for the scaffold
        if(dna.strand(i).type1 /= "scaf") cycle

        ! Loop for bases
        do j = 1, dna.strand(i).n_base

            ! Find the base
            cur_base  = dna.strand(i).base(j)
            up_base   = dna.top(cur_base).up
            down_base = dna.top(cur_base).dn

            ! Find the node for cur_base and up_base
            cur_node = dna.top(cur_base).node

            ! For the unpaired nucleotide of the scaffold
            if(cur_node == -1) cycle

            ! If the upper nt is not the nick
            if(up_base /= -1) then

                ! Find the node for up_base
                up_node = dna.top(up_base).node
            else

                ! Draw the end point of the nick
                if(f_opt == .true.) then
                    write(703, "(a     )"), ".color dark green"
                    write(703, "(a$    )"), ".sphere "
                    write(703, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                    write(703, "(1f9.3 )"), 0.2d0
                end if
                cycle
            end if

            ! Draw starting point with the arrow
            if(down_base == -1 .and. cur_node /= -1 .and. up_node /= -1) then
                pos_1(1:3) = mesh.node(cur_node).pos(1:3)
                pos_2(1:3) = mesh.node(up_node).pos(1:3)
                vec(1:3)   = Normalize(pos_2 - pos_1)

                if(f_opt == .true.) then
                    write(703, "(a     )"), ".color dark green"
                    write(703, "(a$    )"), ".arrow "
                    write(703, "(3f8.2$)"), pos_1(1:3)
                    write(703, "(3f8.2$)"), pos_1(1:3) + vec(1:3) * 1.3d0
                    write(703, "(2f8.2 )"), 0.11d0, 0.35d0
                end if
            end if

            ! For Tn loop, crossover and scaffold
            if(mesh.node(cur_node).up == -1) then

                ! For Tn loop with tan
                do
                    if(dna.top(up_base).node /= -1) exit
                    up_base = dna.top(up_base).up
                end do
                up_node = dna.top(up_base).node

                if(dna.n_scaf == 1) then
                    if(f_opt == .true.)  write(703, "(a)"), ".color tan"
                    if(f_opt == .false.) write(703, "(a)"), ".color steel blue"
                else
                    if(i == 1) write(703, "(a)"), ".color steel blue"
                    if(i == 2) write(703, "(a)"), ".color magenta"
                    if(i == 3) write(703, "(a)"), ".color light sea green"
                    if(i == 4) write(703, "(a)"), ".color sienna"
                    if(i == 5) write(703, "(a)"), ".color orange red"
                    if(i == 6) write(703, "(a)"), ".color dark slate gray"
                    if(i == 7) write(703, "(a)"), ".color dark slate blue"
                end if
            else

                if(dna.top(cur_base).xover == up_base) then
                    ! For crossovers with red
                    if(dna.n_scaf == 1) then
                        write(703, "(a)"), ".color red"
                    else
                        if(i == 1) write(703, "(a)"), ".color steel blue"
                        if(i == 2) write(703, "(a)"), ".color magenta"
                        if(i == 3) write(703, "(a)"), ".color light sea green"
                        if(i == 4) write(703, "(a)"), ".color sienna"
                        if(i == 5) write(703, "(a)"), ".color orange red"
                        if(i == 6) write(703, "(a)"), ".color dark slate gray"
                        if(i == 7) write(703, "(a)"), ".color dark slate blue"
                    end if
                else
                    ! For scaffold with blue
                    if(dna.n_scaf == 1) then
                        write(703, "(a)"), ".color steel blue"
                    else
                        if(i == 1) write(703, "(a)"), ".color steel blue"
                        if(i == 2) write(703, "(a)"), ".color magenta"
                        if(i == 3) write(703, "(a)"), ".color light sea green"
                        if(i == 4) write(703, "(a)"), ".color sienna"
                        if(i == 5) write(703, "(a)"), ".color orange red"
                        if(i == 6) write(703, "(a)"), ".color dark slate gray"
                        if(i == 7) write(703, "(a)"), ".color dark slate blue"
                    end if
                end if
            end if

            ! Draw route path
            if(Is_Same_Vector(mesh.node(cur_node).pos, mesh.node(up_node).pos) == .false.) then
                write(703, "(a$    )"), ".cylinder "
                write(703, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                write(703, "(3f9.3$)"), mesh.node(up_node ).pos(1:3)
                write(703, "(1f9.3 )"), 0.1d0
            end if

            if(para_tecplot == .true.) then

                base_scaf(n_base_scaf + 1, 4) = 0
                base_scaf(n_base_scaf + 2, 4) = 0
                base_scaf(n_base_scaf + 1, 5) = 1
                base_scaf(n_base_scaf + 2, 5) = 1

                base_scaf(n_base_scaf + 1, 1:3) = mesh.node(cur_node).pos(1:3)
                base_scaf(n_base_scaf + 2, 1:3) = mesh.node(up_node ).pos(1:3)
                n_base_scaf = n_base_scaf + 2
            end if
        end do
    end do

    ! --------------------------------------------------
    ! For bases in the staple strand
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! For staple strand
        if(dna.strand(i).type1 /= "stap") cycle

        ! Loop for strand
        do j = 1, dna.strand(i).n_base

            ! Find the base
            cur_base  = dna.strand(i).base(j)
            up_base   = dna.top(cur_base).up
            down_base = dna.top(cur_base).dn

            ! Find the node for cur_base ann up_base
            cur_node = dna.top(cur_base).node

            ! For the unpaired nucleotide of staples
            if(cur_node == -1) cycle

            if(up_base /= -1) then

                ! Find the node for up_base
                up_node = dna.top(up_base).node
            else

                ! Draw end point of the nick
                if(f_opt == .true.) then
                    write(704, "(a     )"), ".color dark green"
                    write(704, "(a$    )"), ".sphere "
                    write(704, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                    write(704, "(1f9.3 )"), 0.2d0
                end if
                cycle
            end if

            ! Draw starting point with the arrow
            if(down_base == -1 .and. cur_node /= -1 .and. up_node /= -1) then
                pos_1(1:3) = mesh.node(cur_node).pos(1:3)
                pos_2(1:3) = mesh.node(up_node).pos(1:3)
                vec(1:3)   = Normalize(pos_2 - pos_1)

                if(f_opt == .true.) then
                    write(704, "(a     )"), ".color dark green"
                    write(704, "(a$    )"), ".arrow "
                    write(704, "(3f8.2$)"), pos_1(1:3)
                    write(704, "(3f8.2$)"), pos_1(1:3) + vec(1:3) * 1.3d0
                    write(704, "(2f8.2 )"), 0.11d0, 0.35d0
                end if
            end if

            ! For Tn loop, crossover and staples
            if(mesh.node(cur_node).dn == -1) then

                ! For Tn loop with tan
                do
                    if(dna.top(up_base).node /= -1) exit
                    up_base = dna.top(up_base).up
                end do
                up_node = dna.top(up_base).node
                if(f_opt == .true.)  write(704, "(a)"), ".color tan"
                if(f_opt == .false.) write(704, "(a)"), ".color orange"
            else

                if(dna.top(cur_base).xover == up_base) then
                    ! For crossovers with red
                    write(704, "(a)"), ".color red"
                else
                    ! For scaffold with blue
                    write(704, "(a)"), ".color orange"
                end if
            end if

            ! Draw route path
            if(Is_Same_Vector(mesh.node(cur_node).pos, mesh.node(up_node).pos) == .false.) then
                write(704, "(a$    )"), ".cylinder "
                write(704, "(3f9.3$)"), mesh.node(cur_node).pos(1:3)
                write(704, "(3f9.3$)"), mesh.node(up_node).pos(1:3)
                write(704, "(1f9.3 )"), 0.1d0
            end if

            if(para_tecplot == .true.) then
                base_stap(n_base_stap + 1, 4) = 0
                base_stap(n_base_stap + 2, 4) = 0
                base_stap(n_base_stap + 1, 5) = 1
                base_stap(n_base_stap + 2, 5) = 1

                base_stap(n_base_stap + 1, 1:3) = mesh.node(cur_node).pos(1:3)
                base_stap(n_base_stap + 2, 1:3) = mesh.node(up_node).pos(1:3)
                n_base_stap = n_base_stap + 2
            end if
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(703)
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(704)
    close(unit = 703)
    close(unit = 704)

    ! --------------------------------------------------
    ! For Tecplot
    ! --------------------------------------------------
    if(para_tecplot == .false.) then
        if(allocated(base_scaf)) deallocate(base_scaf)
        if(allocated(base_stap)) deallocate(base_stap)
        return
    end if

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit = 703, file = trim(path)//"_10_routing_scaf.dat", form = "formatted")
    open(unit = 704, file = trim(path)//"_11_routing_stap.dat", form = "formatted")

    ! For scaffold bases
    write(703, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(703, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(703, "(a$)"), 'ZONE F = FEPOINT'
    write(703, "(a$)"), ', N='//trim(adjustl(Int2Str(n_base_scaf)))
    write(703, "(a$)"), ', E='//trim(adjustl(Int2Str(n_base_scaf/2)))
    write(703, "(a )"), ', ET=LINESEG'

    do i = 1, n_base_scaf
        write(703, "(2f9.3$)"), base_scaf(i, 1:2)
        write(703, "(1f9.3$)"), base_scaf(i, 3) + dble(gap_layer) * (base_scaf(i, 5) - 1)
        write(703, "(1f9.3 )"), base_scaf(i, 4)
    end do
    do i = 1, n_base_scaf/2
        write(703, "(2i7)"), 2*i, 2*i-1
    end do

    ! For crossover for scaffold strand
    write(704, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(704, "(a$)"), 'ZONE F = FEPOINT'
    write(704, "(a$)"), ', N='//trim(adjustl(Int2Str(n_base_stap)))
    write(704, "(a$)"), ', E='//trim(adjustl(Int2Str(n_base_stap/2)))
    write(704, "(a )"), ', ET=LINESEG'

    do i = 1, n_base_stap
        write(704, "(2f9.3$)"), base_stap(i, 1:2)
        write(704, "(1f9.3$)"), base_stap(i, 3) + dble(gap_layer) * (base_stap(i, 5) - 1)
        write(704, "(1f9.3 )"), base_stap(i, 4)
    end do
    do i = 1, n_base_stap/2
        write(704, "(2i7)"), 2*i, 2*i-1
    end do

    if(allocated(base_scaf)) deallocate(base_scaf)
    if(allocated(base_stap)) deallocate(base_stap)

    close(unit = 703)
    close(unit = 704)
end subroutine SeqDesign_Chimera_Route

! -----------------------------------------------------------------------------

! Chimera sequence design
subroutine SeqDesign_Chimera_Sequence_Design(prob, geom, mesh, dna)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(in)    :: dna

    double precision, allocatable :: base_scaf(:,:), base_stap(:,:)
    integer,          allocatable :: col_stap(:)

    double precision :: vec_n(3), vec_jn(3), vec_jt1(3), vec_jt2(3)
    double precision :: pos_1(3), pos_2(3), vec(3), RGB(3)
    integer :: i, j, k, base, up_base, node, up_node
    integer :: croL, iniL, count_scaf, count_stap
    logical :: f_axis, f_unpaired
    character(15)  :: col_list(16)
    character(200) :: path

    if(para_write_705 == .false.) return

    col_list(1:4)   = ["tan",             "salmon",       "orange",        "gold"           ]
    col_list(5:8)   = ["dark green",      "dark cyan",    "medium purple", "rosy brown"     ]
    col_list(9:12)  = ["dark slate gray", "dark magenta", "sea green",     "olive drab"     ]
    col_list(13:16) = ["goldenrod",       "firebrick",    "sienna",        "dark slate blue"]

    ! Set option
    f_axis = para_chimera_axis

    ! File open for sequence design
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 705, file = trim(path)//"_12_routing_all.bild", form = "formatted")

    if(para_tecplot == .true.) then
        allocate(base_scaf(dna.n_base_scaf*2, 3))
        allocate(base_stap(dna.n_base_stap*2, 4))
        allocate(col_stap(dna.n_base_stap*2))
        count_scaf = 0
        count_stap = 0
    end if

    ! --------------------------------------------------
    ! Nucleotides of staples
    ! --------------------------------------------------
    do i = 1, geom.n_iniL
        geom.iniL(i).n_xover = 0
    end do

    do i = 1, dna.n_base_scaf
        if(dna.top(i).xover /= -1) then
            geom.iniL(mesh.node(dna.top(i).node).iniL).n_xover &
                = geom.iniL(mesh.node(dna.top(i).node).iniL).n_xover + 1
        end if
    end do

    ! --------------------------------------------------
    ! For bases of the scaffold strand
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! Only for scaffold strand
        if(dna.strand(i).type1 /= "scaf") cycle

        do j = 1, dna.strand(i).n_base

            ! Find the current and up base
            base    = dna.strand(i).base(j)
            up_base = dna.top(base).up

            ! Skip if there is no upward base
            if(up_base == -1) cycle

            ! Find the node
            node    = dna.top(base).node
            up_node = dna.top(up_base).node

            ! Draw bases in unpaired nucleotides
            f_unpaired = .false.
            if(node == -1 .or. up_node == -1) then
                if(node /= -1 .and. up_node == -1) then
                    do
                        if(up_base == -1) cycle
                        if(dna.top(up_base).node /= -1) exit
                        up_base    = dna.top(up_base).up
                        f_unpaired = .true.
                    end do
                else
                    cycle
                end if

                if(up_base == -1) cycle
                up_node = dna.top(up_base).node
            end if

            ! Get cross-sectional line
            !pos_1(1:3) = 0.0d0
            !pos_2(1:3) = 0.0d0
            !do k = 1, geom.n_sec
            !    croL = mesh.node(node).croL
            !    pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
            !    pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            !end do

            ! Get cross-sectional line
            !iniL  = mesh.node(node).iniL
            !pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            !pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            !vec(1:3) = mesh.node(node).pos(1:3) - pos_2(1:3)
            !vec(1:3) = Normalize(vec(1:3))

            ! Vector projection
            !vec_n(1:3)   = Normalize(pos_2 - pos_1)
            !vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            !vec_jt1(1:3) = Normalize(vec - vec_jn)

            ! Get cross-sectional line #1
            !pos_1(1:3) = 0.0d0
            !pos_2(1:3) = 0.0d0
            !do k = 1, geom.n_sec
            !    croL = mesh.node(up_node).croL
            !    pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
            !    pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            !end do

            ! Get cross-sectional line #2
            !iniL  = mesh.node(up_node).iniL
            !pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            !pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            !vec(1:3) = mesh.node(up_node).pos(1:3) - pos_2(1:3)
            !vec(1:3) = Normalize(vec(1:3))

            ! Vector projection
            !vec_n(1:3)   = Normalize(pos_2 - pos_1)
            !vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            !vec_jt2(1:3) = Normalize(vec - vec_jn)

            ! Draw route path
            if(f_unpaired == .true. ) write(705, "(a)"), ".color red"
            if(f_unpaired == .false.) write(705, "(a)"), ".color steel blue"

            if(Is_Same_Vector(mesh.node(node).pos, mesh.node(up_node).pos) == .false.) then
                write(705, "(a$    )"), ".cylinder "
                write(705, "(3f9.3$)"), mesh.node(node).pos    !+ vec_jt1*0.5d0
                write(705, "(3f9.3$)"), mesh.node(up_node).pos !+ vec_jt2*0.5d0
                write(705, "(1f9.3 )"), 0.1d0
            end if

            ! Draw starting point
            if(dna.top(base).dn == -1) then
                write(705, "(a     )"), ".color red"
                write(705, "(a$    )"), ".sphere "
                write(705, "(3f9.3$)"), mesh.node(node).pos !+ vec_jt1*0.5d0
                write(705, "(1f9.3 )"), 0.4d0
                write(705, "(a     )"), ".color steel blue"
            end if

            if(para_tecplot == .true.) then
                base_scaf(count_scaf + 1, 1:3) = mesh.node(node).pos    !+ vec_jt1*0.5d0
                base_scaf(count_scaf + 2, 1:3) = mesh.node(up_node).pos !+ vec_jt2*0.5d0
                count_scaf = count_scaf + 2
            end if
        end do
    end do

    ! --------------------------------------------------
    ! For bases of the staple strand
    ! --------------------------------------------------
    vec_jt1(1:3) = 0.0d0
    vec_jt2(1:3) = 0.0d0
    do i = 1, dna.n_strand

        ! Only for staple strand
        if(dna.strand(i).type1 /= "stap") cycle

        !RGB(1:3) = [imod(irand(), 256), imod(irand(), 256), imod(irand(), 256)]
        !write(705, "(a, 3f9.4)"), ".color ", dble(RGB(1:3))/255.0d0
        write(705, "(a)"), ".color "//trim(col_list(mod(i-1, 16) + 1))

        !if(i /= 2) cycle
        do j = 1, dna.strand(i).n_base

            ! Find the base
            base    = dna.strand(i).base(j)
            up_base = dna.top(base).up

            ! Ending point as arrow
            if(up_base == -1) then
                pos_1 = mesh.node(dna.top(dna.top(base).dn).node).pos
                pos_2 = mesh.node(dna.top(base).node).pos
                vec   = pos_1 - pos_2
                write(705, "(a$    )"), ".arrow "
                write(705, "(3f9.3$)"), pos_1 + 1.5d0*vec - vec_jt1*0.5d0
                write(705, "(3f9.3$)"), pos_2 - 0.9d0*vec - vec_jt2*0.5d0
                write(705, "(3f9.3 )"), 0.1d0, 0.3d0, 0.3d0
                cycle
            end if

            ! Find node
            node    = dna.top(base).node
            up_node = dna.top(up_base).node

            ! Draw bases in Tn loop
            if(node == -1 .or. up_node == -1) then
                if(node /= -1 .and. up_node == -1) then
                    do
                        if(up_base == -1) cycle
                        if(dna.top(up_base).node /= -1) exit
                        up_base = dna.top(up_base).up
                    end do
                else
                    cycle
                end if

                if(up_base == -1) cycle
                up_node = dna.top(up_base).node
            end if

            ! Get cross-sectional line
            pos_1(1:3) = 0.0d0
            pos_2(1:3) = 0.0d0
            do k = 1, geom.n_sec
                croL = mesh.node(node).croL
                pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
                pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            end do

            ! Get cross-sectional line
            iniL  = mesh.node(node).iniL
            pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            vec(1:3) = mesh.node(node).pos(1:3) - pos_2(1:3)
            vec(1:3) = Normalize(vec(1:3))

            ! Vector projection
            vec_n(1:3)   = Normalize(pos_2 - pos_1)
            vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            vec_jt1(1:3) = Normalize(vec - vec_jn)

            ! Get cross-sectional line #1
            !pos_1(1:3) = 0.0d0
            !pos_2(1:3) = 0.0d0
            !do k = 1, geom.n_sec
            !    croL = mesh.node(up_node).croL
            !    pos_1(1:3) = pos_1 + geom.croP(geom.croL(croL).poi(1)).pos
            !    pos_2(1:3) = pos_2 + geom.croP(geom.croL(croL).poi(2)).pos
            !end do

            ! Get cross-sectional line #2
            iniL  = mesh.node(up_node).iniL
            pos_1 = geom.modP(geom.iniL(iniL).poi(1)).pos(1:3)
            pos_2 = geom.modP(geom.iniL(iniL).poi(2)).pos(1:3)

            ! Find third local vector, t3
            vec(1:3) = mesh.node(up_node).pos(1:3) - pos_2(1:3)
            vec(1:3) = Normalize(vec(1:3))

            ! Vector projection
            vec_n(1:3)   = Normalize(pos_2 - pos_1)
            vec_jn(1:3)  = dot_product(vec, vec_n) * vec_n
            vec_jt2(1:3) = Normalize(vec - vec_jn)

            ! Draw route path
            if(Is_Same_Vector(mesh.node(node).pos, mesh.node(up_node).pos) == .false.) then
                write(705, "(a$    )"), ".cylinder "
                write(705, "(3f9.3$)"), mesh.node(node).pos(1:3)    - vec_jt1*0.5d0
                write(705, "(3f9.3$)"), mesh.node(up_node).pos(1:3) - vec_jt2*0.5d0
                write(705, "(1f9.3 )"), 0.1d0
            end if

            if(para_tecplot == .true.) then
                base_stap(count_stap + 1, 1:3) = mesh.node(node).pos(1:3)    - vec_jt1*0.5d0
                base_stap(count_stap + 2, 1:3) = mesh.node(up_node).pos(1:3) - vec_jt2*0.5d0

                if(mesh.node(node).iniL == mesh.node(up_node).iniL) then
                    col_stap(count_stap + 1) = geom.iniL(mesh.node(node).iniL).n_xover
                    col_stap(count_stap + 2) = geom.iniL(mesh.node(node).iniL).n_xover
                else
                    col_stap(count_stap + 1) = 20
                    col_stap(count_stap + 2) = 20
                end if

                count_stap = count_stap + 2
            end if
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(705)
    close(unit = 705)

    ! ---------------------------------------------
    ! Write the file for Tecplot
    ! ---------------------------------------------
    if(para_tecplot == .false.) return

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit = 705, file = trim(path)//"_12_routing_all.dat", form = "formatted")
    !open(unit = 706, file = trim(path)//"_15_sep_lines.dat", form = "formatted")

    ! For scaffold bases
    write(705, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(705, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(705, "(a$)"), 'ZONE F = FEPOINT'
    write(705, "(a$)"), ', N='//trim(adjustl(Int2Str(count_scaf)))
    write(705, "(a$)"), ', E='//trim(adjustl(Int2Str(count_scaf/2)))
    write(705, "(a )"), ', ET=LINESEG'

    do i = 1, count_scaf
        write(705, "(3f9.3$)"), base_scaf(i, 1:3)
        write(705, "(1f9.3 )"), 1.0d0
    end do
    do i = 1, count_scaf/2
        write(705, "(2i7)"), 2*i, 2*i-1
    end do

    ! For staple bases
    write(705, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(705, "(a$)"), 'ZONE F = FEPOINT'
    write(705, "(a$)"), ', N='//trim(adjustl(Int2Str(count_stap)))
    write(705, "(a$)"), ', E='//trim(adjustl(Int2Str(count_stap/2)))
    write(705, "(a )"), ', ET=LINESEG'

    do i = 1, count_stap
        write(705, "(3f9.3$)"), base_stap(i, 1:3)
        write(705, "(1f9.3 )"), dble(col_stap(i))
    end do
    do i = 1, count_stap/2
        write(705, "(2i7)"), 2*i, 2*i-1
    end do

    ! ---------------------------------------------
    ! New initial geometry
    ! ---------------------------------------------
    !write(706, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    !write(706, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    !write(706, "(a$)"), 'ZONE F = FEPOINT'
    !write(706, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_modP)))
    !write(706, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    !write(706, "(a )"), ', ET=LINESEG'
    !
    !! Write points
    !do i = 1, geom.n_iniL
    !    write(706, "(3f9.3$)"), geom.modP(geom.iniL(i).poi(1)).pos(1:3)
    !    write(706, "(1i7 )"), geom.iniL(i).n_xover
    !    write(706, "(3f9.3$)"), geom.modP(geom.iniL(i).poi(2)).pos(1:3)
    !    write(706, "(1i7 )"), geom.iniL(i).n_xover
    !end do
    !
    !! Write edges
    !do i = 1, geom.n_iniL
    !    write(706, "(1i7$)"), geom.iniL(i).poi(1)
    !    write(706, "(1i7 )"), geom.iniL(i).poi(2)
    !end do

    ! Deallocate memory
    if(allocated(base_scaf)) deallocate(base_scaf)
    if(allocated(base_stap)) deallocate(base_stap)
    if(allocated(col_stap))  deallocate(col_stap)

    close(unit = 705)
    !close(unit = 706)
end subroutine SeqDesign_Chimera_Sequence_Design

! -----------------------------------------------------------------------------

! Write Chimera file for strand and sequence
subroutine SeqDesign_Chimera_Strand(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    double precision :: vec(3), azure(3), tweetybird(3), green(3), carmine(3)
    integer :: i, j, k, base, down, up, across
    logical :: f_axis
    character(15)  :: col_list(16)
    character(200) :: path

    if(para_write_706 == .false.) return

    f_axis = para_chimera_axis

    col_list(1:4)   = ["tan",             "salmon",       "orange",        "gold"           ]
    col_list(5:8)   = ["dark green",      "dark cyan",    "medium purple", "rosy brown"     ]
    col_list(9:12)  = ["dark slate gray", "dark magenta", "sea green",     "olive drab"     ]
    col_list(13:16) = ["goldenrod",       "firebrick",    "sienna",        "dark slate blue"]

    ! Colors for four bases, reference: http://www.umass.edu/molvis/tutorials/dna/atgc.htm
    azure      = [0.000d0, 127.0d0, 255.0d0] / 255.0d0       ! color for A
    tweetybird = [253.0d0, 245.0d0, 1.000d0] / 255.0d0       ! color for T
    green      = [0.000d0, 255.0d0, 0.000d0] / 255.0d0       ! color for G
    carmine    = [150.0d0, 0.000d0, 24.00d0] / 255.0d0       ! color for C

    ! Write the file
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 706, file = trim(path)//"_strand.bild",   form = "formatted")
    open(unit = 707, file = trim(path)//"_sequence.bild", form = "formatted")

    ! Draw the structure using strand data
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base

            ! Find the base
            base = dna.strand(i).base(j)

            ! Color for strand
            write(706, "(a)"), ".color " // trim(col_list(mod(i-1, 16) + 1))

            ! Color for sequence
            if(dna.top(base).seq == "A") then
                write(707, "(a, 3f9.3)"), ".color ", azure(1:3)
            else if(dna.top(base).seq == "T") then
                write(707, "(a, 3f9.3)"), ".color ", tweetybird(1:3)
            else if(dna.top(base).seq == "G") then
                write(707, "(a, 3f9.3)"), ".color ", green(1:3)
            else if(dna.top(base).seq == "C") then
                write(707, "(a, 3f9.3)"), ".color ", carmine(1:3)
            else
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Not assigned sequence                            |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            write(706, "(a$    )"), ".sphere "
            write(706, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(706, "(1f9.3 )"), 0.15d0

            write(707, "(a$    )"), ".sphere "
            write(707, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(707, "(1f9.3 )"), 0.15d0

            do k = 0, 1
                ! Draw the backbones
                down = dna.top(base).dn
                if(down >= 0) then
                    write(706+k, "(a     )"), ".color light gray"
                    write(706+k, "(a$    )"), ".cylinder "
                    write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3)
                    write(706+k, "(3f9.3$)"), dna.top(down).pos(1:3)
                    write(706+k, "(1f9.3 )"), 0.05d0
                end if

                ! Draw the Watson-Crick connections
                across = dna.top(base).across
                if(across >= 0) then
                    write(706+k, "(a     )"), ".color light gray"
                    write(706+k, "(a$    )"), ".cylinder "
                    write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3)
                    write(706+k, "(3f9.3$)"), dna.top(across).pos(1:3)
                    write(706+k, "(1f9.3 )"), 0.025d0
                end if
            end do
        end do

        ! Draw the strand directionality
        if(dna.strand(i).n_base > 1) then
            ! Get last base
            base = dna.strand(i).base(dna.strand(i).n_base)

            if(dna.strand(i).b_circular == .true.) then
                down = dna.strand(i).base(1)
                vec  = Normalize(dna.top(down).pos - dna.top(base).pos)
            else
                up  = dna.strand(i).base(dna.strand(i).n_base - 1)
                vec = Normalize(dna.top(base).pos - dna.top(up).pos)
            end if
        else
            vec = [1.0d0, 0.0d0, 0.0d0]
        end if

        ! Color for strand
        write(706, "(a)"), ".color " // trim(col_list(mod(i-1, 16) + 1))

        ! Color for sequence
        if(dna.top(base).seq == "A") then
            write(707, "(a, 3f9.3)"), ".color ", azure(1:3)
        else if(dna.top(base).seq == "T") then
            write(707, "(a, 3f9.3)"), ".color ", tweetybird(1:3)
        else if(dna.top(base).seq == "G") then
            write(707, "(a, 3f9.3)"), ".color ", green(1:3)
        else if(dna.top(base).seq == "C") then
            write(707, "(a, 3f9.3)"), ".color ", carmine(1:3)
        else
            write(707, "(a)"), ".color light gray"
        end if

        do k = 0, 1
            write(706+k, "(a$    )"), ".cone "
            write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3)
            write(706+k, "(3f9.3$)"), dna.top(base).pos(1:3) + vec(1:3) * 0.36d0
            write(706+k, "(1f9.3 )"), 0.18d0
        end do
    end do

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(706)
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(707)
    close(unit = 706)
    close(unit = 707)
end subroutine SeqDesign_Chimera_Strand

! -----------------------------------------------------------------------------

end module SeqDesign