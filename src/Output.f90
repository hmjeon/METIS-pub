!
! =============================================================================
!
! Module - Output
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
module Output

    use Data_Mesh
    use Data_DNA

    use Mani
    use Para
    use Math

    implicit none

    public  Output_Generation

    private Output_Write_Cylinder_Xover

    private Output_Write_Out_All
    private Output_Write_Out_Sequences
    private Output_Write_Out_Unpaired
    private Output_Write_Out_Graphics
    private Output_Write_Out_Strand_Base
    private Output_Write_Out_Guide_JSON
    private Output_Write_Out_JSON

    private Output_Check_Output
    private Output_Write_Basepair
    private Output_Write_Base
    private Output_Write_CanDo
    private Output_Write_CanDo_New
    private Output_Write_PLY
    private Output_Write_DNA_Info
    private Output_Write_TecPlot
    private Output_Write_ADINA
    private Output_Write_Sequence_CroL

contains

! -----------------------------------------------------------------------------

! Generate output files and run postprocessing tools
subroutine Output_Generation(prob, geom, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: i

    ! Print progress
    write(p_redir, "(a)")
    write(p_redir, "(a)"), " +===================================================+"
    write(p_redir, "(a)"), " | 7. Outputs                                        |"
    write(p_redir, "(a)"), " +===================================================+"
    write(p_redir, "(a)")
    write(p_redir, "(a)"), "  7.1. Write output files"

    ! Check connectivity in dna data
    !call Output_Check_Output(dna)

    ! Write cylinderical model with crossovers
    call Output_Write_Cylinder_Xover(prob, geom, mesh, dna)

    ! Write outputs
    call Output_Write_Out_All(prob, geom, mesh, dna)

    ! Write information related with basepair
    call Output_Write_Basepair(prob, mesh, dna)

    ! Write information related with base
    call Output_Write_Base(prob, dna)

    ! Write cndo file for PDB atom generation and CanDo simulation
    call Output_Write_CanDo(prob, mesh, dna)
    !call Output_Write_CanDo_New(prob, mesh, dna)

    ! Write PLY file
    !call Output_Write_PLY(prob, geom)

    !call Output_Write_DNA_Info(prob, dna)

    ! Write TecPlot input file
    call Output_Write_TecPlot(prob, mesh)

    ! Write ADINA input file
    call Output_Write_ADINA(prob, mesh)

    ! Write sequence based on cross-sectional line
    call Output_Write_Sequence_CroL(prob, mesh, dna)
end subroutine Output_Generation

! -----------------------------------------------------------------------------

! Check connectivity in dna data
subroutine Output_Check_Output(dna)
    type(DNAType), intent(in) :: dna

    integer :: i, j, count, cur_base, across, xover
    integer, allocatable :: base(:)

    ! --------------------------------------------------
    ! Check the # of dnatop and bases in all strands
    ! --------------------------------------------------
    count = 0
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base
            count = count + 1
        end do
    end do

    ! Check the number of bases in all strands
    if(count /= dna.n_top) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error =========================================+"
        write(p_redir, "(a)"), " | # of nts in dna top and strand are not the same.  |"
        write(p_redir, "(a)"), " +===================================================+"
        stop
    end if

    ! --------------------------------------------------
    ! There should be non-circular strand and check # of bases
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! Find starting point
        cur_base = dna.strand(i).base(1)
        do
            if(dna.top(cur_base).dn == -1) exit
            cur_base = dna.top(cur_base).dn
        end do

        count = 1
        do
            if(dna.top(cur_base).up == -1) exit
            cur_base = dna.top(cur_base).up
            count = count + 1
        end do

        ! Check the number of bases in dna strand
        if(count /= dna.strand(i).n_base) then
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error =========================================+"
            write(p_redir, "(a)"), " | # of nts in dna strand are not consistent.        |"
            write(p_redir, "(a)"), " +===================================================+"
            stop
        end if
    end do

    ! --------------------------------------------------
    ! Check across
    ! --------------------------------------------------
    do i = 1, dna.n_top
        if(dna.top(i).across /= -1) then
            across = dna.top(i).across
            if(dna.top(across).across /= dna.top(i).id) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error =========================================+"
                write(p_redir, "(a)"), " | The across ID is not consistent.                  |"
                write(p_redir, "(a)"), " +===================================================+"
                stop
            end if
        end if
    end do

    ! --------------------------------------------------
    ! Check crossover
    ! --------------------------------------------------
    do i = 1, dna.n_top
        if(dna.top(i).xover /= -1) then
            xover = dna.top(i).xover
            if(dna.top(xover).xover /= dna.top(i).id) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error =========================================+"
                write(p_redir, "(a)"), " | The Xover ID is not consistent.                   |"
                write(p_redir, "(a)"), " +===================================================+"
                stop
            end if
        end if
    end do

    ! --------------------------------------------------
    ! Check whether strand has consistent direction (non-circular strand)
    ! --------------------------------------------------
    do i = 1, dna.n_strand

        ! Allocate base data to store entity
        allocate(base(dna.strand(i).n_base))

        ! Find starting point
        cur_base = dna.strand(i).base(1)
        do
            if(dna.top(cur_base).dn == -1) exit
            cur_base = dna.top(cur_base).dn
        end do

        ! Go base up and save base number in base data
        base(1) = cur_base
        do j = 2, dna.strand(i).n_base
            cur_base = dna.top(cur_base).up
            base(j)  = cur_base
        end do

        ! Go base down and check whether the same bases
        do j = 1, dna.strand(i).n_base
            if(base(dna.strand(i).n_base - j + 1) == cur_base) then
                cur_base = dna.top(cur_base).dn
            else
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error =========================================+"
                write(p_redir, "(a)"), " | The strand has not consistent direction.          |"
                write(p_redir, "(a)"), " +===================================================+"
                stop
            end if
        end do

        ! Deallocate base memory
        deallocate(base)
    end do
end subroutine Output_Check_Output

! -----------------------------------------------------------------------------

! Write cylinderical model with crossovers
subroutine Output_Write_Cylinder_Xover(prob, geom, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    double precision :: pos_1(3), pos_2(3), radius, rad_mod
    integer :: i, j, node
    logical :: f_axis
    character(200) :: path

    ! Set flag for drawing option
    f_axis = para_chimera_axis

    path = trim(prob.path_work)//"/"//trim(prob.name_file)//"_13_cylinder_xover"
    open(unit = 701, file = trim(path)//".bild", form = "formatted")

    ! Cylinder radius
    radius = para_rad_helix + para_gap_helix / 2.0d0

    ! Write cylinder model base on edges
    do i = 1, geom.n_croL

        ! Draw cylinder before the mitered design
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).ori_pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).ori_pos(1:3)

        write(701, "(a, 3f9.4)"), ".color ", dble(prob.color(1:3))/255.0d0
        write(701, "(a$    )"), ".cylinder "
        write(701, "(3f9.3$)"), pos_1(1:3)
        write(701, "(3f9.3$)"), pos_2(1:3)
        write(701, "(1f9.3 )"), radius

        ! Draw cylinder for the mitered parts
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(1)).ori_pos(1:3)

        if(Is_Same_Vector(pos_1, pos_2) == .false.) then
            if(Norm(pos_1 - pos_2) > 0.4d0) then
                write(701, "(a)"), ".color dark gray"
                rad_mod = radius * 1.02d0
            else
                write(701, "(a, 3f9.4)"), ".color ", dble(prob.color(1:3))/255.0d0
                rad_mod = radius
            end if
            write(701, "(a$     )"), ".cylinder "
            write(701, "(3f12.5$)"), pos_1(1:3)
            write(701, "(3f12.5$)"), pos_2(1:3)
            write(701, "(1f9.3  )"), rad_mod
        end if

        pos_1(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).ori_pos(1:3)

        if(Is_Same_Vector(pos_1, pos_2) ==.false.) then
            if(Norm(pos_1 - pos_2) > 0.4d0) then
                write(701, "(a)"), ".color dark gray"
                rad_mod = radius * 1.02d0
            else
                write(701, "(a, 3f9.4)"), ".color ", dble(prob.color(1:3))/255.0d0
                rad_mod = radius
            end if
            write(701, "(a$     )"), ".cylinder "
            write(701, "(3f12.5$)"), pos_1(1:3)
            write(701, "(3f12.5$)"), pos_2(1:3)
            write(701, "(1f9.3  )"), rad_mod
        end if
    end do

    ! Write centered scaffold crossovers
    do i = 1, dna.n_strand
        do j = 1, dna.strand(i).n_base
            if(dna.top(dna.strand(i).base(j)).xover /= -1) then

                node = dna.top(dna.strand(i).base(j)).node

                if(mesh.node(node).up == -1) then
                    pos_1(1:3) = mesh.node(node).pos(1:3)
                else
                    pos_1(1:3) = mesh.node(mesh.node(node).up).pos(1:3)
                end if

                if(mesh.node(node).dn == -1) then
                    pos_2(1:3) = mesh.node(node).pos(1:3)
                else
                    pos_2(1:3) = mesh.node(mesh.node(node).dn).pos(1:3)
                end if

                if(dna.strand(i).type1 == "scaf") write(701, "(a)"), ".color medium blue"
                if(dna.strand(i).type1 == "stap") write(701, "(a)"), ".color orange red"

                write(701, "(a$    )"), ".cylinder "
                write(701, "(3f9.3$)"), pos_1(1:3)
                write(701, "(3f9.3$)"), pos_2(1:3)
                write(701, "(1f9.3 )"), radius*1.05d0
            end if
        end do
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write cylindrical model with Xovers"

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(701)
    close(unit = 701)
end subroutine Output_Write_Cylinder_Xover

! -----------------------------------------------------------------------------

! Write output file for sequence and json
subroutine Output_Write_Out_All(prob, geom, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(inout) :: dna

    integer :: max_unpaired

    if(para_write_701 == .false.) return
    open(unit = 701, file = trim(prob.path_work)//"/TXT_Sequence.txt", form = "formatted")

    ! Write staple sequence results
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                      1. Sequence data                          +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    call Output_Write_Out_Sequences(prob, mesh, dna, 701)

    ! Write graphical routing
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                     2. Graphical routing                       +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    write(701, "(a)"), "1. [-], [num], [.], [>], [<] : Nucleotide"
    write(701, "(a)"), "2. [|]                       : Crossover"
    write(701, "(a)"), "3. [A], [T], [G], [C]        : Sequence"
    write(701, "(a)"), "4. [ a: b], c                : From starting base ID(a) to ending base ID(b), total # of bases (c)"
    write(701, "(a)"), "5. [5'], [3']                : Strand direction"
    write(701, "(a)")
    call Output_Write_Out_Graphics(prob, geom, mesh, dna, 701)

    ! Write information on unpaired nucleotides and poly Tn loop
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+    3. Information on unpaired nucleotides and poly Tn loop     +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    max_unpaired = Output_Write_Out_Unpaired(mesh, dna, 701)

    ! Outputs based on strands and nucleotides
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+             4. Strand and nucleotide based outputs             +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    call Output_Write_Out_Strand_Base(mesh, dna, 701)

    ! Output of staple length
    write(701, "(a)")
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|+                       5. Staple length                         +|"
    write(701, "(a)"), "|+                                                                +|"
    write(701, "(a)"), "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|"
    write(701, "(a)"), "--------------------------------------------------------------------"
    write(701, "(a)")
    call Output_Write_Out_Staple_Length(dna, 701)

    ! JSON output
    call Output_Write_Out_Guide_JSON(prob, geom, mesh)
    call Output_Write_Out_JSON(prob, geom, mesh, dna, max_unpaired)

    ! Print progress
    write(p_redir, "(a)"), "   * Write JSON and etc."

    close(unit = 701)
end subroutine Output_Write_Out_All

! -----------------------------------------------------------------------------

! Write output of sequence data
subroutine Output_Write_Out_Sequences(prob, mesh, dna, unit)
    type(ProbType), intent(in)    :: prob
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    integer,        intent(in)    :: unit

    integer :: i, ii, j, base, n_scaf, n_stap, dozen, add
    logical :: seq_info, seq_ordering, on_U, on_S, on_F, on_X
    character(200) :: path

    seq_info     = .true.
    seq_ordering = .true.

    if(seq_info == .true.) then
        write(unit, "(a)"), "   * - unpaired nucleotide"
        write(unit, "(a)"), "   $ - nucleotide in the crossover"
        write(unit, "(a)"), "   # - nucleotide in the 14nt dsDNA domain"
        write(unit, "(a)"), "   @ - nucleotide in the 4nt dsDNA domain"
        write(unit, "(a)")
    end if

    ! Update status for poly Tn loop and xovers
    do i = 1, dna.n_strand
        base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base
            if(dna.top(base).across == -1) dna.top(base).status = "U"
            if(dna.top(base).xover /= -1)  dna.top(base).status = "X"
            base = dna.Top(base).up
        end do
    end do

    n_scaf = 0
    n_stap = 0
    dna.len_min_stap =  10000
    dna.len_max_stap = -10000

    ! Write sequence data
    do i = 1, dna.n_strand
        on_U = .false.
        on_S = .false.
        on_F = .false.
        on_X = .false.

        if(dna.strand(i).type1 == "scaf") then
            n_scaf = n_scaf + 1
            ii     = i
            write(unit, "(i8, a$)"), n_scaf, " - [scaf]"
        else if(dna.strand(i).type1 == "stap") then
            n_stap = n_stap + 1
            if(seq_ordering == .true. ) ii = dna.order_stap(dna.n_stap - n_stap + 1, 1)
            if(seq_ordering == .false.) ii = i
            write(unit, "(i8, a$)"), n_stap, " - (stap)"
            if(dna.strand(ii).n_base < dna.len_min_stap) dna.len_min_stap = dna.strand(ii).n_base
            if(dna.strand(ii).n_base > dna.len_max_stap) dna.len_max_stap = dna.strand(ii).n_base
        end if

        if(0 .and. dna.strand(ii).n_base > 80) then
            do j = 0, unit, unit
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | There are long staples that are not short.       |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end do
        end if

        if(dna.strand(ii).type1 == "stap") then
            write(unit, "(a$)"), "-"//dna.strand(ii).type2
            write(unit, "(a$)"), ", # of 14nt seeds: "//trim(adjustl(Int2Str(dna.strand(ii).n_14nt)))
        end if
        write(unit, "(a$)"), ", # of nts: "//trim(adjustl(Int2Str(dna.strand(ii).n_base)))//" -> "

        base = Mani_Go_Start_Base(dna, ii)
        do j = 1, dna.strand(ii).n_base
            if(seq_info == .true.) then
                if(dna.top(base).status == "U" .and. on_U == .false.) then
                    write(unit, "(a$)"), "*"
                    on_U = .true.
                else if(dna.top(base).status /= "U" .and. on_U == .true.) then
                    write(unit, "(a$)"), "*"
                    on_U = .false.
                end if

                if(dna.top(base).status == "S" .and. on_S == .false.) then
                    write(unit, "(a$)"), "#"
                    on_S = .true.
                else if(dna.top(base).status /= "S" .and. on_S == .true.) then
                    write(unit, "(a$)"), "#"
                    on_S = .false.
                end if

                if(dna.top(base).status == "F" .and. on_F == .false.) then
                    write(unit, "(a$)"), "@"
                    on_F = .true.
                else if(dna.top(base).status /= "F" .and. on_F == .true.) then
                    write(unit, "(a$)"), "@"
                    on_F = .false.
                end if

                if(dna.top(base).status == "X" .and. on_X == .false.) then
                    write(unit, "(a$)"), "$"
                    on_X = .true.
                else if(dna.top(base).status /= "X" .and. on_X == .true.) then
                    write(unit, "(a$)"), "$"
                    on_X = .false.
                end if
            end if

            write(unit, "(a$)"), dna.top(base).seq
            base = dna.Top(base).up

            if(seq_info == .true.) then
                if(j == dna.strand(ii).n_base) then
                    if(on_U == .true.) then
                        write(unit, "(a$)"), "*"; on_U = .false.
                    end if
                    if(on_S == .true.) then
                        write(unit, "(a$)"), "#"; on_S = .false.
                    end if
                    if(on_F == .true.) then
                        write(unit, "(a$)"), "@"; on_F = .false.
                    end if
                    if(on_X == .true.) then
                        write(unit, "(a$)"), "$"; on_X = .false.
                    end if
                end if
            end if
        end do
        write(unit, "(a)")
        if(dna.strand(ii).type1 == "scaf") write(unit, "(a)")
    end do
    write(unit, "(a)")

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 702, file = trim(path)//".csv", form = "formatted")

    ! DNA sequence data according to strands
    write(702, "(a)")
    write(702, "(a)"), "No, Type1, Type2, A12, Strand name, Length, Sequence"
    write(702, "(a)")

    ! Write sequence data
    n_scaf = 0
    n_stap = 0
    add    = 0
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            n_scaf = n_scaf + 1
            write(702, "(a$)"), trim(adjustl(Int2Str(n_scaf)))//", "
            write(702, "(a$)"), "Scaffold, -, -, -, "
        else if(dna.strand(i).type1 == "stap") then
            n_stap = n_stap + 1
            write(702, "(a$)"), trim(adjustl(Int2Str(n_stap)))//", "
            write(702, "(a$)"), "Staple, "//dna.strand(i).type2//", "

            dozen = mod(n_stap, 12)
            if(dozen == 0) dozen = 12
            if(dozen == 1) add = add + 1
            write(702, "(a$)"), achar(65 + add - 1)//trim(adjustl(Int2Str(dozen)))//", "

            write(702, "(a$)"), trim(adjustl(prob.name_prob))//"_"
            if(para_vertex_design == "flat")    write(702, "(a$)"), "F_"
            if(para_vertex_design == "mitered") write(702, "(a$)"), "M_"
            write(702, "(a$)"), trim(adjustl(Int2Str(prob.n_edge_len)))//"bp_"
            write(702, "(a$)"), trim(adjustl(Int2Str(n_stap)))//", "
        end if
        write(702, "(a$)"), trim(adjustl(Int2Str(dna.strand(i).n_base)))//", "

        base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base
            write(702, "(a$)"), dna.top(base).seq
            base = dna.Top(base).up
        end do

        if(dna.strand(i).type1 == "scaf") write(702, "(a)")
        write(702, "(a)")
    end do
    close(unit = 702)
end subroutine Output_Write_Out_Sequences

! -----------------------------------------------------------------------------

! Write graphical routing
subroutine Output_Write_Out_Graphics(prob, geom, mesh, dna, unit)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer,        intent(in) :: unit

    ! Section type data
    type :: SecType
        integer :: start_bp, end_bp
        integer :: start_nei_edge, start_nei_sec
        integer :: end_nei_sec, end_nei_edge

        integer,   allocatable :: bp(:)
        integer,   allocatable :: node(:)
        integer,   allocatable :: strnd_scaf(:), strnd_stap(:)
        integer,   allocatable :: xover_scaf(:), xover_stap(:)
        integer,   allocatable :: nick_scaf(:),  nick_stap(:)
        character, allocatable :: seq_scaf(:),   seq_stap(:)
    end type SecType

    ! Initial line data
    type :: EdgeType
        integer :: n_sec
        type(SecType), allocatable :: sec(:)
    end type EdgeType

    ! Tecplot drawing data type
    type :: TecType
        integer :: types
        double precision :: size
        double precision :: color
        double precision :: pos_cen(3)
        double precision :: pos_node1(3)
        double precision :: pos_node2(3)
    end type TecType

    type(EdgeType), allocatable :: edge(:)
    type(TecType),  allocatable :: tec(:)
    integer,        allocatable :: conn_scaf(:,:), conn_stap(:,:)

    integer :: i, j, k, m, iniL, sec, strnd, factor, pos, type_tec(6), tec_factor
    integer :: node, dn_node, up_node, bp, max_bp, min_bp, mid_bp, mid_sec
    integer :: n_conn_scaf, n_conn_stap, n_edge, n_tec
    logical :: b_conn_scaf, b_conn_stap, b_sec
    character(200) :: path
    character(20)  :: col_list(20)

    col_list(1:4)   = ["magenta",         "orange red",   "deep pink",     "turquoise"      ]
    col_list(5:8)   = ["tan",             "salmon",       "orange",        "gold"           ]
    col_list(9:12)  = ["dark green",      "dark cyan",    "medium purple", "rosy brown"     ]
    col_list(13:16) = ["dark slate gray", "dark magenta", "sea green",     "olive drab"     ]
    col_list(17:20) = ["goldenrod",       "firebrick",    "sienna",        "dark slate blue"]

    ! Find maximum and minimum base pair ID
    max_bp = mesh.node(1).bp
    min_bp = mesh.node(1).bp
    do i = 1, mesh.n_node
        if(mesh.node(i).bp > max_bp) max_bp = mesh.node(i).bp
        if(mesh.node(i).bp < min_bp) min_bp = mesh.node(i).bp
    end do
    min_bp = min_bp - 2
    max_bp = max_bp + 2

    ! Allocate and initialize edge data
    n_edge = geom.n_iniL
    allocate(edge(n_edge))
    do i = 1, n_edge
        edge(i).n_sec = geom.n_sec
        allocate(edge(i).sec(edge(i).n_sec))

        do j = 1, edge(i).n_sec
            allocate(edge(i).sec(j).bp        (min_bp:max_bp))
            allocate(edge(i).sec(j).node      (min_bp:max_bp))
            allocate(edge(i).sec(j).strnd_scaf(min_bp:max_bp))
            allocate(edge(i).sec(j).strnd_stap(min_bp:max_bp))
            allocate(edge(i).sec(j).xover_scaf(min_bp:max_bp))
            allocate(edge(i).sec(j).xover_stap(min_bp:max_bp))
            allocate(edge(i).sec(j).nick_scaf (min_bp:max_bp))
            allocate(edge(i).sec(j).nick_stap (min_bp:max_bp))
            allocate(edge(i).sec(j).seq_scaf  (min_bp:max_bp))
            allocate(edge(i).sec(j).seq_stap  (min_bp:max_bp))

            do k = min_bp, max_bp
                edge(i).sec(j).bp(k)         = k
                edge(i).sec(j).node(k)       = -1
                edge(i).sec(j).strnd_scaf(k) = -1
                edge(i).sec(j).strnd_stap(k) = -1
                edge(i).sec(j).xover_scaf(k) = -1
                edge(i).sec(j).xover_stap(k) = -1
                edge(i).sec(j).nick_scaf(k)  = -1
                edge(i).sec(j).nick_stap(k)  = -1
                edge(i).sec(j).seq_scaf(k)   = "N"
                edge(i).sec(j).seq_stap(k)   = "N"
            end do
        end do
    end do

    ! Allocate conn data
    allocate(conn_scaf(dna.n_xover_scaf*2, 3))
    allocate(conn_stap(dna.n_xover_stap*2, 3))

    ! Build node information based on initial edges with cross-section
    do i = 1, mesh.n_node

        ! Find starting node depending on cross-section ID
        node = mesh.node(i).id
        sec  = mesh.node(i).sec
        iniL = mesh.node(node).iniL

        ! It depends on cross-section ID and direction
        if(mesh.node(i).dn == -1 .and. mod(mesh.node(i).sec, 2) == 0) then

            ! If section ID is even and negative z-direction
            ! Find neiboring edge and section
            if(dna.top(node).dn /= -1) then

                ! To avoid unpaired nucleotide
                dn_node = dna.top(node).dn
                do
                    if(dna.top(dn_node).across /= -1) exit
                    dn_node = dna.top(dn_node).dn
                end do

                edge(iniL).sec(sec+1).start_nei_edge = mesh.node(dn_node).iniL
                edge(iniL).sec(sec+1).start_nei_sec  = mesh.node(dn_node).sec
            else
                edge(iniL).sec(sec+1).start_nei_edge = -1
                edge(iniL).sec(sec+1).start_nei_sec  = -1
            end if

            edge(iniL).sec(sec+1).start_bp = mesh.node(node).bp
            do
                bp = mesh.node(node).bp
                edge(iniL).sec(sec+1).node(bp) = node
                node = mesh.node(node).up
                if(node == -1) exit
            end do
            edge(iniL).sec(sec+1).end_bp = bp

            ! Find neiboring edge and section
            if(dna.top(edge(iniL).sec(sec+1).node(bp)).up /= -1) then

                ! To avoid unpaired nucleotide
                up_node = dna.top(edge(iniL).sec(sec+1).node(bp)).up
                do
                    if(dna.top(up_node).across /= -1) exit
                    up_node = dna.top(up_node).up
                end do

                edge(iniL).sec(sec+1).end_nei_edge = mesh.node(up_node).iniL
                edge(iniL).sec(sec+1).end_nei_sec  = mesh.node(up_node).sec
            else
                edge(iniL).sec(sec+1).end_nei_edge = -1
                edge(iniL).sec(sec+1).end_nei_sec  = -1
            end if
        else if(mesh.node(i).up == -1 .and. mod(mesh.node(i).sec, 2) == 1) then

            ! If section ID is odd and positive z-direction
            ! Find neiboring edge and section
            if(dna.top(node).up /= -1) then

                ! To avoid unpaired nucleotide
                up_node = dna.top(node).up
                do
                    if(dna.top(up_node).across /= -1) exit
                    up_node = dna.top(up_node).up
                end do

                edge(iniL).sec(sec+1).start_nei_edge = mesh.node(up_node).iniL
                edge(iniL).sec(sec+1).start_nei_sec  = mesh.node(up_node).sec
            else
                edge(iniL).sec(sec+1).start_nei_edge = -1
                edge(iniL).sec(sec+1).start_nei_sec  = -1
            end if

            edge(iniL).sec(sec+1).start_bp = mesh.node(node).bp
            do
                bp = mesh.node(node).bp
                edge(iniL).sec(sec+1).node(bp) = node
                node = mesh.node(node).dn
                if(node == -1) exit
            end do
            edge(iniL).sec(sec+1).end_bp = bp

            ! Find neiboring edge and section
            if(dna.top(edge(iniL).sec(sec+1).node(bp)).dn /= -1) then

                ! To avoid unpaired nucleotide
                dn_node = dna.top(edge(iniL).sec(sec+1).node(bp)).dn
                do
                    if(dna.top(dn_node).across /= -1) exit
                    dn_node = dna.top(dn_node).dn
                end do

                edge(iniL).sec(sec+1).end_nei_edge = mesh.node(dn_node).iniL
                edge(iniL).sec(sec+1).end_nei_sec  = mesh.node(dn_node).sec
            else
                edge(iniL).sec(sec+1).end_nei_edge = -1
                edge(iniL).sec(sec+1).end_nei_sec  = -1
            end if
        end if
    end do

    ! Build additional data fields
    do i = 1, dna.n_top
        node = dna.top(i).node

        ! Skip the loop if poly Tn loop and unpaired nucleotides
        if(node == -1) cycle

        strnd = dna.top(i).strand
        sec   = mesh.node(node).sec
        iniL  = mesh.node(node).iniL
        bp    = mesh.node(node).bp

        if(dna.strand(strnd).type1 == "scaf") then

            ! For scaffold strand
            edge(iniL).sec(sec+1).seq_scaf(bp) = dna.top(i).seq
            if(dna.top(i).xover /= -1) then
                edge(iniL).sec(sec+1).xover_scaf(bp) = &
                    mesh.node(dna.top(dna.top(i).xover).node).sec
            end if
            if(dna.top(i).up == -1) then
                if(mod(sec, 2) == 0) then
                    edge(iniL).sec(sec+1).nick_scaf(bp) = 2
                else
                    edge(iniL).sec(sec+1).nick_scaf(bp) = 3
                end if
            end if
            if(dna.top(i).dn == -1) edge(iniL).sec(sec+1).nick_scaf(bp) = 1
            edge(iniL).sec(sec+1).strnd_scaf(bp) = strnd
        else if(dna.strand(strnd).type1 == "stap") then

            ! for staple strand
            edge(iniL).sec(sec+1).seq_stap(bp) = dna.top(i).seq
            if(dna.top(i).xover /= -1) then
                edge(iniL).sec(sec+1).xover_stap(bp) = &
                    mesh.node(dna.top(dna.top(i).xover).node).sec
            end if
            if(dna.top(i).up == -1) then
                if(mod(sec, 2) == 0) then
                    edge(iniL).sec(sec+1).nick_stap(bp) = 3
                else
                    edge(iniL).sec(sec+1).nick_stap(bp) = 2
                end if
            end if
            if(dna.top(i).dn == -1) edge(iniL).sec(sec+1).nick_stap(bp) = 1
            edge(iniL).sec(sec+1).strnd_stap(bp) = strnd
        end if
    end do

    ! Open file stream
    if(para_write_710 == .true.) then
        path = trim(prob.path_work)//"/"//trim(prob.name_file)
        do i = 1, n_edge
            open(unit = 710+i, file = trim(path)//"_design_edge"&
                //trim(adjustl(Int2Str(i)))//".bild", form = "formatted")
        end do

        ! Tecplot drawing
        if(para_tecplot == .true.) then

            ! 1 - normal base, 2 - xover base, 3, xover, 4 - nick, 5 - right arrow, 6- left arrow
            allocate(tec(2*dna.n_top))

            do i = 1, 2*dna.n_top
                tec(i).types          = 0
                tec(i).pos_cen(1:3)   = 0.0d0
                tec(i).pos_node1(1:3) = 0.0d0
                tec(i).pos_node2(1:3) = 0.0d0
                tec(i).size           = 0.0d0
                tec(i).color          = 0.0d0
            end do
            type_tec(1:6) = 0
            n_tec         = 0
        end if
    end if

    ! Print information based on edges
    do i = 1, n_edge

        ! Print base pair ID
        if(i == 1) then
            !do k = min_bp, max_bp
            !    write(unit, "(i5$)"), edge(i).sec(j).bp(k)
            !end do
            !write(unit, "(a)"); write(unit, "(a)")

            !do k = min_bp, max_bp
            !    write(unit, "(a5$)"), "-----"
            !end do
            !write(unit, "(a)"); write(unit, "(a)")
        end if

        ! --------------------------------------------------
        ! Scaffold strand
        ! --------------------------------------------------
        n_conn_scaf = 0
        do j = 1, edge(i).n_sec

            ! Factor : 0,0,1,1,2,2,3,3,...
            mid_sec    = int((6*edge(i).n_sec-5+1)/2)
            factor     = 3*(int(j+1)/2-1) + mid_sec
            tec_factor = (i-1)*(5+6*edge(i).n_sec-5-3*(int(edge(i).n_sec+1)/2-1))

            ! Print edge and point information, Edge 1 - (point 1 -> point 2)
            if(j == 1) then
                write(unit, "(a )"), " ==================================================="
                call Space(unit, 10)
                write(unit, "(a$)"), " [Edge "//trim(adjustl(Int2Str(i)))
                write(unit, "(a$)"), " : point "//trim(adjustl(Int2Str(geom.iniL(i).poi(1))))
                write(unit, "(a )"), " -> point "//trim(adjustl(Int2Str(geom.iniL(i).poi(2))))//"]"
                write(unit, "(a )"), " ==================================================="
            end if

            if(j == 1) then
                write(unit, "(a)")
                write(unit, "(a)"), "       [[ SCAFFOLD STRAND ]]"
                write(unit, "(a)")
            end if

            ! Print section ID and base length, sec xx -->> [xxx:xxx], xxx : - 26 length
            write(unit, "(a5$)"), " sec "
            write(unit, "(a2$)"), trim(adjustl(Int2Str(j-1)))
            write(unit, "(a4$)"), " - ["
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).start_bp + para_start_bp_ID - 1)))
            write(unit, "(a1$)"), ":"
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp + para_start_bp_ID - 1)))
            write(unit, "(a3$)"), "], "
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp - edge(i).sec(j).start_bp + 1)))
            write(unit, "(a2$)"), " :"
            call Space(unit, 4)

            ! Graphic representation, bases with crossover and nick
            do k = min_bp, max_bp

                mid_bp = int((max_bp + min_bp) / 2)

                if(edge(i).sec(j).node(k) == -1) then
                    if(edge(i).sec(j).start_bp - 2 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "5"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "3"
                    else if(edge(i).sec(j).end_bp + 1 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "3"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "5"
                    else if(edge(i).sec(j).start_bp - 1 == k .or. edge(i).sec(j).end_bp + 2 == k) then
                        write(unit, "(a$)"), "'"
                    else
                        write(unit, "(a$)"), " "
                    end if
                else
                    if(edge(i).sec(j).xover_scaf(k) /= -1) then

                        !write(unit, "(a$)"), "+"
                        if(para_output_design == "arrow") then
                            if(edge(i).sec(j).xover_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).xover_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).xover_scaf(k))))
                            end if
                        end if

                        if(para_output_design == "seq") write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)

                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        n_conn_scaf = n_conn_scaf + 1
                        conn_scaf(n_conn_scaf, 1) = edge(i).sec(j).xover_scaf(k)
                        conn_scaf(n_conn_scaf, 2) = k
                        conn_scaf(n_conn_scaf, 3) = edge(i).sec(j).strnd_scaf(k)

                        ! Draw sphere (crossover point)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_tecplot == .true.) then
                                n_tec                 = n_tec + 1
                                type_tec(3)           = type_tec(3) + 1
                                tec(n_tec).types      = 3
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(1)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_scaf(k) == 1) then
                        if(para_output_design == "arrow") write(unit, "(a$)"), "."
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw sphere (nick point)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_tecplot == .true.) then
                                n_tec                 = n_tec + 1
                                type_tec(4)           = type_tec(4) + 1
                                tec(n_tec).types      = 4
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(1)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_scaf(k) == 2) then
                        if(para_output_design == "arrow")  write(unit, "(a$)"), ">"
                        if(para_output_design == "seq")    write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw arrow (right arrow, >)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.4d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_tecplot == .true.) then
                                n_tec                   = n_tec + 1
                                type_tec(5)             = type_tec(5) + 1
                                tec(n_tec).types        = 5
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp) - 0.05d0
                                tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).pos_node1(1) = 0.9d0
                                tec(n_tec).color        = dble(1)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_scaf(k) == 3) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), "<"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw arrow (left arrow, <)
                        if(para_write_710 == .true.) then
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.4d0, -dble(6*j-5-factor), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_tecplot == .true.) then
                                n_tec                   = n_tec + 1
                                type_tec(6)             = type_tec(6) + 1
                                tec(n_tec).types        = 6
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp) + 0.05d0
                                tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                tec(n_tec).pos_node1(1) =-0.9d0
                                tec(n_tec).color        = dble(1)
                            end if
                        end if
                    else
                        if(para_output_design == "arrow") write(unit, "(a$)"), "-"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_scaf(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_scaf(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_scaf(k))
                            else if(edge(i).sec(j).strnd_scaf(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_scaf(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_scaf(k))))
                            end if
                        end if

                        ! Draw cylinder (normal base)
                        strnd = edge(i).sec(j).strnd_scaf(k)

                        if(mesh.node(edge(i).sec(j).node(k)).up == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 0) then
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                                if(para_tecplot == .true.) then
                                    n_tec                   = n_tec + 1
                                    type_tec(5)             = type_tec(5) + 1
                                    tec(n_tec).types        = 5
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                    tec(n_tec).pos_node1(1) = 1.0d0
                                    tec(n_tec).color        = dble(1)
                                end if
                            end if
                        else if(mesh.node(edge(i).sec(j).node(k)).up == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 1) then
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                                if(para_tecplot == .true.) then
                                    n_tec                   = n_tec + 1
                                    type_tec(6)             = type_tec(6) + 1
                                    tec(n_tec).types        = 6
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(6*j-5-factor) - tec_factor
                                    tec(n_tec).pos_node1(1) =-1.0d0
                                    tec(n_tec).color        = dble(1)
                                end if
                            end if
                        else
                            if(para_write_710 == .true.) then
                                write(710+i, "(a     )"), ".color steel blue"
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor), 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(1)          = type_tec(1) + 1
                                    tec(n_tec).types     = 1
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp)-0.5d0, -dble(6*j-5-factor) - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp)+0.5d0, -dble(6*j-5-factor) - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if
                    end if
                end if
            end do

            ! For Neighbor connection status, [ START : 3( 3), END : 5( 0)]
            call Space(unit, 5)
            write(unit, "(a, i2$)"), "[ START : ", edge(i).sec(j).start_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).start_nei_sec
            write(unit, "(a, i2$)"), "), END : ",  edge(i).sec(j).end_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).end_nei_sec
            write(unit, "(a$    )"), ") ]"
            write(unit, "(a     )")

            ! Draw crossovers
            call Space(unit, 30)
            do k = min_bp, max_bp
                if(edge(i).sec(j).node(k) == -1) then
                    write(unit, "(a$)"), " "
                else
                    b_conn_scaf = .false.
                    do m = 1, n_conn_scaf
                        if(conn_scaf(m, 1) >= j .and. conn_scaf(m, 2) == k) then
                            b_conn_scaf = .true.
                            exit
                        end if
                    end do

                    if(b_conn_scaf == .true.) then
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        if(mod(j, 2) == 1) then
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        else
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if
                    else
                        if( edge(i).sec(j).start_nei_edge == i .and. &
                            edge(i).sec(j).start_bp       == k .and. &
                            edge(i).sec(j).start_nei_sec  == j ) then

                        ! Neighbor connection for starting bp
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        if(mod(j, 2) == 1) then
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        else
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if

                        else if( edge(i).sec(j).end_nei_edge == i .and. &
                                 edge(i).sec(j).end_bp       == k .and. &
                                 edge(i).sec(j).end_nei_sec  == j ) then

                        ! Neighbor connection for ending bp
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        if(mod(j, 2) == 1) then
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)-3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        else
                            if(para_write_710 == .true.) then
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-factor)+3.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(1)
                                end if
                            end if
                        end if
                        else
                            write(unit, "(a$)"), " "
                        end if
                    end if
                end if
            end do
            write(unit, "(a)")
        end do

        ! --------------------------------------------------
        ! Staple strand
        ! --------------------------------------------------
        n_conn_stap = 0
        b_sec       = .true.
        do j = 1, edge(i).n_sec

            ! Factor : 0,0,1,1,2,2,3,3,...
            mid_sec    = int((6*edge(i).n_sec-5+1)/2)
            factor     = 3*(int(j + 1)/2 - 1) + mid_sec
            tec_factor = (i-1)*(5+6*edge(i).n_sec-5-3*(int(edge(i).n_sec+1)/2-1))

            if(j == 1) then
                write(unit, "(a)")
                write(unit, "(a)"), "       [[ STAPLE STRAND ]]"
                write(unit, "(a)")
            end if

            ! Print section ID and base length, sec xx -->> [xxx:xxx], xxx : - 26 length
            write(unit, "(a5$)"), " sec "
            write(unit, "(a2$)"), trim(adjustl(Int2Str(j-1)))
            write(unit, "(a4$)"), " - ["
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).start_bp + para_start_bp_ID - 1)))
            write(unit, "(a1$)"), ":"
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp + para_start_bp_ID - 1)))
            write(unit, "(a3$)"), "], "
            write(unit, "(a3$)"), trim(adjustl(Int2Str(edge(i).sec(j).end_bp - edge(i).sec(j).start_bp + 1)))
            write(unit, "(a2$)"), " :"
            call Space(unit, 4)

            ! Graphic representation, bases with crossover and nick
            do k = min_bp, max_bp
                mid_bp = int((min_bp+max_bp)/2)
                if(edge(i).sec(j).node(k) == -1) then
                    if(edge(i).sec(j).start_bp - 2 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "3"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "5"
                    else if(edge(i).sec(j).end_bp + 1 == k) then
                        if(mod(j-1, 2) == 0) write(unit, "(a$)"), "5"
                        if(mod(j-1, 2) == 1) write(unit, "(a$)"), "3"
                    else if(edge(i).sec(j).start_bp - 1 == k .or. edge(i).sec(j).end_bp + 2 == k) then
                        write(unit, "(a$)"), "'"
                    else
                        write(unit, "(a$)"), " "
                    end if
                else

                    if(mod(j, 2) == 1) then
                        pos = 6*int((j+1)/2)-5 + 1
                    else
                        pos = 6*int((j+1)/2)-5 + 2
                    end if
                    pos = 2*pos-1-factor

                    if(edge(i).sec(j).xover_stap(k) /= -1) then

                        !write(unit, "(a$)"), "+"
                        if(para_output_design == "arrow") then
                            if(edge(i).sec(j).xover_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).xover_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).xover_stap(k))))
                            end if
                        end if

                        if(para_output_design == "seq") write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)

                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        n_conn_stap = n_conn_stap + 1
                        conn_stap(n_conn_stap, 1) = edge(i).sec(j).xover_stap(k)
                        conn_stap(n_conn_stap, 2) = k
                        conn_stap(n_conn_stap, 3) = edge(i).sec(j).strnd_stap(k)

                        ! Draw sphere (crossover point)
                        strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                        if(para_write_710 == .true.) then
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(pos), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_tecplot == .true.) then
                                n_tec                 = n_tec + 1
                                type_tec(3)           = type_tec(3) + 1
                                tec(n_tec).types      = 3
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(pos) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(strnd)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_stap(k) == 1) then
                        if(para_output_design == "arrow") write(unit, "(a$)"), "."
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Draw sphere (nick point)
                        strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                        if(para_write_710 == .true.) then
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".sphere "
                            write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(pos), 0.0d0
                            write(710+i, "(1f9.2 )"), 0.25d0

                            if(para_tecplot == .true.) then
                                n_tec                 = n_tec + 1
                                type_tec(4)           = type_tec(4) + 1
                                tec(n_tec).types      = 4
                                tec(n_tec).pos_cen(1) = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2) = -dble(pos) - tec_factor
                                tec(n_tec).size       = 0.25d0
                                tec(n_tec).color      = dble(strnd)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_stap(k) == 2) then
                        if(para_output_design == "arrow") write(unit, "(a$)"), ">"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Draw arrow (right arrow, >)
                        strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                        if(para_write_710 == .true.) then
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_tecplot == .true.) then
                                n_tec                   = n_tec + 1
                                type_tec(5)             = type_tec(5) + 1
                                tec(n_tec).types        = 5
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                tec(n_tec).pos_node1(1) = 0.8d0
                                tec(n_tec).color        = dble(strnd)
                            end if
                        end if
                    else if(edge(i).sec(j).nick_stap(k) == 3) then

                        if(para_output_design == "arrow") write(unit, "(a$)"), "<"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Draw arrow (left arrow, <)
                        strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                        if(para_write_710 == .true.) then
                            write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                            write(710+i, "(a$    )"), ".arrow "
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.4d0, -dble(pos), 0.0d0
                            write(710+i, "(3f8.2 )"), 0.1d0, 0.4d0, 0.55d0

                            if(para_tecplot == .true.) then
                                n_tec                   = n_tec + 1
                                type_tec(6)             = type_tec(6) + 1
                                tec(n_tec).types        = 6
                                tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                tec(n_tec).pos_node1(1) =-0.8d0
                                tec(n_tec).color        = dble(strnd)
                            end if
                        end if
                    else

                        if(para_output_design == "arrow") write(unit, "(a$)"), "-"
                        if(para_output_design == "seq")   write(unit, "(a$)"), edge(i).sec(j).seq_stap(k)
                        if(para_output_design == "strand") then
                            if(edge(i).sec(j).strnd_stap(k) >= 10) then
                                write(unit, "(a$)"), achar(55+edge(i).sec(j).strnd_stap(k))
                            else if(edge(i).sec(j).strnd_stap(k) >= 36) then
                                write(unit, "(a$)"), achar(61+edge(i).sec(j).strnd_stap(k))
                            else
                                write(unit, "(a$)"), trim(adjustl(Int2Str(edge(i).sec(j).strnd_stap(k))))
                            end if
                        end if

                        ! Write Chimera
                        strnd = mod(edge(i).sec(j).strnd_stap(k), 20) + 1
                        if(mesh.node(edge(i).sec(j).node(k)).dn == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 1) then

                            ! Draw arrow (right arrow, >)
                            if(para_write_710 == .true.) then
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.09d0, 0.4d0, 0.55d0

                                if(para_tecplot == .true.) then
                                    n_tec                   = n_tec + 1
                                    type_tec(5)             = type_tec(5) + 1
                                    tec(n_tec).types        = 5
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                    tec(n_tec).pos_node1(1) = 1.0d0
                                    tec(n_tec).color        = dble(strnd)
                                end if
                            end if
                        else if(mesh.node(edge(i).sec(j).node(k)).dn == -1 .and. mod(mesh.node(edge(i).sec(j).node(k)).sec, 2) == 0) then

                            ! Draw arrow (left arrow, <)
                            if(para_write_710 == .true.) then
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".arrow "
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)+0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2$)"), dble(k-mid_bp)-0.5d0, -dble(pos), 0.0d0
                                write(710+i, "(3f8.2 )"), 0.09d0, 0.4d0, 0.55d0

                                if(para_tecplot == .true.) then
                                    n_tec                   = n_tec + 1
                                    type_tec(6)             = type_tec(6) + 1
                                    tec(n_tec).types        = 6
                                    tec(n_tec).pos_cen(1)   = dble(k-mid_bp)
                                    tec(n_tec).pos_cen(2)   = -dble(pos) - tec_factor
                                    tec(n_tec).pos_node1(1) =-1.0d0
                                    tec(n_tec).color        = dble(strnd)
                                end if
                            end if
                        else

                            ! Draw cylinder, standard base
                            if(para_write_710 == .true.) then
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)-0.4d0, -dble(pos), 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp)+0.4d0, -dble(pos), 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(1)          = type_tec(1) + 1
                                    tec(n_tec).types     = 1
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp)-0.4d0, -dble(pos) - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp)+0.4d0, -dble(pos) - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(strnd)
                                end if
                            end if
                        end if
                    end if

                    if(para_write_710 == .true.) then

                        ! Draw sequence
                        if(mod(j, 2) == 0) then
                            write(710+i, "(a, 3f9.3)"), ".cmov ", dble(k-mid_bp-0.3d0), -dble(pos)-1.0d0, 0.0d0
                        else
                            write(710+i, "(a, 3f9.3)"), ".cmov ", dble(k-mid_bp-0.3d0), -dble(pos)+0.4d0, 0.0d0
                        end if
                        write(710+i, "(a)"), ".font arial 20 bold"
                        write(710+i, "(a)"), dna.top(dna.top(edge(i).sec(j).node(k)).across).seq

                        ! Draw section ID
                        if(b_sec == .true.) then
                            write(710+i, "(a)"), ".color red"
                            if(mod(j, 2) == 0) then
                                write(710+i, "(a, 3f9.3)"), ".cmov ", dble(min_bp-mid_bp-5.0d0), -dble(pos)-1.0d0, 0.0d0
                            else
                                write(710+i, "(a, 3f9.3)"), ".cmov ", dble(min_bp-mid_bp-5.0d0), -dble(pos)+0.4d0, 0.0d0
                            end if
                            write(710+i, "(a)"), ".font arial 20 bold"
                            write(710+i, "(a)"), "sec "//trim(adjustl(Int2Str(j-1)))
                            b_sec = .false.
                        end if
                    end if
                end if
            end do
            b_sec = .true.

            ! For Neighbor connection status, [ START : 3( 3), END : 5( 0)]
            call Space(unit, 5)
            do k = 1, edge(i).n_sec
                if(edge(i).sec(k).start_nei_edge == i) then
                    edge(i).sec(k).start_nei_edge = -1
                    edge(i).sec(k).start_nei_sec  = -1
                end if
                if(edge(i).sec(k).end_nei_edge == i) then
                    edge(i).sec(k).end_nei_edge = -1
                    edge(i).sec(k).end_nei_sec  = -1
                end if
            end do
            write(unit, "(a, i2$)"), "[ START : ", edge(i).sec(j).start_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).start_nei_sec
            write(unit, "(a, i2$)"), "), END : ",  edge(i).sec(j).end_nei_edge
            write(unit, "(a, i2$)"), "(",          edge(i).sec(j).end_nei_sec
            write(unit, "(a$    )"), ") ]"
            write(unit, "(a     )")

            ! Draw crossovers
            call Space(unit, 30)
            do k = min_bp, max_bp
                if(edge(i).sec(j).node(k) == -1) then
                    write(unit, "(a$)"), " "
                else
                    b_conn_stap = .false.
                    do m = 1, n_conn_stap
                        if(conn_stap(m, 1) >= j .and. conn_stap(m, 2) == k) then
                            b_conn_stap = .true.
                            exit
                        end if
                    end do

                    if(b_conn_stap == .true.) then
                        write(unit, "(a$)"), "|"

                        ! Draw cylinder
                        strnd = mod(conn_stap(m, 3), 20) + 1
                        if(mod(j, 2) == 0) then
                            if(para_write_710 == .true.) then
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-5+5-factor)+0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-5-factor)-0.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-5+5-factor)+0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-5-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(strnd)
                                end if
                            end if
                        else
                            if(para_write_710 == .true.) then
                                write(710+i, "(a     )"), ".color "//trim(col_list(strnd))
                                write(710+i, "(a$    )"), ".cylinder "
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2+1-factor)+0.2d0, 0.0d0
                                write(710+i, "(3f9.2$)"), dble(k-mid_bp), -dble(6*j-2-1-factor)-0.2d0, 0.0d0
                                write(710+i, "(1f9.2 )"), 0.1d0

                                if(para_tecplot == .true.) then
                                    n_tec                = n_tec + 1
                                    type_tec(2)          = type_tec(2) + 1
                                    tec(n_tec).types     = 2
                                    tec(n_tec).pos_node1 = [dble(k-mid_bp), -dble(6*j-2+1-factor)+0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).pos_node2 = [dble(k-mid_bp), -dble(6*j-2-1-factor)-0.2d0 - tec_factor, 0.0d0]
                                    tec(n_tec).color     = dble(strnd)
                                end if
                            end if
                        end if
                    else
                        write(unit, "(a$)"), " "
                    end if
                end if
            end do
            write(unit, "(a)")
        end do
        write(unit, "(a)")
    end do

    ! Write tecplot output
    if(para_write_710 == .true. .and. para_tecplot == .true.) then

        ! Open file stream
        path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
        open(unit = 709, file = trim(path)//"_design_edge.dat", form = "formatted")
        do i = 1, 6

            if(type_tec(i) == 0) cycle

            write(709, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
            write(709, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3", "c"'
            write(709, "(a$)"), 'ZONE F = FEPOINT'

            if(i == 1 .or. i == 2) then

                write(709, "(a$)"), ', N='//trim(adjustl(Int2Str(2*type_tec(i))))
                write(709, "(a$)"), ', E='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a )"), ', ET=LINESEG'

                ! Set nodal connectivity
                do j = 1, n_tec
                    if(tec(j).types == i) then
                        write(709, "(7f8.2)"), tec(j).pos_node1(1:3), 1.0d0, 1.0d0, 1.0d0, tec(j).color
                        write(709, "(7f8.2)"), tec(j).pos_node2(1:3), 1.0d0, 1.0d0, 1.0d0, tec(j).color
                    end if
                end do

                ! Set connectivity
                do j = 1, type_tec(i)
                    write(709, "(2i8)"), j*2-1, j*2
                end do
            else if(type_tec(i) > 0 .and. (i == 3 .or. i == 4)) then

                write(709, "(a$)"), ', N='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a$)"), ', E='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a )"), ', ET=LINESEG'

                ! Set nodal position
                do j = 1, n_tec
                    if(tec(j).types == i) then
                        write(709, "(7f8.2)"), tec(j).pos_cen(1:3), tec(j).size, 1.0d0, 1.0d0, tec(j).color
                    end if
                end do

                ! Set connectivity
                do j = 1, type_tec(i)
                    write(709, "(2i8)"), j, j
                end do
            else if(type_tec(i) > 0 .and. (i == 5 .or. i == 6)) then
                write(709, "(a$)"), ', N='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a$)"), ', E='//trim(adjustl(Int2Str(type_tec(i))))
                write(709, "(a )"), ', ET=LINESEG'

                ! Set nodal position
                do j = 1, n_tec
                    if(tec(j).types == i) then
                        write(709, "(7f8.2)"), tec(j).pos_cen(1:3), tec(j).pos_node1(1:3), tec(j).color
                    end if
                end do

                ! Set connectivity
                do j = 1, type_tec(i)
                    write(709, "(2i8)"), j, j
                end do
            end if
        end do
    end if

    ! Deallocate memory
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            deallocate(edge(i).sec(j).bp        )
            deallocate(edge(i).sec(j).node      )
            deallocate(edge(i).sec(j).strnd_scaf)
            deallocate(edge(i).sec(j).strnd_stap)
            deallocate(edge(i).sec(j).xover_scaf)
            deallocate(edge(i).sec(j).xover_stap)
            deallocate(edge(i).sec(j).nick_scaf )
            deallocate(edge(i).sec(j).nick_stap )
            deallocate(edge(i).sec(j).seq_scaf  )
            deallocate(edge(i).sec(j).seq_stap  )
        end do
        deallocate(edge(i).sec)
    end do
    deallocate(edge)
    deallocate(conn_scaf)
    deallocate(conn_stap)

    do i = 1, n_edge
        close(unit = 710+i)
    end do

    ! Deallocate memory and close file for Tecplot output
    if(para_write_710 == .true. .and. para_tecplot == .true.) then
        deallocate(tec)
        close(unit = 709)
    end if
end subroutine Output_Write_Out_Graphics

! -----------------------------------------------------------------------------

! Write output of unpaired nucleotides
function Output_Write_Out_Unpaired(mesh, dna, unit) result(max_base)
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    integer,        intent(in)    :: unit

    double precision :: length
    integer :: s_base, e_base, s_iniL, e_iniL, s_sec, e_sec, max_base
    integer :: i, j, n_base, base, count
    character(100) :: seq
    character(4) :: types

    dna.n_nt_unpaired_scaf = 0
    dna.n_nt_unpaired_stap = 0
    dna.n_unpaired_scaf    = 0
    dna.n_unpaired_stap    = 0

    max_base = 0
    n_base   = 0
    do i = 1, dna.n_top

        ! Find the end bases in basepairs
        if(dna.top(i).up /= -1 .and. dna.top(i).dn /= -1) then
            !
            !                i
            ! *--*--*--*--*->*--*--*--*-->
            ! |  |  |  |  |  |
            ! *<-*--*--*--*--*
            if( dna.top(i).across /= -1 .and. &
                dna.top(dna.top(i).up).node == -1 .and. &
                dna.top(dna.top(i).dn).node /= -1 ) then

                base   = dna.top(i).id
                types  = dna.strand(dna.top(i).strand).type1
                s_iniL = mesh.node(dna.top(dna.top(base).dn).node).iniL
                s_sec  = mesh.node(dna.top(dna.top(base).dn).node).sec
                count  = 0
                n_base = n_base + 1
                s_base = base

                if(types == "scaf") dna.n_unpaired_scaf = dna.n_unpaired_scaf + 1
                if(types == "stap") dna.n_unpaired_stap = dna.n_unpaired_stap + 1

                ! Count the number of bases
                do
                    base = dna.top(base).up
                    if(dna.top(base).across /= -1) exit
                    count = count + 1
                    seq(count:count) = dna.top(base).seq

                    if(types == "scaf") dna.n_nt_unpaired_scaf = dna.n_nt_unpaired_scaf + 1
                    if(types == "stap") dna.n_nt_unpaired_stap = dna.n_nt_unpaired_stap + 1
                end do

                e_iniL = mesh.node(dna.top(base).node).iniL
                e_sec  = mesh.node(dna.top(base).node).sec
                e_base = dna.top(base).id
                length = Norm(dna.top(s_base).pos - dna.top(e_base).pos)

                write(unit, "(i10, a$ )"), n_base, " "//trim(types)
                write(unit, "(a$      )"), ", # of unpaired nts : "//trim(adjustl(Int2Str(count)))
                write(unit, "(a, f5.2$)"), " <-- Total length : ", length
                write(unit, "(a$      )"), ", Edge(sec)) : "
                write(unit, "(a$      )"), trim(adjustl(Int2Str(s_iniL)))//"("//trim(adjustl(Int2Str(s_sec)))//") -> "
                write(unit, "(a$      )"), trim(adjustl(Int2Str(e_iniL)))//"("//trim(adjustl(Int2Str(e_sec)))//"), Sequence : "

                ! Check maximum number of unpaired nucleotides
                if(max_base < count) max_base = count

                do j = 1, count
                    write(unit, "(a$)"), seq(j:j)
                end do
                write(unit, "(a)")
            end if
        end if
    end do
    write(unit, "(a)")
end function Output_Write_Out_Unpaired

! -----------------------------------------------------------------------------

! Outputs based on strands and nucleotides
subroutine Output_Write_Out_Strand_Base(mesh, dna, unit)
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna
    integer, intent(in) :: unit

    integer :: i, j

    write(unit, "(a)"), " Output based on strands"
    write(unit, "(a)"), " The total number of strands : "//trim(adjustl(Int2Str(dna.n_strand)))
    write(unit, "(a)")

    ! Strand based output
    do i = 1, dna.n_strand
        write(unit, "(i10,a$)"), i, " "//trim(dna.strand(i).type1)//" ==>"
        write(unit, "(a, i6$)"), " # of nts : ", dna.strand(i).n_base
        write(unit, "(a, l  )"), ", circular : ", dna.strand(i).b_circular

        ! Print bases numbering
        call space(unit, 15)
        write(unit, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(1))))//"->"
        do j = 2, dna.strand(i).n_base - 1
            if(mod(j, 10) == 0) then
                write(unit, "(a9)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//"->"
            else if(mod(j, 10) == 1) then
                call space(unit, 15)
                write(unit, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//"->"
            else
                write(unit, "(a9$)"), trim(adjustl(Int2Str(dna.strand(i).base(j))))//"->"
            end if
        end do
        if(mod(j, 10) == 1) call space(unit, 15)
        write(unit, "(a7)"), trim(adjustl(Int2Str(dna.strand(i).base(dna.strand(i).n_base))))
        write(unit, "(a )")
    end do

    ! DNA base information
    write(unit, "(a)"), " Output based on nucleotides"
    write(unit, "(a)"), " The total number of nucleotides : "//trim(adjustl(Int2Str(dna.n_top)))
    write(unit, "(a)")

    ! node-bp-InitL-sec-up-dn-ac-xo
    do i = 1, dna.n_top
        write(unit, "(i10, a$)"), dna.top(i).id, " nt ==>"
        write(unit, "(a$     )"), trim(dna.strand(dna.top(i).strand).type1)
        write(unit, "(a$     )"), ", seq: "//trim(dna.top(i).seq)
        write(unit, "(a, i6$ )"), ", nde: ", dna.top(i).node

        if(dna.top(i).node /= -1) then
            write(unit, "(a, i6$)"), ", bp: ",   mesh.node(dna.top(i).node).bp
            write(unit, "(a, i3$)"), ", iniL: ", mesh.node(dna.top(i).node).iniL
            write(unit, "(a, i3$)"), ", sec: ",  mesh.node(dna.top(i).node).sec
        else
            write(unit, "(a$)"), ",      UNPAIRED NUCLEOTIDES      "
        end if

        write(unit, "(a, i6$)"), ", up: ",     dna.top(i).up
        write(unit, "(a, i6$)"), ", dn: ",     dna.top(i).dn
        write(unit, "(a, i6$)"), ", xover: ",  dna.top(i).xover
        write(unit, "(a, i6$)"), ", across: ", dna.top(i).across
        write(unit, "(a, i3$)"), ", strand: ", dna.top(i).strand
        write(unit, "(a)")
    end do
    write(unit, "(a)")
end subroutine Output_Write_Out_Strand_Base

! -----------------------------------------------------------------------------

! Output about staple length
subroutine Output_Write_Out_Staple_Length(dna, unit)
    type(DNAType), intent(in) :: dna
    integer,       intent(in) :: unit

    integer, allocatable :: length_stap(:)
    integer :: i

    ! Allocate and initialize memory
    allocate(length_stap(dna.len_min_stap:dna.len_max_stap))
    length_stap(dna.len_min_stap:dna.len_max_stap) = 0

    ! Find staple length
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "stap") then
            length_stap(dna.strand(i).n_base) = length_stap(dna.strand(i).n_base) + 1
        end if
    end do

    ! Write staple length
    write(unit, "(a)"), "   Edge length     Number"
    write(unit, "(a)"), "   ----------------------"
    do i = dna.len_min_stap, dna.len_max_stap
        if(length_stap(i) /= 0) then
            write(unit, "(2i11)"), i, length_stap(i)
        end if
    end do

    ! Deallocate memory
    deallocate(length_stap)
end subroutine Output_Write_Out_Staple_Length

! -----------------------------------------------------------------------------

! Write JSON guide model
subroutine Output_Write_Out_Guide_JSON(prob, geom, mesh)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom
    type(MeshType), intent(in) :: mesh

    double precision :: length, pos_1(3), pos_2(3), pos_c(3), t1(3)
    integer :: i, j, iter, nbp, add, n_conn
    logical :: f_axis, f_info
    character(200) :: path

    ! Set option
    f_axis = para_chimera_axis

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 998, file = trim(path)//"_14_guide.bild", form = "formatted")

    ! Write points for multi lines
    write(998, "(a, 3f6.2)"), ".color ", 0.0d0/255.0d0, 114.0d0/255.0d0, 178.0d0/255.0d0
    do i = 1, geom.n_croP

        write(998, "(a$    )"), ".sphere "
        write(998, "(3f9.3$)"), geom.croP(i).pos(1:3)
        write(998, "(1f9.3 )"), 0.2d0
    end do

    ! Write multi lines
    write(998, "(a, 3f6.2)"), ".color ", 0.0d0/255.0d0, 114.0d0/255.0d0, 178.0d0/255.0d0
    do i = 1, geom.n_croL

        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)

        write(998, "(a$   )"), ".cylinder "
        write(998, "(7f9.3)"), pos_1(1:3), pos_2(1:3), 0.1d0
    end do

    ! Write vertex connection
    n_conn = 0
    write(998, "(a, 3f6.2)"), ".color ", 213.0d0/255.0d0, 94.0d0/255.0d0, 0.0d0/255.0d0
    do i = 1, geom.n_junc
        do j = 1, geom.n_sec*geom.junc(i).n_arm

            pos_1(1:3) = mesh.node(geom.junc(i).conn(j,1)).pos(1:3)
            pos_2(1:3) = mesh.node(geom.junc(i).conn(j,2)).pos(1:3)

            n_conn = n_conn + 1

            if(Is_Same_Vector(pos_1, pos_2) == .false.) then
                write(998, "(a$   )"), ".cylinder "
                write(998, "(7f9.3)"), pos_1(1:3), pos_2(1:3), 0.04d0
            end if
        end do
    end do

    ! Write scaffold direction
    write(998, "(a)"), ".color red"
    do i = 1, geom.n_croL
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        t1(1:3) = geom.croL(i).t(1, 1:3)

        write(998, "(a$    )"), ".arrow "
        write(998, "(3f8.2$)"), pos_c(1:3)
        write(998, "(3f8.2$)"), pos_c(1:3) + t1(1:3) * 1.6d0
        write(998, "(3f8.2 )"), 0.16d0, 0.45d0, 0.5d0
    end do

    ! Write edge number
    add = -1
    write(998, "(a)"), ".color dark green"
    do i = 1, geom.n_croL
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        if(mod(i, geom.n_croL / geom.n_sec) == 1) add = add + 1

        !write(998, "(a$   )"), ".cmov "
        !write(998, "(3f9.3)"), pos_c(1:2) + 0.6d0, pos_c(3) - 0.6d0
        !write(998, "(a    )"), ".font Helvetica 12 bold"
        !write(998, "(i7)"), mod(2*i-1, geom.n_croL) + add - 1

        write(998, "(a$   )"), ".cmov "
        write(998, "(3f9.3)"), pos_c(1:3) + 0.2d0
        write(998, "(a    )"), ".font Helvetica 12 bold"
        write(998, "(i7   )"), mod(geom.n_sec*i-(geom.n_sec-1), geom.n_croL) + add - 1
    end do

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(998)
    close(unit = 998)

    ! --------------------------------------------------
    ! Write the file for Tecplot
    ! --------------------------------------------------
    if(para_tecplot == .false.) return

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit = 998, file = trim(path)//"_14_guide.dat", form = "formatted")

    write(998, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(998, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(998, "(a$)"), 'ZONE F = FEPOINT'
    write(998, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_croP)))
    write(998, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_croL)))
    write(998, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_croP
        write(998, "(3f9.3$)"), geom.croP(i).pos(1:3)
        write(998, "(3f9.3 )"), 1.0d0, 1.0d0, 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_croL
        write(998, "(2i7)"), geom.croL(i).poi(1), geom.croL(i).poi(2)
    end do

    write(998, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(998, "(a$)"), 'ZONE F = FEPOINT'
    write(998, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_croL)))
    write(998, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_croL)))
    write(998, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_croL
        pos_1(1:3) = geom.croP(geom.croL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.croP(geom.croL(i).poi(2)).pos(1:3)
        pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0

        write(998, "(3f9.3$)"), pos_c(1:3)
        write(998, "(3f9.3 )"), geom.croL(i).t(1, 1:3)
    end do

    ! Write edges
    do i = 1, geom.n_croL
        write(998, "(2i7)"), i, i
    end do

    write(998, "(a )"), 'VARIABLES = "X", "Y", "Z", "t1", "t2", "t3"'
    write(998, "(a$)"), 'ZONE F = FEPOINT'
    write(998, "(a$)"), ', N='//trim(adjustl(Int2Str(n_conn*2)))
    write(998, "(a$)"), ', E='//trim(adjustl(Int2Str(n_conn)))
    write(998, "(a )"), ', ET=LINESEG'

    ! Write vertex connection
    do i = 1, geom.n_junc
        do j = 1, geom.n_sec*geom.junc(i).n_arm

            pos_1(1:3) = mesh.node(geom.junc(i).conn(j,1)).pos(1:3)
            pos_2(1:3) = mesh.node(geom.junc(i).conn(j,2)).pos(1:3)

            write(998, "(6f9.3)"), pos_1(1:3), 0.0d0, 0.0d0, 0.0d0
            write(998, "(6f9.3)"), pos_2(1:3), 0.0d0, 0.0d0, 0.0d0
        end do
    end do

    ! Write edges
    do i = 1, n_conn
        write(998, "(2i7)"), 2*i-1, 2*i
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write JSON guide"

    close(unit = 998)
end subroutine Output_Write_Out_Guide_JSON

! -----------------------------------------------------------------------------

! Write JSON output
subroutine Output_Write_Out_JSON(prob, geom, mesh, dna, max_unpaired)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(in)    :: geom
    type(MeshType), intent(in)    :: mesh
    type(DNAType),  intent(inout) :: dna
    integer,        intent(in)    :: max_unpaired

    ! Conn type data
    type :: ConnType
        integer :: conn(4)  ! [ down sec | down bp | up sec | up bp ]
    end type ConnType

    ! Section type data - Global section => (E-1)*iniE + sec
    type :: SecType
        type(ConnType), allocatable :: scaf(:)
        type(ConnType), allocatable :: stap(:)

        integer :: n_stap_col
        integer :: stap_col(40,2)
        integer :: n_xover_stap
        integer :: n_xover_scaf
    end type SecType

    ! Initial line data
    type :: EdgeType
        integer :: n_sec
        type(SecType), allocatable :: sec(:)
    end type EdgeType

    type(EdgeType), allocatable :: edge(:)
    integer :: i, j, k, min_bp, max_bp, n_edge, width, num, shift
    integer :: c_base, c_node, c_edge, c_sec, c_bp, pos_c, pos_r
    integer :: dn_base, dn_node, dn_edge, dn_sec, dn_bp
    integer :: up_base, up_node, up_edge, up_sec, up_bp
    integer :: cup_bp, cdn_bp, col, row, col_shift, row_shift, color(12)
    integer :: n_xover, min_xover, min_xover_edge, min_xover_sec
    logical :: b_stap_json = .true.
    character(200) :: path

    if(dna.n_base_scaf > 20000) return

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 999, file = trim(path)//".json", form = "formatted")

    ! Hex color code
    color( 1) = 13369344    ! #cc0000
    color( 2) = 16204552    ! #f74308
    color( 3) = 16225054    ! #f7931e
    color( 4) = 11184640    ! #aaaa00
    color( 5) = 5749504     ! #57bb00
    color( 6) = 29184       ! #007200
    color( 7) = 243362      ! #03b6a2
    color( 8) = 1507550     ! #1700de
    color( 9) = 7536862     ! #7300de
    color(10) = 12060012    ! #b8056c
    color(11) = 3355443     ! #333333
    color(12) = 8947848     ! #888888

    ! Find maximum and minimum bp ID
    max_bp = mesh.node(1).bp
    min_bp = mesh.node(1).bp
    do i = 2, mesh.n_node
        if(mesh.node(i).bp > max_bp) max_bp = mesh.node(i).bp
        if(mesh.node(i).bp < min_bp) min_bp = mesh.node(i).bp
    end do

    min_bp = min_bp + para_start_bp_ID - 1
    max_bp = max_bp + para_start_bp_ID - 1
    shift  = para_start_bp_ID + 21 + 21 ! + 21

    ! Possible maximum edge length = max_bp - min_bp + 1 + max_unpaired
    width = (max_bp - min_bp + 1 + max_unpaired) + 2 * para_start_bp_ID + 21
    width = (width / 21) * 21 + 21 + 21

    ! Allocate and initialize edge data
    n_edge = geom.n_iniL
    allocate(edge(n_edge))

    do i = 1, n_edge
        edge(i).n_sec = geom.n_sec
        allocate(edge(i).sec(edge(i).n_sec))

        do j = 1, edge(i).n_sec
            allocate(edge(i).sec(j).scaf(width))
            allocate(edge(i).sec(j).stap(width))

            edge(i).sec(j).n_stap_col   = 0
            edge(i).sec(j).stap_col     = 0
            edge(i).sec(j).n_xover_stap = 0
            edge(i).sec(j).n_xover_scaf = 0

            do k = 1, width
                edge(i).sec(j).scaf(k).conn(1:4) = -1
                edge(i).sec(j).stap(k).conn(1:4) = -1
            end do
        end do
    end do

    ! Strand based loop
    do i = 1, dna.n_strand

        c_base = Mani_Go_Start_Base(dna, i)
        do j = 1, dna.strand(i).n_base

            ! Current base
            c_node = dna.top(c_base).node
            if(c_node /= -1) then
                c_edge = mesh.node(c_node).iniL
                c_sec  = mesh.node(c_node).sec
                c_bp   = mesh.node(c_node).bp + shift
            else
                if(mod(c_sec, 2) == 0) then
                    if(dna.strand(i).type1 == "scaf") c_bp = c_bp + 1
                    if(dna.strand(i).type1 == "stap") c_bp = c_bp - 1
                else
                    if(dna.strand(i).type1 == "scaf") c_bp = c_bp - 1
                    if(dna.strand(i).type1 == "stap") c_bp = c_bp + 1
                end if
            end if

            ! Downward base
            dn_base = dna.top(c_base).dn
            if(dn_base /= -1) then
                dn_node = dna.top(dn_base).node
                if(dn_node /= -1) then
                    dn_edge = mesh.node(dn_node).iniL
                    dn_sec  = mesh.node(dn_node).sec
                    dn_bp   = mesh.node(dn_node).bp + shift
                else
                    if(mod(dn_sec, 2) == 0) then
                        if(dna.strand(i).type1 == "scaf") dn_bp = dn_bp + 1
                        if(dna.strand(i).type1 == "stap") dn_bp = dn_bp - 1
                    else
                        if(dna.strand(i).type1 == "scaf") dn_bp = dn_bp - 1
                        if(dna.strand(i).type1 == "stap") dn_bp = dn_bp + 1
                    end if
                end if
                if(dna.strand(i).type1 == "scaf") then
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(1) = (dn_edge - 1) * edge(1).n_sec + dn_sec
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(2) = dn_bp - 1
                else
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(1) = (dn_edge - 1) * edge(1).n_sec + dn_sec
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(2) = dn_bp - 1
                end if
            else if(dn_base == -1 .and. dna.strand(i).type1 == "stap") then

                ! Set color
                edge(c_edge).sec(c_sec+1).n_stap_col = edge(c_edge).sec(c_sec+1).n_stap_col + 1
                num = edge(c_edge).sec(c_sec+1).n_stap_col
                if(num >= 41) then
                    write(p_redir, "(a)")
                    write(p_redir, "(a)"), " +=== error =========================================+"
                    write(p_redir, "(a)"), " | Check # of Xovers per edge                        |"
                    write(p_redir, "(a)"), " +===================================================+"
                    stop
                end if
                edge(c_edge).sec(c_sec+1).stap_col(num,1) = c_bp - 1
                edge(c_edge).sec(c_sec+1).stap_col(num,2) = color(mod(i, 12)+1)
            end if

            ! Upper base
            up_base = dna.top(c_base).up
            if(up_base /= -1) then
                up_node = dna.top(up_base).node
                if(up_node /= -1) then
                    up_edge = mesh.node(up_node).iniL
                    up_sec  = mesh.node(up_node).sec
                    up_bp   = mesh.node(up_node).bp + shift
                else
                    if(mod(c_sec, 2) == 0) then
                        if(dna.strand(i).type1 == "scaf") up_bp = up_bp + 1
                        if(dna.strand(i).type1 == "stap") up_bp = up_bp - 1
                    else
                        if(dna.strand(i).type1 == "scaf") up_bp = up_bp - 1
                        if(dna.strand(i).type1 == "stap") up_bp = up_bp + 1
                    end if
                end if
                if(dna.strand(i).type1 == "scaf") then
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(3) = (up_edge - 1) * edge(1).n_sec + up_sec
                    edge(c_edge).sec(c_sec+1).scaf(c_bp).conn(4) = up_bp - 1
                else
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(3) = (up_edge - 1) * edge(1).n_sec + up_sec
                    edge(c_edge).sec(c_sec+1).stap(c_bp).conn(4) = up_bp - 1
                end if
            end if

            ! Count xover of scaffold and staple
            if(dna.top(c_base).xover /= -1) then
                if(dna.strand(i).type1 == "scaf") then
                    edge(c_edge).sec(c_sec+1).n_xover_scaf = edge(c_edge).sec(c_sec+1).n_xover_scaf + 1
                end if
                if(dna.strand(i).type1 == "stap") then
                    edge(c_edge).sec(c_sec+1).n_xover_stap = edge(c_edge).sec(c_sec+1).n_xover_stap + 1
                end if
            end if

            ! Update base
            c_base = dna.Top(c_base).up
        end do
    end do

    ! Print JSON-style data structure
    write(999, "(a$)"), "{"
    write(999, "(a$)"), '"name":"'//trim(adjustl(prob.name_prob))//'",'
    write(999, "(a$)"), '"vstrands":'
    write(999, "(a$)"), '['

    row_shift = 0
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            write(999, "(a$)"), '{'

            ! 1. Staple colors
            write(999, "(a$)"), '"stap_colors":['
            do k = 1, edge(i).sec(j).n_stap_col
                write(999, "(a$)"), '['//trim(adjustl(Int2Str(edge(i).sec(j).stap_col(k,1))))
                write(999, "(a$)"), ','//trim(adjustl(Int2Str(edge(i).sec(j).stap_col(k,2))))
                
                if(k /= edge(i).sec(j).n_stap_col) write(999, "(a$)"), '],'
                if(k == edge(i).sec(j).n_stap_col) write(999, "(a$)"), ']'
            end do
            write(999, "(a$)"), '],'

            ! Cross-section position
            col_shift = mod(i, 7)
            if(col_shift == 1 .and. j == 1) row_shift = row_shift + 4
            if(col_shift == 0) col_shift = 7
            col_shift = 4*col_shift

            num = (i - 1) * edge(i).n_sec + j

            if(j == 1) then
                pos_c = 1; pos_r = 1
            else if(j == 2) then
                pos_c = 1; pos_r = 2
            else if(j == 3) then
                pos_c = 2; pos_r = 2
            else if(j == 4) then
                pos_c = 3; pos_r = 2
            else if(j == 5) then
                pos_c = 3; pos_r = 1
            else if(j == 6) then
                pos_c = 2; pos_r = 1
            end if

            !pos_c = geom.sec.posR(j)
            !pos_r = geom.sec.posC(j)

            col = col_shift - pos_c
            row = row_shift - pos_r

            ! 2. Num
            write(999, "(a$)"), '"num":'//trim(adjustl(Int2Str(num - 1)))//","

            ! 3. Scaffold loop
            write(999, "(a$)"), '"scafLoop":[],'

            ! 4. Staple
            write(999, "(a$)"), '"stap":['
            do k = 1, width
                if(b_stap_json == .true.) then
                    write(999, "(a$)"), "["//&
                        trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(1))))//","//&
                        trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(2))))//","//&
                        trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(3))))//","//&
                        trim(adjustl(Int2Str(edge(i).sec(j).stap(k).conn(4))))
                else
                    write(999, "(a$)"), "[-1, -1, -1, -1"
                end if

                if(k /= width) write(999, "(a$)"), "],"
                if(k == width) write(999, "(a$)"), "]],"
            end do

            ! 5. Skip
            write(999, "(a$)"), '"skip":['
            do k = 1, width - 1
                write(999, "(a$)"), "0,"
            end do
            write(999, "(a$)"), "0],"

            ! 6. Scaffold
            write(999, "(a$)"), '"scaf":['
            do k = 1, width
                write(999, "(a$)"), "["//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(1))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(2))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(3))))//","//&
                    trim(adjustl(Int2Str(edge(i).sec(j).scaf(k).conn(4))))

                if(k /= width) write(999, "(a$)"), "],"
                if(k == width) write(999, "(a$)"), "]],"
            end do

            ! 7. Staple loop
            write(999, "(a$)"), '"stapLoop":[],'

            ! 8. Col
            write(999, "(a$)"), '"col":'//trim(adjustl(Int2Str(col)))//","

            ! 9. Loop
            write(999, "(a$)"), '"loop":['
            do k = 1, width - 1
                write(999, "(a$)"), "0,"
            end do
            write(999, "(a$)"), "0],"

            ! 10. Row
            write(999, "(a$)"), '"row":'//trim(adjustl(Int2Str(row)))

            if(i == n_edge .and. j == edge(i).n_sec) then
                write(999, "(a$)"), '}'
            else
                write(999, "(a$)"), '},'
            end if
        end do
    end do

    write(999, "(a$)"), ']'
    write(999, "(a$)"), "}"

    ! Count minimum number of xovers of scaffold and staple
    min_xover = 999999
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            n_xover = edge(i).sec(j).n_xover_scaf / 2 + edge(i).sec(j).n_xover_stap / 2
            if(min_xover > n_xover) then
                min_xover      = n_xover
                min_xover_edge = i
                min_xover_sec  = j
            end if
        end do
    end do

    ! Print minimum number of crossovers per edge
    dna.min_xover_scaf = edge(min_xover_edge).sec(min_xover_sec).n_xover_scaf / 2
    dna.min_xover_stap = edge(min_xover_edge).sec(min_xover_sec).n_xover_stap / 2

    ! Deallocate edge data
    do i = 1, n_edge
        do j = 1, edge(i).n_sec
            deallocate(edge(i).sec(j).scaf)
            deallocate(edge(i).sec(j).stap)
        end do
        deallocate(edge(i).sec)
    end do
    deallocate(edge)
end subroutine Output_Write_Out_JSON

! -----------------------------------------------------------------------------

! Write information related with basepair
subroutine Output_Write_Basepair(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    character(200) :: path
    integer :: i

    if(para_write_801 == .false.) return

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 801, file = trim(path)//"_basepair.txt", form = "formatted")

    ! The number of base pairs
    write(801, "(i7)"), mesh.n_node

    ! Write information related with base pair
    do i = 1, mesh.n_node
        write(801, "(i7$   )"), i
        write(801, "(3f8.2$)"), mesh.node(i).pos(1:3)
        write(801, "(3f8.2$)"), mesh.node(i).ori(1, 1:3)
        write(801, "(3f8.2$)"), mesh.node(i).ori(2, 1:3)
        write(801, "(3f8.2$)"), mesh.node(i).ori(3, 1:3)
        write(801, "(i7$   )"), dna.base_scaf(i).xover
        write(801, "(i7    )"), dna.base_stap(i).xover
    end do
    write(801, "(a)")

    ! The number connectivity of base pairs
    write(801, "(i7)"), mesh.n_ele

    ! Write connectivity of basepair
    do i = 1, mesh.n_ele
        write(801, "(i7, 2i8)"), i, mesh.ele(i).cn(1:2)
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write basepair"

    close(unit = 801)
end subroutine Output_Write_Basepair

! -----------------------------------------------------------------------------

! Write information related with base
subroutine Output_Write_Base(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    character(200) :: path
    integer :: i, id, id_up, id_down

    if(para_write_802 == .false.) return

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 802, file = trim(path)//"_base.txt", form = "formatted")

    ! The total number of bases
    write (802, "(i7)"), dna.n_base_scaf + dna.n_base_stap

    ! For scaffold strand
    do i = 1, dna.n_top

        ! ID, up ID, down ID, across ID, position vector
        write(802, "(i7$   )"), dna.top(i).id
        write(802, "(i7$   )"), dna.top(i).up
        write(802, "(i7$   )"), dna.top(i).dn
        write(802, "(i7$   )"), dna.top(i).across
        write(802, "(3f10.4)"), dna.top(i).pos(1:3)
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write nucleotides"

    close(unit = 802)
end subroutine Output_Write_Base

! -----------------------------------------------------------------------------

! Write cndo file for PDB atom generation and CanDo simulation
subroutine Output_Write_CanDo(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    integer, allocatable :: node_nt(:,:)

    double precision :: min_pos(3), pos(3), fscale
    integer :: i, j, base, count, strand
    character(200) :: path

    if(para_write_803 == .false.) return

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 803, file = trim(path)//".cndo", form = "formatted")

    write(803, "(a)"), '"CanDo (.cndo) file format version 1.0"'
    write(803, "(a)")

    ! For dnatop data that is defined by bases
    write(803, "(a)"), 'dnaTop,id,up,down,across,seq'
    do i = 1, dna.n_top

        ! ID, up, down, across and sequence
        write(803, "(a$)"), &
            trim(adjustl(Int2Str(i)))//","//&
            trim(adjustl(Int2Str(dna.top(i).id)))//","//&
            trim(adjustl(Int2Str(dna.top(i).dn)))//","//&   ! cndo convention is opposite
            trim(adjustl(Int2Str(dna.top(i).up)))//","//&   ! cndo convention is opposite
            trim(adjustl(Int2Str(dna.top(i).across)))//","

        write(803, "(a)"), trim(dna.top(i).seq)
    end do
    write(803, "(a)")

    ! Pre-calculation
    min_pos(1:3) = mesh.node(1).pos(1:3)
    do i = 1, mesh.n_node
        if(min_pos(1) > mesh.node(i).pos(1)) min_pos(1) = mesh.node(i).pos(1)
        if(min_pos(2) > mesh.node(i).pos(2)) min_pos(2) = mesh.node(i).pos(2)
        if(min_pos(3) > mesh.node(i).pos(3)) min_pos(3) = mesh.node(i).pos(3)
    end do

    fscale = 10.0d0
    if(min_pos(1) < 0 .and. min_pos(1) < -99.0d0) then
        min_pos(1) = (min_pos(1) + 80.0d0)
    else
        min_pos(1) = 0.0d0
    end if

    if(min_pos(2) < 0 .and. min_pos(2) < -99.0d0) then
        min_pos(2) = (min_pos(2) + 80.0d0)
    else
        min_pos(2) = 0.0d0
    end if

    if(min_pos(3) < 0 .and. min_pos(3) < -99.0d0) then
        min_pos(3) = (min_pos(3) + 80.0d0)
    else
        min_pos(3) = 0.0d0
    end if

    ! The number of base pairs
    write(803, "(a)"), 'dNode,"e0(1)","e0(2)","e0(3)"'
    do i = 1, mesh.n_node

        pos(:) = (mesh.node(i).pos(:)-min_pos(:))*fscale

        write(803, "(a)"), &
            trim(adjustl(Int2Str(i)))//","//&
            trim(adjustl(Dble2Str(pos(1)))) //","//&
            trim(adjustl(Dble2Str(pos(2)))) //","//&
            trim(adjustl(Dble2Str(pos(3))))

        if( (pos(1) > 9999.0d0 .or. pos(1) < -999.0d0) .or. &
            (pos(2) > 9999.0d0 .or. pos(2) < -999.0d0) .or. &
            (pos(3) > 9999.0d0 .or. pos(3) < -999.0d0) ) then

            ! Print error for generating atomic model
            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== error =========================================+"
            write(p_redir, "(a)"), " | Out of bound of PDB                               |"
            write(p_redir, "(a)"), " +===================================================+"
            stop
        end if
    end do
    write(803, "(a)")

    ! Write orientations
    write(803, "(a$)"), 'triad,"e1(1)","e1(2)","e1(3)",'
    write(803, "(a$)"), '"e2(1)","e2(2)","e2(3)",'
    write(803, "(a )"), '"e3(1)","e3(2)","e3(3)"'

    do i = 1, mesh.n_node
        write(803, "(a)"), &
            trim(adjustl(Int2Str(i)))                       //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(1, 1)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(1, 2)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(1, 3)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(2, 1)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(2, 2)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(2, 3)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(3, 1)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(3, 2)))) //","//&
            trim(adjustl(Dble2Str(mesh.node(i).ori(3, 3))))
    end do
    write(803, "(a)")

    ! nt1, and nt2 based on basepair
    allocate(node_nt(mesh.n_node, 2))
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                if(dna.top(base).across /= -1) then
                    node_nt(dna.top(base).node, 1) = dna.top(base).id
                    node_nt(dna.top(base).node, 2) = dna.top(base).across
                end if
            end do
        end if
    end do

    ! Set base ID
    write(803, "(a)"), 'id_nt,id1,id2'
    do i = 1, mesh.n_node
        write(803, "(a)"), &
            trim(adjustl(Int2Str(i))) //","//&
            trim(adjustl(Int2Str(node_nt(i, 1))))//","//&
            trim(adjustl(Int2Str(node_nt(i, 2))))
    end do

    ! Additional information cndo format
    write(803, "(a)")

    ! Set viewpoint and RGB color
    write(803, "(3i)"), prob.color_atom(1:3)
    write(803, "(a)")

    ! For updated cndo with strand type
    do i = 1, dna.n_top
        strand = dna.top(i).strand
        if(dna.strand(strand).type1 == "scaf") then
            write(803, "(a)"), trim("0")
        else if(dna.strand(strand).type1 == "stap") then
            write(803, "(a)"), trim("1")
        end if
    end do
    write(803, "(a)")

    ! Deallocate memory
    deallocate(node_nt)

    ! Print progress
    write(p_redir, "(a)"), "   * Write CanDo file (*.cndo)"

    close(unit = 803)
end subroutine Output_Write_CanDo

! -----------------------------------------------------------------------------

! Write cndo file for PDB atom generation and CanDo simulation
subroutine Output_Write_CanDo_New(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    integer, allocatable :: node_nt(:,:)
    integer :: i, j, base, strand, types
    character(200) :: path

    if(para_write_803 == .false.) return

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 803, file = trim(prob.path_work)//"/cndo_format.cndo", form = "formatted")

    ! For dnatoop data that is defined by bases
    write(803, "(i10)"), dna.n_top
    do i = 1, dna.n_top

        ! ID, up ID, down ID, across ID and sequence
        write(803, "(i10$)"), i
        write(803, "(i10$)"), dna.top(i).id
        write(803, "(i10$)"), dna.top(i).up
        write(803, "(i10$)"), dna.top(i).dn
        write(803, "(i10$)"), dna.top(i).across
        write(803, "(a10$)"), dna.top(i).seq

        strand = dna.top(i).strand
        if(dna.strand(strand).type1 == "scaf") then
            write(803, "(i10)"), 0
        else if(dna.strand(strand).type1 == "stap") then
            write(803, "(i10)"), 1
        end if
    end do
    write(803, "(a)")

    ! The number of base pairs
    write(803, "(i10)"), mesh.n_node
    do i = 1, mesh.n_node
        write(803, "(i10$  )"), i
        write(803, "(f15.6$)"), mesh.node(i).pos(1)*10.0d0
        write(803, "(f15.6$)"), mesh.node(i).pos(2)*10.0d0
        write(803, "(f15.6 )"), mesh.node(i).pos(3)*10.0d0
    end do
    write(803, "(a)")

    ! e1 : the mojor groove
    ! e2 : the preferred nucleotide
    ! e3 : along the duplex axis towards the 3'-direction of
    !      the strand with the preferred nucleotide
    write(803, "(i10)"), mesh.n_node
    do i = 1, mesh.n_node
        write(803, "(i10$  )"), i
        write(803, "(3f15.6$)"), mesh.node(i).ori(1, 1:3)
        write(803, "(3f15.6$)"), mesh.node(i).ori(2, 1:3)
        write(803, "(3f15.6 )"), mesh.node(i).ori(3, 1:3)
    end do
    write(803, "(a)")

    ! nt1, and nt2 based on basepair
    allocate(node_nt(mesh.n_node, 2))
    do i = 1, dna.n_strand
        if(dna.strand(i).type1 == "scaf") then
            do j = 1, dna.strand(i).n_base
                base = dna.strand(i).base(j)
                if(dna.top(base).across /= -1) then
                    node_nt(dna.top(base).node, 1) = dna.top(base).id
                    node_nt(dna.top(base).node, 2) = dna.top(base).across
                end if
            end do
        end if
    end do

    ! Set base ID
    write(803, "(a)"), 'id_nt,id1,id2'
    do i = 1, mesh.n_node
        write(803, "(a)"), &
            trim(adjustl(Int2Str(i))) //","//&
            trim(adjustl(Int2Str(node_nt(i, 1))))//","//&
            trim(adjustl(Int2Str(node_nt(i, 2))))
    end do
    deallocate(node_nt)

    ! Print progress
    write(p_redir, "(a)"), "   * Write CanDo file (*.cndo)"

    close(unit = 803)
end subroutine Output_Write_CanDo_New

! -----------------------------------------------------------------------------

! Write PLY
subroutine Output_Write_PLY(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    integer :: i, j
    character(200) :: path

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 704, file = trim(path)//"_ply.ply", form = "formatted")

    write(704, "(a)"), "ply"
    write(704, "(a)"), "format ascii 1.0"
    write(704, "(a)"), "element vertex "//trim(adjustl(Int2Str(geom.n_iniP)))
    write(704, "(a)"), "property float32 x"
    write(704, "(a)"), "property float32 y"
    write(704, "(a)"), "property float32 z"
    write(704, "(a)"), "element face "//trim(adjustl(Int2Str(geom.n_face)))
    write(704, "(a)"), "property list uint8 int32 vertex_indices"
    write(704, "(a)"), "end_header"

    do i = 1, geom.n_iniP
        write(704, "(3f15.4)"), geom.iniP(i).pos(1:3)
    end do

    do i = 1, geom.n_face
        write(704, "(a$)"), trim(adjustl(Int2Str(geom.face(i).n_poi)))
        do j = 1, geom.face(i).n_poi
            write(704, "(a$)"), " "//trim(adjustl(Int2Str(geom.face(i).poi(j)-1)))
        end do
        write(704, "(a)")
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write PLY"

    close(unit = 704)
end subroutine Output_Write_PLY

! -----------------------------------------------------------------------------

! Write basepair and nucleotide information
subroutine Output_Write_DNA_Info(prob, dna)
    type(ProbType), intent(in) :: prob
    type(DNAType),  intent(in) :: dna

    integer :: i, id, id_up, id_down
    character(200) :: path

    ! open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 100, file = trim(path)//"_dnaInfo.dat", form = "formatted")

    ! the total number of nucleotide
    write (100, "(i8)"), dna.n_top

    ! for scaffold strand
    do i = 1, dna.n_top

        ! id, up id, down id, across id, position vector
        write(100, "(i8$   )"), i
        write(100, "(i8$   )"), dna.top(i).id
        write(100, "(i8$   )"), dna.top(i).up
        write(100, "(i8$   )"), dna.top(i).dn
        write(100, "(i8$   )"), dna.top(i).across
        write(100, "(3f10.4)"), dna.top(i).pos(1:3)
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write DNA Info"

    close(unit = 100)
end subroutine Output_Write_DNA_Info

! -----------------------------------------------------------------------------

! Write TecPlot input file
subroutine Output_Write_TecPlot(prob, mesh)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh

    character(200) :: path
    integer :: i, k

    if(para_write_804 == .false.) return

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 804, file = trim(path)//"_tecplot.dat", form = "formatted")

    ! Loop for 3 orientation vectors
    do k = 1, 3

        ! For Tecplot output
        write(804, "(a     )"), "variables = 'X', 'Y', 'Z', 'e3', 'e2', 'e1' "
        write(804, "(a, i7$)"), "ZONE F=FEPOINT, N=", mesh.n_node
        write(804, "(a, i7$)"), ", E=",               mesh.n_ele
        write(804, "(a     )"), ", ET=QUADRILATERAL"

        ! Write nodal position vector
        do i = 1, mesh.n_node
            write(804, "(3f10.4$)"), mesh.node(i).pos(1:3)
            write(804, "(3f10.4 )"), mesh.node(i).ori(k, 1:3)
        end do

        ! Element connectivity
        do i = 1, mesh.n_ele
            write(804, "(2i6$)"), mesh.ele(i).cn(1:2)
            write(804, "(2i6 )"), mesh.ele(i).cn(1:2)
        end do
    end do

    ! Print progress
    write(p_redir, "(a)"), "   * Write TecPlot input file"

    close(unit = 804)
end subroutine Output_Write_TecPlot

! -----------------------------------------------------------------------------

! Write ADINA input file
subroutine Output_Write_ADINA(prob, mesh)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh

    character(200) :: path
    integer :: i, time(8)

    if(para_write_805 == .false.) return

    call date_and_time(VALUES=time)

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 805, file = trim(path)//"_adina.in", form = "formatted")

    ! Write of ADINA infile
    write(805, "(a     )"), "*"
    write(805, "(a     )"), "* Command file created from METIS for command import"
    write(805, "(a     )"), "*"
    write(805, "(a$    )"), "*--- Command file created "
    write(805, "(i4, a$)"), time(2), "/"
    write(805, "(i4, a$)"), time(3), "/"
    write(805, "(i4, a$)"), time(1), ", "
    write(805, "(i2, a$)"), time(5), ":"
    write(805, "(i2, a$)"), time(6), ":"
    write(805, "(i2, a$)"), time(7), " ---*"
    write(805, "(a     )"), "*--- for ADINA: AUI version 9.0.6 ---*"
    write(805, "(a     )"), "*"
    write(805, "(a     )"), "DATABASE NEW SAVE=NO PROMPT=NO"    ! That creates a new database
    write(805, "(a     )"), "FEPROGRAM ADINA"                   ! For displacement and stress analysis
    write(805, "(a     )"), "COORDINATES POINT SYSTEM=0"        ! Definition of coordinates for geometry points
    write(805, "(a     )"), "@CLEAR"
    write(805, "(a     )")

    ! Write point information
    do i = 1, mesh.n_node
        write(805, "(i, 3f13.5, i)"), i, mesh.node(i).pos(:), 0
    end do
    write(805, "(a)") "@"

    ! Write line information
    do i = 1, mesh.n_ele
        write(805, "(a, i$)"), "LINE STRAIGHT NAME=", i
        write(805, "(a, i$)"), " P1=", mesh.ele(i).cn(1)
        write(805, "(a, i )"), " P2=", mesh.ele(i).cn(2)
        write(805, "(a    )"), "*"
    end do

    ! Set material properties
    write(805, "(a)"), "*"
    write(805, "(a)"), "EGROUP BEAM NAME=1 SUBTYPE=THREE-D DISPLACE=DEFAULT MATERIAL=1 RINT=5,"
    write(805, "(a)"), "     SINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES INITIALS=NONE,"
    write(805, "(a)"), "     CMASS=DEFAULT RIGIDEND=NONE MOMENT-C=NO RIGIDITY=1,"
    write(805, "(a)"), "     MULTIPLY=1000000.00000000 RUPTURE=ADINA OPTION=NONE,"
    write(805, "(a)"), "     BOLT-TOL=0.00000000000000 DESCRIPT='NONE' SECTION=1,"
    write(805, "(a)"), "     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,"
    write(805, "(a)"), "     TDEATH=0.00000000000000 SPOINT=2 BOLTFORC=0.00000000000000,"
    write(805, "(a)"), "     BOLTNCUR=0 TMC-MATE=1 BOLT-NUM=0 BOLT-LOA=0.00000000000000,"
    write(805, "(a)"), "     WARP=NO" ! ENDRELEA=ACCURATE"

    ! Geometry lines are assigned a number of subdivisions
    write(805, "(a)"), "*"
    write(805, "(a)"), "SUBDIVIDE LINE NAME=1 MODE=DIVISIONS NDIV=1 RATIO=1.00000000000000,"
    write(805, "(a)"), "      PROGRESS=GEOMETRIC CBIAS=NO"
    write(805, "(a)"), "@CLEAR"

    ! Element connectivity
    do i = 1, mesh.n_ele
        write(805, "(i7)") i
    end do
    write(805, "(a)") "@"

    ! Mesh generation
    write(805, "(a)"), "*"
    write(805, "(a)"), "GLINE NODES=2 AUXPOINT=0 NCOINCID=ENDS NCENDS=12,"
    write(805, "(a)"), "     NCTOLERA=1.00000000000000E-05 SUBSTRUC=0 GROUP=1 MIDNODES=CURVED,"
    write(805, "(a)"), "     XO=0.00000000000000 YO=1.00000000000000 ZO=0.00000000000000,"
    write(805, "(a)"), "     XYZOSYST=SKEW"
    write(805, "(a)"), "@CLEAR"
    do i = 1, mesh.n_ele
        write(805, "(i7)") i
    end do
    write(805, "(a)") "@"

    ! Print progress
    write(p_redir, "(a)"), "   * Write ADINA input file (*.in)"
    write(p_redir, "(a)")

    ! Print progress
    write(p_redir, "(a)"), "   * Write ADINA input"

    close(unit = 805)
end subroutine Output_Write_ADINA

! -----------------------------------------------------------------------------

! Write sequence based on cross-sectional line
subroutine Output_Write_Sequence_CroL(prob, mesh, dna)
    type(ProbType), intent(in) :: prob
    type(MeshType), intent(in) :: mesh
    type(DNAType),  intent(in) :: dna

    character, allocatable :: seq(:)
    integer :: i, j, base, node
    character(200) :: path

    if(para_write_808 == .false.) return

    ! Allocate memory
    allocate(seq(mesh.n_node))

    ! Find sequence from top data
    do i = 1, dna.n_strand

        ! For staple strand
        if(dna.strand(i).type1 /= "stap") cycle
        do j = 1, dna.strand(i).n_base
            base = dna.strand(i).base(j)
            node = dna.top(base).node

            if(node /= -1) then
                seq(node) = dna.top(base).seq
            end if
        end do
    end do

    ! Open files
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 808, file = trim(path)//"_seq_line.txt", form = "formatted")

    ! Write information on sequence
    do i = 1, mesh.n_node
        write(808, "(i5$)"), mesh.node(i).croL
        write(808, "(a$ )"), " multi line -> node ID : "
        write(808, "(i5$)"), mesh.node(i).id
        write(808, "(a$ )"), ", BP ID : "
        write(808, "(i5$)"), mesh.node(i).bp
        write(808, "(a$ )"), ", sequence for staple : "
        write(808, "(a  )"), seq(i)
    end do

    ! Deallocate memory
    deallocate(seq)

    ! Print progress
    write(p_redir, "(a)"), "   * Write Sequence lines"

    close(unit = 808)
end subroutine Output_Write_Sequence_CroL

! -----------------------------------------------------------------------------

end module Output