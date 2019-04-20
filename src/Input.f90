!
! =============================================================================
!
! Module - Input
! Last Updated : 01/21/2019, by Hyungmin Jun (hyungminjun@outlook.com)
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
module Input

    use Importer

    use Exam_METIS

    use Section

    use Para
    use Math
    use List

    implicit none

    public  Input_Initialize

    private Input_Set_Init
    private Input_Print_Prob_METIS
    private Input_Print_Edge_Len
    private Input_Set_Prob
    private Input_Set_Vertex_Design
    private Input_Sel_Prob_METIS
    private Input_Sel_File
    private Input_Set_Edge_Sec
    private Input_Find_Max_Min_Section
    private Input_Set_Section_Connectivity
    private Input_Set_Edge_Len
    private Input_Convert_Face_To_Line
    private Input_Scale_Init_Geom
    private Input_Set_Work_Dir
    private Input_Write_GEO_File
    private Input_Chimera_Init_Geom
    private Input_Tecplot_Init_Geom
    private Input_Generate_Schlegel_Diagram
    private Input_Print_Para

contains

! -----------------------------------------------------------------------------

! Initialize inputs
subroutine Input_Initialize(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    character(200) :: c_prob, c_seq, c_path_output
    character(30)  :: c_edge_sec, c_v_design, c_edge_ref
    character(30)  :: c_edge_len, c_mesh_spacing, c_run_mode
    integer :: len_char, ppos, ioerr, n_designer
    logical :: results, here

    ! Set initial inputs
    call Input_Set_Init(prob)

    if(iargc() == 0) then

        ! Print pre-defined problems
        call Input_Print_Prob_METIS
        read(*, *), c_prob

        ppos = scan(trim(c_prob), ".", BACK = .true.)
        if(ppos > 0) then
            len_char       = len_trim(c_prob)
            prob.name_file = trim(adjustl(c_prob(1:ppos-1)))
            prob.type_file = trim(adjustl(c_prob(ppos+1:len_char)))

            ! Check file
            inquire(file = trim(prob.path_input)//trim(prob.name_file)//"."//trim(prob.type_file), EXIST=here)
            if(here == .false.) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | The file does not exist.                         |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if
        else
            read(c_prob, *, iostat = ioerr), prob.sel_prob
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Error in input file name                         |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            ! Stop if the the negative value
            if(prob.sel_prob <= 0) stop
        end if

        ! Print and read the edge section
        prob.sel_edge_sec = 3

        ! Print and read the edge length
        call Input_Print_Edge_Len; read(*, *) prob.sel_edge_len; if(prob.sel_edge_len < 0) stop

        ! Reference edge - shortest edge
        prob.sel_edge_ref = 0

        ! Choose specific edge number and take edge-length as an input
        if(prob.sel_edge_len == 0) then

            write(0, "(a)"), "   Type the specific edge ID [Enter] : "
            write(0, "(a)")
            read(*, *) prob.sel_edge_ref
            write(0, "(a)"), "   Type the minimum edge length for the specific edge ID [Enter] : "
            write(0, "(a)")
            read(*, *) prob.sel_edge_len

            ! Stop if the the negative value
            if(prob.sel_edge_len <= 0) stop
        end if
    else if(iargc() == 9) then

        call getarg(1, c_path_output)    ! 1st argument, path for the output
        call getarg(2, c_prob)           ! 2nd argument, geometry file
        call getarg(3, c_seq)            ! 3rd argument, scaf sequence file
        call getarg(4, c_edge_sec)       ! 4th argument, edge section
        call getarg(5, c_v_design)       ! 5th argument, vertex design
        call getarg(6, c_edge_ref)       ! 6th argument, edge
        call getarg(7, c_edge_len)       ! 7th argument, edge length
        call getarg(8, c_mesh_spacing)   ! 8th argument, meshing spacing
        call getarg(9, c_run_mode)       ! 9th argument, run mode

        ppos = scan(trim(c_prob), ".", BACK = .true.)
        if(ppos > 0) then
            ! Already get the name and type in Input_Set_Init
            !len_char       = len_trim(c_prob)
            !prob.name_file = trim(adjustl(c_prob(1:ppos-1)))
            !prob.type_file = trim(adjustl(c_prob(ppos+1:len_char)))
            !
            !! Check file
            !inquire(file = trim(prob.path_input)//trim(prob.name_file)//"."//trim(prob.type_file), EXIST=here)
            !if(here == .false.) then
            !    write(p_redir, "(a)")
            !    write(p_redir, "(a)"), " +=== error ========================================+"
            !    write(p_redir, "(a)"), " | The file does not exist.                         |"
            !    write(p_redir, "(a)"), " +==================================================+"
            !    stop
            !end if
        else
            read(c_prob, *), prob.sel_prob

            ! Stop if the the negative value
            if(prob.sel_prob <= 0) stop
        end if

        ! Edge section
        read(c_edge_sec, *), prob.sel_edge_sec

        ! Vertex design
        read(c_v_design, *), prob.sel_vertex

        ! Edge reference
        read(c_edge_ref, *), prob.sel_edge_ref

        ! Edge length
        read(c_edge_len, *), prob.sel_edge_len

        ! Meshing parameter
        read(c_mesh_spacing, *), prob.p_mesh
    else

        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Not defined input arguments                      |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Set design parameters depending on the designer
    prob.sel_vertex = 2
    if(prob.sel_edge_sec /= 3) stop

    call Input_Set_Vertex_Design(prob)

    ! Set edge section
    call Input_Set_Edge_Sec(prob, geom)

    ! Set edge length
    call Input_Set_Edge_Len(prob, geom)

    ! Set problem
    call Input_Set_Prob(prob, geom)

    ! Convert faces to lines
    call Input_Convert_Face_To_Line(geom)

    ! Set initial geometric scale
    call Input_Scale_Init_Geom(geom, 100.0d0)

    ! Set working directory
    call Input_Set_Work_Dir(prob)

    ! Write *.geo file
    call Input_Write_GEO_File(prob, geom)

    ! Write initial geometry
    call Input_Chimera_Init_Geom(prob, geom)

    ! Write initial geometry for Tecplot
    call Input_Tecplot_Init_Geom(prob, geom)

    ! Generate Schlegel diagram
    call Input_Generate_Schlegel_Diagram(prob, geom)

    ! Print progress
    call Input_Print_Para(prob, geom)
end subroutine Input_Initialize

! -----------------------------------------------------------------------------

! Set initial inputs
subroutine Input_Set_Init(prob)
    type(ProbType), intent(inout) :: prob

    character(200) :: c_prob, c_seq, c_path_output
    logical :: here, results
    integer :: ioerr, ppos, len_char

    character(200) cwd
    integer :: istat

    !istat = getcwd(cwd)
    !print *, "cwd = "//trim(cwd)

    if(p_redir /= 0) then
        open(unit = p_redir, file = "METIS.log",     form = "formatted")
    end if

    ! --------------------------------------------------
    ! Set the input directory
    ! --------------------------------------------------
    if(iargc() == 0) then

        !DEC$ IF DEFINED(_WIN32)
        prob.path_input = "input\"
        !DEC$ ELSE
        prob.path_input = "input/"
        !DEC$ ENDIF

        ! --------------------------------------------------
        ! Set scaffold sequence from the txt file
        ! --------------------------------------------------
        inquire(file = "seq.txt", exist = here)
        if(here == .false.) then

            write(p_redir, "(a)")
            write(p_redir, "(a)"), " +=== warning ======================================+"
            write(p_redir, "(a)"), " | No 'seq.txt' file and M13 seq will be used.      |"
            write(p_redir, "(a)"), " +==================================================+"
            write(p_redir, "(a)")

            para_scaf_seq = "m13"
        else

            ! Open file
            open(unit = 1, file = "seq.txt", form = "formatted")
            read(1, *, iostat = ioerr), para_scaf_seq
            if(ioerr /= 0) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Please check file format in seq.txt.             |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            if(para_scaf_seq /= "m13" .and. para_scaf_seq /= "user" .and. para_scaf_seq /= "rand") then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | Please check the field in the seq.txt file.      |"
                write(p_redir, "(a)"), " +==================================================+"
                stop
            end if

            if(para_scaf_seq == "user") then

                ! Read scaffold sequence
                read(1, "(a)"), prob.scaf_seq

                ! Change upper case if low case
                prob.scaf_seq = Mani_To_Upper(prob.scaf_seq)
            end if
            close(unit = 1)
        end if
    else if(iargc() == 9) then

        ! 1st argument, path for the output
        call getarg(1, c_path_output)

        ! 2nd argument, geometry file
        call getarg(2, c_prob)

        ! 3rd argument, scaffold sequence file
        call getarg(3, c_seq)

        ! Find the path and filename
        !DEC$ IF DEFINED(_WIN32)
        ppos= scan(trim(c_prob), "\", BACK = .true.)
        prob.path_work = trim(c_path_output)
        !DEC$ ELSE
        ppos= scan(trim(c_prob), "/", BACK = .true.)
        prob.path_work = trim(c_path_output)
        !DEC$ ENDIF
        len_char = len_trim(c_prob)
        prob.path_input =trim(adjustl(c_prob(1:ppos)))
        c_prob = trim(adjustl(c_prob(ppos+1:len_char)))

        ppos = scan(trim(c_prob), ".", BACK = .true.)
        len_char       = len_trim(c_prob)
        prob.name_file = trim(adjustl(c_prob(1:ppos-1)))
        prob.type_file = trim(adjustl(c_prob(ppos+1:len_char)))

        ! --------------------------------------------------
        ! Set scaffold sequence from the txt file
        ! --------------------------------------------------
        ! 3rd argument, scaffold sequence
        call getarg(3, c_seq)
        if(trim(c_seq) == "m13" .or. trim(c_seq) == "M13") then

            ! Set m13mp18
            para_scaf_seq = "m13"
        else

            ! Set user scaffold sequence
            para_scaf_seq = "user"

            ! Open file
            inquire(file = trim(c_seq), exist = here)

            if(here == .false.) then
                write(p_redir, "(a)")
                write(p_redir, "(a)"), " +=== error ========================================+"
                write(p_redir, "(a)"), " | No sequence file.                                |"
                write(p_redir, "(a)"), " +==================================================+"
                write(p_redir, "(a)")
                stop
            end if

            open(unit = 1, file = trim(c_seq), form = "formatted")

            ! Read scaffold sequence
            read(1, "(a)"), prob.scaf_seq

            ! Change upper case if low case
            prob.scaf_seq = Mani_To_Upper(prob.scaf_seq)

            close(unit = 1)
        end if
    else

        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | Please check input arguments.                    |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    ! Set command environments
    !DEC$ IF DEFINED(_WIN32)
    !results = systemqq('title Designer')                    ! cmd title
    !results = systemqq('mode con: cols=135 lines=6000')     ! cmd size
    !results = systemqq('color')                             ! convert color, 02, f0, f1, f2
    !results = systemqq('date /t')                           ! display time
    !results = systemqq('hostname')                          ! display hostname of the computer
    !results = systemqq('ver')                               ! display version information
    !DEC$ ENDIF
end subroutine Input_Set_Init

! -----------------------------------------------------------------------------

! Print pre-defined geometries for METIS
subroutine Input_Print_Prob_METIS
    write(0, "(a)")
    write(0, "(a)"), "   +===================================================================+"
    write(0, "(a)"), "   |                                                                   |"
    write(0, "(a)"), "   |          METIS for the 2D DNA wireframe with 6HB edges            |"
    write(0, "(a)"), "   |                                        by Hyungmin Jun            |"
    write(0, "(a)"), "   |                                                                   |"
    write(0, "(a)"), "   +===================================================================+"
    write(0, "(a)")
    write(0, "(a)"), "   A. First input - Pre-defined 2D target geometries"
    write(0, "(a)"), "   ================================================="
    write(0, "(a)")
    write(0, "(a)"), "    [ Triangular-mesh objects ]"
    write(0, "(a)"), "       1. Hexagon,        2. Hexagon Hole,  3. Square,      4. Circle"
    write(0, "(a)")
    write(0, "(a)"), "    [ Quadrilateral-mesh objects ]"
    write(0, "(a)"), "       5. Six Parallelogram,     6. Curved Beam"
    write(0, "(a)"), "       7. Quarter Circle,        8. Annulus"
    write(0, "(a)")
    write(0, "(a)"), "    [ Without the internal mesh - Radius 25nm ]"
    write(0, "(a)"), "       9. Triangle(128), 10. Square(105), 11. Hexagon(74), 12. Octagon(57)"
    write(0, "(a)")
    write(0, "(a)"), "    [ Letters ]"
    write(0, "(a)"), "      13. Letter-A,      14. Letter-T,    15. Letter-G,    16. Letter-C"
    write(0, "(a)")
    write(0, "(a)"), "    [ Variable vertex-number ]"
    write(0, "(a)"), "      17. 4-Side Poly,   18. 5-Side Poly, 19. 6-Side Poly"
    write(0, "(a)")
    write(0, "(a)"), "      20. Controllable Plate"
    write(0, "(a)")
    write(0, "(a)"), "   Select or type geometry file (*.geo, *.igs, *.svg, *.ply) [Enter]: "
    write(0, "(a)")
end subroutine Input_Print_Prob_METIS

! ---------------------------------------------------------------------------------------

! Print pre-defined edge sections
subroutine Input_Print_Edge_Sec
    write(0, "(a)")
    write(0, "(a)"), "   Opt 1. Edge Section - Pre-defined edge section"
    write(0, "(a)"), "   =============================================="
    write(0, "(a)")
    write(0, "(a)"), "    1.  @-@         2.  @-@      "
    write(0, "(a)"), "       @   @          =@   @=    "
    write(0, "(a)"), "       =@-@=            @-@      "
    write(0, "(a)"), "    [6HB-Inner]     [6HB-Middle] "
    write(0, "(a)")
    write(0, "(a)"), "   Select the number [Enter] : "
end subroutine Input_Print_Edge_Sec

! ---------------------------------------------------------------------------------------

! Print pre-defined edge sections
subroutine Input_Print_Edge_Sec_Old
    write(0, "(a)")
    write(0, "(a)"), "   B. Edge Section - Pre-defined edge section"
    write(0, "(a)"), "   =========================================="
    write(0, "(a)")
    write(0, "(a)"), "    1.            2.  @-@         3.  @-@      "
    write(0, "(a)"), "      =@ @=          @   @          =@   @=    "
    write(0, "(a)"), "                     =@-@=            @-@      "
    write(0, "(a)"), "      [DX]        [6HB-Inner]     [6HB-Middle] "
    write(0, "(a)")
    write(0, "(a)")
    write(0, "(a)"), "    4.  @-@      5.               6.    @-@     "
    write(0, "(a)"), "      =@   @=       @-@   @-@        @-@   @-@  "
    write(0, "(a)"), "        @-@       =@   @=@   @=    =@   @=@   @="
    write(0, "(a)"), "      =@   @=       @-@   @-@        @-@   @-@  "
    write(0, "(a)"), "        @-@                             @-@     "
    write(0, "(a)"), "      [10HB]         [12HB]           [16HB]    "
    write(0, "(a)")
    write(0, "(a)"), "   Select the number [Enter] : "
end subroutine Input_Print_Edge_Sec_Old

! -----------------------------------------------------------------------------

! Print pre-defined minimum edge lengths
subroutine Input_Print_Edge_Len

    write(0, "(a)")
    write(0, "(a)"), "   B. Second input - Pre-defined minimum edge lengths"
    write(0, "(a)"), "   =================================================="
    write(0, "(a)")
    write(0, "(a)"), "   * 1.  42 bp =  4 turn * 10.5 bp/turn ->  42 bp * 0.34nm/bp = 14.28nm"
    write(0, "(a)"), "     2.  52 bp =  5 turn * 10.5 bp/turn ->  52 bp * 0.34nm/bp = 17.85nm"
    write(0, "(a)"), "   * 3.  63 bp =  6 turn * 10.5 bp/turn ->  63 bp * 0.34nm/bp = 21.42nm"
    write(0, "(a)"), "     4.  73 bp =  7 turn * 10.5 bp/turn ->  73 bp * 0.34nm/bp = 24.99nm"
    write(0, "(a)"), "   * 5.  84 bp =  8 turn * 10.5 bp/turn ->  84 bp * 0.34nm/bp = 28.56nm"
    write(0, "(a)"), "     6.  94 bp =  9 turn * 10.5 bp/turn ->  94 bp * 0.34nm/bp = 32.13nm"
    write(0, "(a)"), "   * 7. 105 bp = 10 turn * 10.5 bp/turn -> 105 bp * 0.34nm/bp = 35.70nm"
    write(0, "(a)"), "     8. 115 bp = 11 turn * 10.5 bp/turn -> 115 bp * 0.34nm/bp = 39.27nm"
    write(0, "(a)"), "   * 9. 126 bp = 12 turn * 10.5 bp/turn -> 126 bp * 0.34nm/bp = 42.84nm"
    write(0, "(a)")
    write(0, "(a)"), "     0. If needed to choose the specific edge to assign the length"
    write(0, "(a)")
    write(0, "(a)"), "   Select the number or type the min. edge length [Enter] : "
    write(0, "(a)")
end subroutine Input_Print_Edge_Len

! -----------------------------------------------------------------------------

! Set vertex design
subroutine Input_Set_Vertex_Design(prob)
    type(ProbType), intent(inout) :: prob

    if(prob.sel_vertex == 1) then
        para_vertex_design = "flat"
    else
        para_vertex_design = "mitered"
    end if
end subroutine Input_Set_Vertex_Design

! -----------------------------------------------------------------------------

! Set problem
subroutine Input_Set_Prob(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set problem
    if(prob.sel_prob == 0) call Input_Sel_File(prob, geom)
    if(prob.sel_prob /= 0) then
        call Input_Sel_Prob_METIS(prob, geom)
    end if
end subroutine Input_Set_Prob

! -----------------------------------------------------------------------------

! Select geometry file
subroutine Input_Sel_File(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Select file format
    if(prob.type_file == "ply" .or. prob.type_file == "PLY") then
        call Importer_PLY(prob, geom)
    else if(prob.type_file == "geo" .or. prob.type_file == "GEO" .or. &
            prob.type_file == "igs" .or. prob.type_file == "IGS") then
        call Importer_GEO(prob, geom)
    else if(prob.type_file == "svg" .or. prob.type_file == "SVG") then
        call Importer_SVG(prob, geom)
    else
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | No support of the file extension as an input.    |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    prob.name_prob = prob.name_file
    call Mani_Set_Prob(prob, 'blue')

    ! Print filename and type
    write(p_redir, "(a)"), "   * File name: "//trim(prob.name_file)//"."//trim(prob.type_file)
    write(p_redir, "(a)")
end subroutine Input_Sel_File

! -----------------------------------------------------------------------------

! Select the pre-defined geometry for METIS
subroutine Input_Sel_Prob_METIS(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    ! Set file type as primitive
    prob.type_file = "primitive"

    ! Select problem
    select case (prob.sel_prob)

        ! Triangular-mesh objects
        case ( 1); call Exam_METIS_Hexagon_Tri  (prob, geom)
        case ( 2); call Exam_METIS_Hexagon_Hole (prob, geom)
        case ( 3); call Exam_METIS_Square_Tri   (prob, geom)
        case ( 4); call Exam_METIS_Circle       (prob, geom)

        ! Quadrilateral-mesh objects
        case ( 5); call Exam_METIS_Six_Parallelogram (prob, geom)
        case ( 6); call Exam_METIS_Curved_Beam       (prob, geom)
        case ( 7); call Exam_METIS_Quarter_Circle    (prob, geom)
        case ( 8); call Exam_METIS_Annulus           (prob, geom)

        ! Without the intenal mesh
        case ( 9); call Exam_METIS_Triangle (prob, geom)
        case (10); call Exam_METIS_Square   (prob, geom)
        case (11); call Exam_METIS_Hexagon  (prob, geom)
        case (12); call Exam_METIS_Octagon  (prob, geom)

        ! Letters
        case (13); call Exam_METIS_Letter_A (prob, geom)
        case (14); call Exam_METIS_Letter_T (prob, geom)
        case (15); call Exam_METIS_Letter_G (prob, geom)
        case (16); call Exam_METIS_Letter_C (prob, geom)

        ! Variable vertex-number
        case (17); call Exam_METIS_N_Sided_Poly (prob, geom, 4)
        case (18); call Exam_METIS_N_Sided_Poly (prob, geom, 5)
        case (19); call Exam_METIS_N_Sided_Poly (prob, geom, 6)

        case (20); call Exam_METIS_Control_Plate (prob, geom)

        case default
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | The pre-defined geometries does not exist.       |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end select
end subroutine Input_Sel_Prob_METIS

! -----------------------------------------------------------------------------

! Set the edge section
subroutine Input_Set_Edge_Sec(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    integer :: bp_id

    ! The cross-section is defined on the local coordinate, t3-t2
    !        t2
    !        ¡è
    !        |
    !     ---|------¡æ t3
    !        |
    ! The number of columns of the crosssection should be even
    if(prob.sel_edge_sec == 1) then

        ! PERDIX in Science Advances, para_start_bp_ID = 3 + 1
        ! In this case the 63-bp edge has three double crossovers
        if(para_start_bp_ID == -1) para_start_bp_ID = 3 + 1
        !if(para_start_bp_ID == -1) para_start_bp_ID = 13 + 1

        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 3, 13
        !     ¡Ü¡Ü     .00   01.  <-- reference axis
        geom.sec.dir      = 90
        geom.n_sec        = 2
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
    else if(prob.sel_edge_sec == 2) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 13 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 3, 4 / 13, 14       | caDNAno   02    (CW)
        !      ¡Ü¡Ü        04 03              |         03  01
        !     ¡Ü  ¡Ü      05   02             |         04  00
        !      ¡Ü¡Ü       .00 01.   <--- ref  |           05
        geom.sec.dir      = 90
        geom.n_sec        = 6
        geom.sec.ref_row  = 1
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 1; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 2
        geom.sec.id(3) = 2; geom.sec.posR(3) = 2; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 3; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 1
        geom.sec.id(6) = 5; geom.sec.posR(6) = 2; geom.sec.posC(6) = 1

        ! Vertex crash
        para_vertex_crash = "mod2"
    else if(prob.sel_edge_sec == 3) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 11 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 1 / 11              | caDNAno   02    (CW)
        !      ¡Ü¡Ü        05=04              |         03  01
        !     ¡Ü  ¡Ü     .00   03.  <--- ref  |         04  00
        !      ¡Ü¡Ü        01=02              |           05
        geom.sec.dir      = 150
        geom.n_sec        = 6
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id(1) = 0; geom.sec.posR(1) = 2; geom.sec.posC(1) = 1
        geom.sec.id(2) = 1; geom.sec.posR(2) = 1; geom.sec.posC(2) = 1
        geom.sec.id(3) = 2; geom.sec.posR(3) = 1; geom.sec.posC(3) = 2
        geom.sec.id(4) = 3; geom.sec.posR(4) = 2; geom.sec.posC(4) = 2
        geom.sec.id(5) = 4; geom.sec.posR(5) = 3; geom.sec.posC(5) = 2
        geom.sec.id(6) = 5; geom.sec.posR(6) = 3; geom.sec.posC(6) = 1

        ! Vertex crash
        if(para_vertex_crash == "mod") para_vertex_crash = "mod1"
    else if(prob.sel_edge_sec == 4) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 11 + 1
        bp_id = mod(para_start_bp_ID, 21)

        ! Starting BP - 1 / 11
        !      ¡Ü¡Ü        07 06
        !     ¡Ü  ¡Ü      08   05             | caDNAno   04  02    (CW)
        !      ¡Ü¡Ü        09=04              |         05  03  01
        !     ¡Ü  ¡Ü     .00   03.  <--- ref  |         06  08  00
        !      ¡Ü¡Ü        01=02              |           07  09
        geom.sec.dir      = 150
        geom.n_sec        = 10
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 2

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id( 1) = 0; geom.sec.posR( 1) = 2; geom.sec.posC( 1) = 1
        geom.sec.id( 2) = 1; geom.sec.posR( 2) = 1; geom.sec.posC( 2) = 1
        geom.sec.id( 3) = 2; geom.sec.posR( 3) = 1; geom.sec.posC( 3) = 2
        geom.sec.id( 4) = 3; geom.sec.posR( 4) = 2; geom.sec.posC( 4) = 2
        geom.sec.id( 5) = 4; geom.sec.posR( 5) = 3; geom.sec.posC( 5) = 2
        geom.sec.id( 6) = 5; geom.sec.posR( 6) = 4; geom.sec.posC( 6) = 2
        geom.sec.id( 7) = 6; geom.sec.posR( 7) = 5; geom.sec.posC( 7) = 2
        geom.sec.id( 8) = 7; geom.sec.posR( 8) = 5; geom.sec.posC( 8) = 1
        geom.sec.id( 9) = 8; geom.sec.posR( 9) = 4; geom.sec.posC( 9) = 1
        geom.sec.id(10) = 9; geom.sec.posR(10) = 3; geom.sec.posC(10) = 1

        ! Vertex crash
        if(para_vertex_crash == "mod") para_vertex_crash = "mod1"
    else if(prob.sel_edge_sec == 5) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 11 + 1
        bp_id = mod(para_start_bp_ID, 21)

        !                                              |     08
        !                                              |   09  07
        ! Starting BP - 1 / 11                         |   10  06  02    (CW)
        !      ¡Ü¡Ü  ¡Ü¡Ü        05=04   11=10             |     11  03  01
        !     ¡Ü  ¡Ü¡Ü  ¡Ü     .00   03=06   09. <--- ref  |         04  00
        !      ¡Ü¡Ü  ¡Ü¡Ü        01=02   07=08             | caDNAno   05
        geom.sec.dir      = 150
        geom.n_sec        = 12
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id( 1) =  0; geom.sec.posR( 1) = 2; geom.sec.posC( 1) = 1
        geom.sec.id( 2) =  1; geom.sec.posR( 2) = 1; geom.sec.posC( 2) = 1
        geom.sec.id( 3) =  2; geom.sec.posR( 3) = 1; geom.sec.posC( 3) = 2
        geom.sec.id( 4) =  3; geom.sec.posR( 4) = 2; geom.sec.posC( 4) = 2
        geom.sec.id( 5) =  4; geom.sec.posR( 5) = 3; geom.sec.posC( 5) = 2
        geom.sec.id( 6) =  5; geom.sec.posR( 6) = 3; geom.sec.posC( 6) = 1
        geom.sec.id( 7) =  6; geom.sec.posR( 7) = 2; geom.sec.posC( 7) = 1+2
        geom.sec.id( 8) =  7; geom.sec.posR( 8) = 1; geom.sec.posC( 8) = 1+2
        geom.sec.id( 9) =  8; geom.sec.posR( 9) = 1; geom.sec.posC( 9) = 2+2
        geom.sec.id(10) =  9; geom.sec.posR(10) = 2; geom.sec.posC(10) = 2+2
        geom.sec.id(11) = 10; geom.sec.posR(11) = 3; geom.sec.posC(11) = 2+2
        geom.sec.id(12) = 11; geom.sec.posR(12) = 3; geom.sec.posC(12) = 1+2

        ! Vertex crash
        if(para_vertex_crash == "mod") para_vertex_crash = "mod1"
    else if(prob.sel_edge_sec == 6) then

        if(para_start_bp_ID == -1) para_start_bp_ID = 11 + 1
        bp_id = mod(para_start_bp_ID, 21)

        !                                              |       08  12
        ! Starting BP - 1 / 11                         |     09  07  13
        !        ¡Ü¡Ü              15=14                 |     10  06  02    (CW)
        !      ¡Ü¡Ü  ¡Ü¡Ü        05=04   11=10             |       11  03  01
        !     ¡Ü  ¡Ü¡Ü  ¡Ü     .00   03=06   09. <--- ref  |       14  04  00
        !      ¡Ü¡Ü  ¡Ü¡Ü        01=02   07=08             | caDNAno 15   05
        !        ¡Ü¡Ü              13=12
        geom.sec.dir      = 150
        geom.n_sec        = 16
        geom.sec.ref_row  = 2
        geom.sec.ref_minC = 1
        geom.sec.ref_maxC = 4

        call Mani_Allocate_SecType(geom.sec, geom.n_sec)
        call Mani_Init_SecType    (geom.sec, geom.n_sec, "honeycomb")

        geom.sec.id( 1) =  0; geom.sec.posR( 1) = 2+2; geom.sec.posC( 1) = 1
        geom.sec.id( 2) =  1; geom.sec.posR( 2) = 1+2; geom.sec.posC( 2) = 1
        geom.sec.id( 3) =  2; geom.sec.posR( 3) = 1+2; geom.sec.posC( 3) = 2
        geom.sec.id( 4) =  3; geom.sec.posR( 4) = 2+2; geom.sec.posC( 4) = 2
        geom.sec.id( 5) =  4; geom.sec.posR( 5) = 3+2; geom.sec.posC( 5) = 2
        geom.sec.id( 6) =  5; geom.sec.posR( 6) = 3+2; geom.sec.posC( 6) = 1
        geom.sec.id( 7) =  6; geom.sec.posR( 7) = 2+2; geom.sec.posC( 7) = 1+2
        geom.sec.id( 8) =  7; geom.sec.posR( 8) = 1+2; geom.sec.posC( 8) = 1+2
        geom.sec.id( 9) =  8; geom.sec.posR( 9) = 1+2; geom.sec.posC( 9) = 2+2
        geom.sec.id(10) =  9; geom.sec.posR(10) = 2+2; geom.sec.posC(10) = 2+2
        geom.sec.id(11) = 10; geom.sec.posR(11) = 3+2; geom.sec.posC(11) = 2+2
        geom.sec.id(12) = 11; geom.sec.posR(12) = 3+2; geom.sec.posC(12) = 1+2
        geom.sec.id(13) = 12; geom.sec.posR(13) = 0+2; geom.sec.posC(13) = 3
        geom.sec.id(14) = 13; geom.sec.posR(14) = 0+2; geom.sec.posC(14) = 2
        geom.sec.id(15) = 14; geom.sec.posR(15) = 4+2; geom.sec.posC(15) = 3
        geom.sec.id(16) = 15; geom.sec.posR(16) = 4+2; geom.sec.posC(16) = 2

        ! Vertex crash
        if(para_vertex_crash == "mod") para_vertex_crash = "mod1"
    else

        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | The cross-section is not defined.                |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    !! Set section connectivity in the defined initial section
    call Input_Set_Section_Connectivity(prob, geom)

    ! Find maximum and minimum sectional row and column
    call Input_Find_Max_Min_Section(geom)
end subroutine Input_Set_Edge_Sec

! -----------------------------------------------------------------------------

! Find maximum and minimum sectional row and column
subroutine Input_Find_Max_Min_Section(geom)
    type(GeomType), intent(inout) :: geom

    integer :: i

    ! sec.maxR ~ minC were initialized from Mani_Init_SecType
    do i = 1, geom.n_sec
        ! Find max and min of row
        if(geom.sec.posR(i) > geom.sec.maxR) geom.sec.maxR = geom.sec.posR(i)
        if(geom.sec.posR(i) < geom.sec.minR) geom.sec.minR = geom.sec.posR(i)

        ! Find max and min of col
        if(geom.sec.posC(i) > geom.sec.maxC) geom.sec.maxC = geom.sec.posC(i)
        if(geom.sec.posC(i) < geom.sec.minC) geom.sec.minC = geom.sec.posC(i)
    end do

    ! Find the size of columns and rows
    geom.sec.n_row = geom.sec.maxR - geom.sec.minR + 1
    geom.sec.n_col = geom.sec.maxC - geom.sec.minC + 1

    ! Check even number
    if(mod(geom.sec.n_col, 2) /= 0) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | The column number should be even.                |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if
end subroutine Input_Find_Max_Min_Section

! -----------------------------------------------------------------------------

! Set section connectivity in the defined initial section
subroutine Input_Set_Section_Connectivity(prob, geom)
    type(ProbType), intent(in)    :: prob
    type(GeomType), intent(inout) :: geom

    integer :: sec_cur, sec_com, row_cur, row_com
    integer :: i, j, count, id_bp
    logical :: b_connect

    ! If self connection route turns on, find connection to link each other
    ! Section connection was initialized as -1
    ! Loop for current section
    do i = 1, geom.n_sec
        sec_cur = geom.sec.id(i)
        row_cur = geom.sec.posR(i)

        ! Loop for comparing section
        do j = 1, geom.n_sec
            sec_com = geom.sec.id(j)
            row_com = geom.sec.posR(j)

            ! If section numbers are the same
            if(sec_cur == sec_com) cycle

            ! Determine the section connection for scaffold strand
            id_bp     = 1
            b_connect = Section_Connection_Scaf(geom, sec_cur, sec_com, id_bp)

            ! Set section connectivity
            if( (para_vertex_design == "flat"    .and. b_connect == .true.) .or. &
                (para_vertex_design == "mitered" .and. para_vertex_crash == "mod1" .and. &
                 row_cur < geom.sec.ref_row .and. row_cur == row_com ) ) then
                geom.sec.conn(i) = sec_com
                exit
            end if
        end do
    end do

    ! For alternative 6-helice bundle
    if(prob.sel_edge_sec == 2 .and. para_vertex_design == "flat" .and. para_vertex_crash == "mod2") then
        geom.sec.conn(1) = -1       ! Sec 0 -> neighbor
        geom.sec.conn(2) = -1       ! Sec 1 -> neighbor
        geom.sec.conn(3) = 3        ! Sec 2 -> 3
        geom.sec.conn(4) = 2        ! Sec 3 -> 2
        geom.sec.conn(5) = 5        ! Sec 4 -> 5
        geom.sec.conn(6) = 4        ! Sec 5 -> 4
    end if

    ! Print information for self connection route
    count = 0
    write(0, "(a)"), "   --------------------------------------------------"
    write(0, "(a)"), "           Connections for the edge section          "
    write(0, "(a)"), "   --------------------------------------------------"
    do i = 1, geom.n_sec

        write(0, "(i11, a$)"), geom.sec.id(i), " section  ->"

        if(geom.sec.conn(i) /= -1) then
            write(0, "(i3, a)"), geom.sec.conn(i), " section"
        else
            write(0, "(a)"), "  Neighbor connection"
            count = count + 1
        end if
    end do
    write(0, "(a)"), "   --------------------------------------------------"
    write(0, "(a)")

    ! Check neighboring connection
    ! There should be at least two neighboring connection
    ! Also, these number should be even
    if(mod(count, 2) /= 0 .or. count == 0) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | The section connect was wrong.                   |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if
end subroutine Input_Set_Section_Connectivity

! -----------------------------------------------------------------------------

! Set the minimum edge length
subroutine Input_Set_Edge_Len(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(in)    :: geom

    if(prob.sel_edge_len == 1) prob.n_edge_len =  42     ! 10.5bp/turn *  4 turn
    if(prob.sel_edge_len == 2) prob.n_edge_len =  52     ! 10.5bp/turn *  5 turn
    if(prob.sel_edge_len == 3) prob.n_edge_len =  63     ! 10.5bp/turn *  6 turn
    if(prob.sel_edge_len == 4) prob.n_edge_len =  73     ! 10.5bp/turn *  7 turn
    if(prob.sel_edge_len == 5) prob.n_edge_len =  84     ! 10.5bp/turn *  8 turn
    if(prob.sel_edge_len == 6) prob.n_edge_len =  94     ! 10.5bp/turn *  9 turn
    if(prob.sel_edge_len == 7) prob.n_edge_len = 105     ! 10.5bp/turn * 10 turn
    if(prob.sel_edge_len == 8) prob.n_edge_len = 115     ! 10.5bp/turn * 11 turn
    if(prob.sel_edge_len == 9) prob.n_edge_len = 126     ! 10.5bp/turn * 12 turn

    if(prob.sel_edge_len >= 10 .and. prob.sel_edge_len <= 37) then
        write(p_redir, "(a)")
        write(p_redir, "(a)"), " +=== error ========================================+"
        write(p_redir, "(a)"), " | The minimum edge length should be over 38-bp.    |"
        write(p_redir, "(a)"), " +==================================================+"
        stop
    end if

    if(prob.sel_edge_len >= 21) then
        prob.n_edge_len = prob.sel_edge_len
    end if
end subroutine Input_Set_Edge_Len

! -----------------------------------------------------------------------------

! Convert surface to line connectivity
subroutine Input_Convert_Face_To_Line(geom)
    type(GeomType), intent(inout) :: geom
    
    integer :: i, j, k, f_zero, flag
    integer :: point_1, point_2, point_1_com, point_2_com

    ! Mesh data structure
    type :: MeshType
        integer :: cn(100)   ! Maximum connectivity
    end type MeshType

    type(MeshType), allocatable :: Basepair_con(:)  ! 1st: # of meshes, 2nd: points
    type(ListConn), pointer :: straight_con
    type(ListConn), pointer :: ptr, ptr1

    allocate(Basepair_con(geom.n_face))

    ! Nullify the linked list for the junction data
    nullify(straight_con)
    nullify(ptr)
    nullify(ptr1)

    ! Initialze variable
    f_zero      = 0
    geom.n_iniL = 0

    ! Read mesh
    do i = 1, geom.n_face

        ! Read the number of vectices in the mesh and connectivity
        Basepair_con(i).cn(1:geom.face(i).n_poi) = geom.face(i).poi(1:geom.face(i).n_poi)

        ! If there is zero value in connectivity
        do j = 1, geom.face(i).n_poi
            if(Basepair_con(i).cn(j) == 0) f_zero = 1
        end do

        do j = 1, geom.face(i).n_poi

            flag = 1    ! Flag 1 means that there is no connectivity in the existing array

            if(j == geom.face(i).n_poi) then
                point_1 = Basepair_con(i).cn(j)
                point_2 = Basepair_con(i).cn(1)
            else
                point_1 = Basepair_con(i).cn(j)
                point_2 = Basepair_con(i).cn(j+1)
            end if

            if(i == 1) then     ! First connectivities are always saved
                flag = 1
            else
                allocate(ptr1)
                ptr1 => straight_con

                do k = 1, geom.n_iniL

                    point_1_com = ptr1%point(1)
                    point_2_com = ptr1%point(2)

                    ! Check where there is the same connectivity in the existing array
                    if((point_1_com == point_1 .and. point_2_com == point_2) ) then
                        flag = 0
                        exit
                    else if((point_1_com == point_2 .and. point_2_com == point_1) ) then
                        flag = 0
                        exit
                    end if

                    ! Pointer to move next list
                    ptr1 => ptr1%next
                end do
            end if

            ! Set connectivity adding this into list
            if(flag == 1) then

                allocate(ptr)
                geom.n_iniL  = geom.n_iniL + 1
                ptr%point(1) = point_1
                ptr%point(2) = point_2
                straight_con => List_Insert_Conn(straight_con, ptr)

                if(i == 1 .and. geom.face(i).n_poi == 2) exit
            end if

        end do
    end do

    ! Allocate stragiht line data structure
    allocate(geom.iniL(geom.n_iniL))

    ! Initialize line data type
    call Mani_Init_LineType(geom.iniL, geom.n_iniL)

    ! Set straight line information
    ptr1 => straight_con
    do i = 1, geom.n_iniL

        ! If there is zero in connectivity, add one
        if(f_zero == 1) then
            geom.iniL(geom.n_iniL+1-i).poi(1:2) = ptr1%point(1:2) + 1
        else
            geom.iniL(geom.n_iniL+1-i).poi(1:2) = ptr1%point(1:2)
        end if

        ! Pointer to move next linked list
        ptr1 => ptr1%next
    end do

    if(f_zero == 1) then
        do i = 1, geom.n_face
            do j = 1, geom.face(i).n_poi
                geom.face(i).poi(j) = geom.face(i).poi(j) + 1
            end do
        end do
    end if

    ! Copy from poi to iniP that always save initial points
    do i = 1, geom.n_iniL
        geom.iniL(i).iniP(1:2) = geom.iniL(i).poi(1:2)
    end do

    ! Delete linked list allocated
    call List_Delete_Conn(straight_con)
    !call List_Delete_Conn(ptr)
    !call List_Delete_Conn(ptr1)

    ! Deallocate mesh connectivity data
    deallocate(Basepair_con)
end subroutine Input_Convert_Face_To_Line

! -----------------------------------------------------------------------------

! Set initial geometric scale
subroutine Input_Scale_Init_Geom(geom, scale)
    type(GeomType),   intent(inout) :: geom
    double precision, intent(in)    :: scale

    double precision :: pos_c(3), length, min_len
    integer :: i, poi_1, poi_2

    ! Find center position of the structure
    pos_c(1:3) = 0.0d0
    do i = 1, geom.n_iniP
        pos_c(1:3) = pos_c + geom.iniP(i).pos
    end do
    pos_c(1:3) = pos_c / dble(geom.n_iniP)

    ! Move the geometry to the centered position
    do i = 1, geom.n_iniP
        geom.iniP(i).pos(1:3) = geom.iniP(i).pos(1:3) - pos_c(1:3)
    end do

    ! Find the edge that has the minimum length
    poi_1 = geom.iniL(1).poi(1)
    poi_2 = geom.iniL(1).poi(2)
    min_len = Norm(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)

    do i = 2, geom.n_iniL

        poi_1  = geom.iniL(i).poi(1)
        poi_2  = geom.iniL(i).poi(2)
        length = Norm(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)

        ! If the short edge exists
        if(length < min_len) then
            min_len = length
        end if
    end do

    ! Recalucate the edge length
    do i = 1, geom.n_iniP
        geom.iniP(i).pos(1:3) = geom.iniP(i).pos / min_len * scale
    end do

    ! Find edge that has minimum length
    !do i = 1, geom.n_iniL
    !    poi_1  = geom.iniL(i).poi(1)
    !    poi_2  = geom.iniL(i).poi(2)
    !    length = Norm(geom.iniP(poi_1).pos - geom.iniP(poi_2).pos)
    !    write(0, "(i, f15.5)"), i, length
    !end do
end subroutine Input_Scale_Init_Geom

! -----------------------------------------------------------------------------

! Set working directory
subroutine Input_Set_Work_Dir(prob)
    type(ProbType), intent(inout) :: prob

    character :: run_mode
    logical   :: results

    if(iargc() == 0) then

        ! Set the working directory
        prob.path_work = "outputs"
    else if(iargc() == 9) then

        ! 9th argument, run_mode
        call getarg(9, run_mode)

        if(run_mode == "m") then
            !DEC$ IF DEFINED(_WIN32)
            prob.path_work = trim(prob.path_work)//"\"//trim(prob.name_file)
            !DEC$ ELSE
            prob.path_work = trim(prob.path_work)//"/"//trim(prob.name_file)
            !DEC$ ENDIF
        end if
    end if

    !DEC$ IF DEFINED(_WIN32)
    results = systemqq("rd "//'"'//trim(prob.path_work)//'"'//' /s /q')
    results = systemqq("md "//'"'//trim(prob.path_work)//'"')
    if(para_tecplot == .true.) then
        results = systemqq("md "//'"'//trim(prob.path_work)//"\tecplot"//'"')
    end if
    !DEC$ ELSE
    results = systemqq("rm -rf "//'"'//trim(prob.path_work)//'"')
    results = systemqq("mkdir -p "//'"'//trim(prob.path_work)//'"')
    if(para_tecplot == .true.) then
        results = systemqq("mkdir -p "//'"'//trim(prob.path_work)//'"'//"/tecplot")
    end if
    !DEC$ ENDIF

    write(0, "(a)"), "   ...Removed the existing working directory"
    write(0, "(a)")
    if(p_redir /= 0) write(0, "(a)"), "   ...Please wait..."
end subroutine Input_Set_Work_Dir

! -----------------------------------------------------------------------------

! Write geo file file
subroutine Input_Write_GEO_File(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    character(200) :: path
    integer :: i, j

    ! Exception
    if(para_write_101 == .false.) return

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 101, file = trim(path)//".geo", form = "formatted")

    ! Write points
    write(101, "(i7)"), geom.n_iniP
    do i = 1, geom.n_iniP
        write(101,"(i7, 3f10.4)"), i, geom.iniP(i).pos(1:3)
    end do
    write(101, "(a)")

    ! Write faces
    write(101, "(i7)"), geom.n_face
    do i = 1, geom.n_face
        write(101,"(i7, i16$)"), geom.face(i).n_poi, geom.face(i).poi(1)
        do j = 2, geom.face(i).n_poi - 1
            write(101,"(i7$)"), geom.face(i).poi(j)
        end do
        write(101,"(i7)"), geom.face(i).poi(geom.face(i).n_poi)
    end do
    write(101, "(a)")

    ! Write edges
    write(101, "(i7)"), geom.n_iniL
    do i = 1, geom.n_iniL
        write(101, "(3i8)"), 2, geom.iniL(i).poi(1:2)
    end do
    write(101, "(a)")

    ! Write geometric for input module
    write(101, "(a)")
    write(101, "(a)"), "--------------------------------------------------"
    write(101, "(a)")
    call Space(101, 4)
    write(101, "(a)") "! The number of points and faces"
    call Space(101, 4)
    write(101, "(a, i4)") "geom.n_iniP = ", geom.n_iniP
    call Space(101, 4)
    write(101, "(a, i4)") "geom.n_face = ", geom.n_face
    write(101, "(a)")
    call Space(101, 4)
    write(101, "(a)") "allocate(geom.iniP(geom.n_iniP))"
    call Space(101, 4)
    write(101, "(a)") "allocate(geom.face(geom.n_face))"
    write(101, "(a)")

    ! For points
    call Space(101, 4)
    write(101, "(a)") "! Set point position vectors"
    do i = 1, geom.n_iniP
        call Space(101, 4)
        write(101, "(a, i7, a$)"), "geom.iniP(", i, ").pos(1:3) = [ "
        write(101, "(f10.4, a$)"), geom.iniP(i).pos(1), "d0, "
        write(101, "(f10.4, a$)"), geom.iniP(i).pos(2), "d0, "
        write(101, "(f10.4, a )"), geom.iniP(i).pos(3), "d0 ]"
    end do
    write(101, "(a)")

    ! For faces
    call Space(101, 4)
    write(101, "(a)") "! Set face connectivity"
    do i = 1, geom.n_face
        call Space(101, 4)
        write(101, "(a, i7$)"), "geom.face(",            i
        write(101, "(a, i7$)"), ").n_poi = ",            geom.face(i).n_poi
        write(101, "(a, i7$)"), "; allocate(geom.face(", i
        write(101, "(a, i7$)"), ").poi(",                geom.face(i).n_poi
        write(101, "(a, i7$)"), ")); geom.face(",        i
        write(101, "(a, i7$)"), ").poi(1:",              geom.face(i).n_poi
        write(101, "(a$    )"), ") = ["

        do j = 1, geom.face(i).n_poi - 1
            write(101, "(i7, a$)"), geom.face(i).poi(j), ", "
        end do
        write(101, "(i7, a)"), geom.face(i).poi(geom.face(i).n_poi), " ]"
    end do

    write(101, "(a)")
    write(101, "(a)"), "--------------------------------------------------"
    write(101, "(a)")

    ! For edges
    write(101, "(i7)"), geom.n_iniL
    do i = 1, geom.n_iniL
        write(101, "(a, i7, a$)"), "line(", i, ", 1:2) = [ "
        write(101, "(i7, a$   )"), geom.iniL(i).poi(1), ", "
        write(101, "(i7, a    )"), geom.iniL(i).poi(2), " ]"
    end do
    write(101, "(a)")

    close(unit = 101)
end subroutine Input_Write_GEO_File

! -----------------------------------------------------------------------------

! Write initial geometry for Chimera
subroutine Input_Chimera_Init_Geom(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: length, pos_1(3), pos_2(3), pos_c(3)
    integer :: i, j
    logical :: f_info, f_axis
    character(200) :: path

    if(para_write_102 == .false.) return

    ! Set option
    f_axis = para_chimera_axis
    f_info = para_chimera_102_info

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 102, file = trim(path)//"_01_target_geometry.bild", form = "formatted")

    ! Write initial points
    write(102, "(a)"), ".color red"
    do i = 1, geom.n_iniP
        write(102, "(a$    )"), ".sphere "
        write(102, "(3f9.2$)"), geom.iniP(i).pos(1:3)
        write(102, "(1f9.2 )"), 0.75d0
    end do

    ! Write initial edges
    write(102, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL

        pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
        pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)

        write(102,"(a$    )"), ".cylinder "
        write(102,"(3f9.2$)"), pos_1(1:3)
        write(102,"(3f9.2$)"), pos_2(1:3)
        write(102,"(1f9.2 )"), 0.3d0
    end do

    ! Information on initial geometry
    if(f_info == .true.) then

        ! For points
        do i = 1, geom.n_iniP
            write(102, "(a$   )"), ".cmov "
            write(102, "(3f9.2)"), geom.iniP(i).pos + 1.0d0
            write(102, "(a    )"), ".color red"
            write(102, "(a    )"), ".font Helvetica 12 bold"
            write(102, "(i    )"), i
        end do

        ! For edges
        do i = 1, geom.n_iniL
            pos_1(1:3) = geom.iniP(geom.iniL(i).poi(1)).pos(1:3)
            pos_2(1:3) = geom.iniP(geom.iniL(i).poi(2)).pos(1:3)
            pos_c(1:3) = (pos_1(1:3) + pos_2(1:3)) / 2.0d0
            length     = Norm(pos_2 - pos_1)

            write(102, "(a$   )"), ".cmov "
            write(102, "(3f9.2)"), pos_c(1:3) + 0.5d0
            write(102, "(a    )"), ".color dark green"
            write(102, "(a    )"), ".font Helvetica 12 bold"
            write(102, "(i,    a$)"), i, "("
            write(102, "(f5.1, a )"), length, ")"
        end do

        ! For faces
        do i = 1, geom.n_face

            ! Find center position in mesh
            pos_c(1:3) = 0.0d0
            do j = 1, geom.face(i).n_poi
                pos_c(1:3) = pos_c + geom.iniP(geom.face(i).poi(j)).pos
            end do
            pos_c(1:3) = pos_c / dble(geom.face(i).n_poi)

            ! Write face number
            write(102, "(a$   )"), ".cmov "
            write(102, "(3f9.2)"), pos_c(1:3) + 1.0d0
            write(102, "(a    )"), ".color black"
            write(102, "(a    )"), ".font Helvetica 12 bold"
            write(102, "(i7   )"), i
        end do
    end if

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(102)
    close(unit = 102)

    ! --------------------------------------------------
    ! Write the FE format output
    ! --------------------------------------------------
    !open(unit = 102, file = trim(path)//"_19_FE_Format.txt", form = "formatted")
    !write(102, "(i)"), geom.n_iniP
    !do i = 1, geom.n_iniP
    !    write(102, "(4f10.3, 2i10)"), geom.iniP(i).pos(1:2), 0.0d0, 0.0d0, 1, 1
    !end do

    !write(102, "(i)"), geom.n_iniL
    !do i = 1, geom.n_iniL
    !    write(102, "(2i10)"), geom.iniL(i).poi(1), geom.iniL(i).poi(2)
    !end do
    !close(unit = 102)

    ! --------------------------------------------------
    ! Write the file for Tecplot
    ! --------------------------------------------------
    if(para_tecplot == .false.) return

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit = 102, file = trim(path)//"_01_target_geometry.dat", form = "formatted")

    write(102, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(102, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(102, "(a$)"), 'ZONE F = FEPOINT'
    write(102, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniP)))
    write(102, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(102, "(a )"), ', ET=LINESEG'

    ! Write points
    do i = 1, geom.n_iniP
        write(102, "(3f9.3$)"), geom.iniP(i).pos(1:3)
        write(102, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(102, "(1i7$)"), geom.iniL(i).poi(1)
        write(102, "(1i7 )"), geom.iniL(i).poi(2)
    end do

    close(unit = 102)
end subroutine Input_Chimera_Init_Geom

! -----------------------------------------------------------------------------

! Write initial geometry for Tecplot
subroutine Input_Tecplot_Init_Geom(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision :: length, pos_1(3), pos_2(3), pos_c(3)
    logical :: f_info, f_axis
    integer :: i, j, nline
    character(200) :: path

    if(para_write_103 == .false.) return

    ! Open file
    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 102, file = trim(path)//"_init_geo_face.dat", form = "formatted")

    ! Find the number of lines
    nline = 0
    do i = 1, geom.n_face
        nline = nline + geom.face(i).n_poi
    end do

    ! For Tecplot output
    write(102, "(a)"), 'VARIABLES = "X", "Y", "Z"'
    write(102, "(a)"), 'ZONE T    = "'//trim(prob.name_file)//'"'
    write(102, "(a)"), 'ZONETYPE  = FEPOLYGON'
    write(102, "(a)"), 'NODES     = '//trim(adjustl(Int2Str(geom.n_iniP)))
    write(102, "(a)"), 'ELEMENTS  = '//trim(adjustl(Int2Str(geom.n_face)))
    write(102, "(a)"), 'FACES     = '//trim(adjustl(Int2Str(nline)))
    write(102, "(a)"), 'NumConnectedBoundaryFaces   = 0'
    write(102, "(a)"), 'TotalNumBoundaryConnections = 0'
    write(102, "(a)")

    ! Write x-direction position
    do i = 1, geom.n_iniP
        write(102, "(f8.2$)"), geom.iniP(i).pos(1)
    end do
    write(102, "(a)")

    ! Write y-direction position
    do i = 1, geom.n_iniP
        write(102, "(f8.2$)"), geom.iniP(i).pos(2)
    end do
    write(102, "(a)")

    ! Write z-direction position
    do i = 1, geom.n_iniP
        write(102, "(f8.2$)"), geom.iniP(i).pos(3)
    end do
    write(102, "(a)"); write(102, "(a)")

    ! Write line connectivity
    write(102, "(a)"), "# Face Nodes List"
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi
            if(j == geom.face(i).n_poi) then
                write(102, "(2i8)"), geom.face(i).poi(j), geom.face(i).poi(1)
            else
                write(102, "(2i8)"), geom.face(i).poi(j), geom.face(i).poi(j+1)
            end if
        end do
    end do
    write(102, "(a)"); write(102, "(a)")

    ! # Face Left Elements (In this case, they all point to the single element in this zone)
    write(102, "(a)"), "# Face Left Elements"
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi
            write(102, "(i7$)"), i
        end do
    end do
    write(102, "(a)"); write(102, "(a)")

    ! # Face Right elements (0 means no boundary, -n means use the nth boundary connection)
    write(102, "(a)"), "# Face Right elements"
    do i = 1, geom.n_face
        do j = 1, geom.face(i).n_poi
            write(102, "(i7$)"), 0
        end do
    end do
    write(102, "(a)"); write(102, "(a)")

    close(unit = 102)
end subroutine Input_Tecplot_Init_Geom

! -----------------------------------------------------------------------------

! Generate Schlegel diagram
subroutine Input_Generate_Schlegel_Diagram(prob, geom)
    type(ProbType), intent(in) :: prob
    type(GeomType), intent(in) :: geom

    double precision, allocatable :: pos_xy(:,:), pos_nei(:,:)
    integer, allocatable :: face(:), vert_row(:), vert_col(:), nei(:)

    double precision :: angle, pos_mid(2)
    integer :: i, j, k, n_vert, n_rept, n_nei, max_n_vert, max_face, nbr
    logical :: b_continue

    ! Number of vertices and iterations to calculate pos_xy
    n_vert = geom.n_iniP
    n_rept = 10*n_vert

    ! Choose face that has the maximum number of vertices
    max_n_vert = geom.face(1).n_poi
    max_face   = 1
    do i = 2, geom.n_face
        if(geom.face(i).n_poi > max_n_vert) then
            max_n_vert = geom.face(i).n_poi
            max_face   = i
        end if
    end do

    ! Identify vertices associated with biggest face
    allocate(face(max_n_vert))
    do i = 1, max_n_vert
        face(i) = geom.face(max_face).poi(i)
    end do

    ! Initialize pos_xy
    allocate(pos_xy(n_vert, 2))
    pos_xy(1:n_vert, 1:2) = 0.0d0

    ! Set big face on unit circle
    angle = 2.0d0*pi / dble(max_n_vert)

    ! For each vertex in big face
    do i = 1, max_n_vert
        pos_xy(face(i), 1) = dcos(angle*dble(i-1))
        pos_xy(face(i), 2) = dsin(angle*dble(i-1))
    end do

    ! Calculate positions of other vertices through iterative process
    do i = 1, n_rept
        do j = 1, n_vert

            b_continue = .true.
            do k = 1, max_n_vert
                if(face(k) == j) then
                    b_continue = .false.
                    exit
                end if
            end do

            if(b_continue == .true.) then

                ! Find j th vert
                allocate(vert_row(geom.n_iniL))
                allocate(vert_col(geom.n_iniL))

                n_nei = 0
                do k = 1, geom.n_iniL
                    if(geom.iniL(k).poi(1) == j) then
                        n_nei           = n_nei + 1
                        vert_row(n_nei) = k
                        vert_col(n_nei) = 1
                    else if(geom.iniL(k).poi(2) == j) then
                        n_nei           = n_nei + 1
                        vert_row(n_nei) = k
                        vert_col(n_nei) = 2
                    end if
                end do

                ! Find neighbors
                allocate(nei(n_nei))
                nei(1:n_nei) = 0

                do k = 1, n_nei
                    nbr    = 3 - vert_col(k)
                    nei(k) = geom.iniL(vert_row(k)).poi(nbr)
                end do

                ! Get neighbor position
                allocate(pos_nei(n_nei,2))
                pos_nei(:, 1:2) = pos_xy(nei, 1:2)

                ! Find mid position
                pos_mid(1:2) = 0.0d0
                do k = 1, n_nei
                    pos_mid(1:2) = pos_mid(1:2) + pos_nei(k, 1:2)
                end do
                pos_mid(1:2) = pos_mid(1:2)/dble(n_nei)

                ! Store as new pos_xy for jth vertex
                pos_xy(j, 1:2) = pos_mid(1:2)

                ! Deallocate memory
                deallocate(vert_row, vert_col, pos_nei, nei)
            end if
        end do
    end do

    ! Write Schlegel diagram
    call Input_Chimera_Schlegel_Diagram(prob, geom, pos_xy)

    ! Deallocate memory
    deallocate(face)
    deallocate(pos_xy)
end subroutine Input_Generate_Schlegel_Diagram

! -----------------------------------------------------------------------------

! Write Schlegel diagram
subroutine Input_Chimera_Schlegel_Diagram(prob, geom, pos_xy)
    type(ProbType),   intent(in) :: prob
    type(GeomType),   intent(in) :: geom
    double precision, intent(in) :: pos_xy(:,:)

    double precision :: pos_1(3), pos_2(3)
    integer :: i, j
    logical :: f_axis
    character(200) :: path

    if(para_write_104 == .false.) return

    ! Set option
    f_axis = para_chimera_axis

    path = trim(prob.path_work)//"/"//trim(prob.name_file)
    open(unit = 104, file = trim(path)//"_schlegel.bild", form = "formatted")

    ! Write initial points
    write(104, "(a)"), ".color red"
    do i = 1, geom.n_iniP
        write(104, "(a$    )"), ".sphere "
        write(104, "(3f9.2$)"), pos_xy(i,1:2) * 30.0d0, 0.0d0
        write(104, "(1f9.2 )"), 0.75d0
    end do

    ! Write initial edges
    write(104, "(a)"), ".color dark green"
    do i = 1, geom.n_iniL
        pos_1(3)   = 0.0d0
        pos_2(3)   = 0.0d0
        pos_1(1:2) = pos_xy(geom.iniL(i).poi(1),1:2) * 30.0d0
        pos_2(1:2) = pos_xy(geom.iniL(i).poi(2),1:2) * 30.0d0

        write(104,"(a$    )"), ".cylinder "
        write(104,"(3f9.2$)"), pos_1(1:3)
        write(104,"(3f9.2$)"), pos_2(1:3)
        write(104,"(1f9.2 )"), 0.3d0
    end do

    ! Write global axis
    if(f_axis == .true.) call Mani_Set_Chimera_Axis(104)
    close(unit = 104)

    ! ---------------------------------------------
    ! Write the file for Tecplot
    ! ---------------------------------------------
    if(para_tecplot == .false.) return

    path = trim(prob.path_work)//"/tecplot/"//trim(prob.name_file)
    open(unit = 104, file = trim(path)//"_schlegel.dat", form = "formatted")

    write(104, "(a )"), 'TITLE = "'//trim(prob.name_file)//'"'
    write(104, "(a )"), 'VARIABLES = "X", "Y", "Z", "weight"'
    write(104, "(a$)"), 'ZONE F = FEPOINT'
    write(104, "(a$)"), ', N='//trim(adjustl(Int2Str(geom.n_iniP)))
    write(104, "(a$)"), ', E='//trim(adjustl(Int2Str(geom.n_iniL)))
    write(104, "(a )"), ', ET=LINESEG'

    ! Write vertices
    do i = 1, geom.n_iniP
        write(104, "(3f9.3$)"), pos_xy(i, 1:2), 0.0d0
        write(104, "(1f9.3 )"), 1.0d0
    end do

    ! Write edges
    do i = 1, geom.n_iniL
        write(104, "(1i7$)"), geom.iniL(i).poi(1)
        write(104, "(1i7 )"), geom.iniL(i).poi(2)
    end do

    close(unit = 104)
end subroutine Input_Chimera_Schlegel_Diagram

! -----------------------------------------------------------------------------

! Print progress
subroutine Input_Print_Para(prob, geom)
    type(ProbType), intent(inout) :: prob
    type(GeomType), intent(inout) :: geom

    write(p_redir, "(a)"), " +==================================================+"
    write(p_redir, "(a)"), " | 1. Inputs - geometry                             |"
    write(p_redir, "(a)"), " +==================================================+"
    write(p_redir, "(a)")
    write(p_redir, "(a)"), "  1.1. Geometry"
    write(p_redir, "(a)"), "   * Geometric name       : "//trim(prob.name_prob)
    write(p_redir, "(a)"), "   * Geometric file type  : "//trim(prob.type_file)
    write(p_redir, "(a)"), "   * The number of faces  : "//trim(adjustl(Int2Str(geom.n_face)))
    write(p_redir, "(a)"), "   * The number of points : "//trim(adjustl(Int2Str(geom.n_iniP)))
    write(p_redir, "(a)"), "   * The number of edges  : "//trim(adjustl(Int2Str(geom.n_iniL)))
    write(p_redir, "(a)"), "   * Input path           : "//trim(adjustl(prob.path_input))
    write(p_redir, "(a)"), "   * Output path          : "//trim(adjustl(prob.path_work))
    write(p_redir, "(a)")
    write(p_redir, "(a)"), "  1.2. Edge section"
    write(p_redir, "(a)"), "   * Section type           : "//trim(geom.sec.types)//" lattice"
    write(p_redir, "(a)"), "   * The number of duplexes : "//trim(adjustl(Int2Str(geom.n_sec)))
    write(p_redir, "(a)"), "   * The number of rows     : "//trim(adjustl(Int2Str(geom.sec.maxR-geom.sec.minR+1)))
    write(p_redir, "(a)"), "   * The number of columns  : "//trim(adjustl(Int2Str(geom.sec.maxC-geom.sec.minC+1)))
    write(p_redir, "(a)"), "   * Reference row          : "//trim(adjustl(Int2Str(geom.sec.ref_row)))
    write(p_redir, "(a)"), "   * Reference min/max col  : "//trim(adjustl(Int2Str(geom.sec.ref_minC)))//" / "//trim(adjustl(Int2Str(geom.sec.ref_maxC)))
    write(p_redir, "(a)")
    write(p_redir, "(a)"), "  1.3. Edge length"
    write(p_redir, "(a)"), "   * The minimum edge length : "//trim(adjustl(Int2Str(prob.n_edge_len)))
    write(p_redir, "(a)")
end subroutine Input_Print_Para

! -----------------------------------------------------------------------------

end module Input