!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
 MODULE generateparameters
!**************************************************************************************************
!==================================================================================================
! This module gets all the parameters relevant to the calculation from a file.
! It is one of the first call of the program.
! Author: Steve Ndengue (Haverford College)
! Date: 07/24/2020
! Last Change: 01/25/2026
!==================================================================================================
 INTEGER       :: unitlog, unitdvr, unitoutput               ! Units name identifiers
 INTEGER       :: unitinput, unitoperator, unitpsi           ! Units name identifiers
 CHARACTER(15) :: ca_type                                    ! Run calculation type
 INTEGER       :: ca_spec                                    ! Run calculation specification 
 CHARACTER(15) :: ca_sym                                     ! Run calculation symmetry
 INTEGER       :: ca_nmd                                     ! Run calculation Number of modes
 INTEGER       :: ca_nest                                    ! Run calculation Number of electronic states
 REAL(8)       :: r0                                         ! Run calculation r0, like CAP position
 REAL(8)       :: en_start, en_final, en_step                ! Run calculation energy: start, end, step
 INTEGER       :: Jmin, Jmax                                 ! Run calculation total J: min and max
 !INTEGER       :: pbdim                                     ! Full dimension of the primitive basis 
 TYPE pbastype                                               ! Primitive basis structure type
 character(5)  :: pb_nme                                     ! Name of primitive basis
 character(10) :: pb_typ                                     ! Type of Primitive basis
 integer       :: pb_nbr                                     ! Number of primitive basis functions
 integer       :: pb_pa1                                     ! Parameter 1 for primitive basis - k_x
 real(8)       :: pb_min                                     ! Minimum range for primitive basis
 real(8)       :: pb_max                                     ! Maximum range for primitive basis
 real(8)       :: pb_pa2                                     ! Parameter 2 for primitive basis
 END TYPE pbastype                                           ! End Primitive basis structure type
 TYPE(pbastype), ALLOCATABLE :: pbasst(:)                    ! Primitive basis structure vector
 real(8)       :: Arot, Brot, Crot                           ! System parameters - Rotational constants
 real(8)       :: mu_R                                       ! System parameters - Reduced mass
!
 LOGICAL       :: atom_atom, atom_diatom, atom_triatom, &    ! Calculation type - type of systems
                  diatom_diatom                              ! Calculation type - type of systems
 LOGICAL       :: closed_shell, open_shell                   ! Calculation type - molecule shell type
 LOGICAL       :: rigid_rotor, full_dimension                ! Calculation type - dimensionality
 NAMELIST  /Run_Parameters/ ca_type, ca_spec, ca_sym, &      ! Namelist Run Parameters
                            ca_nmd, ca_nest, r0, en_start, & ! Namelist Run Parameters
                            en_final, en_step, Jmin, Jmax    ! Namelist Run Parameters
 NAMELIST  /Primitive_Basis/ pbasst                          ! Namelist Primitive Basis structure
 NAMELIST  /System_Parameters/ mu_R, Arot, Brot, Crot        ! Namelist System parameters 
!
 CONTAINS
!========
!
!**************************************************************************************************
!**************************************************************************************************
 SUBROUTINE read_input
!**************************************************************************************************
!==================================================================================================
! This subroutine reads the 'input.nml' file to get all the relevant parameters
! that will be used for the SKVP calculation.
!
!       Date            :       11.21.2024
!       Author          :       Steve Ndengue (Haverford College)
!       Last Change     :       01.26.2026
!
!==================================================================================================
!USE AtomDiatomskvp
!
 IMPLICIT NONE
 INTEGER :: i, j, ier
! 
!Print parameters of the calculation
!-----------------------------------
 write(*,*) 'Parameters of the Computation:'
 write(*,*) ' '
 open(newunit=unitinput,file='input.nml',status='unknown')
 read(unitinput,nml=Run_Parameters,iostat=ier)
 write(*,Run_Parameters)
 allocate(pbasst(1:ca_nmd), stat=ier)
 read(unitinput,nml=Primitive_Basis,iostat=ier)
 write(*,Primitive_Basis)
 read(unitinput,nml=System_Parameters,iostat=ier)
 write(*,System_Parameters)
 !read(unitinput,nml=Calculation_Parameters,iostat=ier) 
 !read(unitinput,nml=SPF_Basis,iostat=ier)
 !write(*,SPF_Basis)
 !allocate(ifunst(1:ca_nmd), stat=ier)
 close(unitinput)

 atom_atom = .FALSE.
 atom_diatom = .FALSE.
 atom_triatom = .FALSE.
 diatom_diatom = .FALSE.
 closed_shell = .FALSE.
 open_shell = .FALSE.
 rigid_rotor = .FALSE.
 full_dimension = .FALSE.

 select case (ca_spec)

    case (1102)
      if (ca_nmd /= 5) then
         print*, " For diatom-diatom, ca_nmd should be 5."
         stop
      endif

      print*, " Rigid rotor Diatom-Diatom closed-shell system "
      diatom_diatom = .TRUE.
      closed_shell = .TRUE.
      rigid_rotor = .TRUE.

    case (1103)
    print*," Rigid rotor Diatom-Triatom closed-shell system " 
    atom_triatom = .TRUE.
    closed_shell = .TRUE.
    rigid_rotor = .TRUE.

    case default
    print*, " Provide a valid calculation specification *ca_spec* "

 end select
!==================================================================================================
!**************************************************************************************************
 END SUBROUTINE read_input
!**************************************************************************************************
!**************************************************************************************************
!
!**************************************************************************************************
!
 END MODULE generateparameters
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
