! Author: Ishaan Desai (desaii)
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

      include 'vaba_param.inc'
#include <SMAAspUserSubroutines.hdr>

      ! Contents of i_Array
      parameter(i_int_nTotalNodes=1,
     * i_int_nTotalElements = 2,
     * i_int_kStep = 3,
     * i_int_kInc = 4,
     * i_int_iStatus = 5,
     * i_int_lWriteRestart = 6,
     * i_int_ExtraOutputFrame = 7)

      ! Possible values for the lOp argument
      parameter(j_int_StartAnalysis=0,
     * j_int_StartStep = 1,
     * j_int_SetupIncrement = 2,
     * j_int_StartIncrement = 3,
     * j_int_EndIncrement = 4,
     * j_int_EndStep = 5,
     * j_int_EndAnalysis = 6)

      ! Possible values for i_Array(i_int_iStatus)
      parameter(j_int_Continue=0,
     * j_int_TerminateStep = 1,
     * j_int_TerminateAnalysis = 2)

      ! Contents of r_Array
      parameter(i_flt_TotalTime=1,
     * i_flt_StepTime = 2,
     * i_flt_dTime = 3)

      ! preCICE variables
      integer :: rank, size, ongoing, dimensions, bool, numberOfVertices
      double precision :: preCICE_dt
      double precision, dimension(:), allocatable :: strainsToWrite,
     * stressesToRead, couplingVertices
      integer, dimension(:), allocatable :: vertexIDs

      ! For defining shared arrays
      ! max3DSize = nblock * ndir, max6DSize = nblock * (ndir + nshr)
      parameter(max3DSize=1000, max6DSize=1000)

      integer :: nblock, ndir, nshr

      integer, dimension(3) :: intsFromVUMATArray
      pointer(ptr_intsFromVUMATArray, intsFromVUMATArray)

      double precision, dimension(max3DSize) :: macroVertices
      pointer(ptr_macroVertices, macroVertices)

      double precision, dimension(max6DSize) :: stresses, strains
      pointer(ptr_stresses, stresses)
      pointer(ptr_strains, strains)

      dimension i_Array(niArray), r_Array(nrArray)

      kStep = i_Array(i_int_kStep)
      kInc = i_Array(i_int_kInc)

      ! Start of the analysis
      if (lOp .eq. j_int_StartAnalysis) then
         ! Continuation from a previous analysis (restart)
         if (kStep .ne. 0) then
         end if

      ! Start of the step
      else if (lOp .eq. j_int_StartStep) then

         !  The initial values may need to match those at the point of restart.
         if (kInc .ne. 0) then
         end if

         write(*,*) "VEXTERNALDB: Start of the step."
         ! Get MPI rank and size (total number of MPI processors in this job)
         call vgetnumcpus(size)
         call vgetrank(rank)

         ! Create preCICE participant
         call precicef_create("Laminate-3D-ply",
     *    "/home/desaii/composite-multiscale/precice-config.xml", rank, size)

         write(*,*) "VEXTERNALDB: After precicef_create"

         ! Get problem dimensions from preCICE
         call precicef_get_mesh_dimensions("laminate-macro-mesh", dimensions)

         write(*,*) "VEXTERNALDB: After precicef_get_mesh_dimensions"

         ! Get number of vertices from VUMAT global array
         !ptrIntsFromVUMATArray = SMALocalIntArrayAccess(1000)
         !nblock = intsFromVUMATArray(1)
         !ndir = intsFromVUMATArray(2)
         !nshr = intsFromVUMATArray(3)

         write(*,*) "VEXTERNALDB: After first shared array access"

         ! Change from VUMAT terminology to preCICE terminology
         numberOfVertices = nblock

         allocate (couplingVertices(numberOfVertices*dimensions))
         allocate (vertexIDs(numberOfVertices))
         allocate (stressesToRead(numberOfVertices*dimensions))
         allocate (strainsToWrite(numberOfVertices*dimensions))

         ! Get coordinates of vertices from VUMAT global array
         ptr_macroVertices = SMALocalFloatArrayAccess(1001)

         write(*,*) "Macro coordinates: ", couplingVertices

         ! Set coupling mesh vertices in preCICE
         call precicef_set_vertices("laminate-macro-mesh",
     *    numberOfVertices, couplingVertices, vertexIDs)
         deallocate (couplingVertices)

         ! Set up or exchange (import and export) initial values with external programs.
         call precicef_requires_initial_data(bool)

         if (bool .eq. 1) then
            ptr_strains = SMALocalFloatArrayAccess(1002)
            call precicef_write_data("laminate-macro-mesh",
     *       "strains", numberOfVertices, vertexIDs, strainsToWrite)
         end if

         ! Initialize preCICE
         call precicef_initialize()
         call precicef_get_max_time_step_size(preCICE_dt)

         ! Create the stress array to share with VUMAT
         ptr_stresses = SMALocalFloatArrayCreate(1003,
     *    max6DSize, 0.0)

      ! Setup the increment
      else if (lOp .eq. j_int_SetupIncrement) then
      !    Change i_Array(i_int_lWriteRestart) and i_Array(i_int_iStatus) if desired.
      !    Change r_Array(i_flt_dTime) if desired.

      ! Start of the increment
      else if (lOp .eq. j_int_StartIncrement) then

         ! Check if coupling is still going on
         call precicef_is_coupling_ongoing(ongoing)

         if (ongoing .ne. 0) then
            ! Read stresses from preCICE
            call precicef_read_data("laminate-macro-mesh",
     *       "stresses", numberOfVertices, vertexIDs, dt, stressesToRead)
         end if

      ! End of the increment
      else if (lOp .eq. j_int_EndIncrement) then

         ! Check if coupling is still going on
         call precicef_is_coupling_ongoing(ongoing)

         if (ongoing .ne. 0) then
            ! Write strains to preCICE
            call precicef_write_data("laminate-macro-mesh", "strains",
     *       numberOfVertices, vertexIDs, strainsToWrite)

            ! Get preCICE time step
            call precicef_get_max_time_step_size(preCICE_dt)

            ! Get Abaqus time step
            dt = r_Array(i_flt_dTime)

            ! Reset Abaqus time step (origial value or preCICE_dt, whichever is smaller)
            dt = min(dt, preCICE_dt)
            r_Array(i_flt_dTime) = dt ! Set the dt in Abaqus

            ! Advance the coupling
            call precicef_advance(dt)
         end if

      ! Change i_Array(i_int_iStatus) if desired.
      ! Gather and export data  from the configuration at the end of the current increment
      ! to external programs.

      ! End of the step
      else if (lOp .eq. j_int_EndStep) then
         ! In the case of multiple steps, prepare the transition to the next step.
         ! For example, these data can serve as initial values for the next step.

      ! End of the analysis
      else if (lOp .eq. j_int_EndAnalysis) then
         ! Finalize the coupling
         call precicef_finalize()

      end if

      return
      end ! Subroutine





! Authors: Minh Hoang Nguyen (mhoangn) and Ishaan Desai (desaii)
      subroutine vumat(
      ! Read only (unmodifiable)variables -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
   ! Write only (modifiable) variables -
     7     stressNew, stateNew, enerInternNew, enerInelasNew)
   !
      include 'vaba_param.inc'
#include <SMAAspUserSubroutines.hdr>

      dimension props(nprops), density(nblock), coordMp(nblock, *),
     1     charLength(nblock), strainInc(nblock, ndir + nshr),
     2     relSpinInc(nblock, nshr), tempOld(nblock),
     3     stretchOld(nblock, ndir + nshr),
     4     defgradOld(nblock, ndir + nshr + nshr),
     5     fieldOld(nblock, nfieldv), stressOld(nblock, ndir + nshr),
     6     stateOld(nblock, nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock, ndir + nshr),
     8     defgradNew(nblock, ndir + nshr + nshr),
     9     fieldNew(nblock, nfieldv),
     1     stressNew(nblock, ndir + nshr), stateNew(nblock, nstatev),
     2     enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname

      ! = User-defined code ========================== [11, 22, 33, 12, 23, 13]

      parameter(zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, three=3.d0)
      parameter(pi=3.14159265358979323846264)

      parameter(i_prp_E11=1,
     * i_prp_E22 = 2,
     * i_prp_G12 = 3,
     * i_prp_nu12 = 4,
     * i_prp_nu23 = 5)

      parameter(i_sdv_eps11=1,
     * i_sdv_eps22 = 2,
     * i_sdv_eps33 = 3,
     * i_sdv_gamma12 = 4,
     * i_sdv_gamma23 = 5,
     * i_sdv_gamma13 = 6,
     * i_sdv_active = 7,
     * i_sdv_t_f = 8)

      ! For defining shared arrays
      ! directComponents == ndir, indirectComponents == nshr
      ! maxTensorComponents == ndir + nshr, max3DSize = nblock * ndir,
      ! max6DSize = nblock * (ndir + nshr)
      parameter(directComponents=3,
     * indirectComponents=3,
     * maxTensorComponents=6,
     * max3DSize=1000,
     * max6DSize=1000)

      real*8 :: state(nstatev), strains_total(6), stresses(6), Q(6,6)
      real*8 :: E11, E220, E33, nu12, nu23, nu13, G13, G23, nu31, nu21
      integer :: d, counter

      double precision, dimension(max3DSize) :: coordsToShare
      pointer(ptr_coordsToShare, coordsToShare)

      double precision, dimension(max6DSize) :: strainsToWrite, stressesToRead

      pointer(ptr_strainsToWrite, strainsToWrite)
      pointer(ptr_stressesToRead, stressesToRead)

      integer, dimension(3) :: intsToShare
      pointer(ptr_intsToShare, intsToShare)

      ! Initializations for stress calculation
      E11  = props(i_prp_E11)
      E22  = props(i_prp_E22)
      E33  = props(i_prp_E22)
      nu12 = props(i_prp_nu12)
      nu23 = props(i_prp_nu23)
      nu13 = props(i_prp_nu12)
      G12  = props(i_prp_G12)
      G13  = props(i_prp_G12)
      G23  = props(i_prp_E22)/2./(1.+props(i_prp_nu23))
      nu31 = nu13*E33/E11
      nu21 = nu12*E22/E11
      nu31 = nu13*E33/E11
      nu32 = nu23*E33/E22

      ! Create shared arrays
      ! Share integers
      ptr_intsToShare = SMALocalIntArrayCreate(1000, 3, 0)
      ! Share material point coordinates
      ptr_coordsToShare = SMALocalFloatArrayCreate(1001,
     * max3DSize, 0.0)
      ! Share strains
      ptr_strainsToWrite = SMALocalFloatArrayCreate(1002,
     * max6DSize, 0.0)

      ! Initialize -- run only once
      if (totalTime < 2.*dt) then

         state(i_sdv_eps11) = zero
         state(i_sdv_eps22) = zero
         state(i_sdv_eps33) = zero
         state(i_sdv_gamma12) = zero
         state(i_sdv_gamma23) = zero
         state(i_sdv_gamma13) = zero
         state(i_sdv_active) = one
         state(i_sdv_t_f) = zero

         ! Share integers nblock, ndir, nshr
         intsToShare(1) = nblock
         intsToShare(2) = ndir
         intsToShare(3) = nshr

         ! Put coordinates in the shared array
         counter = 1
         do k = 1, nblock
            do d = 1, directComponents
               coordsToShare(counter) = coordMp(k, d)
               counter = counter + 1
            end do
         end do

         write(*,*) "(t = ", totalTime, ") VUMAT: Initialization done."
      end if

      ! Loop through material points to collect strains
      counter = 1
      do k = 1, nblock

         state = stateOld(k, :)

         strains_total(1) = state(i_sdv_eps11) + strainInc(k, 1)
         strains_total(2) = state(i_sdv_eps22) + strainInc(k, 2)
         strains_total(3) = state(i_sdv_eps33) + strainInc(k, 3)
         strains_total(4) = state(i_sdv_gamma12) + two*strainInc(k, 4)
         strains_total(5) = state(i_sdv_gamma23) + two*strainInc(k, 5)
         strains_total(6) = state(i_sdv_gamma13) + two*strainInc(k, 6)
         state(i_sdv_eps11:i_sdv_gamma13) = strains_total

         do d = 1, maxTensorComponents
            strainsToWrite(counter) = strains_total(d)
            counter = counter + 1

            stateNew(k, d) = state(d)

         end do ! maxTensorComponents

      end do ! nblock

      write(*,*) "(t = ", totalTime, ") VUMAT: Strains collected."

      if (totalTime > 2.*dt) then ! Only after VEXTERNALDB is run

         write(*,*) "(t = ", totalTime, ") VUMAT: Trying to read stresses."

         ! Get stresses from VEXTERNALDB via shared array
         ptr_stressesToRead = SMALocalFloatArrayAccess(1003)

         ! Loop through material points to apply stresses
         counter = 1
         do k = 1, nblock
            do d = 1, maxTensorComponents

               stressNew(k, d) = stressesToRead(counter)
               counter = counter + 1

            end do ! maxTensorComponents
         end do ! nblock

         write(*,*) "(t = ", totalTime, ") VUMAT: Stresses applied."

      else ! Before VEXTERNALDB is run

         ! Calculate stresses for Abaqus preparation step
         temp   = nu12*nu21 + nu23*nu32 + nu13*nu31 + 2.*nu21*nu32*nu13
         Q      = zero
         Q(1,1) = (1.-nu23*nu32)*E11 / (1.-temp)
         Q(2,2) = (1.-nu13*nu31)*E22 / (1.-temp)
         Q(3,3) = (1.-nu12*nu21)*E33 / (1.-temp)
         Q(1,2) = (nu21+nu31*nu23)*E11/(1.-temp)
         Q(2,1) = Q(1,2)
         Q(1,3) = (nu31+nu21*nu32)*E11/(1.-temp)
         Q(3,1) = Q(1,3)
         Q(2,3) = (nu32+nu12*nu31)*E22/(1.-temp)
         Q(3,2) = Q(2,3)
         Q(4,4) = G12
         Q(5,5) = G23
         Q(6,6) = G13
         stresses = MATMUL(Q,strains_total)

         do k = 1, nblock
            stressNew(k, :) = stresses
         end do

         write(*,*) "(t = ", totalTime, ") VUMAT: Material model solved."

      end if

      write(*,*) "(t = ", totalTime, ") VUMAT: Complete."
      return
      end ! Subroutine
