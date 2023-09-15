! Ishaan Desai (desaii)
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
         include 'vaba_param.inc'
#include <SMAAspUserSubroutines.hdr>

         !   Contents of i_Array
         parameter(i_int_nTotalNodes=1,
     *             i_int_nTotalElements = 2,
     *             i_int_kStep = 3,
     *             i_int_kInc = 4,
     *             i_int_iStatus = 5,
     *             i_int_lWriteRestart = 6,
     *             i_int_ExtraOutputFrame = 7)

         !   Possible values for the lOp argument
         parameter(j_int_StartAnalysis=0,
     *             j_int_StartStep = 1,
     *             j_int_SetupIncrement = 2,
     *             j_int_StartIncrement = 3,
     *             j_int_EndIncrement = 4,
     *             j_int_EndStep = 5,
     *             j_int_EndAnalysis = 6)

         !    Possible values for i_Array(i_int_iStatus)
         parameter(j_int_Continue=0,
     *             j_int_TerminateStep = 1,
     *             j_int_TerminateAnalysis = 2)

         !    Contents of r_Array
         parameter(i_flt_TotalTime=1,
     *             i_flt_StepTime = 2,
     *             i_flt_dTime = 3)

         dimension i_Array(niArray), r_Array(nrArray)

         kStep = i_Array(i_int_kStep)
         kInc = i_Array(i_int_kInc)

         ! preCICE variables
         integer :: rank, size, ongoing, dimensions, bool, numberOfVertices
         double precision :: preCICE_dt
         integer, dimension(:), allocatable :: vertexIDs

         ! Variables acquired from VUMAT
         integer :: nblock, ndir, nshr
         integer, dimension(3) :: intsFromVUMATArray
         pointer(ptr_intsFromVUMATArray, intsFromVUMATArray)
         double precision, dimension(:), allocatable :: couplingVerticesCoords
         double precision, dimension(:), allocatable :: stresses, strains
         pointer(ptr_couplingVertices, couplingVertices)
         pointer(ptr_stresses, stresses)
         pointer(ptr_strains, strains)

         ! For defining shared arrays (same as VUMAT)
         parameter(maxMaterialPts=5000,
     *             directComponents=3,
     *             indirectComponents=3)

         ! Start of the analysis
         if (lOp .eq. j_int_StartAnalysis) then
            ! Get MPI rank and size (total number of MPI processors in this job)
            call getnumcpus(size)
            call getrank(rank)

            ! Create preCICE participant
            call precicef_create("laminate_3ply", "../precice-config.xml", rank, size)

            ! Get problem dimensions from preCICE
            call precicef_get_mesh_dimensions("laminate-macro-mesh", dimensions)

            ! Get number of vertices from VUMAT global array
            ptrIntsFromVUMATArray = SMAIntArrayAccess(1000)
            nblock = intsFromVUMATArray(1)
            ndir = intsFromVUMATArray(2)
            nshr = intsFromVUMATArray(3)

            ! Change from VUMAT terminology to preCICE terminology
            numberOfVertices = nblock

            allocate (vertices(numberOfVertices*dimensions))
            allocate (vertexIDs(numberOfVertices))
            allocate (readData(numberOfVertices*dimensions))
            allocate (writeData(numberOfVertices*dimensions))
            allocate (couplingVertices(numberOfVertices, dimensions))

            ! Get coordinates of vertices from VUMAT global array
            ptrcouplingVertices = SMAFloatArrayAccess(1001)

            ! Set coupling mesh vertices in preCICE
            call precicef_set_vertices("laminate-macro-mesh", numberOfVertices, couplingVertices, vertexIDs)
            deallocate (couplingVertices)

            ! Continuation from a previous analysis (restart)
            if (kStep .ne. 0) then
            end if

            ! Start of the step
         else if (lOp .eq. j_int_StartStep) then
            ! Set up or exchange (import and export) initial values with external programs.
            call precicef_requires_initial_data(bool)
            if (bool .eq. 1) then
               ptr_strains = SMALocalFloatArrayAccess(1002)
               call precicef_write_data("laminate-macro-mesh",
     *              "strains", numberOfVertices, vertexIDs, strains)
            end if

            ! Initialize preCICE
            call precicef_initialize()
            call precicef_get_max_time_step_size(preCICE_dt)

            ! Create the stress array to share with VUMAT
            ptr_stresses = SMALocalFloatArrayCreate(1003,
     *      (maxMaterialPts, directComponents + indirectComponents), 0.0)

            !  The initial values may need to match those at the point of restart.
            if (kInc .ne. 0) then
            end if

            !   Setup the increment
         else if (lOp .eq. j_int_SetupIncrement) then
            !    Change i_Array(i_int_lWriteRestart) and i_Array(i_int_iStatus) if desired.
            !    Change r_Array(i_flt_dTime) if desired.

            !     Start of the increment
         else if (lOp .eq. j_int_StartIncrement) then

            ! Check if coupling is still going on
            call precicef_is_coupling_ongoing(ongoing)

            if (ongoing .ne. 0) then
               ! Read stresses from preCICE
               call precicef_read_data("laminate-macro-mesh",
     *              "stresses", numberOfVertices, vertexIDs, dt, stresses)
            end if

            !    End of the increment
         else if (lOp .eq. j_int_EndIncrement) then

            ! Check if coupling is still going on
            call precicef_is_coupling_ongoing(ongoing)

            if (ongoing .ne. 0) then
               ! Write strains to preCICE
               call precicef_write_data("laminate-macro-mesh", "strains",
               *numberOfVertices, vertexIDs, writeData)

               ! Get preCICE time step
               call precicef_get_max_time_step_size(preCICE_dt)

               ! Get Abaqus time step
               dt = r_Array(i_flt_dTime)

               ! Reset Abaqus time step (origial value or preCICE_dt, whichever is smaller)
               dt = min(dt, preCICE_dt)
               r_Array(i_flt_dTime) = dt

               ! Advance the coupling
               call precicef_advance(dt)
            end if

            !    Change i_Array(i_int_iStatus) if desired.
            !    Gather and export data  from the configuration at the end of the current increment
            !    to external programs.

            !    End of the step
         else if (lOp .eq. j_int_EndStep) then

            !    In the case of multiple steps, prepare the transition to the next step.
            !    For example, these data can serve as initial values for the next step.

            !     End of the analysis
         else if (lOp .eq. j_int_EndAnalysis) then
            ! Finalize the coupling
            call precicef_finalize()

         end if

         return
      end
