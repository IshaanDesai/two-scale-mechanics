! Ishaan Desai (desaii)
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
      include 'vaba_param.inc'
      #include <SMAAspUserSubroutines.hdr>

      !   Contents of i_Array
      parameter( i_int_nTotalNodes     = 1,
      *          i_int_nTotalElements  = 2,
      *          i_int_kStep           = 3,
      *          i_int_kInc            = 4,
      *          i_int_iStatus         = 5,
      *          i_int_lWriteRestart   = 6,  
      *          i_int_ExtraOutputFrame= 7 )

      !   Possible values for the lOp argument
      parameter( j_int_StartAnalysis    = 0,      
      *          j_int_StartStep        = 1,      
      *          j_int_SetupIncrement   = 2,      
      *          j_int_StartIncrement   = 3,      
      *          j_int_EndIncrement     = 4,      
      *          j_int_EndStep          = 5,      
      *          j_int_EndAnalysis      = 6 )     
    
      !    Possible values for i_Array(i_int_iStatus)
      parameter( j_int_Continue          = 0,      
      *          j_int_TerminateStep     = 1,      
      *          j_int_TerminateAnalysis = 2)      
    
      !    Contents of r_Array
      parameter( i_flt_TotalTime   = 1,
         *       i_flt_StepTime    = 2,
         *       i_flt_dTime       = 3 )

      dimension i_Array(niArray), r_Array(nrArray)
    
      kStep = i_Array(i_int_kStep)
      kInc  = i_Array(i_int_kInc)

      ! preCICE variables
      character*50 meshName, readDataName, writeDataName
      integer :: rank, size, ongoing, dimensions, bool, numberOfVertices
      double precision :: preCICE_dt
      integer, dimension(:), allocatable :: vertexIDs

      ! Variables acquired from VUMAT
      integer :: nblock, ndir, nshr
      integer, dimension(3) :: intsFromVUMATArray
      pointer(ptr_intsFromVUMATArray, intsFromVUMATArray)
      double precision, dimension(:), allocatable :: couplingVertices, stresses, strains
      pointer(ptr_couplingVertices, couplingVertices)
      pointer(ptr_stresses)
      pointer(ptr_strains, strains)
  
      ! Start of the analysis
      if (lOp .eq. j_int_StartAnalysis) then
            ! Get MPI rank and size (total number of MPI processors in this job)
            call GETNUMCPUS(size)
            call GETRANK(rank)

            ! Define names of mesh, read data and write data
            meshName = "laminate-macro-mesh"
            readDataName = "stresses"
            writeDataName = "strains"

            ! Create preCICE participant
            call precicef_create("laminate_3ply", "../precice-config.xml", rank, size)

            ! Get problem dimensions from preCICE
            call precicef_get_mesh_dimensions(meshName, dimensions)

            ! Get number of vertices from VUMAT global array
            ptrIntsFromVUMATArray =  SMAIntArrayAccess(1000)
            nblock = intsFromVUMATArray(1)
            ndir = intsFromVUMATArray(2)
            nshr = intsFromVUMATArray(3)

            ! Change from VUMAT terminology to preCICE terminology
            numberOfVertices = nblock

            allocate(vertices(numberOfVertices*dimensions))
            allocate(vertexIDs(numberOfVertices))
            allocate(readData(numberOfVertices*dimensions))
            allocate(writeData(numberOfVertices*dimensions))
            allocate(couplingVertices(numberOfVertices,dimensions))
            
            ! Get coordinates of vertices from VUMAT global array
            ptrcouplingVertices = SMAFloatArrayAccess(1001)

            ! Set coupling mesh vertices in preCICE
            call precicef_set_vertices(meshName, numberOfVertices, couplingVertices, vertexIDs)
            deallocate(couplingVertices)

      !     continuation from a previous analysis (restart)
            if (kStep .ne. 0) then 
            end if 
    
      ! Start of the step
        else if (lOp .eq. j_int_StartStep) then
            ! Set up or exchange (import and export) initial values with external programs.
            call precicef_requires_initial_data(bool)
            if (bool.eq.1) then
                ptr_stresses = SMAFloatArrayAccess(1002)
                call precicef_write_data(meshName, writeDataName, numberOfVertices, vertexIDs, stresses)
            end if
            
      !    The initial values may need to match those at the point of restart.
              if ( kInc .ne. 0) then
              end if 
        
      !     Setup the increment       
          else if (lOp .eq. j_int_SetupIncrement) then        
      !    Change i_Array(i_int_lWriteRestart) and i_Array(i_int_iStatus) if desired.      
      !    Change r_Array(i_flt_dTime) if desired.      
          
      !     Start of the increment
          else if (lOp .eq. j_int_StartIncrement) then
    
      !    The time increment is finalized.  Use r_Array(i_flt_dTime) if desired.
      !    If needed, gather and export data  from the configuration at the end of the previous                               increment to external programs.      
      !    Import and scatter data from external program to influence the current Abaqus increment.      
    
      !    End of the increment
          else if (lOp .eq. j_int_EndIncrement) then
    
      !    Change i_Array(i_int_iStatus) if desired.       
      !    Gather and export data  from the configuration at the end of the current increment 
      !    to external programs.      
    
      !    End of the step
          else if (lOp .eq. j_int_EndStep) then
    
      !    In the case of multiple steps, prepare the transition to the next step.
      !    For example, these data can serve as initial values for the next step.
    
      !     End of the analysis
          else if (lOp .eq. j_int_EndAnalysis) then
    
      !    User coding to close  files and disconnect any external programs, etc.      
    
          end if 
    
          return
          end
