! Ishaan Desai (desaii)
#include <SMAAspUserSubroutines.hdr>

      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
      include 'vaba_param.inc'

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
      double precision, dimension(:), allocatable :: couplingVertices, writeData, readData
      integer, dimension(:), allocatable :: vertexIDs
      integer, dimension(1) :: numberOfVerticesArray

      double precision, dimension(:), allocatable :: coords_local, coords_global
    
      ! Note that you  can use the MPI communication between parallel Abaqus processes to gather 
      ! and scatter the data.
    
      ! Start of the analysis
      if (lOp .eq. j_int_StartAnalysis) then
            ! Get MPI rank and size (total number of MPI processors in this job)
            call GETNUMCPUS(size)
            call GETRANK(rank)

            ! Create preCICE participant
            call precicef_create("laminate_3ply", "../precice-config.xml", rank, size)

            ! Get problem dimensions from preCICE
            call precicef_get_mesh_dimensions(meshName, dimensions)

            ! Get number of vertices from VUMAT global array
            pointer(ptrNumberOfvertices, numberOfVerticesArray)
            ptrNumberOfvertices =  SMAIntArrayAccess(1)
            numberOfVertices = numberOfVerticesArray(1)

            allocate(vertices(numberOfVertices*dimensions))
            allocate(vertexIDs(numberOfVertices))
            allocate(readData(numberOfVertices*dimensions))
            allocate(writeData(numberOfVertices*dimensions))
            allocate(couplingVertices(numberOfVertices,dimensions))
            
            ! Get coordinates of vertices from VUMAT global array
            pointer(ptrcouplingVertices, couplingVertices)
            ptrcouplingVertices = SMAFloatArrayAccess(1001)

      !      continuation from a previous analysis (restart)
             if (kStep .ne. 0) then 
             end if 
    
      !   Start of the step
          else if (lOp .eq. j_int_StartStep) then
    
      !    Set up or exchange (import and export) initial values with external programs.      
          
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
