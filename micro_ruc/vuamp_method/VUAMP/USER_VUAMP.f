

      SUBROUTINE VUAMP(
     1     ampName, time, ampValueOld, dt, nprops, props, nSvars, 
     2     svars, lFlagsInfo, nSensor, sensorValues, sensorNames,	
     3     jSensorLookUpTable,
     4     ampValueNew,
     5     lFlagsDefine,
     6     AmpDerivative, AmpSecDerivative, AmpIncIntegral)

      INCLUDE 'VABA_PARAM.INC'
      
!     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
!     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           ikStep            = 3,
     *           nFlagsInfo        = 3)
!     optional flags to be defined
      parameter (iComputeDeriv     = 1,
     *           iComputeSecDeriv  = 2,
     *           iComputeInteg     = 3,
     *           iStopAnalysis     = 4,
     *           iConcludeStep     = 5,
     *           nFlagsDefine      = 5)
      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine),
     *          sensorValues(nSensor),
     *          props(nprops),
     *          svars(nSvars)
      
      character*80 sensorNames(nSensor)
      character*80 ampName
      dimension jSensorLookUpTable(*)
      
      ! user coding to define AmpValueNew, and 
      ! optionally lFlagsDefine, AmpDerivative, AmpSecDerivative, AmpIncIntegral
      
      ! H: Our codes ==========================================================================
      
      REAL*8 :: RF11, RF22, RF33 ! from Abaqus via sensorValues
      REAL*8 :: sig11, sig22, sig33 ! to give to Micro Manager
      
      REAL*8 :: Lx, Ly, Lz ! can be defined as properties in the input file
      REAL*8 :: eps11, eps22, eps33 ! from Micro Manager
      
      ! print *, 'time(iTotalTime) = ', time(iTotalTime)
      
      Lx =  2.6E-03   ! props(1)
      Ly =  6.447E-03 ! props(2)
      Lz = 11.167E-03 ! props(3)
      
      ! --- Get stresses -------------------------------------------------------------------
      
      RF11 = vGetSensorValue('SENSOR_11',jSensorLookUpTable, sensorValues)
      RF22 = vGetSensorValue('SENSOR_22',jSensorLookUpTable, sensorValues)
      RF33 = vGetSensorValue('SENSOR_33',jSensorLookUpTable, sensorValues)
      
      sig11 = RF11/(Ly*Lz) ! to give to Micro Manager
      sig22 = RF22/(Lx*Lz) ! to give to Micro Manager
      sig33 = RF33/(Lx*Ly) ! to give to Micro Manager
      tau12 = 0.           ! to give to Micro Manager
      tau23 = 0.           ! to give to Micro Manager
      tau13 = 0.           ! to give to Micro Manager
      
      ! --- Apply strains ------------------------------------------------------------------
      
      eps11 = time(iTotalTime) * 0.01 ! example, Should get from Micro Manager
      eps22 = time(iTotalTime) * 0.01 ! example, Should get from Micro Manager
      eps33 = time(iTotalTime) * 0.01 ! example, Should get from Micro Manager
      gamma12 = 0.                    ! example, Should get from Micro Manager
      gamma23 = 0.                    ! example, Should get from Micro Manager
      gamma13 = 0.                    ! example, Should get from Micro Manager
      
      U11 = eps11*Lx
      U22 = eps22*Ly
      U33 = eps33*Lz
      
      IF     (ampName(1:6) .eq. 'AMP_11') THEN
         ampValueNew = U11
      ELSEIF (ampName(1:6) .eq. 'AMP_22') THEN
         ampValueNew = U22
      ELSEIF (ampName(1:6) .eq. 'AMP_33') THEN
         ampValueNew = U33
      END IF
      
      

      RETURN
      END
      
      
      
      