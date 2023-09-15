! mhoangn and desaii
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
!
      character*80 cmname

      ! = User-defined code ======================================================= [11, 22, 33, 12, 23, 13]

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
      parameter(maxMaterialPts=1000,
     * maxTensorComps=6,
     * directComponents=3,
     * indirectComponents=3)

! - Variables --------------------------------------------------------- [11, 22, 33, 12, 23, 13]

      real*8 :: state(nstatev)
      real*8 :: strains_total(6), stresses(6), Q(6, 6)

      real*8 :: E11, E220, E33, nu12, nu23, nu13, G13, G23, nu31, nu21

      integer :: i, j

      double precision, dimension(maxMaterialPts, directComponents) :: coordsToShare
      pointer(ptr_coordsToShare, coordsToShare)
      double precision, dimension(maxTensorComps) :: strainsToWrite
      pointer(ptr_strainsToWrite, strainsToWrite)
      double precision, dimension(maxMaterialPts, directComponents + indirectComponents) :: stressesToRead
      pointer(ptr_stressesToRead, stressesToRead)
      integer, dimension(3) :: intsToShare
      pointer(ptr_intsToShare, intsToShare)

      ! Create shared array to share strains with VEXTERNALDB
      ptr_strainsToWrite = SMALocalFloatArrayCreate(1002,
     * (maxMaterialPts, directComponents + indirectComponents), 0.0)

! = Loop through points =========================================================================================================
      do k = 1, nblock

         ! --- Initialize ----------------------------------------------------------------------------------------------
         state = stateOld(k, :)

         if (totalTime < 2.*dt) then ! First increment (run only once)

            state(i_sdv_eps11) = zero
            state(i_sdv_eps22) = zero
            state(i_sdv_eps33) = zero
            state(i_sdv_gamma12) = zero
            state(i_sdv_gamma23) = zero
            state(i_sdv_gamma13) = zero

            state(i_sdv_active) = one
            state(i_sdv_t_f) = zero

            ! Share integers nblock, ndir, nshr
            ptr_intsToShare = SMALocalIntArrayCreate(1000, 3, 0)
            intsToShare(1) = nblock
            intsToShare(2) = ndir
            intsToShare(3) = nshr
      
            ! Share coordinates of material points
            ptr_coordsToShare = SMALocalFloatArrayCreate(1001,
     *      (maxMaterialPts, directComponents), 0.0)
            coordsToShare = coordMp

         end if

         strains_total(1) = state(i_sdv_eps11) + strainInc(k, 1)
         strains_total(2) = state(i_sdv_eps22) + strainInc(k, 2)
         strains_total(3) = state(i_sdv_eps33) + strainInc(k, 3)
         strains_total(4) = state(i_sdv_gamma12) + two*strainInc(k, 4)
         strains_total(5) = state(i_sdv_gamma23) + two*strainInc(k, 5)
         strains_total(6) = state(i_sdv_gamma13) + two*strainInc(k, 6)
         state(i_sdv_eps11:i_sdv_gamma13) = strains_total

         do i = 1, 6
            strainsToWrite(i) = strains_total(i)
         end do

         ! Get stresses from VEXTERNALDB via shared array
         ptr_stressesToRead = SMALocalFloatArrayAccess(1003)
         stresses(k, :) = stressesToRead(k, :)

         stressNew(k, :) = stresses
         stateNew(k, :) = state

      end do ! nBlock

      return
      end ! SUBROUTINE
