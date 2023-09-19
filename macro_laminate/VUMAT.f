! Author: Ishaan Desai (desaii)
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

! - Variables --------------------------------------------------------- [11, 22, 33, 12, 23, 13]

      ! For defining shared arrays
      ! maxMaterialPts == nblock, directComponents == ndir, indirectComponents == nshr
      ! maxTensorComponents == ndir + nshr
      parameter(maxMaterialPts=1000,
     * directComponents=3,
     * indirectComponents=3,
     * maxTensorComponents=6,
     * max3Dsize=3000,
     * max6Dsize=6000)

      real*8 :: state(nstatev), strains_total(6)
      real*8 :: E11, E220, E33, nu12, nu23, nu13, G13, G23, nu31, nu21
      integer :: d, ArraySize2D

      double precision, dimension(max3Dsize) :: coordsToShare
      pointer(ptr_coordsToShare, coordsToShare)

      double precision, dimension(max6Dsize) :: strainsToWrite, stressesToRead

      pointer(ptr_strainsToWrite, strainsToWrite)
      pointer(ptr_stressesToRead, stressesToRead)

      integer, dimension(3) :: intsToShare
      pointer(ptr_intsToShare, intsToShare)

      ! Create shared array to share strains with VEXTERNALDB
      ptr_strainsToWrite = SMALocalFloatArrayCreate(1002,
     * max6Dsize, 0.0)

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
         ptr_intsToShare = SMALocalIntArrayCreate(1000, 3, 0)
         intsToShare(1) = nblock
         intsToShare(2) = ndir
         intsToShare(3) = nshr
      
         ! Create shared array for material point coordinates
         ptr_coordsToShare = SMALocalFloatArrayCreate(1001,
     *    max3Dsize, 0.0)      

         ! Put coordinates in the shared array
         do k = 1, nblock
            do d = 1, directComponents
               coordsToShare(k, d) = coordMp(k, d)
            end do
         end do

      end if

      ! Loop through material points to collect strains
      do k = 1, nblock

         state = stateOld(k, :)

         strains_total(1) = state(i_sdv_eps11) + strainInc(k, 1)
         strains_total(2) = state(i_sdv_eps22) + strainInc(k, 2)
         strains_total(3) = state(i_sdv_eps33) + strainInc(k, 3)
         strains_total(4) = state(i_sdv_gamma12) + two*strainInc(k, 4)
         strains_total(5) = state(i_sdv_gamma23) + two*strainInc(k, 5)
         strains_total(6) = state(i_sdv_gamma13) + two*strainInc(k, 6)
         state(i_sdv_eps11:i_sdv_gamma13) = strains_total

         do i = 1, maxTensorComponents
            strainsToWrite(k, i) = strains_total(i)
         end do ! maxTensorComponents

      end do ! nblock

      ! Get stresses from VEXTERNALDB via shared array
      ptr_stressesToRead = SMALocalFloatArrayAccess(1003)

      ! Loop through material points to apply stresses
      do k = 1, nblock
         do d = 1, maxTensorComponents

            stressNew(k, d) = stressesToRead(k, d)
            stateNew(k, d) = state(d)

         end do ! maxTensorComponents
      end do ! nblock

      return
      end ! SUBROUTINE
