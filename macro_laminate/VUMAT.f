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

      real*8, dimension(nblock, nstatev) :: state
      real*8, dimension(nblock, ndir + nshr) :: strains, stresses
      integer :: d, counter

      ! preCICE variables
      integer :: rank, size, ongoing, dimensions, bool
      double precision :: preCICE_dt
      double precision, dimension(nblock*ndir) :: strains1to3, strains4to6,
     * stresses1to3, stresses4to6, couplingVertices
      integer, dimension(nblock) :: vertexIDs

      write(*,*) "t = ", totalTime, ",", " dt = ", dt, " VUMAT: Entering VUMAT."

      if (dt /= 1) then
         if (totalTime < 2.*dt) then ! only the first step

            ! Get MPI rank and size (total number of MPI processors in this job)
            call vgetnumcpus(size)
            call vgetrank(rank)

            ! Create preCICE participant
            call precicef_create("Laminate-3D-ply",
     *       "/home/desaii/composite-multiscale/precice-config.xml",
     *       rank, size)

            ! Get problem dimensions from preCICE
            call precicef_get_mesh_dimensions("laminate-macro-mesh", dimensions)

            counter = 1
            do k = 1, nblock
               do d = 1, ndir
                  couplingVertices(counter) = coordMp(k, d)
                  counter = counter + 1
               end do
            end do

            ! Set coupling mesh vertices in preCICE
            call precicef_set_vertices("laminate-macro-mesh",
     *       nblock, couplingVertices, vertexIDs)

            call get_strains(nblock, ndir, nshr, nstatev,
     *       strainInc, stateOld, state, strains)

            ! Apply the calculated stresses
            do k = 1, nblock
               stateNew(k, :) = state(k, :)
            end do

            counter = 1
            do k = 1, nblock

               do d = 1, ndir
                  strains1to3(counter) = strains(k, d)
                  strains4to6(counter) = strains(k, d+ndir)
                  counter = counter + 1

               end do ! ndir

               do d = 1, nstatev

                  stateNew(k, d) = state(k, d)
               
               end do ! nstatev

            end do ! nblock

            call precicef_requires_initial_data(bool)

            if (bool == 1) then
               call precicef_write_data("laminate-macro-mesh",
     *          "strains1to3", nblock, vertexIDs, strains1to3)
               call precicef_write_data("laminate-macro-mesh",
     *          "strains4to6", nblock, vertexIDs, strains4to6)
            end if

            write(*,*) "VUMAT: Initial data written to preCICE."

            ! Initialize preCICE
            call precicef_initialize()
            call precicef_get_max_time_step_size(preCICE_dt)

            write(*,*) "VUMAT: preCICE initialization complete."

         end if ! if (totalTime < 2.*dt)

         ! When not in initialization, vertexIDs of preCICE are defined
         ! manually, because VUMAT cannot hold global data.
         ! TODO: Need to find a way to store global data in VUMAT
         do k = 0, nblock - 1
            vertexIDs(k + 1) = k
         end do

         ! Check if coupling is still going on
         call precicef_is_coupling_ongoing(bool)
         !call assert(bool.eq.1)

! Get strains and write them to preCICE ===============================
         call get_strains(nblock, ndir, nshr, nstatev,
     *    strainInc, stateOld, state, strains)
 
         counter = 1
         do k = 1, nblock

            do d = 1, ndir

               strains1to3(counter) = strains(k, d)
               strains4to6(counter) = strains(k, d+ndir)
               counter = counter + 1

            end do ! ndir

            do d = 1, nstatev
               stateNew(k, d) = state(k, d)
            end do ! nstatev

         end do ! nblock

         call precicef_write_data("laminate-macro-mesh",
     *       "strains1to3", nblock, vertexIDs, strains1to3)
         call precicef_write_data("laminate-macro-mesh",
     *       "strains4to6", nblock, vertexIDs, strains4to6)

         write(*,*) "(t = ", totalTime, ") VUMAT: Strains written to preCICE."

! ==========================================================================

         call precicef_advance(dt)
         write(*,*) "(t = ", totalTime, ") VUMAT: Coupling has been advanced."

! Read stresses from preCICE and apply them ===============================

         if (totalTime > dt) then ! from the second step onward

            call precicef_read_data("laminate-macro-mesh",
     *       "stresses1to3", nblock, vertexIDs, dt, stresses1to3)
            call precicef_read_data("laminate-macro-mesh",
     *       "stresses4to6", nblock, vertexIDs, dt, stresses4to6)

            ! Loop through material points to apply stresses
            counter = 1
            do k = 1, nblock
               do d = 1, ndir
                  stressNew(k, d) = stresses1to3(counter)
                  stressNew(k, d+ndir) = stresses4to6(counter)
                  counter = counter + 1

               end do ! ndir
            end do ! nblock

            write(*,*) "(t = ", totalTime, ") VUMAT: Stresses applied."

         else ! only in the first step
            call solve_material_model(nblock, ndir, nshr,
     *       nprops, props, strains, stresses) 

            ! Apply the calculated stresses
            do k = 1, nblock
               stressNew(k, :) = stresses(k, :)
            end do

            write(*,*) "(t = ", totalTime, ") VUMAT: stresses not read from preCICE, but instead the material model solved."

         end if
! ==========================================================================

         if (totalTime == 1.0) then ! Hardcode the end time because we cannot access it
            call precicef_finalize()
         end if

         write(*,*) "(t = ", totalTime, ") VUMAT: run complete."

      else ! if (dt == 1)

         call get_strains(nblock, ndir, nshr, nstatev,
     *    strainInc, stateOld, state, strains)

         ! Apply the calculated stresses
         do k = 1, nblock
            stateNew(k, :) = state(k, :)
         end do

         call solve_material_model(nblock, ndir, nshr,
     *    nprops, props, strains, stresses)

         ! Apply the calculated stresses
         do k = 1, nblock
            stressNew(k, :) = stresses(k, :)
         end do

         write(*,*) "(t = ", totalTime, ") VUMAT: Material model solved at dt = 1"

      end if ! if (dt /= 1)

      return
      end ! Subroutine vumat



      subroutine get_strains(
     *   nblock, ndir, nshr, nstatev,
     *   strainInc, stateOld, state, strains)
         ! Orientation of strains: [11, 22, 33, 12, 23, 13]

         integer, intent(in) :: nblock, ndir, nshr, nstatev
         real*8, dimension(nblock, nstatev), intent(in) :: stateOld
         real*8, dimension(nblock, ndir + nshr), intent(in) :: strainInc

         real*8, dimension(nblock, nstatev), intent(out) :: state
         real*8, dimension(nblock, ndir + nshr), intent(out) :: strains

         real*8, dimension(ndir + nshr) :: strain
         integer :: k, d

         parameter(i_sdv_eps11=1,
     *    i_sdv_eps22 = 2,
     *    i_sdv_eps33 = 3,
     *    i_sdv_gamma12 = 4,
     *    i_sdv_gamma23 = 5,
     *    i_sdv_gamma13 = 6,
     *    i_sdv_active = 7,
     *    i_sdv_t_f = 8)

         state(:, i_sdv_eps11) = zero
         state(:, i_sdv_eps22) = zero
         state(:, i_sdv_eps33) = zero
         state(:, i_sdv_gamma12) = zero
         state(:, i_sdv_gamma23) = zero
         state(:, i_sdv_gamma13) = zero
         state(:, i_sdv_active) = one
         state(:, i_sdv_t_f) = zero

         ! Loop through material points to collect strains
         do k = 1, nblock
            state(k, :) = stateOld(k, :)

            strain(1) = state(k, i_sdv_eps11) + strainInc(k, 1)
            strain(2) = state(k, i_sdv_eps22) + strainInc(k, 2)
            strain(3) = state(k, i_sdv_eps33) + strainInc(k, 3)
            strain(4) = state(k, i_sdv_gamma12) + two*strainInc(k, 4)
            strain(5) = state(k, i_sdv_gamma23) + two*strainInc(k, 5)
            strain(6) = state(k, i_sdv_gamma13) + two*strainInc(k, 6)
            state(k, i_sdv_eps11:i_sdv_gamma13) = strain

            do d = 1, ndir + nshr

               strains(k, d) = strain(d)

            end do ! ndir + nshr

         end do ! nblock

      return
      end ! subroutine get_strains



      subroutine solve_material_model(
     *   nblock, ndir, nshr, nprops,
     *   props, strains, stresses)
      ! Orientation of stresses: [11, 22, 33, 12, 23, 13]

      parameter(zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, three=3.d0)

      parameter(i_prp_E11=1,
     * i_prp_E22 = 2,
     * i_prp_G12 = 3,
     * i_prp_nu12 = 4,
     * i_prp_nu23 = 5)

      integer, intent(in) :: nblock, ndir, nshr, nprops
      real*8, dimension(nprops), intent(in) :: props
      real*8, dimension(nblock, ndir + nshr), intent(in) :: strains

      real*8, dimension(nblock, ndir + nshr), intent(out) :: stresses

      real*8, dimension(ndir + nshr, ndir + nshr) :: Q
      real*8 :: E11, E220, E33, nu12, nu23, nu13, G13, G23, nu31, nu21
      real*8 :: temp

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

      ! Calculate stresses using a material model
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

      do k = 1, nblock
         stresses(k, :) = MATMUL(Q, strains(k, :))
      end do

      return
      end ! Subroutine solve_material_model
