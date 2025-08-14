!==================================================================================================
 module mpas_atmphys_fdda
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars

 use diffusion_filter_coefficients
 
   !
   ! Abstract interface for routine used to communicate halos of fields
   ! in a named group
   !
 
 
 implicit none
 private
 public:: apply_fdda
 
 abstract interface
    subroutine halo_exchange_routine(domain, halo_group, ierr)
 
       use mpas_derived_types, only : domain_type
 
       type (domain_type), intent(inout) :: domain
       character(len=*), intent(in) :: halo_group
       integer, intent(out), optional :: ierr
 
    end subroutine halo_exchange_routine
 end interface

 integer, private:: i,k

!>\brief MPAS grid nudging.
!>\author O. Russell Bullock Jr. (send comments to bullock.russell@epa.gov).
!>\date 2016-03-30.
!>
!>\details
!>
!> subroutines in mpas_atmphys_fdda:
!> ----------------------------------------------
!> apply_fdda      : main fdda subroutine (called from subroutine physics_driver).
!>

!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>      Original sourcecode by O. Russell Bullock Jr. (bullock.russell@epa.gov) / 2016-03-30. 
!>
!>    * Added various controls for length-scaled FDDA application 
!>      O. Russell Bullock Jr. (bullock.russell@epa.gov) / 2016-08-11. 
!>    * Modified to only require reads of "fdda_new" data arrays.  "fdda_old" is initially 
!>      set to initial conditions.  At end of FDDA data interval, "fdda_old" is set to 
!>      "fdda_new" and updated "fdda_new" data comes from FDDA input stream.
!>      O. Russell Bullock Jr. (bullock.russell@epa.gov) / 2016-08-30.


 contains

!==================================================================================================
 subroutine apply_fdda(domain,exchange_halo_group,xtime_s,configs,mesh,state,time_lev,fdda,diag,diag_physics,tend_physics,its,ite)
!==================================================================================================
 

 !input variables:
 procedure (halo_exchange_routine) :: exchange_halo_group
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics
 real(kind=RKIND),intent(in):: xtime_s
 integer,intent(in):: time_lev
 integer,intent(in):: its,ite
 
 !inout variables:
 type(domain_type),intent(inout):: domain
 type(mpas_pool_type),intent(inout):: fdda 
 type(mpas_pool_type),intent(inout):: tend_physics
 
 !local pointers:
 integer,pointer:: index_qv,fdda_t_min_layer,fdda_q_min_layer,fdda_uv_min_layer,n_steps
 integer,dimension(:),pointer :: kpbl
 integer,dimension(:),pointer :: nEdgesOnCell
 integer,dimension(:,:),pointer :: edgesOnCell,cellsOnEdge
 real(kind=RKIND),pointer:: fdda_t_coef,fdda_q_coef,fdda_uv_coef,areacellmax
 real(kind=RKIND),pointer:: fdda_int,fdda_scale_min,fdda_scale_max,sigma_gaussianfilter,filter_scale,dx_min
 real(kind=RKIND),pointer:: config_dt
 real(kind=RKIND),dimension(:),pointer    :: dcEdge,xCell,yCell,zCell,invAreaCell,dvEdge,invDcEdge
 real(kind=RKIND),dimension(:),pointer    :: p_dfilter
 real(kind=RKIND),dimension(:,:),pointer  :: rthfddaten,rqvfddaten
 real(kind=RKIND),dimension(:,:),pointer  :: rufddaten,rvfddaten
 real(kind=RKIND),dimension(:,:),pointer  :: qv_fdda_old,th_fdda_old,u_fdda_old,v_fdda_old
 real(kind=RKIND),dimension(:,:),pointer  :: qv_fdda_new,th_fdda_new,u_fdda_new,v_fdda_new
 real(kind=RKIND),dimension(:,:),pointer  :: edgesOnCell_sign
 real(kind=RKIND),dimension(:,:),pointer  :: filtervar_prev,filtervar_curr,filtervar_next,laplacian_term,filtered
 real(kind=RKIND),dimension(:,:),pointer  :: u,v
 real(kind=RKIND),dimension(:,:),pointer  :: qv
 real(kind=RKIND),dimension(:,:),pointer  :: theta_m
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 logical,pointer:: fdda_t, fdda_t_in_pbl
 logical,pointer:: fdda_q, fdda_q_in_pbl
 logical,pointer:: fdda_uv,fdda_uv_in_pbl
 logical,pointer:: config_do_restart
 logical,pointer:: nudging_with_spatial_filter,gaussian_filter,diffusion_filter

 character (len=StrKIND), pointer:: fdda_scheme,filter_shape

!local variables:
 integer, pointer  :: nCells, nEdges, nVertLevels, maxEdges
 real(kind=RKIND)  :: tfrac, next_tfrac
 real(kind=RKIND)  :: t_target,q_target,u_target,v_target
 real(kind=RKIND)  :: rufddaten_temp,rvfddaten_temp,sum_weights,weights_insuccess,weights_success,distance_2cells,weight
 real(kind=RKIND)  :: sfac
 real(kind=RKIND)  :: theta_p
 real(kind=RKIND)  :: srange
 real(kind=RKIND)  :: edge_sign,inv_s_max,kappacell,kappacell1,kappacell2

 integer :: i,k,filter_iCell,diffusion_step,filter_iEdge,cell1,cell2

 logical :: first_call = .true.

!local parameters:
 real(kind=RKIND),parameter:: rvrd   = R_v/R_d

!==================================================================================================
!  write(0,*) '--- enter apply_fdda:  at xtime_s = ',xtime_s
 call mpas_log_write('--- apply_fdda at the time step = $r', realArgs=(/xtime_s/))

 call mpas_pool_get_config(configs,'config_fdda_scheme' ,fdda_scheme )

 call mpas_pool_get_config(configs,'config_do_restart' ,config_do_restart )
 call mpas_pool_get_config(configs,'config_dt' ,config_dt )
 call mpas_pool_get_config(configs,'config_fdda_int',fdda_int)
 call mpas_pool_get_config(configs,'config_fdda_t',fdda_t)
 call mpas_pool_get_config(configs,'config_fdda_q',fdda_q)
 call mpas_pool_get_config(configs,'config_fdda_uv',fdda_uv)
 call mpas_pool_get_config(configs,'config_fdda_t_in_pbl',fdda_t_in_pbl)
 call mpas_pool_get_config(configs,'config_fdda_q_in_pbl',fdda_q_in_pbl)
 call mpas_pool_get_config(configs,'config_fdda_uv_in_pbl',fdda_uv_in_pbl)
 call mpas_pool_get_config(configs,'config_fdda_t_min_layer',fdda_t_min_layer)
 call mpas_pool_get_config(configs,'config_fdda_q_min_layer',fdda_q_min_layer)
 call mpas_pool_get_config(configs,'config_fdda_uv_min_layer',fdda_uv_min_layer)
 call mpas_pool_get_config(configs,'config_fdda_t_coef',fdda_t_coef)
 call mpas_pool_get_config(configs,'config_fdda_q_coef',fdda_q_coef)
 call mpas_pool_get_config(configs,'config_fdda_uv_coef',fdda_uv_coef)
 call mpas_pool_get_config(configs,'config_fdda_scale_min',fdda_scale_min)
 call mpas_pool_get_config(configs,'config_fdda_scale_max',fdda_scale_max)
 call mpas_pool_get_config(configs,'config_nudging_with_spatial_filter',nudging_with_spatial_filter)
 call mpas_pool_get_config(configs,'config_gaussian_filter',gaussian_filter)
 call mpas_pool_get_config(configs,'config_diffusion_filter',diffusion_filter)
 call mpas_pool_get_config(configs,'sigma_gaussianfilter',sigma_gaussianfilter)
 call mpas_pool_get_config(configs,'config_filter_scale',filter_scale)
 call mpas_pool_get_config(configs,'config_dx_min',dx_min)
 call mpas_pool_get_config(configs,'config_n_steps',n_steps)
 call mpas_pool_get_config(configs,'config_filter_shape',filter_shape)
 call mpas_pool_get_config(configs,'config_areacellmax',areacellmax)

 call mpas_pool_get_array(fdda,'qv_fdda_old',qv_fdda_old)
 call mpas_pool_get_array(fdda,'th_fdda_old',th_fdda_old)
 call mpas_pool_get_array(fdda,'u_fdda_old', u_fdda_old )
 call mpas_pool_get_array(fdda,'v_fdda_old', v_fdda_old )
 call mpas_pool_get_array(fdda,'qv_fdda_new',qv_fdda_new)
 call mpas_pool_get_array(fdda,'th_fdda_new',th_fdda_new)
 call mpas_pool_get_array(fdda,'u_fdda_new', u_fdda_new )
 call mpas_pool_get_array(fdda,'v_fdda_new', v_fdda_new )
 call mpas_pool_get_array(fdda,'filtervar_prev', filtervar_prev )
 call mpas_pool_get_array(fdda,'filtervar_curr', filtervar_curr )
 call mpas_pool_get_array(fdda,'filtervar_next', filtervar_next )
 call mpas_pool_get_array(fdda,'laplacian_term', laplacian_term )
 call mpas_pool_get_array(fdda,'filtered', filtered )
 call mpas_pool_get_array(tend_physics,'rthfddaten',rthfddaten)
 call mpas_pool_get_array(tend_physics,'rqvfddaten',rqvfddaten)
 call mpas_pool_get_array(tend_physics,'rufddaten' ,rufddaten )
 call mpas_pool_get_array(tend_physics,'rvfddaten' ,rvfddaten )
 call mpas_pool_get_array(mesh,'dcEdge',dcEdge)
 call mpas_pool_get_array(mesh,'xCell',xCell)
 call mpas_pool_get_array(mesh,'yCell',yCell)
 call mpas_pool_get_array(mesh,'zCell',zCell)
 call mpas_pool_get_array(mesh,'nEdgesOnCell',nEdgesOnCell)
 call mpas_pool_get_array(mesh,'edgesOnCell',edgesOnCell)
 call mpas_pool_get_array(mesh,'invAreaCell',invAreaCell)
 call mpas_pool_get_array(mesh,'edgesOnCell_sign',edgesOnCell_sign)
 call mpas_pool_get_array(mesh,'dvEdge',dvEdge)
 call mpas_pool_get_array(mesh,'invDcEdge',invDcEdge)
 call mpas_pool_get_array(mesh,'cellsOnEdge',cellsOnEdge)
 call mpas_pool_get_array(diag_physics,'kpbl',kpbl)
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u)
 call mpas_pool_get_array(diag,'uReconstructMeridional',v)
 
 call mpas_pool_get_dimension(state,'index_qv',index_qv)

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
 call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
 call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)

 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)

!check to see if 'old' FDDA target values need to be defined from initial conditions
 if(first_call .and. .not. config_do_restart) then
    qv_fdda_old = qv
    th_fdda_old = theta_m / (1._RKIND + rvrd * qv)
    !  th_fdda_old = theta_m
    u_fdda_old  = u
    v_fdda_old  = v
    first_call  = .false.
 endif


 !calculate time fraction within FDDA target data interval
 tfrac = mod(xtime_s,fdda_int)/fdda_int
 next_tfrac = mod(xtime_s+config_dt,fdda_int)/fdda_int
 call mpas_log_write('----- check in tfrac = $r', realArgs=(/tfrac/))
 call mpas_log_write('----- check in next tfrac = $r', realArgs=(/next_tfrac/))
 fdda_select: select case(fdda_scheme)

    case ("scaled")

       srange = fdda_scale_max-fdda_scale_min

       DO i=its,ite

          if (dcEdge(i) .lt. fdda_scale_max) then
             if (dcEdge(i) .lt. fdda_scale_min) then
                sfac = 0.0
             else
                sfac = (dcEdge(i)-fdda_scale_min)/srange
             endif
          else
             sfac = 1.0
          endif

       DO k=kts,kte
         !  theta_p = theta_m(k,i) / (1._RKIND + rvrd * qv(k,i))
         ! Theta_m is already the moisture potential temperature, don't add rvrd to theta_p or it will cause instability
         theta_p = theta_m(k,i)
         
         if (fdda_t .and. k.gt.fdda_t_min_layer .and. (fdda_t_in_pbl .or. k.gt.kpbl(i))) then
            t_target = (1.0-tfrac)*th_fdda_old(k,i)+tfrac*th_fdda_new(k,i)
            rthfddaten(k,i) = fdda_t_coef * ( t_target - theta_p ) * sfac
         else
             rthfddaten(k,i) = 0.
            endif
            
            if (fdda_q .and. k.gt.fdda_q_min_layer .and. (fdda_q_in_pbl .or. k.gt.kpbl(i))) then
               q_target = (1.0-tfrac)*qv_fdda_old(k,i)+tfrac*qv_fdda_new(k,i)
               rqvfddaten(k,i) = fdda_q_coef * ( q_target - qv(k,i) ) * sfac
            else
               rqvfddaten(k,i) = 0.
            endif
            
            if (fdda_uv .and. k.gt.fdda_uv_min_layer .and. (fdda_uv_in_pbl .or. k.gt.kpbl(i))) then
               u_target = (1.0-tfrac)*u_fdda_old(k,i)+tfrac*u_fdda_new(k,i)
               rufddaten(k,i) = fdda_uv_coef * ( u_target - u(k,i) ) * sfac
               v_target = (1.0-tfrac)*v_fdda_old(k,i)+tfrac*v_fdda_new(k,i)
               rvfddaten(k,i) = fdda_uv_coef * ( v_target - v(k,i) ) * sfac
            else
               rufddaten(k,i) = 0.
               rvfddaten(k,i) = 0.
            endif
            
         ENDDO
       ENDDO

      case("analysis")
         
         DO i=its,ite
            DO k=kts,kte
               
          ! Theta_m is already the moisture potential temperature, don't add rvrd to theta_p or it will cause instability
          !  theta_p = theta_m(k,i) / (1._RKIND + rvrd * qv(k,i))
          theta_p = theta_m(k,i)

          if (fdda_t .and. k.gt.fdda_t_min_layer .and. (fdda_t_in_pbl .or. k.gt.kpbl(i))) then
             t_target = (1.0-tfrac)*th_fdda_old(k,i)+tfrac*th_fdda_new(k,i)
             rthfddaten(k,i) = fdda_t_coef * ( t_target - theta_p )
          else
             rthfddaten(k,i) = 0.
          endif

          if (fdda_q .and. k.gt.fdda_q_min_layer .and. (fdda_q_in_pbl .or. k.gt.kpbl(i))) then
             q_target = (1.0-tfrac)*qv_fdda_old(k,i)+tfrac*qv_fdda_new(k,i)
             rqvfddaten(k,i) = fdda_q_coef * ( q_target - qv(k,i) )
          else
             rqvfddaten(k,i) = 0.
          endif

          if (fdda_uv .and. k.gt.fdda_uv_min_layer .and. (fdda_uv_in_pbl .or. k.gt.kpbl(i))) then
             u_target = (1.0-tfrac)*u_fdda_old(k,i)+tfrac*u_fdda_new(k,i)
             rufddaten(k,i) = fdda_uv_coef * ( u_target - u(k,i) )
             v_target = (1.0-tfrac)*v_fdda_old(k,i)+tfrac*v_fdda_new(k,i)
             rvfddaten(k,i) = fdda_uv_coef * ( v_target - v(k,i) )
          else
             rufddaten(k,i) = 0.
             rvfddaten(k,i) = 0.
          endif

            ENDDO
         ENDDO

    case default

 end select fdda_select

if (nudging_with_spatial_filter) then
   ! spatial Filter wih low pass gaussian filter based on distance added by yycheng 2025.02

   if (gaussian_filter) then
      ! gaussian filter based on distance over the block, for test the following diffusion-based type
      call mpas_log_write('--- apply gaussian filter for u & v')
      weights_insuccess = 0.
      weights_success   = 0.

      DO k=kts,kte
         DO i=its,ite
            
            ! apply over rufddaten rvfddaten
            sum_weights    = 0.
            rufddaten_temp = 0.
            rvfddaten_temp = 0.
            DO filter_iCell=its,ite
               ! filter weight based on distance
               distance_2cells = sqrt( (xCell(i) - xCell(filter_iCell))**2 + &
               (yCell(i) - yCell(filter_iCell))**2 + &
               (zCell(i) - zCell(filter_iCell))**2 )
               weight = exp(-0.5 * distance_2cells ** 2 / sigma_gaussianfilter **2)
               rufddaten_temp = rufddaten_temp + rufddaten(k,filter_iCell) * weight
               rvfddaten_temp = rvfddaten_temp + rvfddaten(k,filter_iCell) * weight
               sum_weights = sum_weights + weight
               sum_weights = sum_weights + 1
               ! call mpas_log_write('> loop filter cells')
            ENDDO
            if (sum_weights > 0.0) then
               rufddaten(k,i) =  rufddaten_temp / sum_weights
               rvfddaten(k,i)  = rvfddaten_temp / sum_weights
               weights_success = weights_success + 1.
            else
               weights_insuccess = weights_insuccess + 1.
            endif

         ENDDO
      ENDDO

   endif
   
   
   if (diffusion_filter) then
      ! diffusion-based filter from Grooms et al. 2021
      call mpas_log_write('> start diffusion-based filter')
      
      ! The anisotropic coefficient, defined on the basis of Cell area
      ! takes into account the square approximation of ratio and resolution to achieve similar filter scale at different locations
      ! dx_min is defined in terms of the largest grid spacing and considers dx_min * * 2 to be the largest area cell
      

      allocate(p_dfilter(n_steps+1))

      call compute_filter_coefficients(n_steps, filter_scale, dx_min, filter_shape, p_dfilter)
      
      if (fdda_uv) then
         !----------- apply diffusion-based filter : rufddaten -----------!
         filtered(:,:)       = 0._RKIND
         filtervar_prev(:,:) = 0._RKIND
         filtervar_curr(:,:) = 0._RKIND
         filtervar_next(:,:) = 0._RKIND
         laplacian_term(:,:) = 0._RKIND
         call mpas_log_write('Max rufddaten before filter : $r', realArgs=(/maxval(rufddaten)/))
         call apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        rufddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q)                           
                        
         !----------- apply diffusion-based filter : rvfddaten -----------!
         filtered(:,:)       = 0._RKIND
         filtervar_prev(:,:) = 0._RKIND
         filtervar_curr(:,:) = 0._RKIND
         filtervar_next(:,:) = 0._RKIND
         laplacian_term(:,:) = 0._RKIND
         call mpas_log_write('Max rvfddaten before filter : $r', realArgs=(/maxval(rvfddaten)/))
         call apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        rvfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q)                           

      endif
                     
      if (fdda_t) then
         !----------- apply diffusion-based filter : rvfddaten -----------!
         filtered(:,:)       = 0._RKIND
         filtervar_prev(:,:) = 0._RKIND
         filtervar_curr(:,:) = 0._RKIND
         filtervar_next(:,:) = 0._RKIND
         laplacian_term(:,:) = 0._RKIND
         call mpas_log_write('Max rthfddaten before filter : $r', realArgs=(/maxval(rthfddaten)/))
         call apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        rthfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q)                           
      endif
      
      if (fdda_q) then
         !----------- apply diffusion-based filter : rqvfddaten -----------!
         filtered(:,:)       = 0._RKIND
         filtervar_prev(:,:) = 0._RKIND
         filtervar_curr(:,:) = 0._RKIND
         filtervar_next(:,:) = 0._RKIND
         laplacian_term(:,:) = 0._RKIND
         call mpas_log_write('Max rqvfddaten before filter : $r', realArgs=(/maxval(rqvfddaten)/))
         call apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        rqvfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q)                           
      endif

      deallocate(p_dfilter)
   endif
endif

!check to see if this is the last time step within the current FDDA target data interval and 
!if so, save the current "new" targets as the "old" targets for the next time step
 if(mod(xtime_s+config_dt,fdda_int)/fdda_int .lt. tfrac) then
    th_fdda_old = th_fdda_new
    qv_fdda_old = qv_fdda_new
    u_fdda_old  = u_fdda_new
    v_fdda_old  = v_fdda_new
   !  write(0,*) '--- exit apply_fdda'
    call mpas_log_write('***** apply_fdda at the last time step, get new targets next!')
    
    !  deallocate(th_fdda_new)
    !  deallocate(qv_fdda_new)
    !  deallocate(u_fdda_new)
    !  deallocate(v_fdda_new)
   endif
   
   ! write(0,*) '--- exit apply_fdda'
   call mpas_log_write('--- exit apply_fdda')
   
 end subroutine apply_fdda



subroutine apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        varfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q)         

!=======================================================================
!  This subroutine applies diffusion-based spatial filtering 
!  beyond grid scale (supergrid smoothing).
!  Developed by Y. Y. Cheng (yycheng), 2025.06
!  References:
!    - Grooms et al. (2021), Diffusion-Based Smoothers for Spatial Filtering
!    - Loose et al. (2022), GCM-Filters: A Python Package for Diffusion-based Filtering
!    - Bullock Jr. et al. (2018), Analysis Nudging in MPAS-A (v4.0)
!=======================================================================

  implicit none

  ! Subroutine Inputs
  procedure (halo_exchange_routine) :: exchange_halo_group
  integer, intent(in) :: its, ite, kts, kte, n_steps
  real(RKIND), intent(in) :: dx_min
  real(RKIND), intent(in) :: p_dfilter(:)
  real(RKIND), intent(in) :: areacellmax
  real(RKIND), intent(in) :: invAreaCell(:)
  integer, intent(in) :: nEdgesOnCell(:)
  integer, intent(in) :: edgesOnCell(:,:)
  real, intent(in) :: edgesOnCell_sign(:,:)
  real(RKIND), intent(in) :: dvEdge(:), invDcEdge(:)
  integer, intent(in) :: cellsOnEdge(:,:)
  logical, intent(in) :: fdda_uv, fdda_t, fdda_q
  
  ! Input/output
  type(domain_type),intent(inout):: domain
  real(RKIND), intent(inout) :: varfddaten(:,:),filtered(:,:), filtervar_prev(:,:), filtervar_curr(:,:), filtervar_next(:,:), laplacian_term(:,:)
  
  ! Local variables
  integer :: i, k, filter_iCell, filter_iEdge, cell1, cell2, diffusion_step
  real(RKIND) :: inv_s_max, kappacell, kappacell1, kappacell2, edge_sign

  ! Initialization
  filtered(:,:) = 0._RKIND
  filtervar_prev(:,:) = 0._RKIND
  filtervar_curr(:,:) = 0._RKIND
  filtervar_next(:,:) = 0._RKIND
  laplacian_term(:,:) = 0._RKIND

  inv_s_max = (dx_min / 2.0_RKIND)**2.0 / 2.0_RKIND

  ! Store initial varfddaten into filtervar_prev
  do i = its, ite
    do k = kts, kte
      filtervar_prev(k,i) = varfddaten(k,i)
    end do
  end do

  ! Compute the 2nd-order laplacian filter: filtervar_prev
  laplacian_term(:,:) = 0._RKIND
  do filter_iCell = its, ite
    kappacell = 1.0_RKIND / areacellmax / invAreaCell(filter_iCell)
    do i = 1, nEdgesOnCell(filter_iCell)
      filter_iEdge = edgesOnCell(i, filter_iCell)
      edge_sign = invAreaCell(filter_iCell) * edgesOnCell_sign(i, filter_iCell) * dvEdge(filter_iEdge) * invDcEdge(filter_iEdge) * kappacell
      cell1 = cellsOnEdge(1, filter_iEdge)
      cell2 = cellsOnEdge(2, filter_iEdge)
      kappacell1 = 1.0_RKIND / areacellmax / invAreaCell(cell1)
      kappacell2 = 1.0_RKIND / areacellmax / invAreaCell(cell2)
      do k = kts, kte
        laplacian_term(k, filter_iCell) = laplacian_term(k, filter_iCell) + edge_sign * (filtervar_prev(k, cell2) - filtervar_prev(k, cell1)) + &
                                         (filtervar_prev(k, cell2) - filtervar_prev(k, cell1)) * invDcEdge(filter_iEdge) * &
                                         (kappacell2 - kappacell1) * invDcEdge(filter_iEdge)
      end do
    end do
  end do
  
  ! Update filtervar_curr and filtered
  do i = its, ite
   do k = kts, kte
      filtervar_curr(k, i) = (-2.0_RKIND * laplacian_term(k, i) * inv_s_max - filtervar_prev(k, i))
   end do
end do

do i = its, ite
   do k = kts, kte
      filtered(k, i) = p_dfilter(1) * filtervar_prev(k, i) + p_dfilter(2) * filtervar_curr(k, i)
   end do
end do

! Diffusion steps
do diffusion_step = 3, n_steps + 1
   
   laplacian_term(:,:) = 0._RKIND

   do filter_iCell = its, ite
      kappacell = 1.0_RKIND / areacellmax / invAreaCell(filter_iCell)
      do i = 1, nEdgesOnCell(filter_iCell)
         filter_iEdge = edgesOnCell(i, filter_iCell)
         edge_sign = invAreaCell(filter_iCell) * edgesOnCell_sign(i, filter_iCell) * dvEdge(filter_iEdge) * invDcEdge(filter_iEdge) * kappacell
         cell1 = cellsOnEdge(1, filter_iEdge)
         cell2 = cellsOnEdge(2, filter_iEdge)
         kappacell1 = 1.0_RKIND / areacellmax / invAreaCell(cell1)
         kappacell2 = 1.0_RKIND / areacellmax / invAreaCell(cell2)
         do k = kts, kte
            laplacian_term(k, filter_iCell) = laplacian_term(k, filter_iCell) + edge_sign * (filtervar_curr(k, cell2) - filtervar_curr(k, cell1)) + &
            (filtervar_curr(k, cell2) - filtervar_curr(k, cell1)) * invDcEdge(filter_iEdge) * &
            (kappacell2 - kappacell1) * invDcEdge(filter_iEdge)
         end do
      end do
   end do



   do i = its, ite
      do k = kts, kte
         filtervar_next(k, i) = 2.0_RKIND * (-2.0_RKIND * laplacian_term(k, i) * inv_s_max - filtervar_curr(k, i)) - filtervar_prev(k, i)
         filtered(k, i) = filtered(k, i) + p_dfilter(diffusion_step) * filtervar_next(k, i)
         filtervar_prev(k, i) = filtervar_curr(k, i)
         filtervar_curr(k, i) = filtervar_next(k, i)
      end do
   end do   

     ! 预先存放在fddaten中，用于后续直接进行halo_group的交换
     ! Store final results back into varfddaten
     do i = its, ite
       do k = kts, kte
         varfddaten(k, i) = filtered(k, i)
       end do
     end do


   !$OMP BARRIER
      if(fdda_uv) then
         call exchange_halo_group(domain, 'physics:uvfddaten')
      endif
   !$OMP BARRIER
   
   !$OMP BARRIER
      if(fdda_t) then
         call exchange_halo_group(domain, 'physics:thfddaten')
      endif
   !$OMP BARRIER
   
   !$OMP BARRIER
      if(fdda_q) then
         call exchange_halo_group(domain, 'physics:qvfddaten')
      endif
   !$OMP BARRIER

end do

  ! Deallocate arrays
!   deallocate(filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term)

end subroutine apply_diffusion_based_filter



!==================================================================================================
 end module mpas_atmphys_fdda
!==================================================================================================
