!==================================================================================================
 module mpas_atmphys_fdda
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars
 use mpas_geometry_utils, only: mpas_distance

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
 integer,pointer:: fdda_t_max_layer,fdda_q_max_layer,fdda_uv_max_layer
 integer,dimension(:),pointer :: kpbl
 integer,dimension(:),pointer :: nEdgesOnCell
 integer,dimension(:,:),pointer :: edgesOnCell,cellsOnEdge
 real(kind=RKIND),pointer:: fdda_t_coef,fdda_q_coef,fdda_uv_coef,areacellmax
 real(kind=RKIND),pointer:: fdda_int,fdda_scale_min,fdda_scale_max,sigma_gaussianfilter,filter_scale,dx_min
 real(kind=RKIND),pointer:: config_dt,sphere_radius
 real(kind=RKIND),pointer:: config_fdda_mask_lon0,config_fdda_mask_lat0
 real(kind=RKIND),pointer:: config_fdda_mask_radius,config_fdda_mask_transition
 real(kind=RKIND),pointer:: config_fdda_mask_radius_x,config_fdda_mask_radius_y
 real(kind=RKIND),pointer:: config_fdda_mask_rotation
 real(kind=RKIND),pointer:: config_fdda_mask_lon_min,config_fdda_mask_lon_max
 real(kind=RKIND),pointer:: config_fdda_mask_lat_min,config_fdda_mask_lat_max
 real(kind=RKIND),dimension(:),pointer    :: dcEdge,xCell,yCell,zCell,latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: invAreaCell,dvEdge,invDcEdge
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
 logical,pointer:: config_do_restart,on_a_sphere,config_fdda_region_mask
 logical,pointer:: nudging_with_spatial_filter,gaussian_filter,diffusion_filter
 logical,pointer:: diffusion_filter_anisotropic

 character (len=StrKIND), pointer:: fdda_scheme,filter_shape,config_fdda_mask_type

!local variables:
 integer, pointer  :: nCells, nEdges, nVertLevels, maxEdges
 real(kind=RKIND)  :: tfrac, next_tfrac
 real(kind=RKIND)  :: t_target,q_target,u_target,v_target
 real(kind=RKIND)  :: rufddaten_temp,rvfddaten_temp,sum_weights,weights_insuccess,weights_success,distance_2cells,weight
 real(kind=RKIND)  :: sfac,dist_center,pi_rkind,lon0_rad,lat0_rad,lon_diff
 real(kind=RKIND)  :: east_dist,north_dist,ellipse_radius_x,ellipse_radius_y
 real(kind=RKIND)  :: box_lon_width,box_lon_center,box_lat_min,box_lat_max
 real(kind=RKIND)  :: box_lat_center,box_half_width_x,box_half_width_y
 real(kind=RKIND)  :: theta_p
 real(kind=RKIND)  :: edge_sign,inv_s_max,kappacell,kappacell1,kappacell2
 real(kind=RKIND), dimension(3) :: center_xyz, cell_xyz
 real(kind=RKIND), dimension(:), allocatable :: fdda_region_mask

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
 call mpas_pool_get_config(configs,'config_fdda_t_max_layer',fdda_t_max_layer)
 call mpas_pool_get_config(configs,'config_fdda_q_max_layer',fdda_q_max_layer)
 call mpas_pool_get_config(configs,'config_fdda_uv_max_layer',fdda_uv_max_layer)
 call mpas_pool_get_config(configs,'config_fdda_t_coef',fdda_t_coef)
 call mpas_pool_get_config(configs,'config_fdda_q_coef',fdda_q_coef)
 call mpas_pool_get_config(configs,'config_fdda_uv_coef',fdda_uv_coef)
 call mpas_pool_get_config(configs,'config_fdda_scale_min',fdda_scale_min)
 call mpas_pool_get_config(configs,'config_fdda_scale_max',fdda_scale_max)
 call mpas_pool_get_config(configs,'config_nudging_with_spatial_filter',nudging_with_spatial_filter)
 call mpas_pool_get_config(configs,'config_gaussian_filter',gaussian_filter)
 call mpas_pool_get_config(configs,'config_diffusion_filter',diffusion_filter)
 call mpas_pool_get_config(configs,'config_fdda_diffusion_anisotropic',diffusion_filter_anisotropic)
 call mpas_pool_get_config(configs,'sigma_gaussianfilter',sigma_gaussianfilter)
 call mpas_pool_get_config(configs,'config_filter_scale',filter_scale)
 call mpas_pool_get_config(configs,'config_dx_min',dx_min)
 call mpas_pool_get_config(configs,'config_n_steps',n_steps)
 call mpas_pool_get_config(configs,'config_filter_shape',filter_shape)
 call mpas_pool_get_config(configs,'config_areacellmax',areacellmax)
 call mpas_pool_get_config(configs,'config_fdda_region_mask',config_fdda_region_mask)
 call mpas_pool_get_config(configs,'config_fdda_mask_type',config_fdda_mask_type)
 call mpas_pool_get_config(configs,'config_fdda_mask_lon0',config_fdda_mask_lon0)
 call mpas_pool_get_config(configs,'config_fdda_mask_lat0',config_fdda_mask_lat0)
 call mpas_pool_get_config(configs,'config_fdda_mask_radius',config_fdda_mask_radius)
 call mpas_pool_get_config(configs,'config_fdda_mask_radius_x',config_fdda_mask_radius_x)
 call mpas_pool_get_config(configs,'config_fdda_mask_radius_y',config_fdda_mask_radius_y)
 call mpas_pool_get_config(configs,'config_fdda_mask_rotation',config_fdda_mask_rotation)
 call mpas_pool_get_config(configs,'config_fdda_mask_lon_min',config_fdda_mask_lon_min)
 call mpas_pool_get_config(configs,'config_fdda_mask_lon_max',config_fdda_mask_lon_max)
 call mpas_pool_get_config(configs,'config_fdda_mask_lat_min',config_fdda_mask_lat_min)
 call mpas_pool_get_config(configs,'config_fdda_mask_lat_max',config_fdda_mask_lat_max)
 call mpas_pool_get_config(configs,'config_fdda_mask_transition',config_fdda_mask_transition)

 call mpas_pool_get_array(fdda,'qv_fdda_old',qv_fdda_old)
 call mpas_pool_get_array(fdda,'th_fdda_old',th_fdda_old)
 call mpas_pool_get_array(fdda,'u_fdda_old', u_fdda_old )
 call mpas_pool_get_array(fdda,'v_fdda_old', v_fdda_old )
 call mpas_pool_get_array(fdda,'qv_fdda_new',qv_fdda_new)
 call mpas_pool_get_array(fdda,'th_fdda_new',th_fdda_new)
 call mpas_pool_get_array(fdda,'u_fdda_new', u_fdda_new )
 call mpas_pool_get_array(fdda,'v_fdda_new', v_fdda_new )
 call mpas_pool_get_array(tend_physics,'filtervar_prev', filtervar_prev )
 call mpas_pool_get_array(tend_physics,'filtervar_curr', filtervar_curr )
 call mpas_pool_get_array(tend_physics,'filtervar_next', filtervar_next )
 call mpas_pool_get_array(tend_physics,'laplacian_term', laplacian_term )
 call mpas_pool_get_array(tend_physics,'filtered', filtered )
 call mpas_pool_get_array(tend_physics,'rthfddaten',rthfddaten)
 call mpas_pool_get_array(tend_physics,'rqvfddaten',rqvfddaten)
 call mpas_pool_get_array(tend_physics,'rufddaten' ,rufddaten )
 call mpas_pool_get_array(tend_physics,'rvfddaten' ,rvfddaten )
 call mpas_pool_get_array(mesh,'dcEdge',dcEdge)
 call mpas_pool_get_array(mesh,'xCell',xCell)
 call mpas_pool_get_array(mesh,'yCell',yCell)
 call mpas_pool_get_array(mesh,'zCell',zCell)
 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
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

 call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
 call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)

! build the horizontal FDDA mask, using 1.0 inside the selected region and 0.0 outside
 if (config_fdda_region_mask) then
   if (.not. allocated(fdda_region_mask)) allocate(fdda_region_mask(nCells))
      pi_rkind = acos(-1.0_RKIND)
      lon0_rad = config_fdda_mask_lon0 * pi_rkind / 180.0_RKIND
      lat0_rad = config_fdda_mask_lat0 * pi_rkind / 180.0_RKIND

      select case (trim(config_fdda_mask_type))
         case ("scaled")
            do i = 1, nCells
               fdda_region_mask(i) = resolution_taper_weight(dcEdge(i), fdda_scale_min, fdda_scale_max)
            enddo

         case ("circle")
            center_xyz(1) = sphere_radius * cos(lat0_rad) * cos(lon0_rad)
            center_xyz(2) = sphere_radius * cos(lat0_rad) * sin(lon0_rad)
            center_xyz(3) = sphere_radius * sin(lat0_rad)

            do i = 1, nCells
               cell_xyz(1) = xCell(i)
               cell_xyz(2) = yCell(i)
               cell_xyz(3) = zCell(i)
               dist_center = mpas_distance(cell_xyz, center_xyz, on_a_sphere)
               fdda_region_mask(i) = radial_taper_weight(dist_center, config_fdda_mask_radius, &
                                                       config_fdda_mask_transition)
            enddo

         case ("ellipse")
            ellipse_radius_x = config_fdda_mask_radius_x
            ellipse_radius_y = config_fdda_mask_radius_y
            if (ellipse_radius_x <= 0.0_RKIND) ellipse_radius_x = config_fdda_mask_radius
            if (ellipse_radius_y <= 0.0_RKIND) ellipse_radius_y = config_fdda_mask_radius

            do i = 1, nCells
               lon_diff = shortest_lon_diff_rad(lonCell(i), lon0_rad)
               east_dist = sphere_radius * cos(lat0_rad) * lon_diff
               north_dist = sphere_radius * (latCell(i) - lat0_rad)
               fdda_region_mask(i) = ellipse_taper_weight(east_dist, north_dist, &
                    ellipse_radius_x, ellipse_radius_y, config_fdda_mask_rotation, &
                    config_fdda_mask_transition)
            enddo

         case ("box")
            box_lon_width = config_fdda_mask_lon_max - config_fdda_mask_lon_min
            if (box_lon_width < 0.0_RKIND) box_lon_width = box_lon_width + 360.0_RKIND
            if (box_lon_width <= 0.0_RKIND) box_lon_width = 360.0_RKIND
            box_lon_center = wrap_lon_rad((config_fdda_mask_lon_min + 0.5_RKIND * box_lon_width) * &
                 pi_rkind / 180.0_RKIND)
            box_lon_width = box_lon_width * pi_rkind / 180.0_RKIND

            box_lat_min = min(config_fdda_mask_lat_min, config_fdda_mask_lat_max) * pi_rkind / 180.0_RKIND
            box_lat_max = max(config_fdda_mask_lat_min, config_fdda_mask_lat_max) * pi_rkind / 180.0_RKIND
            box_lat_center = 0.5_RKIND * (box_lat_min + box_lat_max)
            box_half_width_x = 0.5_RKIND * sphere_radius * abs(cos(box_lat_center)) * box_lon_width
            box_half_width_y = 0.5_RKIND * sphere_radius * max(0.0_RKIND, box_lat_max - box_lat_min)

            do i = 1, nCells
               lon_diff = shortest_lon_diff_rad(lonCell(i), box_lon_center)
               east_dist = sphere_radius * cos(box_lat_center) * lon_diff
               north_dist = sphere_radius * (latCell(i) - box_lat_center)
               fdda_region_mask(i) = box_taper_weight(east_dist, north_dist, &
                    box_half_width_x, box_half_width_y, config_fdda_mask_transition)
            enddo

         case default
            call mpas_log_write('WARNING: unknown config_fdda_mask_type; FDDA region mask is set to one.')
            fdda_region_mask(:) = 1.0_RKIND
      end select
 endif

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

       DO i=its,ite

          if (config_fdda_region_mask .and. trim(config_fdda_mask_type) == "scaled") then
             sfac = 1.0_RKIND
          else
             sfac = resolution_taper_weight(dcEdge(i), fdda_scale_min, fdda_scale_max)
          endif

       DO k=kts,kte
         !  theta_p = theta_m(k,i) / (1._RKIND + rvrd * qv(k,i))
         ! Theta_m is already the moisture potential temperature, don't add rvrd to theta_p or it will cause instability
         theta_p = theta_m(k,i)
         
         if (fdda_t .and. fdda_layer_enabled(k, fdda_t_min_layer, fdda_t_max_layer) .and. &
             (fdda_t_in_pbl .or. k.gt.kpbl(i))) then
            t_target = (1.0-tfrac)*th_fdda_old(k,i)+tfrac*th_fdda_new(k,i)
            rthfddaten(k,i) = fdda_t_coef * ( t_target - theta_p ) * sfac
         else
             rthfddaten(k,i) = 0.
            endif
            
            if (fdda_q .and. fdda_layer_enabled(k, fdda_q_min_layer, fdda_q_max_layer) .and. &
                (fdda_q_in_pbl .or. k.gt.kpbl(i))) then
               q_target = (1.0-tfrac)*qv_fdda_old(k,i)+tfrac*qv_fdda_new(k,i)
               rqvfddaten(k,i) = fdda_q_coef * ( q_target - qv(k,i) ) * sfac
            else
               rqvfddaten(k,i) = 0.
            endif
            
            if (fdda_uv .and. fdda_layer_enabled(k, fdda_uv_min_layer, fdda_uv_max_layer) .and. &
                (fdda_uv_in_pbl .or. k.gt.kpbl(i))) then
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

          if (fdda_t .and. fdda_layer_enabled(k, fdda_t_min_layer, fdda_t_max_layer) .and. &
              (fdda_t_in_pbl .or. k.gt.kpbl(i))) then
             t_target = (1.0-tfrac)*th_fdda_old(k,i)+tfrac*th_fdda_new(k,i)
             rthfddaten(k,i) = fdda_t_coef * ( t_target - theta_p )
          else
             rthfddaten(k,i) = 0.
          endif

          if (fdda_q .and. fdda_layer_enabled(k, fdda_q_min_layer, fdda_q_max_layer) .and. &
              (fdda_q_in_pbl .or. k.gt.kpbl(i))) then
             q_target = (1.0-tfrac)*qv_fdda_old(k,i)+tfrac*qv_fdda_new(k,i)
             rqvfddaten(k,i) = fdda_q_coef * ( q_target - qv(k,i) )
          else
             rqvfddaten(k,i) = 0.
          endif

          if (fdda_uv .and. fdda_layer_enabled(k, fdda_uv_min_layer, fdda_uv_max_layer) .and. &
              (fdda_uv_in_pbl .or. k.gt.kpbl(i))) then
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

      ! call mpas_log_write('> check block filter application rate = $r', realArgs=(/weights_success / (weights_insuccess + weights_success)/))
      ! call mpas_log_write('> check block filter application success number = $r', realArgs=(/weights_success/))
      ! call mpas_log_write('> check block filter application insuccess number = $r', realArgs=(/weights_insuccess/))
   endif
   
   
   if (diffusion_filter) then
      ! diffusion-based filter from Grooms et al. 2021
      call mpas_log_write('> start diffusion-based filter')
      
      ! The default Laplacian is isotropic. When requested, the anisotropic coefficient
      ! is scaled by Cell area to keep a similar filter scale across variable-resolution meshes.
      

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
                        rufddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q, &
                        diffusion_filter_anisotropic)
                        
         !----------- apply diffusion-based filter : rvfddaten -----------!
         filtered(:,:)       = 0._RKIND
         filtervar_prev(:,:) = 0._RKIND
         filtervar_curr(:,:) = 0._RKIND
         filtervar_next(:,:) = 0._RKIND
         laplacian_term(:,:) = 0._RKIND
         call mpas_log_write('Max rvfddaten before filter : $r', realArgs=(/maxval(rvfddaten)/))
         call apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        rvfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q, &
                        diffusion_filter_anisotropic)

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
                        rthfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q, &
                        diffusion_filter_anisotropic)
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
                        rqvfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q, &
                        diffusion_filter_anisotropic)
      endif

      deallocate(p_dfilter)
   endif
endif

 ! re-apply region mask after spatial filtering so outside-mask tendency remains zero
if (config_fdda_region_mask) then
   if (fdda_uv) then
      call apply_region_mask_to_tendency(rufddaten, fdda_region_mask, its, ite, kts, kte)
      call apply_region_mask_to_tendency(rvfddaten, fdda_region_mask, its, ite, kts, kte)
   endif
   if (fdda_t) then
      call apply_region_mask_to_tendency(rthfddaten, fdda_region_mask, its, ite, kts, kte)
   endif
   if (fdda_q) then
      call apply_region_mask_to_tendency(rqvfddaten, fdda_region_mask, its, ite, kts, kte)
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



pure logical function fdda_layer_enabled(k, min_layer, max_layer)
   implicit none
   integer, intent(in) :: k, min_layer, max_layer

   fdda_layer_enabled = k > min_layer
   if (max_layer > 0) fdda_layer_enabled = fdda_layer_enabled .and. k <= max_layer
end function fdda_layer_enabled


pure real(kind=RKIND) function resolution_taper_weight(dx, scale_min, scale_max)
   implicit none
   real(kind=RKIND), intent(in) :: dx, scale_min, scale_max

   if (scale_max <= scale_min) then
      if (dx >= scale_max) then
         resolution_taper_weight = 1.0_RKIND
      else
         resolution_taper_weight = 0.0_RKIND
      endif
   else if (dx <= scale_min) then
      resolution_taper_weight = 0.0_RKIND
   else if (dx >= scale_max) then
      resolution_taper_weight = 1.0_RKIND
   else
      resolution_taper_weight = (dx - scale_min) / (scale_max - scale_min)
   endif
end function resolution_taper_weight


pure real(kind=RKIND) function radial_taper_weight(dist, radius_outer, transition_width)
   implicit none
   real(kind=RKIND), intent(in) :: dist, radius_outer, transition_width
   real(kind=RKIND) :: radius_inner, pi_rkind

   pi_rkind = acos(-1.0_RKIND)

   if (transition_width <= 0.0_RKIND) then
      if (dist <= radius_outer) then
         radial_taper_weight = 1.0_RKIND
      else
         radial_taper_weight = 0.0_RKIND
      endif
      return
   endif

   radius_inner = max(0.0_RKIND, radius_outer - transition_width)

   if (dist <= radius_inner) then
      radial_taper_weight = 1.0_RKIND
   else if (dist >= radius_outer) then
      radial_taper_weight = 0.0_RKIND
   else
      radial_taper_weight = 0.5_RKIND * (1.0_RKIND + &
           cos(pi_rkind * (dist - radius_inner) / transition_width))
   endif
end function radial_taper_weight


pure real(kind=RKIND) function wrap_lon_rad(lon)
   implicit none
   real(kind=RKIND), intent(in) :: lon
   real(kind=RKIND) :: pi_rkind

   pi_rkind = acos(-1.0_RKIND)
   wrap_lon_rad = modulo(lon + pi_rkind, 2.0_RKIND * pi_rkind) - pi_rkind
end function wrap_lon_rad


pure real(kind=RKIND) function shortest_lon_diff_rad(lon, lon0)
   implicit none
   real(kind=RKIND), intent(in) :: lon, lon0

   shortest_lon_diff_rad = wrap_lon_rad(lon - lon0)
end function shortest_lon_diff_rad


pure real(kind=RKIND) function ellipse_taper_weight(east_dist, north_dist, radius_x, radius_y, &
      rotation_deg, transition_width)
   implicit none
   real(kind=RKIND), intent(in) :: east_dist, north_dist, radius_x, radius_y
   real(kind=RKIND), intent(in) :: rotation_deg, transition_width
   real(kind=RKIND) :: angle, x_rot, y_rot, dist_center, unit_x, unit_y, radius_outer, pi_rkind

   ellipse_taper_weight = 0.0_RKIND
   if (radius_x <= 0.0_RKIND .or. radius_y <= 0.0_RKIND) return

   pi_rkind = acos(-1.0_RKIND)
   angle = rotation_deg * pi_rkind / 180.0_RKIND
   x_rot = east_dist * cos(angle) + north_dist * sin(angle)
   y_rot = -east_dist * sin(angle) + north_dist * cos(angle)
   dist_center = sqrt(x_rot * x_rot + y_rot * y_rot)

   if (dist_center <= 0.0_RKIND) then
      ellipse_taper_weight = 1.0_RKIND
      return
   endif

   unit_x = x_rot / dist_center
   unit_y = y_rot / dist_center
   radius_outer = 1.0_RKIND / sqrt((unit_x / radius_x)**2 + (unit_y / radius_y)**2)
   ellipse_taper_weight = radial_taper_weight(dist_center, radius_outer, transition_width)
end function ellipse_taper_weight


pure real(kind=RKIND) function box_taper_weight(east_dist, north_dist, half_width_x, half_width_y, &
      transition_width)
   implicit none
   real(kind=RKIND), intent(in) :: east_dist, north_dist, half_width_x, half_width_y, transition_width
   real(kind=RKIND) :: dist_to_edge, transition, pi_rkind

   box_taper_weight = 0.0_RKIND
   if (half_width_x <= 0.0_RKIND .or. half_width_y <= 0.0_RKIND) return

   dist_to_edge = min(half_width_x - abs(east_dist), half_width_y - abs(north_dist))
   if (dist_to_edge < 0.0_RKIND) return
   if (transition_width <= 0.0_RKIND) then
      box_taper_weight = 1.0_RKIND
      return
   endif

   transition = min(transition_width, min(half_width_x, half_width_y))
   if (transition <= 0.0_RKIND .or. dist_to_edge >= transition) then
      box_taper_weight = 1.0_RKIND
   else
      pi_rkind = acos(-1.0_RKIND)
      box_taper_weight = 0.5_RKIND * (1.0_RKIND - cos(pi_rkind * dist_to_edge / transition))
   endif
end function box_taper_weight


subroutine apply_region_mask_to_tendency(var2d, cell_mask, its, ite, kts, kte)
   implicit none
   integer, intent(in) :: its, ite, kts, kte
   real(kind=RKIND), intent(inout) :: var2d(:,:)
   real(kind=RKIND), intent(in)    :: cell_mask(:)
   integer :: i, k

   do i = its, ite
      do k = kts, kte
         var2d(k,i) = var2d(k,i) * cell_mask(i)
      end do
   end do
end subroutine apply_region_mask_to_tendency



subroutine apply_diffusion_based_filter(domain, exchange_halo_group, its, ite, kts, kte, dx_min, p_dfilter, n_steps, areacellmax, invAreaCell, nEdgesOnCell, &
                        edgesOnCell, edgesOnCell_sign, dvEdge, invDcEdge, cellsOnEdge, &
                        varfddaten, filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term, fdda_uv, fdda_t, fdda_q, &
                        diffusion_filter_anisotropic)


   ! use mpas_atm_halos, only: exchange_halo_group

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
  logical, intent(in) :: fdda_uv, fdda_t, fdda_q, diffusion_filter_anisotropic
  
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

   !   ! Store initial varfddaten into filtervar_prev
   !   do i = its, ite
   !     do k = kts, kte
   !       filtervar_prev(k,i) = varfddaten(k,i)
   !     end do
   !   end do

  filtervar_prev(:,:) = varfddaten(:,:)

  ! Compute the 2nd-order laplacian filter: filtervar_prev
  laplacian_term(:,:) = 0._RKIND
  if (diffusion_filter_anisotropic) then
     do filter_iCell = its, ite
        kappacell = 1.0_RKIND / areacellmax / invAreaCell(filter_iCell)
        do i = 1, nEdgesOnCell(filter_iCell)
           filter_iEdge = edgesOnCell(i, filter_iCell)
           edge_sign = invAreaCell(filter_iCell) * edgesOnCell_sign(i, filter_iCell) * &
              dvEdge(filter_iEdge) * invDcEdge(filter_iEdge) * kappacell
           cell1 = cellsOnEdge(1, filter_iEdge)
           cell2 = cellsOnEdge(2, filter_iEdge)
           kappacell1 = 1.0_RKIND / areacellmax / invAreaCell(cell1)
           kappacell2 = 1.0_RKIND / areacellmax / invAreaCell(cell2)
           do k = kts, kte
              laplacian_term(k, filter_iCell) = laplacian_term(k, filter_iCell) + &
                 edge_sign * (filtervar_prev(k, cell2) - filtervar_prev(k, cell1)) + &
                 (filtervar_prev(k, cell2) - filtervar_prev(k, cell1)) * invDcEdge(filter_iEdge) * &
                 (kappacell2 - kappacell1) * invDcEdge(filter_iEdge)
           end do
        end do
     end do
  else
     do filter_iCell = its, ite
        do i = 1, nEdgesOnCell(filter_iCell)
           filter_iEdge = edgesOnCell(i, filter_iCell)
           edge_sign = invAreaCell(filter_iCell) * edgesOnCell_sign(i, filter_iCell) * &
              dvEdge(filter_iEdge) * invDcEdge(filter_iEdge)
           cell1 = cellsOnEdge(1, filter_iEdge)
           cell2 = cellsOnEdge(2, filter_iEdge)
           do k = kts, kte
              laplacian_term(k, filter_iCell) = laplacian_term(k, filter_iCell) + &
                 edge_sign * (filtervar_prev(k, cell2) - filtervar_prev(k, cell1))
           end do
        end do
     end do
  endif
  
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


! filtervar_prev(:,:) = varfddaten(:,:)
!$OMP BARRIER
   call exchange_halo_group(domain, 'physics:filtervar_curr')
!$OMP BARRIER


! Diffusion steps
do diffusion_step = 3, n_steps + 1
   
   laplacian_term(:,:) = 0._RKIND

   if (diffusion_filter_anisotropic) then
      do filter_iCell = its, ite
         kappacell = 1.0_RKIND / areacellmax / invAreaCell(filter_iCell)
         do i = 1, nEdgesOnCell(filter_iCell)
            filter_iEdge = edgesOnCell(i, filter_iCell)
            edge_sign = invAreaCell(filter_iCell) * edgesOnCell_sign(i, filter_iCell) * &
               dvEdge(filter_iEdge) * invDcEdge(filter_iEdge) * kappacell
            cell1 = cellsOnEdge(1, filter_iEdge)
            cell2 = cellsOnEdge(2, filter_iEdge)
            kappacell1 = 1.0_RKIND / areacellmax / invAreaCell(cell1)
            kappacell2 = 1.0_RKIND / areacellmax / invAreaCell(cell2)
            do k = kts, kte
               laplacian_term(k, filter_iCell) = laplacian_term(k, filter_iCell) + &
                  edge_sign * (filtervar_curr(k, cell2) - filtervar_curr(k, cell1)) + &
                  (filtervar_curr(k, cell2) - filtervar_curr(k, cell1)) * invDcEdge(filter_iEdge) * &
                  (kappacell2 - kappacell1) * invDcEdge(filter_iEdge)
            end do
         end do
      end do
   else
      do filter_iCell = its, ite
         do i = 1, nEdgesOnCell(filter_iCell)
            filter_iEdge = edgesOnCell(i, filter_iCell)
            edge_sign = invAreaCell(filter_iCell) * edgesOnCell_sign(i, filter_iCell) * &
               dvEdge(filter_iEdge) * invDcEdge(filter_iEdge)
            cell1 = cellsOnEdge(1, filter_iEdge)
            cell2 = cellsOnEdge(2, filter_iEdge)
            do k = kts, kte
               laplacian_term(k, filter_iCell) = laplacian_term(k, filter_iCell) + &
                  edge_sign * (filtervar_curr(k, cell2) - filtervar_curr(k, cell1))
            end do
         end do
      end do
   endif



   do i = its, ite
      do k = kts, kte
         filtervar_next(k, i) = 2.0_RKIND * (-2.0_RKIND * laplacian_term(k, i) * inv_s_max - filtervar_curr(k, i)) - filtervar_prev(k, i)
         filtered(k, i) = filtered(k, i) + p_dfilter(diffusion_step) * filtervar_next(k, i)
         ! filtervar_prev(k, i) = filtervar_curr(k, i)
         ! filtervar_curr(k, i) = filtervar_next(k, i)
      end do
   end do   

   filtervar_prev(:, :) = filtervar_curr(:, :)
   filtervar_curr(:, :) = filtervar_next(:, :)
   
   !$OMP BARRIER
      call exchange_halo_group(domain, 'physics:filtervar_curr')
   !$OMP BARRIER
   

   
end do



  ! Deallocate arrays
!   deallocate(filtered, filtervar_prev, filtervar_curr, filtervar_next, laplacian_term)

end subroutine apply_diffusion_based_filter



!==================================================================================================
 end module mpas_atmphys_fdda
!==================================================================================================
