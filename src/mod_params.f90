module mod_params
!-----------------------------------------------------------------------------
   use mod_prec
   implicit none
!-----------------------------------------------------------------------------
   complex(rp), parameter :: ZI = (0.0_rp,1.0_rp) 
!-----------------------------------------------------------------------------
   real(rp), protected :: black_hole_mass
   real(rp), protected :: black_hole_spin
   real(rp), protected :: compactification_length
   real(rp), protected :: evolve_time
   real(rp), protected :: start_multiple

   real(rp), protected :: dt
   real(rp), protected :: horizon
   real(rp), protected :: R_max 
   real(rp), protected :: constraint_damping
   real(rp), protected :: integrate_psi4_start_time
   real(rp), protected :: scd_order_start_time

   logical, protected :: sparse_save
   logical, protected :: metric_recon
   logical, protected :: scd_order
   logical, protected :: constrained_evo
   logical, protected :: write_indep_res
   logical, protected :: write_metric_recon_fields
   logical, protected :: write_scd_order_source
   logical, protected :: write_coefs
   logical, protected :: write_coefs_swal

   integer(ip), protected :: nx
   integer(ip), protected :: nl
   integer(ip), protected :: max_l
   integer(ip), protected :: ny
   integer(ip), protected :: nt
   integer(ip), protected :: t_step_save

   integer(ip), protected :: max_m, min_m

   integer(ip), protected :: psi_spin
   integer(ip), protected :: psi_boost

   integer(ip), parameter :: max_s =  3_ip
   integer(ip), parameter :: min_s = -3_ip

   integer(ip), protected :: len_lin_m
   integer(ip), protected :: len_scd_m
   integer(ip), protected :: len_lin_pos_m
   integer(ip), protected :: len_scd_pos_m
   integer(ip), protected :: len_write_lin_m
   integer(ip), protected :: len_write_scd_m

   integer(ip), allocatable, protected :: lin_m(:)
   integer(ip), allocatable, protected :: scd_m(:)
   integer(ip), allocatable, protected :: lin_pos_m(:)
   integer(ip), allocatable, protected :: scd_pos_m(:)
   integer(ip), allocatable, protected :: write_lin_m(:)
   integer(ip), allocatable, protected :: write_scd_m(:)

   character(:), allocatable, protected :: initial_data_direction
   real(rp),     allocatable, protected :: amp_re(:)
   real(rp),     allocatable, protected :: amp_im(:)
   real(rp),     allocatable, protected :: rl(:)
   real(rp),     allocatable, protected :: ru(:)
   integer(ip),  allocatable, protected :: l_ang(:)

   character(:), allocatable, protected :: output_dir
   character(:), allocatable, protected :: tables_dir
   
   !--------------------------------------------------------------------------

   public :: read_params

   interface read_param
      module procedure &
         read_param_bool, &
         read_param_int, &
         read_param_real, &
         read_param_complex, &
         read_param_int_vec, &
         read_param_real_vec, &
         read_param_char_vec
   end interface
!==============================================================================
contains
!==============================================================================
   subroutine read_param_bool(label,buffer,val)
      character(*), intent(in)    :: label
      character(*), intent(in)    :: buffer
      logical,      intent(out)   :: val

      integer(ip) :: ios
      read(buffer, *, iostat=ios) val
      print *, trim(label)//': ', val 
   end subroutine read_param_bool
!------------------------------------------------------------------------------
   subroutine read_param_int(label,buffer,val)
      character(*), intent(in)    :: label
      character(*), intent(in)    :: buffer
      integer(ip),  intent(out)   :: val

      integer(ip) :: ios
      read(buffer, *, iostat=ios) val
      print *, trim(label)//': ', val 
   end subroutine read_param_int
!------------------------------------------------------------------------------
   subroutine read_param_real(label,buffer,val)
      character(*), intent(in)    :: label
      character(*), intent(in)    :: buffer
      real(rp),     intent(out)   :: val

      integer(ip) :: ios
      read(buffer, *, iostat=ios) val
      print *, trim(label)//': ', val 
   end subroutine read_param_real
!------------------------------------------------------------------------------
   subroutine read_param_complex(label,buffer,val)
      character(*), intent(in)    :: label
      character(*), intent(in)    :: buffer
      complex(rp),  intent(out)   :: val

      integer(ip) :: ios
      read(buffer, *, iostat=ios) val
      print *, trim(label)//': ', val 
   end subroutine read_param_complex
!------------------------------------------------------------------------------
   subroutine read_param_int_vec(label,buffer,n,val)
      character(*),              intent(in)    :: label
      character(*),              intent(in)    :: buffer
      integer(ip),               intent(in)    :: n
      integer(ip),  allocatable, intent(inout) :: val(:)

      integer(ip) :: ios

      allocate(val(n))

      read(buffer, *, iostat=ios) val
      print *, trim(label)//': ', val
   end subroutine read_param_int_vec
!------------------------------------------------------------------------------
   subroutine read_param_real_vec(label,buffer,n,val)
      character(*),              intent(in)    :: label
      character(*),              intent(in)    :: buffer
      integer(ip),               intent(in)    :: n
      real(rp),     allocatable, intent(inout) :: val(:)

      integer(ip) :: ios

      allocate(val(n))

      read(buffer, *,iostat=ios) val
      print *, trim(label)//': ', val
   end subroutine read_param_real_vec
!------------------------------------------------------------------------------
   subroutine read_param_char_vec(label,buffer,val)
      character(*),              intent(in)    :: label
      character(*),              intent(in)    :: buffer
      character(:), allocatable, intent(inout) :: val

      val = trim(buffer)

      print *, trim(label)//': ', val
   end subroutine read_param_char_vec
!==============================================================================
   subroutine read_params()

      ! Input related variables
      integer(ip), parameter :: cl = 2000
      character(cl) :: buffer, label
      integer(ip)   :: pos
      integer(ip)   :: fh
      integer(ip)   :: ios = 0
      integer(ip)   :: line = 0
      character(cl) :: output_dir_inter

      call get_command_argument(1,output_dir_inter)
      output_dir = trim(output_dir_inter)
     
      print *, 'output_dir: ', output_dir
 
      open(newunit=fh, file=output_dir//'/sim_params.txt', action='read')

      ! ios is negative if an end of record condition is encountered or if
      ! an endfile condition was detected.  It is positive if an error was
      ! detected.  ios is zero otherwise.

      read_lines: do while (ios == 0)
         read(fh, '(A)', iostat=ios) buffer
         ios_ok: if (ios == 0) then
            line = line + 1

            ! Find the first instance of whitespace.  Split label and data.
            pos = scan(buffer, ' ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            param_select: select case (label)
               !--------------------------------------------------------------
               case ('black_hole_mass')
                  call read_param(label,buffer,black_hole_mass)
               !--------------------------------------------------------------
               case ('black_hole_spin')
                  call read_param(label,buffer,black_hole_spin)
               !--------------------------------------------------------------
               case ('compactification_length')
                  call read_param(label,buffer,compactification_length)
               !--------------------------------------------------------------
               case ('evolve_time')
                  call read_param(label,buffer,evolve_time)
               !--------------------------------------------------------------
               case ('start_multiple')
                  call read_param(label,buffer,start_multiple)

               !--------------------------------------------------------------
               case ('dt')
                  call read_param(label,buffer,dt)
               !--------------------------------------------------------------
               case ('horizon')
                  call read_param(label,buffer,horizon)
               !--------------------------------------------------------------
               case ('R_max')
                  call read_param(label,buffer,R_max)
               !--------------------------------------------------------------
               case ('constraint_damping')
                  call read_param(label,buffer,constraint_damping)
               !--------------------------------------------------------------
               case ('integrate_psi4_start_time')
                  call read_param(label,buffer,integrate_psi4_start_time)
               !--------------------------------------------------------------
               case ('scd_order_start_time')
                  call read_param(label,buffer,scd_order_start_time)

               !--------------------------------------------------------------
               case ('sparse_save')
                  call read_param(label,buffer,sparse_save)
               !--------------------------------------------------------------
               case ('metric_recon')
                  call read_param(label,buffer,metric_recon)
               !--------------------------------------------------------------
               case ('scd_order')
                  call read_param(label,buffer,scd_order)
               !--------------------------------------------------------------
               case ('constrained_evo')
                  call read_param(label,buffer,constrained_evo)
               !--------------------------------------------------------------
               case ('write_indep_res')
                  call read_param(label,buffer,write_indep_res)
               !--------------------------------------------------------------
               case ('write_metric_recon_fields')
                  call read_param(label,buffer,write_metric_recon_fields)
               !--------------------------------------------------------------
               case ('write_scd_order_source')
                  call read_param(label,buffer,write_scd_order_source)
               !--------------------------------------------------------------
               case ('write_coefs')
                  call read_param(label,buffer,write_coefs)
               !--------------------------------------------------------------
               case ('write_coefs_swal')
                  call read_param(label,buffer,write_coefs_swal)

               !--------------------------------------------------------------
               case ('nx')
                  call read_param(label,buffer,nx)
               !--------------------------------------------------------------
               case ('nl')
                  call read_param(label,buffer,nl)
               !--------------------------------------------------------------
               case ('max_l')
                  call read_param(label,buffer,max_l)
               !--------------------------------------------------------------
               case ('ny')
                  call read_param(label,buffer,ny)
               !--------------------------------------------------------------
               case ('nt')
                  call read_param(label,buffer,nt)
               !--------------------------------------------------------------
               case ('t_step_save')
                  call read_param(label,buffer,t_step_save)

               !--------------------------------------------------------------
               case ('psi_spin')
                  call read_param(label,buffer,psi_spin)
               !--------------------------------------------------------------
               case ('psi_boost')
                  call read_param(label,buffer,psi_boost)

               !--------------------------------------------------------------
               case ('len_lin_m')
                  call read_param(label,buffer,len_lin_m)
               !--------------------------------------------------------------
               case ('len_scd_m')
                  call read_param(label,buffer,len_scd_m)
               !--------------------------------------------------------------
               case ('len_lin_pos_m')
                  call read_param(label,buffer,len_lin_pos_m)
               !--------------------------------------------------------------
               case ('len_scd_pos_m')
                  call read_param(label,buffer,len_scd_pos_m)
               !--------------------------------------------------------------
               case ('len_write_lin_m')
                  call read_param(label,buffer,len_write_lin_m)
               !--------------------------------------------------------------
               case ('len_write_scd_m')
                  call read_param(label,buffer,len_write_scd_m)

               !--------------------------------------------------------------
               case ('lin_m')
                  call read_param(label,buffer,len_lin_m,lin_m)
               !--------------------------------------------------------------
               case ('scd_m')
                  call read_param(label,buffer,len_scd_m,scd_m)
               !--------------------------------------------------------------
               case ('lin_pos_m')
                  call read_param(label,buffer,len_lin_pos_m,lin_pos_m)
               !--------------------------------------------------------------
               case ('scd_pos_m')
                  call read_param(label,buffer,len_scd_pos_m,scd_pos_m)
               !--------------------------------------------------------------
               case ('write_lin_m')
                  call read_param(label,buffer,len_write_lin_m,write_lin_m)
               !--------------------------------------------------------------
               case ('write_scd_m')
                  call read_param(label,buffer,len_write_scd_m,write_scd_m)

               !--------------------------------------------------------------
               case ('initial_data_direction')
                  call read_param(label,buffer,initial_data_direction)
               !--------------------------------------------------------------
               case ('amp_re')
                  call read_param(label,buffer,len_lin_m,amp_re)
               !--------------------------------------------------------------
               case ('amp_im')
                  call read_param(label,buffer,len_lin_m,amp_im)
               !--------------------------------------------------------------
               case ('rl')
                  call read_param(label,buffer,len_lin_m,rl)
               !--------------------------------------------------------------
               case ('ru')
                  call read_param(label,buffer,len_lin_m,ru)
               !--------------------------------------------------------------
               case ('l_ang')
                  call read_param(label,buffer,len_lin_m,l_ang)

               !--------------------------------------------------------------
               case ('tables_dir')
                  call read_param(label,buffer,tables_dir)

               !--------------------------------------------------------------
               case ('min_m')
                  call read_param(label,buffer,min_m)
               case ('max_m')
                  call read_param(label,buffer,max_m)

               !--------------------------------------------------------------
               case default
                  print *, 'Skipping invalid label: ', trim(label)
            end select param_select
         end if ios_ok
      end do read_lines

   end subroutine read_params
!==============================================================================
end module mod_params
