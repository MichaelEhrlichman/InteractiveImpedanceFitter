module phase_advance_mod

use, intrinsic :: iso_c_binding, only: c_int, c_double
use bmad
use bmad_parser_mod, only: bp_com

implicit none

integer, parameter :: MAX_BPMS = 200
integer, parameter :: MAX_PARAMS = 1001

type(lat_struct) lat
integer bpm_ixs(MAX_BPMS)
integer nparams
integer nbpms

!parameters
character(200) lat_file
character(20) bpm_list(MAX_BPMS)
character(20) params_list(MAX_PARAMS)
integer params_ixs(MAX_PARAMS)

contains

subroutine parse_in_file(c_in_file, in_file_len, c_nbpms, c_nparams, c_param_names, upper_bounds, lower_bounds, initial_params, c_datafile, current_low, current_high) bind(c)
  implicit none

  integer(c_int), intent(in), value :: in_file_len
  character, dimension(in_file_len), intent(in) :: c_in_file
  integer(c_int), intent(out) :: c_nbpms, c_nparams
  character, dimension(20000), intent(out) :: c_param_names
  real(c_double), dimension(MAX_PARAMS), intent(out) :: upper_bounds, lower_bounds, initial_params
  character, dimension(200), intent(out) :: c_datafile
  real(c_double), intent(out) :: current_high, current_low

  character(20) one_param_name
  character(len=:), allocatable :: in_file
  character(200) datafile

  integer i,j,k

  namelist / parameters / lat_file, datafile, current_low, current_high, params_list, bpm_list, upper_bounds, lower_bounds, initial_params

  allocate(character(len=in_file_len) :: in_file)
  in_file = transfer(c_in_file(1:in_file_len), in_file)

  bpm_list = ''
  params_list = ''
  datafile = ''
  upper_bounds = 0.0
  lower_bounds = 0.0
  initial_params = 0.0

  open (unit = 10, file = in_file)
  read (10, nml = parameters)
  close (10)

  nbpms = 0
  do i=1,MAX_BPMS
    if(bpm_list(i) == '') exit
    nbpms = nbpms + 1 
  enddo

  nparams = 0
  do i=1,MAX_PARAMS
    if(params_list(i) == '') exit
    nparams = nparams + 1
  enddo

  k=0
  c_param_names=''
  do i=1,nparams
    one_param_name = params_list(i)
    do j=1,len(trim(adjustl(one_param_name)))
      k=k+1
      c_param_names(k) = one_param_name(j:j)
    enddo
    k=k+1
    c_param_names(k) = ","
  enddo

  do i=1,200
    c_datafile(i) = datafile(i:i)
  enddo

  c_nbpms = nbpms
  c_nparams = nparams
end subroutine

subroutine init() bind(c)
  use output_mod
  use bmad_routine_interface

  integer i,j
  integer n_loc
  type(ele_pointer_struct), allocatable :: eles(:)
  logical found

  params_ixs = -1

  !call output_direct(do_print=.false.,min_level=s_blank$,max_level=s_warn$)
  call output_direct(print_and_capture=.false.,min_level=s_blank$,max_level=s_warn$)

  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  bpm_ixs = 0

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, lat)

  do i=1,nbpms
    found = .false.
    do j=0,lat%n_ele_track
      if( bpm_list(i) == trim(adjustl(lat%ele(j)%name(5:7))) ) then
        bpm_ixs(i) = j
        found = .true.
        exit
      endif
    enddo
    if(.not. found) then
      write(*,*) "Initialization failed: BPM ", bpm_list(i), " not found!"
      call err_exit
    endif
  enddo

  do i=1,nparams
    if(trim(adjustl(params_list(i))) == '') EXIT
    if(params_list(i) == '__OFFSET') then
      params_ixs(i) = -2
    else
      if(allocated(eles)) deallocate(eles)
      call lat_ele_locator(params_list(i), lat, eles, n_loc)
      if(n_loc .ne. 1) then
        write(*,*) "Error: for parameter ", trim(adjustl(params_list(i))), ", n_loc is not 1!"
        call err_exit
      endif
      params_ixs(i) = eles(1)%ele%ix_ele
    endif
  enddo
end subroutine

subroutine phase_vertical_calc(p,low_current,high_current,delta_phi) bind(c)
  use bmad
  use output_mod

  implicit none

  real(c_double), dimension(nparams), intent(in) :: p
  real(c_double), intent(in), value :: low_current, high_current
  real(c_double), dimension(nbpms), intent(out) :: delta_phi

  real(c_double) phi_low(nbpms), phi_high(nbpms)

  type(coord_struct), allocatable :: orb(:)

  logical error

  integer i,j,k
  integer status

  character(18) var_str
  character(30) set_str

  real(rp) offset

  error = .false.
  do k=1,2
    do i=1,nparams
      if(params_list(i) == '__OFFSET') then
        offset = p(i)
      else
        if(k==1) write(var_str,'(e18.8)') p(i)*low_current
        if(k==2) write(var_str,'(e18.8)') p(i)*high_current
        set_str = 'k1l='//trim(adjustl(var_str))
        call set_ele_attribute(lat%ele(params_ixs(i)), set_str, error)
        if(error) then
          write(*,'(3a)') "Set ele attribute error.  Terminating.", params_list(i), set_str
          call err_exit
        endif
      endif
    enddo
    call lattice_bookkeeper(lat)

    call twiss_and_track(lat,orb,status)
    if(status .ne. ok$) then
      error = .true.
    endif

    do i=1,nbpms
      if(k==1) phi_low(i) = lat%ele(bpm_ixs(i))%b%phi
      if(k==2) phi_high(i) = lat%ele(bpm_ixs(i))%b%phi
    enddo
  enddo

  if(error) then
    delta_phi = 999.9
  else
    delta_phi = phi_high - phi_low + offset
  endif
end subroutine

end module












