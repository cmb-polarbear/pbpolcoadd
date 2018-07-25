!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Julien Peloton, Dominic Beck
! FORTRAN routines to compute low level products
! Main purposes is interfacing with python (using f2py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module polcoadd
use omp_lib

contains

	subroutine tod2map_simple(signal_map, hits_map, pointing, array, npix, nt, mask, nskypix)
		implicit none

		integer, parameter		 :: I4B = 4
		integer, parameter		 :: DP = 8
		real(DP), parameter		 :: pi = 3.141592

		integer(I4B), intent(in) :: npix, nt, nskypix
		integer(I4B), intent(in) :: pointing(0:nt - 1), mask(0:npix*nt - 1)
		real(DP), intent(in)	 :: array(0:npix*nt - 1)

		real(DP), intent(inout)	 :: signal_map(0:npix*nskypix -1)
		integer(I4B), intent(inout) :: hits_map(0:npix*nskypix -1)

		integer(I4B)			 :: i, j, ipix
		integer(I4B)			 :: pixel

		!$omp parallel do private(i,j,ipix,pixel) shared(hits_map, signal_map)
		do j=0, npix - 1
			do i=0, nt - 1
				ipix = i + j*nt
				if (mask(ipix) .gt. 0 .and. pointing(i) .gt. 0) then
					pixel = pointing(i)+ j*nskypix
					hits_map(pixel) = hits_map(pixel) + 1
					signal_map(pixel) = signal_map(pixel) + array(ipix)
				endif
			enddo
		enddo
		
		!$omp parallel do private(i,j,pixel) shared(signal_map)
		do j=0, npix - 1
			do i=0, nskypix - 1
				pixel=i+j*nskypix
				if (hits_map(pixel) .gt. 0) then
					signal_map(pixel)=signal_map(pixel)/hits_map(pixel)
				endif
			enddo
		enddo
					
		
	end subroutine
	
	subroutine tod2map_temp_preweighted(d0, w0, nhit, waferi1d, &
	waferts, weight0, npix, nt, &
	wafermask_pixel, nskypix)
		implicit none

		integer, parameter		 :: I4B = 4
		integer, parameter		 :: DP = 8
		real(DP), parameter		 :: pi = 3.141592

		integer(I4B), intent(in) :: npix, nt, nskypix
		integer(I4B), intent(in) :: waferi1d(0:npix*nt - 1), wafermask_pixel(0:npix*nt - 1)
		real(DP), intent(in)	 :: waferts(0:npix*nt - 1)
		real(DP), intent(in)	 :: weight0(0:npix - 1)

		real(DP), intent(inout)	 :: d0(0:nskypix - 1)
		real(DP), intent(inout)	 :: w0(0:nskypix - 1)
		integer(I4B), intent(inout) :: nhit(0:nskypix - 1)

		integer(I4B)			 :: i, j, ipix
		integer(I4B)			 :: pixel
		integer(I4B)			 :: if0

		!$omp parallel do private(i,j,ipix,if0,pixel) shared(nhit, w0, d0)
		do j=0, npix - 1
			do i=0, nt - 1
				ipix = i + j*nt
				if (wafermask_pixel(ipix) .gt. 0 .and. waferi1d(ipix) .gt. 0) then
					if0 = i + j*nt

					pixel = waferi1d(ipix)


					nhit(pixel) = nhit(pixel) + 1

					w0(pixel) = w0(pixel) + weight0(j)
					d0(pixel) = d0(pixel) + waferts(if0)
				endif
			enddo
		enddo
	end subroutine
	
	subroutine tod2map_temp(d0, w0, nhit, waferi1d, &
	waferts, weight0, npix, nt, &
	wafermask_pixel, nskypix)
		implicit none

		integer, parameter		 :: I4B = 4
		integer, parameter		 :: DP = 8
		real(DP), parameter		 :: pi = 3.141592

		integer(I4B), intent(in) :: npix, nt, nskypix
		integer(I4B), intent(in) :: waferi1d(0:npix*nt - 1), wafermask_pixel(0:npix*nt - 1)
		real(DP), intent(in)	 :: waferts(0:npix*nt - 1)
		real(DP), intent(in)	 :: weight0(0:npix - 1)

		real(DP), intent(inout)	 :: d0(0:nskypix - 1)
		real(DP), intent(inout)	 :: w0(0:nskypix - 1)
		integer(I4B), intent(inout) :: nhit(0:nskypix - 1)

		integer(I4B)			 :: i, j, ipix
		integer(I4B)			 :: pixel
		integer(I4B)			 :: if0

		!$omp parallel do private(i,j,ipix,if0,pixel) shared(nhit, w0, d0)
		do j=0, npix - 1
			do i=0, nt - 1
				ipix = i + j*nt
				if (wafermask_pixel(ipix) .gt. 0 .and. waferi1d(ipix) .gt. 0) then
					if0 = i + j*nt

					pixel = waferi1d(ipix)


					nhit(pixel) = nhit(pixel) + 1

					w0(pixel) = w0(pixel) + weight0(j)
					d0(pixel) = d0(pixel) + waferts(if0) * weight0(j)
				endif
			enddo
		enddo
	end subroutine

	subroutine tod2map_pol(d0, d4r, d4i, w0, w4, nhit, waferi1d, &
	waferpa, waferts, weight4, weight0, npix, nt, &
	wafermask_pixel, nskypix)
		implicit none

		integer, parameter		 :: I4B = 4
		integer, parameter		 :: DP = 8
		real(DP), parameter		 :: pi = 3.141592

		integer(I4B), intent(in) :: npix, nt, nskypix
		integer(I4B), intent(in) :: waferi1d(0:npix*nt - 1), wafermask_pixel(0:npix*nt - 1)
		real(DP), intent(in)	 :: waferpa(0:npix*nt - 1), waferts(0:npix*nt*3 - 1)
		real(DP), intent(in)	 :: weight0(0:npix - 1), weight4(0:npix - 1)

		real(DP), intent(inout)	 :: d0(0:nskypix - 1), d4r(0:nskypix - 1), d4i(0:nskypix - 1)
		real(DP), intent(inout)	 :: w0(0:nskypix - 1), w4(0:nskypix - 1)
		integer(I4B), intent(inout) :: nhit(0:nskypix - 1)

		integer(I4B)			 :: i, j, ipix
		integer(I4B)			 :: pixel
		integer(I4B)			 :: if0, i4r, i4i
		real(DP)				 :: c, s

		!$omp parallel do private(i,j,ipix,if0,i4r,i4i,c,s,pixel) shared(nhit, w0, w4, d0, d4r, d4i)
		do j=0, npix - 1
			do i=0, nt - 1
				ipix = i + j*nt
				if (wafermask_pixel(ipix) .gt. 0 .and. waferi1d(ipix) .gt. 0) then
					if0 = i + j*3*nt
					i4r = i + nt + j*3*nt
					i4i = i + nt*2 + j*3*nt

					pixel = waferi1d(ipix)

					c = cos(2.0*waferpa(ipix))
					s = sin(2.0*waferpa(ipix))

					nhit(pixel) = nhit(pixel) + 1

					w0(pixel) = w0(pixel) + weight0(j)
					w4(pixel) = w4(pixel) + weight4(j)
					d0(pixel) = d0(pixel) + waferts(if0) * weight0(j)
					d4r(pixel) = d4r(pixel) + (c*waferts(i4r)+s*waferts(i4i)) * weight4(j)
					d4i(pixel) = d4i(pixel) + (s*waferts(i4r)-c*waferts(i4i)) * weight4(j)
				endif
			enddo
		enddo
	end subroutine

end module