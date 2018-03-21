module utils
contains

subroutine mts(wn, pts, qn, M, T, S)

	implicit none
	integer :: i
	integer, intent(in) :: wn, qn
	real(kind=kind(0.0d0)), dimension(0:100), intent(in) :: pts
	real(kind=kind(0.0d0)), intent(out) :: M, T, S

	M = 0; T = 0; S = 0
	do i=1,qn-1
		M = M + fn(wn,(pts(i)+pts(i-1))/2.0d0)*(pts(i)-pts(i-1))
		T = T + ((fn(wn,pts(i-1))+fn(wn,pts(i)))*(pts(i)-pts(i-1)))/2.0d0
	enddo
	S = (2.0d0/3.0d0)*M + T/3.0d0


endsubroutine mts

real(kind=kind(0.0d0)) function fn(kw, point)

	implicit none
	integer, intent(in) :: kw
	real(kind=kind(0.0d0)), intent(in) :: point

	fn = cos(100.0d0*point-kw*sin(point))

endfunction fn

endmodule utils
