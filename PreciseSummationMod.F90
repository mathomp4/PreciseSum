module precise_sum

   implicit none

   private

   public kahan_sum
   public pair_sum
   public dcomp_sum

   interface kahan_sum
      module procedure kahan_summation_R4_1
      module procedure kahan_summation_R8_1
      module procedure kahan_summation_R4_2
      module procedure kahan_summation_R8_2
      module procedure kahan_summation_R4_3
      module procedure kahan_summation_R8_3
   end interface

   interface pair_sum
      module procedure pairwise_summation_R4_1
      module procedure pairwise_summation_R8_1
   end interface

   interface shell_sort
      module procedure shell_sort_R4_1
      module procedure shell_sort_R8_1
   end interface

   interface dcomp_sum
      module procedure doubly_compensated_summation_R4_1
      module procedure doubly_compensated_summation_R8_1
   end interface

   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: sp = kind(1.0e0)

   contains

      function kahan_summation_R4_1(input) result(msum)

         implicit none

         real(kind=sp), intent(in ) :: input(:)
         real(kind=sp) :: msum

         real(kind=sp) :: comp, y, tmp

         integer :: i

         msum = 0.0_sp
         comp = 0.0_sp

         do i = 1, size(input,dim=1)
            y = comp + input(i)
            tmp = msum + y
            comp = (msum - tmp) + y
            msum = tmp
         end do

         msum = msum + comp

      end function kahan_summation_R4_1

      function kahan_summation_R8_1(input) result(msum)

         implicit none

         real(kind=dp), intent(in ) :: input(:)
         real(kind=dp) :: msum

         real(kind=dp) :: comp, y, tmp

         integer :: i

         msum = 0.0_dp
         comp = 0.0_dp

         do i = 1, size(input,dim=1)
            y = comp + input(i)
            tmp = msum + y
            comp = (msum - tmp) + y
            msum = tmp
         end do

         msum = msum + comp

      end function kahan_summation_R8_1

      function kahan_summation_R4_2(input) result(msum)

         implicit none

         real(kind=sp), intent(in ) :: input(:,:)
         real(kind=sp) :: msum

         real(kind=sp) :: comp, y, tmp

         integer :: i, j

         msum = 0.0_sp
         comp = 0.0_sp

         do i = 1, size(input,dim=1)
            do j = 1, size(input,dim=2)
               y = comp + input(i,j)
               tmp = msum + y
               comp = (msum - tmp) + y
               msum = tmp
            end do
         end do

         msum = msum + comp

      end function kahan_summation_R4_2

      function kahan_summation_R8_2(input) result(msum)

         implicit none

         real(kind=dp), intent(in ) :: input(:,:)
         real(kind=dp) :: msum

         real(kind=dp) :: comp, y, tmp

         integer :: i, j

         msum = 0.0_dp
         comp = 0.0_dp

         do i = 1, size(input,dim=1)
            do j = 1, size(input,dim=2)
               y = comp + input(i,j)
               tmp = msum + y
               comp = (msum - tmp) + y
               msum = tmp
            end do
         end do

         msum = msum + comp

      end function kahan_summation_R8_2

      function kahan_summation_R4_3(input) result(msum)

         implicit none

         real(kind=sp), intent(in ) :: input(:,:,:)
         real(kind=sp) :: msum

         real(kind=sp) :: comp, y, tmp

         integer :: i, j, k

         msum = 0.0_sp
         comp = 0.0_sp

         do i = 1, size(input,dim=1)
            do j = 1, size(input,dim=2)
               do k = 1, size(input,dim=3)
                  y = comp + input(i,j,k)
                  tmp = msum + y
                  comp = (msum - tmp) + y
                  msum = tmp
               end do
            end do
         end do

         msum = msum + comp

      end function kahan_summation_R4_3

      function kahan_summation_R8_3(input) result(msum)

         implicit none

         real(kind=dp), intent(in ) :: input(:,:,:)
         real(kind=dp) :: msum

         real(kind=dp) :: comp, y, tmp

         integer :: i, j, k

         msum = 0.0_dp
         comp = 0.0_dp

         do i = 1, size(input,dim=1)
            do j = 1, size(input,dim=2)
               do k = 1, size(input,dim=3)
                  y = comp + input(i,j,k)
                  tmp = msum + y
                  comp = (msum - tmp) + y
                  msum = tmp
               end do
            end do
         end do
         
         msum = msum + comp

      end function kahan_summation_R8_3

      function pairwise_summation_R4_1(input) result(psum)

         implicit none

         real(kind=sp), intent(in ) :: input(:)
         real(kind=sp) :: psum

         real(kind=sp), dimension(size(input)) :: working

         integer :: i, inext, icurrent, N

         psum = 0.0_sp

         N = size(input)

         working = input

         icurrent = N
         inext = ceiling(real(N)/2)

         do while (inext > 1)
            do i = 1, inext
               if ( 2*i <= icurrent ) working(i)=working(i)+working(i+inext)
            end do
            icurrent = inext
            inext = ceiling(real(inext)/2)
         end do

         psum = working(1)+working(2)

      end function pairwise_summation_R4_1

      function pairwise_summation_R8_1(input) result(psum)

         implicit none

         real(kind=dp), intent(in ) :: input(:)
         real(kind=dp) :: psum

         real(kind=dp), dimension(size(input)) :: working

         integer :: i, inext, icurrent, N

         psum = 0.0_dp

         N = size(input)

         working = input

         icurrent = N
         inext = ceiling(real(N)/2)

         do while (inext > 1)
            do i = 1, inext
               if ( 2*i <= icurrent ) working(i)=working(i)+working(i+inext)
            end do
            icurrent = inext
            inext = ceiling(real(inext)/2)
         end do

         psum = working(1)+working(2)

      end function pairwise_summation_R8_1

      subroutine shell_sort_R4_1(a)
 
         implicit none
         integer :: i, j, increment
         real(kind=sp) :: temp
         real(kind=sp), intent(inout) :: a(:)

         increment = size(a) / 2
         do while (increment > 0)
            do i = increment+1, size(a)
               j = i
               temp = a(i)
               do while (j >= increment+1 .and. a(j-increment) > temp)
                  a(j) = a(j-increment)
                  j = j - increment
               end do
               a(j) = temp
            end do
            if (increment == 2) then
               increment = 1
            else
               increment = increment * 5 / 11
            end if      
         end do

      end subroutine shell_sort_R4_1

      subroutine shell_sort_R8_1(a)
 
         implicit none
         integer :: i, j, increment
         real(kind=dp) :: temp
         real(kind=dp), intent(inout) :: a(:)

         increment = size(a) / 2
         do while (increment > 0)
            do i = increment+1, size(a)
               j = i
               temp = a(i)
               do while (j >= increment+1 .and. a(j-increment) > temp)
                  a(j) = a(j-increment)
                  j = j - increment
               end do
               a(j) = temp
            end do
            if (increment == 2) then
               increment = 1
            else
               increment = increment * 5 / 11
            end if      
         end do

      end subroutine shell_sort_R8_1

      function doubly_compensated_summation_R4_1(input) result(s)

         implicit none

         real(kind=sp), intent(in) :: input(:)
         real(kind=sp) :: s

         real(kind=sp) :: c, y, u, t, v, z
         real(kind=sp), dimension(size(input)) :: x
         integer :: k, N

         s = 0.0_sp
         c = 0.0_sp
         y = 0.0_sp
         u = 0.0_sp
         t = 0.0_sp
         v = 0.0_sp
         z = 0.0_sp

         N = size(input)

         ! Sort input into x
         x = input
         call shell_sort(x)

         s = x(1)
         
         do k = 2, N
            y = c + x(k)
            u = x(k) - (y - c)
            t = y + s
            v = y - (t - s)
            z = u + v
            s = t + z
            c = z - (s - t)
         end do

      end function doubly_compensated_summation_R4_1

      function doubly_compensated_summation_R8_1(input) result(s)

         implicit none

         real(kind=dp), intent(in) :: input(:)
         real(kind=dp) :: s

         real(kind=dp) :: c, y, u, t, v, z
         real(kind=dp), dimension(size(input)) :: x
         integer :: k, N

         s = 0.0_dp
         c = 0.0_dp
         y = 0.0_dp
         u = 0.0_dp
         t = 0.0_dp
         v = 0.0_dp
         z = 0.0_dp

         N = size(input)

         ! Sort input into x
         x = input
         call shell_sort(x)

         s = x(1)
         
         do k = 2, N
            y = c + x(k)
            u = x(k) - (y - c)
            t = y + s
            v = y - (t - s)
            z = u + v
            s = t + z
            c = z - (s - t)
         end do

      end function doubly_compensated_summation_R8_1

end module precise_sum
