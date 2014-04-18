
      subroutine froutine
      integer i, j
      common /vars/ i, j
      print *, "hello from fortran", i
      end

      subroutine finit
      integer i, j
      common /vars/ i, j
      i=555 ;
      end
