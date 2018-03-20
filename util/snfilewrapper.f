*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     snFileOpenRead   snFileOpenAppend   snFileClose
*
*  This file is here for the Python/C/C++ interface.
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFileOpenRead
     &     ( unitno, filename )
      implicit
     &     none
      integer
     &     unitno
      character*(*)
     &     filename
!     ==================================================================
!     Intel compiler doesn't like opening a file in one library and
!     passing it to another.
!
!     Wrapper to open files for reading inside SNOPT library.
!
!     10 Jan 2017: First version.
!     ==================================================================
      integer len

      call s1trim(filename,len)
      open(unit=unitno, file=filename, status='old')

      end ! subroutine snFileOpenRead

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFileOpenAppend
     &     ( unitno, filename )
      implicit
     &     none
      integer
     &     unitno
      character*(*)
     &     filename
!     ==================================================================
!     Intel compiler doesn't like opening a file in one library and
!     passing it to another.
!
!     Wrapper to open files for append inside SNOPT library.
!
!     10 Jan 2017: First version.
!     ==================================================================
      integer len

      call s1trim(filename,len)
      open(unit=unitno, file=filename,
     &     status='unknown', position='append')

      end ! subroutine snFileOpenAppend

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snFileClose
     &     ( unitno )
      implicit
     &     none
      integer
     &     unitno
!     ==================================================================
!     Close the given unit file.
!
!     10 Jan 2017: First version.
!     ==================================================================

      close(unitno)

      end ! subroutine snFileClose

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
