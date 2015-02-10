!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! snopt7_python.f90
!
! To create sig file (*.pyf):
!     f2py snopt7_python.f90 -m snopt7_python -h snopt7_python.pyf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module snopt_interfaces
  public

  interface
     subroutine ifuncon ( modeC, nnCon, nnJac, negCon, &
                          x, fCon, gCon, statusUser, &
                          cu, lencu, iu, leniu, ru, lenru )
       implicit none
       integer,          intent(in)    :: nnCon, nnJac, negCon, statusUser, &
                                          lencu, leniu, lenru
       double precision, intent(in)    :: x(nnJac)
       character(8),     intent(inout) :: cu(lencu)
       integer,          intent(inout) :: modeC, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       double precision, intent(out)   :: fCon(nnCon), gCon(negCon)
     end subroutine ifuncon

     subroutine ifunobj ( modeF, nnObjU, x, fObj, gObj, statusUser, &
                          cu, lencu, iu, leniu, ru, lenru )
       implicit none
       integer,          intent(in)    :: nnObjU, statusUser, lencu, leniu, lenru
       double precision, intent(in)    :: x(nnObjU)
       character(8),     intent(inout) :: cu(lencu)
       integer,          intent(inout) :: modeF, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       double precision, intent(out)   :: fObj, gObj(nnObjU)
     end subroutine ifunobj

     subroutine iusrfun ( Status, n, xN, needF, nF, F, needG, lenG, G, &
                          cu, lencu, iu, leniu, ru, lenru )
       implicit none
       integer,          intent(in)    :: n, nF, needF, needG, lenG, &
                                          lencu, leniu, lenru
       double precision, intent(in)    :: xN(n)
       character(8),     intent(inout) :: cu(lencu)
       integer,          intent(inout) :: Status, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       double precision, intent(out)   :: F(nF), G(lenG)
     end subroutine iusrfun

  end interface
end module snopt_interfaces

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine sntitle_wrap ( iw, leniw )
  implicit none
  integer,          intent(in)    :: leniw
  integer,          intent(inout) :: iw(leniw)

  !=============================================================================
  ! Print title for python interface.
  !   snInit should have been called prior to this.
  !
  !=============================================================================
  character(30) :: title
  character(30), parameter :: dashes = '=============================='

  call snTitle( title )

  call snPRNT (11, '         '//dashes, iw, leniw )
  call snPRNT ( 1, '         '//title , iw, leniw )
  call snPRNT ( 1, '         '//dashes, iw, leniw )

  call snPRNT (12, ' '//dashes, iw, leniw )
  call snPRNT ( 2, ' '//title , iw, leniw )
  call snPRNT ( 2, ' '//dashes, iw, leniw )

end subroutine sntitle_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snopen_wrap ( iPrint, prtfile, iSumm, sumfile )
  implicit none
  character(*),     intent(in)    :: prtfile, sumfile
  integer,          intent(in)    :: iPrint, iSumm
  !=============================================================================
  !=============================================================================
  if ( iPrint > 0 )                  &
       open ( iPrint, file=prtfile, status='unknown', position='append')
  if ( iSumm  > 0 .and. iSumm /= 6 ) &
       open ( iSumm,  file=sumfile, status='unknown', position='append')

end subroutine snopen_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine sninit_wrap ( cw, lencw, iw, leniw, rw, lenrw )
  implicit none
  integer,          intent(in)    :: lencw, leniw, lenrw
  character(8),     intent(inout) :: cw(lencw)
  integer,          intent(inout) :: iw(leniw)
  double precision, intent(inout) :: rw(lenrw)
  !=============================================================================
  !=============================================================================

  call snInit ( 0, 0, cw, lencw, iw, leniw, rw, lenrw )

end subroutine sninit_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snspec_wrap ( info, iSpecs, spcfile, cw, lencw, iw, leniw, rw, lenrw )
  implicit none
  character(*),     intent(in)    :: spcfile
  integer,          intent(in)    :: iSpecs, lencw, leniw, lenrw
  character(8),     intent(inout) :: cw(lencw)
  integer,          intent(inout) :: iw(leniw)
  double precision, intent(inout) :: rw(lenrw)
  integer,          intent(out)   :: info
  !=============================================================================
  !=============================================================================


  if ( iSpecs > 0 ) then
     open ( iSpecs, file=spcfile, status='unknown')
     call snSpec ( iSpecs, info, cw, lencw, iw, leniw, rw, lenrw )
     close (iSpecs )
     if ( info /= 101 .and. info /= 107) return
  else
     info = 131
  end if

end subroutine snspec_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmema_wrap ( info, nF, n, nxname, nFname, neA, neG, &
                         mincw, miniw, minrw, &
                         cw, lencw, iw, leniw, rw, lenrw )
  implicit none
  integer,          intent(in)    :: nF, n, nxname, nFname, neA, neG, &
                                     lencw, leniw, lenrw
  character(8),     intent(inout) :: cw(lencw)
  integer,          intent(inout) :: iw(leniw)
  double precision, intent(inout) :: rw(lenrw)
  integer,          intent(out)   :: info, mincw, miniw, minrw
  !=============================================================================
  !=============================================================================

  call snMemA ( info, nF, n, nxname, nFname, neA, neG, &
                mincw, miniw, minrw, cw, lencw, iw, leniw, rw, lenrw )

end subroutine snmema_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snmemb_wrap ( info, m, n, ne, neG, nnCon, nnJac, nnObj, &
                         mincw, miniw, minrw, &
                         cw, lencw, iw, leniw, rw, lenrw )
  implicit none
  integer,          intent(in)    :: m, n, ne, neG, nnCon, nnJac, nnObj, &
                                     lencw, leniw, lenrw
  character(8),     intent(inout) :: cw(lencw)
  integer,          intent(inout) :: iw(leniw)
  double precision, intent(inout) :: rw(lenrw)
  integer,          intent(out)   :: info, mincw, miniw, minrw
  !=============================================================================
  !=============================================================================

  call snMemB ( info, m, n, ne, neG, nnCon, nnJac, nnObj, &
                mincw, miniw, minrw, cw, lencw, iw, leniw, rw, lenrw )

end subroutine snmemb_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snjac_wrap ( INFO, nF, n, usrfun,                &
                        iAfun, jAvar, lenA, neA, A,         &
                        iGfun, jGvar, lenG, neG,            &
                        x, xlow, xupp, mincw, miniw, minrw, &
                        cu, lencu, iu, leniu, ru, lenru,    &
                        cw, lencw, iw, leniw, rw, lenrw )
  use snopt_interfaces
  implicit none

  procedure(iusrfun)              :: usrfun
  integer,          intent(in)    :: nF, n, lenA, lenG, &
                                     lencu, leniu, lenru, &
                                     lencw, leniw, lenrw
  double precision, intent(in)    :: xlow(n), xupp(n)

  character(8),     intent(inout) :: cu(lencu), cw(lencw)
  integer,          intent(inout) :: iu(leniu), iw(leniw)
  double precision, intent(inout) :: x(n), ru(lenru), rw(lenrw)

  integer,          intent(out)   :: info, mincw, miniw, minrw, neA, neG, &
                                     iAfun(lenA), jAvar(lenA), &
                                     iGfun(lenG), jGvar(lenG)
  double precision, intent(out)   :: A(lenA)

  !=============================================================================
  !=============================================================================

  call snJac ( INFO, nF, n, usrfun,                &
               iAfun, jAvar, lenA, neA, A,         &
               iGfun, jGvar, lenG, neG,            &
               x, xlow, xupp, mincw, miniw, minrw, &
               cu, lencu, iu, leniu, ru, lenru,    &
               cw, lencw, iw, leniw, rw, lenrw )

end subroutine snjac_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snopta_wrap ( start, nF, n, nxname, nFname,           &
                         objUAdd, objRow, prob, usrfun,          &
                         iAfun, jAvar, lenA, neA, A,             &
                         iGfun, jGvar, lenG, neG,                &
                         xlow, xupp, xnames, Flow, Fupp, Fnames, &
                         x, xstate, xmul, F, Fstate, Fmul,       &
                         INFO, itn, mjritn,                      &
                         mincw, miniw, minrw,                    &
                         nS, nInf, sInf, Obj,                    &
                         cu, lencu, iu, leniu, ru, lenru,        &
                         cw, lencw, iw, leniw, rw, lenrw )
  use snopt_interfaces
  implicit none

  procedure(iusrfun)              :: usrfun

  character(8),     intent(in)    :: Prob
  integer,          intent(in)    :: start, nF, n, nxname, nFname, &
                                     ObjRow, neA, lenA, lenG, neG, &
                                     iAfun(lenA), jAvar(lenA), &
                                     iGfun(lenG), jGvar(lenG), &
                                     lencu, leniu, lenru, &
                                     lencw, leniw, lenrw
  double precision, intent(in)    :: ObjUAdd, A(lenA), xlow(n), xupp(n), &
                                     Flow(nF), Fupp(nF)
  character(8),     intent(in)    :: xNames(nxname), Fnames(nFname)

  character(8),     intent(inout) :: cu(lencu), cw(lencw)
  integer,          intent(inout) :: xstate(n), Fstate(nF), iu(leniu), iw(leniw)
  double precision, intent(inout) :: x(n), xmul(n), F(nF), Fmul(nF), &
                                     ru(lenru), rw(lenrw)

  integer,          intent(out)   :: info, itn, mjritn, nS, nInf, &
                                     mincw, miniw, minrw
  double precision, intent(out)   :: sInf, Obj

  !=============================================================================
  !=============================================================================

  call snopta ( start, nF, n, nxname, nFname,           &
                objUAdd, objRow, trim(Prob), usrfun,    &
                iAfun, jAvar, lenA, neA, A,             &
                iGfun, jGvar, lenG, neG,                &
                xlow, xupp, xnames, Flow, Fupp, Fnames, &
                x, xstate, xmul, F, Fstate, Fmul,       &
                INFO, mincw, miniw, minrw,              &
                nS, nInf, sInf,                         &
                cu, lencu, iu, leniu, ru, lenru,        &
                cw, lencw, iw, leniw, rw, lenrw )

  itn    = iw(421)
  mjritn = iw(422)
  Obj    = 0.0
  if ( ObjRow > 0 ) Obj = F(objRow)

end subroutine snopta_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine snoptb_wrap ( Start, m, n, neJ, nNames,        &
                         nnCon, nnObjU, nnJac,            &
                         iObjU, objUAdd, Prob,            &
                         funcon, funobj,                  &
                         valJ, indJ, locJ, bl, bu, Names, &
                         hs, x, pi, rc,                   &
                         INFO, itn, mjritn,               &
                         mincw, miniw, minrw,             &
                         nS, nInf, sInf, obj,             &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
  use snopt_interfaces
  implicit none

  procedure(ifuncon)              :: funcon
  procedure(ifunobj)              :: funobj

  character(*),     intent(in)    :: Start
  character(8),     intent(in)    :: Prob
  integer,          intent(in)    :: m, n, neJ, nNames, nnCon, nnJac, nnObjU, &
                                     iObjU, indJ(neJ), locJ(n+1), &
                                     lencu, leniu, lenru, &
                                     lencw, leniw, lenrw
  double precision, intent(in)    :: ObjUAdd, valJ(neJ), bl(n+m), bu(n+m)
  character(8),     intent(in)    :: Names(nNames)

  character(8),     intent(inout) :: cu(lencu), cw(lencw)
  integer,          intent(inout) :: hs(n+m), iu(leniu), iw(leniw)
  double precision, intent(inout) :: x(n+m), pi(m), rc(n+m), ru(lenru), rw(lenrw)

  integer,          intent(out)   :: info, itn, mjritn, nS, nInf, &
                                     mincw, miniw, minrw
  double precision, intent(out)   :: sInf, Obj

  !=============================================================================
  !=============================================================================

  call snoptb ( Start, m, n, neJ, nNames,        &
                nnCon, nnObjU, nnJac,            &
                iObjU, objUAdd, trim(Prob),      &
                funcon, funobj,                  &
                valJ, indJ, locJ, bl, bu, Names, &
                hs, x, pi, rc,                   &
                INFO, mincw, miniw, minrw,       &
                nS, nInf, sInf, obj,             &
                cu, lencu, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw )

  itn    = iw(421)
  mjritn = iw(422)

end subroutine snoptb_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine close_wrap ( iPrint, iSumm )
  integer, intent(in) :: iPrint, iSumm
  !=============================================================================
  !=============================================================================

  if ( iPrint > 0 )                  close ( iPrint )
  if ( iSumm  > 0 .and. iSumm /= 6 ) close ( iSumm  )

end subroutine close_wrap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
