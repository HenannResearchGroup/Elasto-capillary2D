!************************************************************************
! User element (UEL) for large-deformation viscoelastic deformation 
!  in two-dimensions with surface tension -- plane-strain and axisymmetric
!
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-2).
!
! Material behavior is Neo-Hookean or Gent rubber elasticity 
!  with a finite-deformation Maxwell element in parallel.
!
!                 |
!           -------------
!           |           |  Maxwell element:
! Neo-Hooke |           |
!  or Gent  |           \
!  spring:  \           /    Hencky spring:Gneq
!   G0      /           |    
!   Kbulk   \           |
!  (Im)     /           |    Linearly viscous
!           |         | - |   Mises dashpot:eta
!           |         |___|
!           |           |
!           -------------
!                 |
! 
! This subroutine is for a two-dimensional 4-node isoparametric
!  quadrilateral element as shown below with 4pt (full) or 1pt 
!  (reduced) integration, as well as plane-strain or
!  axisymmetric settings.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
! Surface tension contributions are included in this element.  
!  Based on our convention, the face on which the surface tension
!  acts is the "label", i.e.
!  - U1,U2,U3,U4 refer to surface tension acting
!    on faces 1,2,3,4, respectively,
!  Mechanical traction- and pressure-type boundary conditions 
!  may be applied to the dummy mesh using the Abaqus built-in 
!  commands *Dload or *Dsload.
!     
!
!              A xi_2
!  4-node      |
!   quad       | Face 3
!        4-----------3
!        |     |     |
!        |     |     |
! Face 4 |     ------|---> xi_1
!        |           |
!        |           |  Face 2
!        1-----------2
!            Face 1
!
! David L. Henann, July 2013
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=4,Type=U1,Iproperties=2,Properties=5,Coordinates=2,Variables=??,Unsymm
!  1,2
!
!     State Variables
!       jj = 0
!       do k = 1,nIntPt
!          svars(1+jj) = Fv(1,1) ---- Fv(1,1) at integ pt k
!          svars(2+jj) = Fv(1,2) ---- Fv(1,2) at integ pt k
!          svars(3+jj) = Fv(1,3) ---- Fv(1,3) at integ pt k
!          svars(4+jj) = Fv(2,1) ---- Fv(2,1) at integ pt k
!          svars(5+jj) = Fv(2,2) ---- Fv(2,2) at integ pt k
!          svars(6+jj) = Fv(2,3) ---- Fv(2,3) at integ pt k
!          svars(7+jj) = Fv(3,1) ---- Fv(3,1) at integ pt k
!          svars(8+jj) = Fv(3,2) ---- Fv(3,2) at integ pt k
!          svars(9+jj) = Fv(3,3) ---- Fv(3,3) at integ pt k
!          jj = jj + nlSdv
!       end loop over k
!
! A heuristic displacement step restriction is included to aid in the simulation
!  of instabilities.  In the subroutine UEL, set maxDisp to the maximum allowable
!  displacement increment in terms of a fraction of the element size.
!
! In the subroutine UEL, set 'nInt' = number of integration points
!  Options are nInt=4 (full integration) or nInt=1 (reduced integration).
!
! In the input file, in the *User Element command, set Variables=(nInt*9).
!
! In the input file, set the parameter pe=1 for plane-strain or 
!  pe=0 for axisymmetric.
!
! In the input file, set the parameter matflag=1 for Neo-Hookean material 
!  behavior or matflag=2 for Gent.
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     G0      = props(1)  ! Ground-state shear modulus
!     Kbulk   = props(2)  ! Bulk modulus
!     Im      = props(3)  ! Limited chain extensibility parameter (Gent only)
!     Gneq    = props(4)  ! Nonequilibrium (Maxwell element) stiffness
!     eta     = props(5)  ! Maxwell element viscosity
!     pe      = jprops(1) ! Plane strain=1, axisymmetric=0
!     matflag = jprops(2) ! Neo-Hookean=1, Gent=2
!
!************************************************************************

      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     +     props,nprops,coords,mcrd,nnode,uall,duall,vel,accn,jtype,
     +     time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     +     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,
     +     njprop,period)

      implicit none
      !
      ! variables defined in uel, passed back to Abaqus
      !
      real*8 rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
      !
      ! variables passed into UEL
      !
      integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      real*8 props(*),coords(mcrd,nnode),uall(ndofel),duall(mlvarx,*),
     +  vel(ndofel),accn(ndofel),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      !
      ! variables defined and used in the UEL
      !
      integer i,j,k,jj,A11,B11,A12,B12,nInt,nIntPt,intpt,nDim,nSdv,stat,
     +  pe,matflag,face
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=2) ! number of spatial dimensions, do not change
      parameter(nSdv=9) ! number of state variables per integ pt, do not change
      parameter(nInt=4) ! number of integration points
      !
      ! nInt=4: fully-integrated, nInt=1: reduced integration
      ! 
      ! When reduced integration is selected, be sure to
      !  specify an hourglass stiffness of about 0.005*G0
      !  for the dummy mesh in the input file.  Also, change
      !  the dummy mesh element to reduced.
      !
      real*8 maxDisp
      parameter(maxDisp=0.20d0) ! maximum allowed displacement
                                !  increment in terms of a fraction 
                                !  of the element size
      !
      ! For the cavitation simulation maxDisp was taken as 0.20.
      ! For the filament/tube simulations maxDisp was taken as 0.25.
      !  Since simulating instabilities can be tricky, these parameters
      !  were chosen by trial and error to yield a robust simulation
      !  for a given geometry.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 u(nNode,2),uOld(nNode,ndofel),du(nNode,ndofel),
     +  coordsC(mcrd,nNode),Ru(2*nNode,1),Kuu(2*nNode,2*nNode),
     +  Iden(3,3),xi(nInt,2),w(nInt),sh0(nNode),sh(nNode),
     +  dsh0(nNode,2),dshC0(nNode,2),dsh(nNode,2),dshC(nNode,2),
     +  dshxi(nNode,2),detMapJ0,detMapJ0C,detMapJ,detMapJC,
     +  Fc_tau(3,3),detFc_tau,Fc_t(3,3),detFc_t,F_tau(3,3),detF_tau,
     +  F_t(3,3),detF_t,T_tau(3,3),Fv_t(3,3),Fv_tau(3,3),
     +  SpTanMod(3,3,3,3),Smat(3,1),Bmat(3,2*nNode),Gmat(4,2*nNode),
     +  G0mat(4,2*nNode),Amat(4,4),Qmat(4,4),body(3),
     +  BodyForceRes(2*nNode,1),SmatAx(4,1),BmatAx(4,2*nNode),
     +  GmatAx(5,2*nNode),G0MatAx(5,2*nNode),AmatAx(5,5),QmatAx(5,5),
     +  BodyForceResAx(2*nNode,1),AR,AR0,ARc,AR_t,Le,gamma
      !
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     +     Pi=3.141592653d0,three=3.d0,third=1.d0/3.d0)


      ! Check the procedure type; this should be a 
      !  *Static step, which is either 1 or 2
      !
      if((lflags(1).eq.1).or.(lflags(1).eq.2)) then
         !
         ! correct procedure specified
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear perturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear perturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear perturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return


      ! Get flag for plane strain or axisymmetric
      !
      pe = jprops(1)


      ! Get flag for material behavior
      !
      matflag = jprops(2)


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices and energy to zero.
      !
      Ru  = zero
      Kuu = zero
      Energy = zero


      ! Obtain nodal displacements
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         end do
      end do


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose the heuristic restriction that the displacement 
      !  increment is less than some fraction of the element diagonal,
      !  which is taken as an approximate element size.
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +                    ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.maxDisp*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the 'F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33, 3277-3296.
      !
      ! Obtain shape functions and their local gradients at the element
      !  centroid, that means xi_1=xi_2=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      AR0  = one ! for plane-strain you could put the depth here
      ARc  = one !
      AR_t = one !
      if(pe.ne.1) then ! this is axisymmetric
         AR0 = zero
         ARc = zero
         AR_t = zero
         do i=1,nNode
            ! radial coord in ref config at centroid
            AR0 = AR0 + sh0(i)*coords(1,i)
            ! radial coord in current config at centroid
            ARc = ARc + sh0(i)*(coords(1,i) + u(i,1))
            ! radial coord at the beginning of the step at centroid
            AR_t = AR_t + sh0(i)*(coords(1,i) + uOld(i,1))
         enddo
      endif


      ! Calculate the deformation gradient at the element centroid
      !  at the beginning and end of the increment for use in the 
      !  'F-bar' method.  The subscript tau denotes the time at the 
      !  end of the increment, while t denotes the time at the 
      !  beginning of the increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify F(3,3) for plane-strain or axisymmetric
      !
      if(pe.ne.1) then
         !
         ! modify for axisymmetric
         !
         Fc_tau(3,3) = ARc/AR0 ! 'hoop' stretch
         Fc_t(3,3) = AR_t/AR0
         !
         ! axisymmetric implementation of detF
         !
         call mdet(Fc_tau,detFc_tau)
         call mdet(Fc_t,detFc_t)
      else
         !
         ! modify for plane-strain
         !
         Fc_tau(3,3) = one
         Fc_t(3,3) = one
         !
         ! 2D plane-strain implementation detF
         !
         detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
         detFc_t = Fc_t(1,1)*Fc_t(2,2) - Fc_t(1,2)*Fc_t(2,1)
      endif
      !
      ! With the deformation gradient known at the element centroid
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! This is the first increment of the first step.
            !  Give initial conditions.
            !
            Fv_t = Iden
            !
         else
            !
            ! This is not the first increment; read old values.
            !
            Fv_t(1,1) = svars(1+jj)
            Fv_t(1,2) = svars(2+jj)
            Fv_t(1,3) = svars(3+jj)
            Fv_t(2,1) = svars(4+jj)
            Fv_t(2,2) = svars(5+jj)
            Fv_t(2,3) = svars(6+jj)
            Fv_t(3,1) = svars(7+jj)
            Fv_t(3,2) = svars(8+jj)
            Fv_t(3,3) = svars(9+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            call xit 
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            call xit 
         endif


         AR0  = one ! for plane-strain you could put the depth here
         AR   = one !
         AR_t = one !
         if(pe.ne.1) then ! this is axisymmetric
            AR0  = zero
            AR   = zero
            AR_t = zero
            do i=1,nNode
               ! radial coord in reference config
               AR0 = AR0 + sh(i)*coords(1,i)
               ! radial coord in current config
               AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
               ! radial coord at beginning of step
               AR_t = AR_t + sh(i)*(coords(1,i) + uOld(i,1))
            enddo
            AR0  = two*Pi*AR0
            AR   = two*Pi*AR
            AR_t = two*Pi*AR_t
         endif


         ! Obtain the deformation gradient at this integration point.
         ! The subscript tau denotes the time at the end of the 
         !  increment.
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for plane-strain or axisymmetric
         !
         if(pe.ne.1) then
            !
            ! this is axisymmetric
            !
            F_tau(3,3) = AR/AR0
            F_t(3,3) = AR_t/AR0
         else
            !
            ! this is plane strain
            !
            F_tau(3,3) = one
            F_t(3,3) = one
         endif
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            if(pe.eq.1) then
               !
               !  2D plane-strain implementation
               !
               detF_tau = F_tau(1,1)*F_tau(2,2) - F_tau(1,2)*F_tau(2,1)
               detF_t = F_t(1,1)*F_t(2,2) - F_t(1,2)*F_t(2,1)
               do i=1,nDim
                  do j=1,nDim
                     F_tau(i,j) =((detFc_tau/detF_tau)**half)*F_tau(i,j)
                     F_t(i,j) = ((detFc_t/detF_t)**half)*F_t(i,j)
                  enddo
               enddo
            else
               !
               ! 2D axisymmetric implementation
               !
               call mdet(F_tau,detF_tau)
               call mdet(F_t,detF_t)
               F_tau = ((detFc_tau/detF_tau)**third)*F_tau
               F_t = ((detFc_tau/detF_tau)**third)*F_t
            endif
         endif


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive update at this integ. point
         !
         if (matflag.eq.1) then
            call NeoHookean(props,nprops,dtime,F_tau,Fv_t,
     +           T_tau,Fv_tau,SpTanMod,stat)
         elseif (matflag.eq.2) then
            call Gent(props,nprops,dtime,F_tau,Fv_t,
     +           T_tau,Fv_tau,SpTanMod,stat)
         else
            write(*,*) 'Invalid matflag: matflag.ne.1 or 2'
            call xit
         endif
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment.
         !
         svars(1+jj)  = Fv_tau(1,1)
         svars(2+jj)  = Fv_tau(1,2)
         svars(3+jj)  = Fv_tau(1,3)
         svars(4+jj)  = Fv_tau(2,1)
         svars(5+jj)  = Fv_tau(2,2)
         svars(6+jj)  = Fv_tau(2,3)
         svars(7+jj)  = Fv_tau(3,1)
         svars(8+jj)  = Fv_tau(3,2)
         svars(9+jj)  = Fv_tau(3,3)
         jj = jj + nSdv ! setup for the next intPt


         ! Compute/update the displacement residual vector
         !
         if(pe.eq.1) then
            !
            ! this is plane strain
            !
            Smat(1,1) = T_tau(1,1)
            Smat(2,1) = T_tau(2,2)
            Smat(3,1) = T_tau(1,2)
            !
            Bmat = zero
            do k=1,nNode
               Bmat(1,1+nDim*(k-1)) = dshC(k,1)
               Bmat(2,2+nDim*(k-1)) = dshC(k,2)
               Bmat(3,1+nDim*(k-1)) = dshC(k,2)
               Bmat(3,2+nDim*(k-1)) = dshC(k,1)
            enddo
            !
            body = zero ! The body force vector may be specified here
            !
            BodyForceRes = zero
            do k=1,nNode
               BodyForceRes(1+nDim*(k-1),1) = sh(k)*body(1)
               BodyForceRes(2+nDim*(k-1),1) = sh(k)*body(2)
            enddo
            !
            Ru = Ru 
     +           - matmul(transpose(Bmat),Smat)*detmapJC*w(intpt)*AR
     +           + BodyForceRes*detmapJC*w(intpt)*AR
         else
            !
            ! this is axisymetric
            !
            SmatAx(1,1) = T_tau(1,1)
            SmatAx(2,1) = T_tau(2,2)
            SmatAx(3,1) = T_tau(1,2)
            SmatAx(4,1) = T_tau(3,3)
            !
            BmatAx = zero
            do k=1,nNode
               BmatAx(1,1+nDim*(k-1)) = dshC(k,1)
               BmatAx(2,2+nDim*(k-1)) = dshC(k,2)
               BmatAx(3,1+nDim*(k-1)) = dshC(k,2)
               BmatAx(3,2+nDim*(k-1)) = dshC(k,1)
               BmatAx(4,1+nDim*(k-1)) = sh(k)/(AR/(two*Pi))
            enddo
            !
            body = zero ! The body force vector may be specified here
            !
            BodyForceResAx = zero
            do k=1,nNode
               BodyForceResAx(1+nDim*(k-1),1) = sh(k)*body(1)
               BodyForceResAx(2+nDim*(k-1),1) = sh(k)*body(2)
            enddo
            !
            Ru = Ru
     +          - matmul(transpose(BmatAx),SmatAx)*detmapJC*w(intpt)*AR
     +          + BodyForceResAx*detmapJC*w(intpt)*AR
         endif


         ! Compute/update the displacement tangent matrix
         !
         if(pe.eq.1) then
            !
            ! this is plane strain
            !
            Gmat = zero
            do k=1,nNode
               Gmat(1,1+nDim*(k-1)) = dshC(k,1)
               Gmat(2,2+nDim*(k-1)) = dshC(k,1)
               Gmat(3,1+nDim*(k-1)) = dshC(k,2)
               Gmat(4,2+nDim*(k-1)) = dshC(k,2)
            enddo

            G0mat = zero
            do k=1,nNode
               G0mat(1,1+nDim*(k-1)) = dshC0(k,1)
               G0mat(2,2+nDim*(k-1)) = dshC0(k,1)
               G0mat(3,1+nDim*(k-1)) = dshC0(k,2)
               G0mat(4,2+nDim*(k-1)) = dshC0(k,2)
            enddo
            !
            Amat = zero
            Amat(1,1) = SpTanMod(1,1,1,1)
            Amat(1,2) = SpTanMod(1,1,2,1)
            Amat(1,3) = SpTanMod(1,1,1,2)
            Amat(1,4) = SpTanMod(1,1,2,2)
            Amat(2,1) = SpTanMod(2,1,1,1)
            Amat(2,2) = SpTanMod(2,1,2,1)
            Amat(2,3) = SpTanMod(2,1,1,2)
            Amat(2,4) = SpTanMod(2,1,2,2)
            Amat(3,1) = SpTanMod(1,2,1,1)
            Amat(3,2) = SpTanMod(1,2,2,1)
            Amat(3,3) = SpTanMod(1,2,1,2)
            Amat(3,4) = SpTanMod(1,2,2,2)
            Amat(4,1) = SpTanMod(2,2,1,1)
            Amat(4,2) = SpTanMod(2,2,2,1)
            Amat(4,3) = SpTanMod(2,2,1,2)
            Amat(4,4) = SpTanMod(2,2,2,2)
            !
            Qmat = zero
            Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
            Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
            Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
            Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
            Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
            Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
            Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            !
            if((nNode.eq.4).and.(nInt.eq.4)) then
               !
               ! This is the tangent using the F-bar method with the
               !  4 node fully integrated linear element
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(Gmat),Amat),Gmat)
     +              *detMapJC*w(intpt)*AR
     +              +matmul(transpose(Gmat),matmul(Qmat,(G0mat - Gmat)))
     +              *detMapJC*w(intpt)*AR
            else
               !
               ! This is the tangent not using the F-bar method with all
               !  other elements
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(Gmat),Amat),Gmat)
     +              *detMapJC*w(intpt)*AR
            endif
            !
         else
            !
            ! this is axisymmetric
            !
            GmatAx = zero
            do k=1,nNode
               GmatAx(1,1+nDim*(k-1)) = dshC(k,1)
               GmatAx(2,2+nDim*(k-1)) = dshC(k,1)
               GmatAx(3,1+nDim*(k-1)) = dshC(k,2)
               GmatAx(4,2+nDim*(k-1)) = dshC(k,2)
               GmatAx(5,1+nDim*(k-1)) = sh(k)/(AR/(two*Pi))
            enddo
            !
            G0matAx = zero
            do k=1,nNode
               G0matAx(1,1+nDim*(k-1)) = dshC0(k,1)
               G0matAx(2,2+nDim*(k-1)) = dshC0(k,1)
               G0matAx(3,1+nDim*(k-1)) = dshC0(k,2)
               G0matAx(4,2+nDim*(k-1)) = dshC0(k,2)
               G0matAx(5,1+nDim*(k-1)) = sh0(k)/ARc
            enddo
            !
            AmatAx = zero
            AmatAx(1,1) = SpTanMod(1,1,1,1)
            AmatAx(1,2) = SpTanMod(1,1,2,1)
            AmatAx(1,3) = SpTanMod(1,1,1,2)
            AmatAx(1,4) = SpTanMod(1,1,2,2)
            AmatAx(1,5) = SpTanMod(1,1,3,3)
            AmatAx(2,1) = SpTanMod(2,1,1,1)
            AmatAx(2,2) = SpTanMod(2,1,2,1)
            AmatAx(2,3) = SpTanMod(2,1,1,2)
            AmatAx(2,4) = SpTanMod(2,1,2,2)
            AmatAx(2,5) = SpTanMod(2,1,3,3)
            AmatAx(3,1) = SpTanMod(1,2,1,1)
            AmatAx(3,2) = SpTanMod(1,2,2,1)
            AmatAx(3,3) = SpTanMod(1,2,1,2)
            AmatAx(3,4) = SpTanMod(1,2,2,2)
            AmatAx(3,5) = SpTanMod(1,2,3,3)
            AmatAx(4,1) = SpTanMod(2,2,1,1)
            AmatAx(4,2) = SpTanMod(2,2,2,1)
            AmatAx(4,3) = SpTanMod(2,2,1,2)
            AmatAx(4,4) = SpTanMod(2,2,2,2)
            AmatAx(4,5) = SpTanMod(2,2,3,3)
            AmatAx(5,1) = SpTanMod(3,3,1,1)
            AmatAx(5,2) = SpTanMod(3,3,2,1)
            AmatAx(5,3) = SpTanMod(3,3,1,2)
            AmatAx(5,4) = SpTanMod(3,3,2,2)
            AmatAx(5,5) = SpTanMod(3,3,3,3)
            !
            QmatAx = zero
            QmatAx(1,1) = third*(AmatAx(1,1)+AmatAx(1,4)+AmatAx(1,5)) 
     +           - (two/three)*T_tau(1,1)
            QmatAx(2,1) = third*(AmatAx(2,1)+AmatAx(2,4)+AmatAx(2,5))
     +           - (two/three)*T_tau(1,2)
            QmatAx(3,1) = third*(AmatAx(3,1)+AmatAx(3,4)+AmatAx(3,5))
     +           - (two/three)*T_tau(1,2)
            QmatAx(4,1) = third*(AmatAx(4,1)+AmatAx(4,4)+AmatAx(4,5))
     +           - (two/three)*T_tau(2,2)
            QmatAx(5,1) = third*(AmatAx(5,1)+AmatAx(5,4)+AmatAx(5,5))
     +           - (two/three)*T_tau(3,3)
            QmatAx(1,4) = QmatAx(1,1)
            QmatAx(2,4) = QmatAx(2,1)
            QmatAx(3,4) = QmatAx(3,1)
            QmatAx(4,4) = QmatAx(4,1)
            QmatAx(5,4) = QmatAx(5,1)
            QmatAx(1,5) = QmatAx(1,1)
            QmatAx(2,5) = QmatAx(2,1)
            QmatAx(3,5) = QmatAx(3,1)
            QmatAx(4,5) = QmatAx(4,1)
            QmatAx(5,5) = QmatAx(5,1)
            !
            if((nNode.eq.4).and.(nInt.eq.4)) then
               !
               ! This is the tangent using the F-bar method with the
               !  4 node fully integrated linear element
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +              *detMapJC*w(intpt)*AR
     +              + matmul(transpose(GmatAx),matmul(QmatAx,
     +              (G0matAx-GmatAx)))*detMapJC*w(intpt)*AR
            else
               !
               ! This is the tangent NOT using the F-bar method with all
               !  other elements
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +              *detMapJC*w(intpt)*AR
            endif
            !
         endif

      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface tension terms
      !
      ! Based on our convention, the face on which the surface tension
      !  acts is the "label", i.e.
      !  - U1,U2,U3,U4 refer to surface tension acting
      !    on faces 1,2,3,4, respectively,
      !  Mechanical traction- and pressure-type boundary conditions 
      !  may be applied to the dummy mesh using the Abaqus built-in 
      !  commands *Dload or *Dsload.
      !
      if (pe.eq.1) then
        !
        ! This is plane strain
        !
        if(ndload.gt.0) then
          !
          ! Loop over faces and make proper modifications to
          !  residuals and tangents as needed.
          !
          do i=1,ndload
            !
            face = jdltyp(i,1)  ! face label
            gamma = adlmag(i,1) ! surface tension
            !
            if(face.eq.1) then
               !
               ! surface tension on face 1 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                    ((coordsC(2,1)-coordsC(2,2))**two))
               Ru(1,1) = Ru(1,1) - 
     +                        gamma*AR*(coordsC(1,1)-coordsC(1,2))/Le
               Ru(2,1) = Ru(2,1) - 
     +                        gamma*AR*(coordsC(2,1)-coordsC(2,2))/Le
               Ru(3,1) = Ru(3,1) - 
     +                        gamma*AR*(coordsC(1,2)-coordsC(1,1))/Le
               Ru(4,1) = Ru(4,1) - 
     +                        gamma*AR*(coordsC(2,2)-coordsC(2,1))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(1,1) = Kuu(1,1) + gamma*AR/Le - 
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(1,1)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(1,2) = Kuu(1,2) - 
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(2,1)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(1,3) = Kuu(1,3) - gamma*AR/Le - 
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(1,4) = Kuu(1,4) - 
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(2,1) = Kuu(2,1) - 
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(1,1)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(2,2) = Kuu(2,2) + gamma*AR/Le -
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(2,1)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(2,3) = Kuu(2,3) -
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(1,2)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(2,4) = Kuu(2,4) - gamma*AR/Le -
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(2,2)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(3,1) = Kuu(3,1) - gamma*AR/Le - 
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(3,2) = Kuu(3,2) - 
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(3,3) = Kuu(3,3) + gamma*AR/Le - 
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(1,2)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(3,4) = Kuu(3,4) - 
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(2,2)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(4,1) = Kuu(4,1) -
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(1,1)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(4,2) = Kuu(4,2) - gamma*AR/Le -
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(2,1)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(4,3) = Kuu(4,3) -
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(1,2)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(4,4) = Kuu(4,4) + gamma*AR/Le -
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(2,2)-coordsC(2,1))*gamma*AR/(Le**three)
               !
            elseif(face.eq.2) then
               !
               ! surface tension on face 2 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                    ((coordsC(2,2)-coordsC(2,3))**two))
               Ru(3,1) = Ru(3,1) - 
     +                        gamma*AR*(coordsC(1,2)-coordsC(1,3))/Le
               Ru(4,1) = Ru(4,1) - 
     +                        gamma*AR*(coordsC(2,2)-coordsC(2,3))/Le
               Ru(5,1) = Ru(5,1) - 
     +                        gamma*AR*(coordsC(1,3)-coordsC(1,2))/Le
               Ru(6,1) = Ru(6,1) - 
     +                        gamma*AR*(coordsC(2,3)-coordsC(2,2))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(3,3) = Kuu(3,3) + gamma*AR/Le - 
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(1,2)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(3,4) = Kuu(3,4) - 
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(2,2)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(3,5) = Kuu(3,5) - gamma*AR/Le - 
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(3,6) = Kuu(3,6) - 
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(4,3) = Kuu(4,3) -
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(1,2)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(4,4) = Kuu(4,4) + gamma*AR/Le -
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(2,2)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(4,5) = Kuu(4,5) -
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(1,3)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(4,6) = Kuu(4,6) - gamma*AR/Le -
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(2,3)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(5,3) = Kuu(5,3) - gamma*AR/Le - 
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(5,4) = Kuu(5,4) - 
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(5,5) = Kuu(5,5) + gamma*AR/Le - 
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(1,3)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(5,6) = Kuu(5,6) - 
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(2,3)-coordsC(2,2))*gamma*AR/(Le**three)
               Kuu(6,3) = Kuu(6,3) -
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(1,2)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(6,4) = Kuu(6,4) - gamma*AR/Le -
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(2,2)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(6,5) = Kuu(6,5) -
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(1,3)-coordsC(1,2))*gamma*AR/(Le**three)
               Kuu(6,6) = Kuu(6,6) + gamma*AR/Le -
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(2,3)-coordsC(2,2))*gamma*AR/(Le**three)
               !
            elseif(face.eq.3) then
               !
               ! surface tension on face 3 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                    ((coordsC(2,3)-coordsC(2,4))**two))
               Ru(5,1) = Ru(5,1) - 
     +                        gamma*AR*(coordsC(1,3)-coordsC(1,4))/Le
               Ru(6,1) = Ru(6,1) - 
     +                        gamma*AR*(coordsC(2,3)-coordsC(2,4))/Le
               Ru(7,1) = Ru(7,1) - 
     +                        gamma*AR*(coordsC(1,4)-coordsC(1,3))/Le
               Ru(8,1) = Ru(8,1) - 
     +                        gamma*AR*(coordsC(2,4)-coordsC(2,3))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(5,5) = Kuu(5,5) + gamma*AR/Le - 
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(1,3)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(5,6) = Kuu(5,6) - 
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(2,3)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(5,7) = Kuu(5,7) - gamma*AR/Le - 
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(5,8) = Kuu(5,8) - 
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(6,5) = Kuu(6,5) -
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(1,3)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(6,6) = Kuu(6,6) + gamma*AR/Le -
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(2,3)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(6,7) = Kuu(6,7) -
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(1,4)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(6,8) = Kuu(6,8) - gamma*AR/Le -
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(2,4)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(7,5) = Kuu(7,5) - gamma*AR/Le - 
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(7,6) = Kuu(7,6) - 
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(7,7) = Kuu(7,7) + gamma*AR/Le - 
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(1,4)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(7,8) = Kuu(7,8) - 
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(2,4)-coordsC(2,3))*gamma*AR/(Le**three)
               Kuu(8,5) = Kuu(8,5) -
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(1,3)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(8,6) = Kuu(8,6) - gamma*AR/Le -
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(2,3)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(8,7) = Kuu(8,7) -
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(1,4)-coordsC(1,3))*gamma*AR/(Le**three)
               Kuu(8,8) = Kuu(8,8) + gamma*AR/Le -
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(2,4)-coordsC(2,3))*gamma*AR/(Le**three)
               !
            elseif(face.eq.4) then
               !
               ! surface tension on face 4 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                    ((coordsC(2,4)-coordsC(2,1))**two))
               Ru(7,1) = Ru(7,1) - 
     +                        gamma*AR*(coordsC(1,4)-coordsC(1,1))/Le
               Ru(8,1) = Ru(8,1) - 
     +                        gamma*AR*(coordsC(2,4)-coordsC(2,1))/Le
               Ru(1,1) = Ru(1,1) - 
     +                        gamma*AR*(coordsC(1,1)-coordsC(1,4))/Le
               Ru(2,1) = Ru(2,1) - 
     +                        gamma*AR*(coordsC(2,1)-coordsC(2,4))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(7,7) = Kuu(7,7) + gamma*AR/Le - 
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(1,4)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(7,8) = Kuu(7,8) - 
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(2,4)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(7,1) = Kuu(7,1) - gamma*AR/Le - 
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(7,2) = Kuu(7,2) - 
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(8,7) = Kuu(8,7) -
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(1,4)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(8,8) = Kuu(8,8) + gamma*AR/Le -
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(2,4)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(8,1) = Kuu(8,1) -
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(1,1)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(8,2) = Kuu(8,2) - gamma*AR/Le -
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(2,1)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(1,7) = Kuu(1,7) - gamma*AR/Le - 
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(1,8) = Kuu(1,8) - 
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(1,1) = Kuu(1,1) + gamma*AR/Le - 
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(1,1)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(1,2) = Kuu(1,2) - 
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(2,1)-coordsC(2,4))*gamma*AR/(Le**three)
               Kuu(2,7) = Kuu(2,7) -
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(1,4)-coordsC(1,1))*gamma*AR/(Le**three)
               Kuu(2,8) = Kuu(2,8) - gamma*AR/Le -
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(2,4)-coordsC(2,1))*gamma*AR/(Le**three)
               Kuu(2,1) = Kuu(2,1) -
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(1,1)-coordsC(1,4))*gamma*AR/(Le**three)
               Kuu(2,2) = Kuu(2,2) + gamma*AR/Le -
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(2,1)-coordsC(2,4))*gamma*AR/(Le**three)
               !
            else
               write(*,*) 'Incorrect dload type',face
               call xit
            endif
            !
          end do
          !
        endif
        !
      else
        !
        ! This is axisymmetric
        !
        if(ndload.gt.0) then
          !
          ! Loop over faces and make proper modifications to
          !  residuals and tangents as needed.
          !
          do i=1,ndload
            !
            face = jdltyp(i,1)  ! face label
            gamma = adlmag(i,1) ! surface tension
            !
            if(face.eq.1) then
               !
               ! surface tension on face 1 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                    ((coordsC(2,1)-coordsC(2,2))**two))
               Ru(1,1) = Ru(1,1) - pi*gamma*(coordsC(1,1)+coordsC(1,2))
     +                       *(coordsC(1,1)-coordsC(1,2))/Le -
     +                       pi*gamma*Le
               Ru(2,1) = Ru(2,1) - pi*gamma*(coordsC(1,1)+coordsC(1,2))
     +                       *(coordsC(2,1)-coordsC(2,2))/Le
               Ru(3,1) = Ru(3,1) - pi*gamma*(coordsC(1,1)+coordsC(1,2))
     +                       *(coordsC(1,2)-coordsC(1,1))/Le -
     +                       pi*gamma*Le
               Ru(4,1) = Ru(4,1) - pi*gamma*(coordsC(1,1)+coordsC(1,2))
     +                       *(coordsC(2,2)-coordsC(2,1))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(1,1) = Kuu(1,1) + two*coordsC(1,1)*pi*gamma/Le - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(1,1)-coordsC(1,2))*pi*gamma/(Le**three) +
     +               (coordsC(1,1)-coordsC(1,2))*pi*gamma/Le
               Kuu(1,2) = Kuu(1,2) - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/(Le**three) +
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/Le
               Kuu(1,3) = Kuu(1,3) - two*coordsC(1,2)*pi*gamma/Le - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,1))*pi*gamma/(Le**three) +
     +               (coordsC(1,2)-coordsC(1,1))*pi*gamma/Le
               Kuu(1,4) = Kuu(1,4) - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,1)-coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,1))*pi*gamma/(Le**three) +
     +               (coordsC(2,2)-coordsC(2,1))*pi*gamma/Le
               Kuu(2,1) = Kuu(2,1) + 
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(1,1)-coordsC(1,2))*pi*gamma/(Le**three)
               Kuu(2,2) = Kuu(2,2) + 
     +               (coordsC(1,1)+coordsC(1,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/(Le**three)
               Kuu(2,3) = Kuu(2,3) + 
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(1,2)-coordsC(1,1))*pi*gamma/(Le**three)
               Kuu(2,4) = Kuu(2,4) - 
     +               (coordsC(1,1)+coordsC(1,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,1)-coordsC(2,2))*
     +               (coordsC(2,2)-coordsC(2,1))*pi*gamma/(Le**three)
               Kuu(3,1) = Kuu(3,1) - two*coordsC(1,1)*pi*gamma/Le - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,2))*pi*gamma/(Le**three) +
     +               (coordsC(1,1)-coordsC(1,2))*pi*gamma/Le
               Kuu(3,2) = Kuu(3,2) - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/(Le**three) +
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/Le
               Kuu(3,3) = Kuu(3,3) + two*coordsC(1,2)*pi*gamma/Le - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(1,2)-coordsC(1,1))*pi*gamma/(Le**three) +
     +               (coordsC(1,2)-coordsC(1,1))*pi*gamma/Le
               Kuu(3,4) = Kuu(3,4) - 
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,1))*
     +               (coordsC(2,2)-coordsC(2,1))*pi*gamma/(Le**three) +
     +               (coordsC(2,2)-coordsC(2,1))*pi*gamma/Le
               Kuu(4,1) = Kuu(4,1) - 
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(1,1)-coordsC(1,2))*pi*gamma/(Le**three)
               Kuu(4,2) = Kuu(4,2) - 
     +               (coordsC(1,1)+coordsC(1,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/(Le**three)
               Kuu(4,3) = Kuu(4,3) - 
     +               (coordsC(2,1)-coordsC(2,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(1,2)-coordsC(1,1))*pi*gamma/(Le**three)
               Kuu(4,4) = Kuu(4,4) + 
     +               (coordsC(1,1)+coordsC(1,2))*pi*gamma/Le -
     +               (coordsC(1,1)+coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,1))*
     +               (coordsC(2,2)-coordsC(2,1))*pi*gamma/(Le**three)
               !
            elseif(face.eq.2) then
               !
               ! surface tension on face 2 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                    ((coordsC(2,2)-coordsC(2,3))**two))
               Ru(3,1) = Ru(3,1) - pi*gamma*(coordsC(1,2)+coordsC(1,3))
     +                       *(coordsC(1,2)-coordsC(1,3))/Le -
     +                       pi*gamma*Le
               Ru(4,1) = Ru(4,1) - pi*gamma*(coordsC(1,2)+coordsC(1,3))
     +                       *(coordsC(2,2)-coordsC(2,3))/Le
               Ru(5,1) = Ru(5,1) - pi*gamma*(coordsC(1,2)+coordsC(1,3))
     +                       *(coordsC(1,3)-coordsC(1,2))/Le -
     +                       pi*gamma*Le
               Ru(6,1) = Ru(6,1) - pi*gamma*(coordsC(1,2)+coordsC(1,3))
     +                       *(coordsC(2,3)-coordsC(2,2))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(3,3) = Kuu(3,3) + two*coordsC(1,2)*pi*gamma/Le - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(1,2)-coordsC(1,3))*pi*gamma/(Le**three) +
     +               (coordsC(1,2)-coordsC(1,3))*pi*gamma/Le
               Kuu(3,4) = Kuu(3,4) - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/(Le**three) +
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/Le
               Kuu(3,5) = Kuu(3,5) - two*coordsC(1,3)*pi*gamma/Le - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,2))*pi*gamma/(Le**three) +
     +               (coordsC(1,3)-coordsC(1,2))*pi*gamma/Le
               Kuu(3,6) = Kuu(3,6) - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,2)-coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,2))*pi*gamma/(Le**three) +
     +               (coordsC(2,3)-coordsC(2,2))*pi*gamma/Le
               Kuu(4,3) = Kuu(4,3) + 
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(1,2)-coordsC(1,3))*pi*gamma/(Le**three)
               Kuu(4,4) = Kuu(4,4) + 
     +               (coordsC(1,2)+coordsC(1,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/(Le**three)
               Kuu(4,5) = Kuu(4,5) + 
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(1,3)-coordsC(1,2))*pi*gamma/(Le**three)
               Kuu(4,6) = Kuu(4,6) - 
     +               (coordsC(1,2)+coordsC(1,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,2)-coordsC(2,3))*
     +               (coordsC(2,3)-coordsC(2,2))*pi*gamma/(Le**three)
               Kuu(5,3) = Kuu(5,3) - two*coordsC(1,2)*pi*gamma/Le - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(1,2)-coordsC(1,3))*pi*gamma/(Le**three) +
     +               (coordsC(1,2)-coordsC(1,3))*pi*gamma/Le
               Kuu(5,4) = Kuu(5,4) - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/(Le**three) +
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/Le
               Kuu(5,5) = Kuu(5,5) + two*coordsC(1,3)*pi*gamma/Le - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(1,3)-coordsC(1,2))*pi*gamma/(Le**three) +
     +               (coordsC(1,3)-coordsC(1,2))*pi*gamma/Le
               Kuu(5,6) = Kuu(5,6) - 
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,2))*
     +               (coordsC(2,3)-coordsC(2,2))*pi*gamma/(Le**three) +
     +               (coordsC(2,3)-coordsC(2,2))*pi*gamma/Le
               Kuu(6,3) = Kuu(6,3) - 
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(1,2)-coordsC(1,3))*pi*gamma/(Le**three)
               Kuu(6,4) = Kuu(6,4) - 
     +               (coordsC(1,2)+coordsC(1,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/(Le**three)
               Kuu(6,5) = Kuu(6,5) - 
     +               (coordsC(2,2)-coordsC(2,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(1,3)-coordsC(1,2))*pi*gamma/(Le**three)
               Kuu(6,6) = Kuu(6,6) + 
     +               (coordsC(1,2)+coordsC(1,3))*pi*gamma/Le -
     +               (coordsC(1,2)+coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,2))*
     +               (coordsC(2,3)-coordsC(2,2))*pi*gamma/(Le**three)
               !
            elseif(face.eq.3) then
               !
               ! surface tension on face 3 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                    ((coordsC(2,3)-coordsC(2,4))**two))
               Ru(5,1) = Ru(5,1) - pi*gamma*(coordsC(1,3)+coordsC(1,4))
     +                       *(coordsC(1,3)-coordsC(1,4))/Le -
     +                       pi*gamma*Le
               Ru(6,1) = Ru(6,1) - pi*gamma*(coordsC(1,3)+coordsC(1,4))
     +                       *(coordsC(2,3)-coordsC(2,4))/Le
               Ru(7,1) = Ru(7,1) - pi*gamma*(coordsC(1,3)+coordsC(1,4))
     +                       *(coordsC(1,4)-coordsC(1,3))/Le -
     +                       pi*gamma*Le
               Ru(8,1) = Ru(8,1) - pi*gamma*(coordsC(1,3)+coordsC(1,4))
     +                       *(coordsC(2,4)-coordsC(2,3))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(5,5) = Kuu(5,5) + two*coordsC(1,3)*pi*gamma/Le - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(1,3)-coordsC(1,4))*pi*gamma/(Le**three) +
     +               (coordsC(1,3)-coordsC(1,4))*pi*gamma/Le
               Kuu(5,6) = Kuu(5,6) - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/(Le**three) +
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/Le
               Kuu(5,7) = Kuu(5,7) - two*coordsC(1,4)*pi*gamma/Le - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,3))*pi*gamma/(Le**three) +
     +               (coordsC(1,4)-coordsC(1,3))*pi*gamma/Le
               Kuu(5,8) = Kuu(5,8) - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,3)-coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,3))*pi*gamma/(Le**three) +
     +               (coordsC(2,4)-coordsC(2,3))*pi*gamma/Le
               Kuu(6,5) = Kuu(6,5) + 
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(1,3)-coordsC(1,4))*pi*gamma/(Le**three)
               Kuu(6,6) = Kuu(6,6) + 
     +               (coordsC(1,3)+coordsC(1,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/(Le**three)
               Kuu(6,7) = Kuu(6,7) + 
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(1,4)-coordsC(1,3))*pi*gamma/(Le**three)
               Kuu(6,8) = Kuu(6,8) - 
     +               (coordsC(1,3)+coordsC(1,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,3)-coordsC(2,4))*
     +               (coordsC(2,4)-coordsC(2,3))*pi*gamma/(Le**three)
               Kuu(7,5) = Kuu(7,5) - two*coordsC(1,3)*pi*gamma/Le - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(1,3)-coordsC(1,4))*pi*gamma/(Le**three) +
     +               (coordsC(1,3)-coordsC(1,4))*pi*gamma/Le
               Kuu(7,6) = Kuu(7,6) - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/(Le**three) +
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/Le
               Kuu(7,7) = Kuu(7,7) + two*coordsC(1,4)*pi*gamma/Le - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(1,4)-coordsC(1,3))*pi*gamma/(Le**three) +
     +               (coordsC(1,4)-coordsC(1,3))*pi*gamma/Le
               Kuu(7,8) = Kuu(7,8) - 
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,3))*
     +               (coordsC(2,4)-coordsC(2,3))*pi*gamma/(Le**three) +
     +               (coordsC(2,4)-coordsC(2,3))*pi*gamma/Le
               Kuu(8,5) = Kuu(8,5) - 
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(1,3)-coordsC(1,4))*pi*gamma/(Le**three)
               Kuu(8,6) = Kuu(8,6) - 
     +               (coordsC(1,3)+coordsC(1,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/(Le**three)
               Kuu(8,7) = Kuu(8,7) - 
     +               (coordsC(2,3)-coordsC(2,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(1,4)-coordsC(1,3))*pi*gamma/(Le**three)
               Kuu(8,8) = Kuu(8,8) + 
     +               (coordsC(1,3)+coordsC(1,4))*pi*gamma/Le -
     +               (coordsC(1,3)+coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,3))*
     +               (coordsC(2,4)-coordsC(2,3))*pi*gamma/(Le**three)
               
            elseif(face.eq.4) then
               !
               ! surface tension on face 4 of the element
               !
               ! Modify the displacement residual
               !
               Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                    ((coordsC(2,4)-coordsC(2,1))**two))
               Ru(7,1) = Ru(7,1) - pi*gamma*(coordsC(1,4)+coordsC(1,1))
     +                       *(coordsC(1,4)-coordsC(1,1))/Le -
     +                       pi*gamma*Le
               Ru(8,1) = Ru(8,1) - pi*gamma*(coordsC(1,4)+coordsC(1,1))
     +                       *(coordsC(2,4)-coordsC(2,1))/Le
               Ru(1,1) = Ru(1,1) - pi*gamma*(coordsC(1,4)+coordsC(1,1))
     +                       *(coordsC(1,1)-coordsC(1,4))/Le -
     +                       pi*gamma*Le
               Ru(2,1) = Ru(2,1) - pi*gamma*(coordsC(1,4)+coordsC(1,1))
     +                       *(coordsC(2,1)-coordsC(2,4))/Le
               !
               ! Modify the tangent matrix
               !
               Kuu(7,7) = Kuu(7,7) + two*coordsC(1,4)*pi*gamma/Le - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(1,4)-coordsC(1,1))*pi*gamma/(Le**three) +
     +               (coordsC(1,4)-coordsC(1,1))*pi*gamma/Le
               Kuu(7,8) = Kuu(7,8) - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/(Le**three) +
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/Le
               Kuu(7,1) = Kuu(7,1) - two*coordsC(1,1)*pi*gamma/Le - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,4))*pi*gamma/(Le**three) +
     +               (coordsC(1,1)-coordsC(1,4))*pi*gamma/Le
               Kuu(7,2) = Kuu(7,2) - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,4)-coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,4))*pi*gamma/(Le**three) +
     +               (coordsC(2,1)-coordsC(2,4))*pi*gamma/Le
               Kuu(8,7) = Kuu(8,7) + 
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(1,4)-coordsC(1,1))*pi*gamma/(Le**three)
               Kuu(8,8) = Kuu(8,8) + 
     +               (coordsC(1,4)+coordsC(1,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/(Le**three)
               Kuu(8,1) = Kuu(8,1) + 
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(1,1)-coordsC(1,4))*pi*gamma/(Le**three)
               Kuu(8,2) = Kuu(8,2) - 
     +               (coordsC(1,4)+coordsC(1,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,4)-coordsC(2,1))*
     +               (coordsC(2,1)-coordsC(2,4))*pi*gamma/(Le**three)
               Kuu(1,7) = Kuu(1,7) - two*coordsC(1,4)*pi*gamma/Le - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(1,4)-coordsC(1,1))*pi*gamma/(Le**three) +
     +               (coordsC(1,4)-coordsC(1,1))*pi*gamma/Le
               Kuu(1,8) = Kuu(1,8) - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/(Le**three) +
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/Le
               Kuu(1,1) = Kuu(1,1) + two*coordsC(1,1)*pi*gamma/Le - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(1,1)-coordsC(1,4))*pi*gamma/(Le**three) +
     +               (coordsC(1,1)-coordsC(1,4))*pi*gamma/Le
               Kuu(1,2) = Kuu(1,2) - 
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(1,1)-coordsC(1,4))*
     +               (coordsC(2,1)-coordsC(2,4))*pi*gamma/(Le**three) +
     +               (coordsC(2,1)-coordsC(2,4))*pi*gamma/Le
               Kuu(2,7) = Kuu(2,7) - 
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(1,4)-coordsC(1,1))*pi*gamma/(Le**three)
               Kuu(2,8) = Kuu(2,8) - 
     +               (coordsC(1,4)+coordsC(1,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/(Le**three)
               Kuu(2,1) = Kuu(2,1) - 
     +               (coordsC(2,4)-coordsC(2,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(1,1)-coordsC(1,4))*pi*gamma/(Le**three)
               Kuu(2,2) = Kuu(2,2) + 
     +               (coordsC(1,4)+coordsC(1,1))*pi*gamma/Le -
     +               (coordsC(1,4)+coordsC(1,1))*
     +               (coordsC(2,1)-coordsC(2,4))*
     +               (coordsC(2,1)-coordsC(2,4))*pi*gamma/(Le**three)
               !
            else
               write(*,*) 'Incorrect dload type',face
               call xit
            endif
            !
          end do
          !
        endif
        !
      endif
      !
      ! End loop over surface tension terms
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return Abaqus the right hand side vector
      !
      do i=1,nNode
         A11 = (nDim+2)*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A12,1) = Ru(A12,1)
         rhs(A12+1,1) = Ru(A12+1,1)
      enddo
      !
      ! Return Abaqus the tangent matrix
      !
      amatrx = Kuu
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine uel

!************************************************************************
!     Material subroutines
!************************************************************************

      subroutine NeoHookean(props,nprops,dtime,F_tau,Fv_t,
     +        T_tau,Fv_tau,SpTanMod,stat)

      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fv_t(3,3),T_tau(3,3),
     +  Fv_tau(3,3),SpTanMod(3,3,3,3),Iden(3,3),G0,Kbulk,eta,
     +  Gneq,detF,Finv(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),
     +  I1_tau,I1bar,TR_tau(3,3),dTRdF(3,3,3,3),Fv_t_inv(3,3),
     +  det_Fv_t,Fe_tr(3,3),Re_tr(3,3),Ue_tr(3,3),Ee_tr(3,3),
     +  trEe_tr,Ee0_tr(3,3),Me_tr(3,3),tauBar_tr,Nv_tau(3,3),
     +  nuv_tau,tauBar_tau,Dv_tau(3,3),Dv_eig(3),Dv_vec(3,3),
     +  expdtDv(3,3),tmp,Me_tau(3,3),Gshear_tilde,Lambda_tilde,
     +  fac,c3
      !
      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0)
 

      ! Identity matrix
      !
      call onem(Iden)


      ! Obtain material properties
      !
      G0     = props(1) ! Ground-state shear modulus
      Kbulk  = props(2) ! Bulk modulus
      Gneq   = props(4) ! Nonequilibrium (Maxwell element) stiffness
      eta    = props(5) ! Maxwell element viscosity
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Neo-Hookean spring
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT = transpose(F_tau)
      FinvT = transpose(Finv)


      ! Compute the right Cauchy-Green tensor
      !  and the first stretch invariant
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      I1_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 

      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*I1_tau
 

      ! Compute the 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*G0*(F_tau-third*I1bar*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
 

      ! Compute the Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))


      ! Calculate the material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*G0*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*I1bar*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*I1bar*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the spatial tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) +
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Maxwell element
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the trial elastic deformation gradient
      !
      call matInv3D(Fv_t,Fv_t_inv,det_Fv_t,stat)
      Fe_tr = matmul(F_tau,Fv_t_inv)


      ! Compute the trial kinematics
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr,stat)
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3)
      Ee0_tr = Ee_tr - (one/three)*trEe_tr*Iden


      ! Compute the trial Mandel stress, which is deviatoric
      !
      Me_tr = two*Gneq*Ee0_tr


      ! Compute the trial equiv. tensile stress
      !
      tauBar_tr = dsqrt(one/two)*dsqrt(sum(Me_tr*Me_tr))


      ! Compute the direction of viscous flow
      !
      if(tauBar_tr.le.zero) then
         Nv_tau = zero
      else
         Nv_tau = dsqrt(one/two)*(Me_tr/tauBar_tr)
      endif


      ! Compute the equivalent shear viscous strain rate
      !
      nuv_tau = tauBar_tr/(eta + Gneq*dtime)


      ! Compute the equivalent shear stress
      !
      tauBar_tau = tauBar_tr - Gneq*dtime*nuv_tau


      ! Compute the viscous stretching
      !
      Dv_tau = dsqrt(one/two)*nuv_tau*Nv_tau


      ! Compute the viscous deformation gradient at the
      !  end of the increment using the exponential map
      !
      if(nuv_tau.le.zero) then
         Fv_tau = Fv_t
      else
         call spectral(dtime*Dv_tau,Dv_eig,Dv_vec,stat)
         expdtDv = zero
         expdtDv(1,1) = dexp(Dv_eig(1))
         expdtDv(2,2) = dexp(Dv_eig(2))
         expdtDv(3,3) = dexp(Dv_eig(3))
         expdtDv = matmul(matmul(Dv_vec,expdtDv),transpose(Dv_vec))
         Fv_tau = matmul(expdtDv,Fv_t)
      endif


      ! Check to make sure that det(Fv_tau)>0
      !
      call mdet(Fv_tau,tmp)
      if(tmp.le.zero) then
         write(*,*) 'det(Fv_tau).le.zero in INTEG'
        stat = 0
        return
      endif


      ! Compute the Mandel stres at the end of the increment
      !
      Me_tau = Me_tr - two*Gneq*dtime*Dv_tau


      ! Compute the contribution to the Cauchy stress due to
      !  viscoelasticity
      !
      T_tau = T_tau + 
     +        matmul(Re_tr,matmul(Me_tau,transpose(Re_tr)))/detF


      ! Compute contribution to the spatial tangent due to
      !  viscoelasticity
      !
      Nv_tau = matmul(Re_tr,matmul(Nv_tau,transpose(Re_tr)))
      if (tauBar_tr.gt.zero) then
	   Gshear_tilde = (tauBar_tau/tauBar_tr)*Gneq
      else
	   Gshear_tilde = Gneq
      end if
      fac = (one + Gneq*dtime/eta)**(-one)
      Lambda_tilde =  - Gshear_tilde*two/three
      if (tauBar_tr.gt.zero) then
         c3 = -two*Gneq*((tauBar_tau/tauBar_tr) - fac)
      else
         c3 = -two*Gneq*(one - fac)
      end if
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) 
     +          +Gshear_tilde*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k))
     +          +Lambda_tilde*Iden(i,j)*Iden(k,l)
     +          +c3*Nv_tau(i,j)*Nv_tau(k,l)
               enddo
            enddo
         enddo
      enddo


      return
      end subroutine NeoHookean 
      
!****************************************************************************

      subroutine Gent(props,nprops,dtime,F_tau,Fv_t,
     +        T_tau,Fv_tau,SpTanMod,stat)
      !
      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fv_t(3,3),T_tau(3,3),
     +  Fv_tau(3,3),SpTanMod(3,3,3,3),Iden(3,3),G0,Kbulk,eta,Im,
     +  Gneq,detF,Finv(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),Cinv(3,3),
     +  detC,trC,I1bar,fac,GShearGent,dGdF(3,3),TR_tau(3,3),
     +  dTRdF(3,3,3,3),Fv_t_inv(3,3),det_Fv_t,Fe_tr(3,3),Re_tr(3,3),
     +  Ue_tr(3,3),Ee_tr(3,3),trEe_tr,Ee0_tr(3,3),Me_tr(3,3),
     +  tauBar_tr,Nv_tau(3,3),nuv_tau,tauBar_tau,Dv_tau(3,3),Dv_eig(3),
     +  Dv_vec(3,3),expdtDv(3,3),tmp,Me_tau(3,3),Gshear_tilde,
     +  Lambda_tilde,fac2,c3
      !
      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
      !
      G0     = props(1) ! Ground-state shear modulus
      Kbulk  = props(2) ! Bulk modulus
      Im     = props(3) ! Limiting chain extensibility parameter
      Gneq   = props(4) ! Nonequilibrium (Maxwell element) stiffness
      eta    = props(5) ! Maxwell element viscosity
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Gent spring
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT = transpose(F_tau)
      FinvT = transpose(Finv)


      ! Compute the right Cauchy-Green tensor, its inverse, 
      !  its determinant, and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      call matInv3D(C_tau,Cinv,detC,stat)
      trC = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 

      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*trC
 

      ! Compute the ``I1-3'' factor appearing in the Gent model
      !
      fac = (I1bar - three)/Im
      if(fac.gt.0.95d0) fac = 0.95d0
      fac = one/(one - fac)
 

      ! Compute the ``shear'' modulus. Note: this is not really the shear
      !  modulus, but it will help when computing the material tangent later
      !
      GShearGent = G0*fac
 

      ! Compute the derivative of the ``shear modulus'' with respect
      !  to the deformation gradient for use in the material tangent
      !
      dGdF = two*(G0/Im)*(detF**(-two/three))*
     +     fac*fac*(F_tau - third*trC*FinvT)
 

      ! Compute the 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*GShearGent*(F_tau-third*trC*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
 

      ! Compute the Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))


      ! Compute the material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3                  
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*dGdF(k,l)*
     +                 (
     +                 F_tau(i,j) - third*trC*Finv(j,i)
     +                 )
     +                 + (detF**(-two/three))*GshearGent*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*trC*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*trC*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the spatial tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Maxwell element
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the trial elastic deformation gradient
      !
      call matInv3D(Fv_t,Fv_t_inv,det_Fv_t,stat)
      Fe_tr = matmul(F_tau,Fv_t_inv)


      ! Compute the trial kinematics
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr,stat)
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3)
      Ee0_tr = Ee_tr - (one/three)*trEe_tr*Iden


      ! Compute the trial Mandel stress, which is deviatoric
      !
      Me_tr = two*Gneq*Ee0_tr


      ! Compute the trial equiv. tensile stress
      !
      tauBar_tr = dsqrt(one/two)*dsqrt(sum(Me_tr*Me_tr))


      ! Compute the direction of viscous flow
      !
      if(tauBar_tr.le.zero) then
         Nv_tau = zero
      else
         Nv_tau = dsqrt(one/two)*(Me_tr/tauBar_tr)
      endif


      ! Compute the equivalent shear viscous strain rate
      !
      nuv_tau = tauBar_tr/(eta + Gneq*dtime)


      ! Compute the equivalent shear stress
      !
      tauBar_tau = tauBar_tr - Gneq*dtime*nuv_tau


      ! Compute the viscous stretching
      !
      Dv_tau = dsqrt(one/two)*nuv_tau*Nv_tau


      ! Compute the viscous deformation gradient at the
      !  end of the increment using the exponential map
      !
      if(nuv_tau.le.zero) then
         Fv_tau = Fv_t
      else
         call spectral(dtime*Dv_tau,Dv_eig,Dv_vec,stat)
         expdtDv = zero
         expdtDv(1,1) = dexp(Dv_eig(1))
         expdtDv(2,2) = dexp(Dv_eig(2))
         expdtDv(3,3) = dexp(Dv_eig(3))
         expdtDv = matmul(matmul(Dv_vec,expdtDv),transpose(Dv_vec))
         Fv_tau = matmul(expdtDv,Fv_t)
      endif


      ! Check to make sure that det(Fv_tau)>0
      !
      call mdet(Fv_tau,tmp)
      if(tmp.le.zero) then
         write(*,*) 'det(Fv_tau).le.zero in INTEG'
        stat = 0
        return
      endif


      ! Compute the Mandel stres at the end of the increment
      !
      Me_tau = Me_tr - two*Gneq*dtime*Dv_tau


      ! Compute the contribution to the Cauchy stress due to
      !  viscoelasticity
      !
      T_tau = T_tau + 
     +        matmul(Re_tr,matmul(Me_tau,transpose(Re_tr)))/detF


      ! Compute contribution to the spatial tangent due to
      !  viscoelasticity
      !
      Nv_tau = matmul(Re_tr,matmul(Nv_tau,transpose(Re_tr)))
      if (tauBar_tr.gt.zero) then
	   Gshear_tilde = (tauBar_tau/tauBar_tr)*Gneq
      else
	   Gshear_tilde = Gneq
      end if
      fac2 = (one + Gneq*dtime/eta)**(-one)
      Lambda_tilde =  - Gshear_tilde*two/three
      if (tauBar_tr.gt.zero) then
         c3 = -two*Gneq*((tauBar_tau/tauBar_tr) - fac2)
      else
         c3 = -two*Gneq*(one - fac2)
      end if
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) 
     +          +Gshear_tilde*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k))
     +          +Lambda_tilde*Iden(i,j)*Iden(k,l)
     +          +c3*Nv_tau(i,j)*Nv_tau(k,l)
               enddo
            enddo
         enddo
      enddo
      

      return
      end subroutine Gent
      
!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi_1,xi_2 coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2),w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi_1,xi_2 coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2),w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt
      
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i):      shape function of node i at the intpt.
      ! dshxi(i,j): derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

!************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi_1-xi_2 domain
      !  to x-y domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi_1-xi_2 domain
      !  to x-y domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da
      
!****************************************************************************
!     The next subroutine calculates various kinematical quantities
!      associated with the deformation gradient
!****************************************************************************

      subroutine skinem(F,R,U,E,istat)
      !
      ! This subroutine performs the right polar decomposition
      !  F = RU of the deformation gradient F into a rotation
      !  R and the right stretch tensor U.  The logarithmic 
      !  strain E = ln(U) is also returned.
      !
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	U(3,3):       the right stretch tensor; output
      !	Uinv(3,3):    the inverse of U
      !	C(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Ueigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of U
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
      implicit none
      !
      integer istat
      !
      real*8 F(3,3),C(3,3),omega(3),Ueigval(3),eigvec(3,3),
     +  U(3,3),E(3,3),Uinv(3,3),R(3,3),detF
     

      !	Store the identity matrix in R, U, and Uinv
      !
      call onem(R)
      call onem(U)
      call onem(Uinv)
      

      ! Store the zero matrix in E
      !
      call zerom(E)
      

      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the 
      !  time increment.
      !
      call mdet(F,detF)
      if (detF.le.0.d0) then
        write(*,'(/5X,A/)') '--problem in kinematics-- the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
      

      ! Calculate the right Cauchy-Green tensor C
      !
      C = matmul(transpose(F),F)
      
 
      ! Calculate the eigenvalues and eigenvectors of C
      !
      call spectral(C,omega,eigvec,istat)
      

      ! Calculate the principal values of U and E
      !
      Ueigval(1) = dsqrt(omega(1))
      Ueigval(2) = dsqrt(omega(2))
      Ueigval(3) = dsqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
      !
      E(1,1) = dlog(Ueigval(1))
      E(2,2) = dlog(Ueigval(2))
      E(3,3) = dlog(Ueigval(3))
      

      ! Calculate the complete tensors U and E
      !
      U = matmul(matmul(eigvec,U),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
      

      ! Calculate Uinv
      !
      call matInv3D(U,Uinv,detF,istat)
      

      ! calculate R
      !
      R = matmul(F,Uinv)
      

      return
      end subroutine skinem

!****************************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!****************************************************************************

      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)


      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      call eigsrt(D,V,3,np)
	

      return
      end subroutine spectral
	
!****************************************************************************

      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
     
      
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
        B(ip) = A(ip,ip)
        D(ip) = B(ip)
        Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
        !
        ! Sum off-diagonal elements
        !
        sm = 0.d0
        do ip=1,n-1
          do iq=ip+1,n
            sm = sm + dabs(A(ip,iq))
          end do
        end do
        !
        ! If sm = 0., then return.  This is the normal return,
        !  which relies on quadratic convergence to machine
        !  underflow.
        !
        if (sm.eq.0.d0) return
        !
        ! In the first three sweeps carry out the PQ rotation only if
        !  |A_PQ| > tresh, where tresh is some threshold value,
        !  see equation (11.1.25).  Thereafter tresh = 0.
        !
        if (i.lt.4) then
          tresh = 0.2d0*sm/n**2
        else
          tresh = 0.d0
        end if
        !
        do ip=1,n-1
          do iq=ip+1,n
            G = 100.d0*dabs(A(ip,iq))
            !
            ! After four sweeps, skip the rotation if the 
            !  off-diagonal element is small.
            !
            if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +          .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
              A(ip,iq) = 0.d0
            else if (dabs(A(ip,iq)).gt.tresh) then
              H = D(iq) - D(ip)
              if (dabs(H)+G.eq.dabs(H)) then
                !
                ! T = 1./(2.*theta), equation (11.1.10)
                !
                T =A(ip,iq)/H
              else
                theta = 0.5d0*H/A(ip,iq)
                T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
                if (theta.lt.0.d0) T = -T
              end if
              C = 1.d0/dsqrt(1.d0 + T**2.d0)
              S = T*C
              tau = S/(1.d0 + C)
              H = T*A(ip,iq)
              Z(ip) = Z(ip) - H
              Z(iq) = Z(iq) + H
              D(ip) = D(ip) - H
              D(iq) = D(iq) + H
              A(ip,iq) = 0.d0
              !
              ! Case of rotations 1 <= J < P
              !		
              do j=1,ip-1
                G = A(j,ip)
                H = A(j,iq)
                A(j,ip) = G - S*(H + G*tau)
                A(j,iq) = H + S*(G - H*tau)
              end do
              !
              ! Case of rotations P < J < Q
              !
              do j=ip+1,iq-1
                G = A(ip,j)
                H = A(j,iq)
                A(ip,j) = G - S*(H + G*tau)
                A(j,iq) = H + S*(G - H*tau)
              end do
              !
              ! Case of rotations Q < J <= N
              !
              do j=iq+1,n
                G = A(ip,j)
                H = A(iq,j)
                A(ip,j) = G - S*(H + G*tau)
                A(iq,j) = H + S*(G - H*tau)
              end do
              do j = 1,n
                G = V(j,ip)
                H = V(j,iq)
                V(j,ip) = G - S*(H + G*tau)
                V(j,iq) = H + S*(G - H*tau)
              end do
              nrot = nrot + 1
            end if
          end do
        end do
        !
        ! Update D with the sum of T*A_PQ, and reinitialize Z
        !
        do ip=1,n
          B(ip) = B(ip) + Z(ip)
          D(ip) = B(ip)
          Z(ip) = 0.d0
        end do
      end do


      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
      

      return
      end subroutine jacobi
	
!****************************************************************************

      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  jacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      

      do i=1,n-1
        k = i
        P = D(i)
        do j=i+1,n
          if (D(j).ge.P) then
            k = j
            P = D(j)
          end if
        end do
        if (k.ne.i) then
          D(k) = D(i)
          D(i) = P
          do j=1,n
            P = V(j,i)
            V(j,i) = V(j,k)
            V(j,k) = P
          end do
        end if
      end do
      

      return
      end subroutine eigsrt

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet

!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
        do J=1,3
          if (i .eq. j) then
            A(i,j) = 1.d0
          else
            A(i,j) = 0.d0
          end if
        end do
      end do


      return
      end subroutine onem

!****************************************************************************

      subroutine zerom(A)
      !
      ! This subroutine sets all entries of a 3 by 3 matrix to zero
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)
      
      
      do i=1,3
        do j=1,3
          A(i,j) = 0.d0
        end do
      end do
      
      
      return
      end subroutine zerom

!****************************************************************************