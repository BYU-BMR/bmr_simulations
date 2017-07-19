!********************************************************************
!**  Program Keffective                                            **
!**  Hybrid calculation of effective conductivity in x,y,z         **
!**  directions for a microstructure contained on a 3D rectangular ** 
!**  grid or lattice. The structure is previously segmented such   **
!**  that each voxel contains one of three distinct materials or   **
!**  domains. The domain identities are read into this code from   **
!**  a text file. Separate calculations are then done for ionic    **
!**  and electronic conductivity. Intrinsic ionic and electronic   **
!**  conductivity values for each domain are hardcoded here, but   **
!**  may be changed.                                               **
!**                                                                **
!**  The "hybrid" designation is because initially the grid can be **
!**  coarse-grained to a smaller size before finite volume (FV)    **
!**  calculations are performed. A renormalization scheme is used  **
!**  to do the coarse-graining (see subroutine k_coarse). Once     **
!**  that is done, the grid is of a more manageable size and local **
!**  variations in conductivity have been partially smoothed.      **
!**                                                                **
!**  For the FV calculations, the user can choose the boundary     **
!**  conditions: either mixed (part insulating and part fixed-     **
!**  potential) or periodic BCs are used. User-adjustable          **
!**  parameters are found below and in subroutine k_FVM.           **
!**                                                                **
!**  Program memory requirements scale with NX*NY*NZ, so don't     **
!**  make the grid too large. Use of coarse-graining will          **
!**  significantly speed up the FV calculation, and will generate  **
!**  more consistent results, particularly when the conductivity   **
!**  of the original grid varies rapidly or differences between    **
!**  domain conductivities are large. Sensitivity of the final     **
!**  results to the amount of coarse-graining should be checked.   **
!**                                                                **
!**  Fortran90 code by Dean Wheeler, Brigham Young University,     **
!**  2013-2016. Known to work with g95 compiler in DOS/Windows.    **
!**  e.g. compile with command: g95 -O3 Keffective8.f90            **
!**                                                                **
!********************************************************************

 module nodes !storage for common data

!********************************************************************
!       ADJUST THESE VALUES
!       *******************

!       First make sure NX,NY,NZ are correct for the size and
!       and number of input text files for the structure
   integer, parameter:: NX=74, NY=74, NZ=44


!       Variable CG is the minimal degree of coarse-graining. The 
!       original grid is reduced to by factor 2**CG in each direc- 
!       tion, rounded up to the next integer. Therefore, if CG=0 
!       no coarse-graining is used. If CG=1 then grid is cut in  
!       half in each direction, etc.
   integer, parameter :: CG = 1


!       Lengths of original nodes/voxels in the x,y,z directions
!       (arbitrary units, though microns is a convenient choice)
   real, parameter,dimension(3):: L = (/ 0.4, 0.4, 0.4 /)


!       Base name of geometry input files 
!       (a number increment and .txt get appended to this)
   character*14:: filebase = 'toda523'


!      Variable BC says whether mixed or periodic boundary
!      conditions are used on outer surfaces of the cell.
   integer,parameter:: BC = 1

!      'Mixed' (BC=1) means insulating on four outer surfaces  
!      parallel to current flow and fixed potential on two outer
!      surfaces orthogonal to current flow. The difference between 
!      these two constant-potential surfaces drives the current flow.
!
!      'Periodic' (BC=2) means periodic BCs are used in all three 
!      directions. An external homogeneous field is used to drive
!      current flow in one direction.
!
!      Generally, BC=2 converges to a solution faster, but you 
!      should compare both to see what effect BC has. If the
!      difference is significant, this suggests that your
!      sampling volume is not sufficiently large or 
!      representative to generate a homogeneous transport 
!      property.


!       Intrinsic conductivities of the domains:
   real,parameter,dimension(0:2,2) :: kdomain = reshape((/ &

            0.0,    0.188,   1.0,    &! IONIC

            3.0,  3000.,    0.01    &! ELECTRONIC

                          /),(/3,2/)) !(don't change this line)

!       The first 3 numbers are IONIC conductivities of Active, 
!       Carbon, and Pore domains (corresponds to carrier = 1); 
!       the last 3 numbers are ELECTRONIC conductivities of Active,
!       Carbon, and Pore domains (corresponds to carrier = 2).
!       In each case the conductivities are in arbitrary units or
!       are relative, because an effective conductivity will be 
!       produced with the same units. Leave the (/3,2/) part alone.
!
!********************************************************************
!         common data that specifies structure
   integer, parameter:: ngrid=NX*NY*NZ
   integer, parameter::  maxpass =ceiling(log(dfloat(max(NX,NY,NZ)))&
                                         /log(2.d0))
   integer (kind=1) ::  domain(NX,NY,NZ)
   integer          ::  nnode(3) = (/NX,NY,NZ/)
   character(len=10) :: ctype(2) 
   data ctype/'Ionic','Electronic'/

 end module nodes

!********************************************************************

 program Keffective
   use nodes
   implicit none

!          variables used for file input
   character*21:: filename
   character*3 :: increm

!          other variables
   integer ::  pass,skip2,cnode(3)
   integer ::  steps(3,0:maxpass)
   integer ::  i,j,k,i1,j1,k1,dir,carrier,dd
   integer ::  IFN
   integer ::  d0,d1,d2,d3,d4,d5,d6,d7,ndomain(0:2)
   integer ::  NNL(0:2,0:2),NND(0:2,0:2),NNS(0:2,0:2)
   real    ::  di,dj,dk
   real    ::  PNNL(0:2,0:2),PNND(0:2,0:2),PNNS(0:2,0:2),frac(0:2)
   real    ::  dvoxel         ! average length of voxel
   real    ::  eps = 1.e-16   ! parameter to ensure num stability
   real    ::  k1avg(0:maxpass),k2avg(0:maxpass)
   real    ::  keff(3,0:maxpass),kerr(3,0:maxpass)
   real    ::  kmax(2) = 0., kmin(2) = 0., kmed(2) = 1.
   real    ::  condmax,condmed
   real,allocatable :: kstor(:,:,:)
   character(len=8) :: bctype(2) 
   data bctype/'mixed','periodic'/

!******************************************************************

!       write out diagnostic data on input files

!       output to screen
   write(6,*)
   write(6,*) 'Reading in domain identities from image file:'
   write(6,*) '    Need at least ',NX*NY*NZ,' lines'
   write(6,*) '    Each should have (x,y,z,voxel) entries'


!       output to text file
   if (BC == 1) then  
      open(9,file='Cond_results_BC1.txt')
   else
      open(9,file='Cond_results_BC2.txt')
   endif                           
   rewind (9)
   write(9,*) 'Calculation of effective conductivities'
   write(9,*) '***************************************'
   write(9,*)
   write(9,*) 'Actual domain size (arbitrary units) is ', &
               L(1)*NX,'x',L(2)*NY,'x',L(3)*NZ


!*******************************************************************
!             Input domain identity for each node

!          open input file
      filename = trim(filebase)//'.txt'
      IFN = 100
      open(IFN,file=filename,status='old',err=100)

!           read row by row through file 
   do k = 1, NZ
     do i = 1, NX
       do j = 1, NY
         read(IFN,*,err=105,end=110) di,dj,dk,dd
         if (dd == 255) then
            domain(i,j,k) = 2
         elseif (dd == 140) then
            domain(i,j,k) = 0
         else
            domain(i,j,k) = 1
         endif  
      enddo 
     enddo
   enddo
      close(IFN) 
      goto 200

!        branches that get used if file has problem
  100 write(6,*) 'problem with ',filename
      write(6,*) 'could not open'
      stop
  105 write(6,*) 'problem with ',filename
      write(6,*) 'data wrong format'
      stop
  110 write(6,*) 'problem with ',filename
      write(6,*) 'not enough data'
      stop
  200 continue 
 

!*************************************************************
!          Validate domain input and accumulate number of
!          voxels of each domain type.

   do d0 = 0, 2
      ndomain(d0) = 0
   enddo
   do k = 1, NZ
     do j = 1, NY
        do i = 1, NX
          d0 = domain(i,j,k)
          if ((d0 < 0).or.(d0 > 2)) then
            write(6,*) 'an input domain is ',d0
            write(6,*) 'but must be 0,1,2'
            stop
          endif
          ndomain(d0) = ndomain(d0) + 1
       enddo
     enddo
   enddo

!     Output domain fractions and conductivities

   write(9,*)
   write(9,*) 'Domain fractions and conductivities (arbitrary units)'
   write(9,*) '----------------------------------------------------'

   do d0 = 0, 2
      frac(d0) = float(ndomain(d0))/float(NX*NY*NZ)
      write(9,fmt='(A,I1,A,F6.3,A,I1,A,F8.3,A,F8.3)') &
        '   frac',d0,' =',frac(d0),'      k',d0,    &
        ' =',kdomain(d0,1),', ',kdomain(d0,2)

!     Calculate theoretical max and min effective conductivities 

      kmax(1) = kmax(1) + frac(d0)*kdomain(d0,1)
      kmax(2) = kmax(2) + frac(d0)*kdomain(d0,2)
      kmin(1) = kmin(1) + frac(d0)/(kdomain(d0,1)+eps)
      kmin(2) = kmin(2) + frac(d0)/(kdomain(d0,2)+eps)
      kmed(1) = kmed(1) * kdomain(d0,1)**frac(d0)
      kmed(2) = kmed(2) * kdomain(d0,2)**frac(d0)
   enddo
   kmin(1) = 1./(kmin(1)+eps)
   kmin(2) = 1./(kmin(2)+eps)

!******************************************************************
!    Make conductivity maps for different levels of coarse-grain:
!    Compute maximum domain conductivity and approximate
!    median conductivity for both types of conduction

   do carrier = 1, 2
      condmax = 0.
      do d0 = 0, 2
         if (kdomain(d0,carrier) > condmax) &
            condmax = kdomain(d0,carrier)
      enddo
      condmed = max(sqrt(kmax(carrier)*kmin(carrier)), &
                             kmax(carrier)/16)
 
      write(6,*)
      write(6,*) 'making conductivity maps for carrier=',carrier
      call k_maps(carrier,condmax,condmed)	
   enddo

!******************************************************************
!           Calculate nearest-neighbor probabilities
!
!       periodic boundaries are currently used.
!       to eliminate this, make loops over NX-1 in place of NX,
!       NY-1 in place of NY, etc. Also change normalization of
!       PNNL, etc.

   NNL = 0
   NND = 0
   NNS = 0
   do k = 1, NZ
     k1 = mod(k,NZ)+1
     do j = 1, NY
        j1 = mod(j,NY)+1
        do i = 1, NX
          i1 = mod(i,NX)+1
!           Identify voxel identities for a 2x2x2 subgrid
          d0 = domain(i, j, k )
          d1 = domain(i1,j, k )
          d2 = domain(i, j1,k )
          d3 = domain(i1,j1,k )
          d4 = domain(i, j, k1)
          d5 = domain(i1,j, k1)
          d6 = domain(i, j1,k1)
          d7 = domain(i1,j1,k1)
!           Find unique pairs of lateral neighbors
!           and increment counters
          NNL(d0,d1) = NNL(d0,d1) + 1
          NNL(d0,d2) = NNL(d0,d2) + 1
          NNL(d0,d4) = NNL(d0,d4) + 1
!           Find unique pairs of diagonal neighbors
!           and increment counters
          NND(d0,d3) = NND(d0,d3) + 1
          NND(d0,d5) = NND(d0,d5) + 1
          NND(d0,d6) = NND(d0,d6) + 1
          NND(d1,d2) = NND(d1,d2) + 1
          NND(d1,d4) = NND(d1,d4) + 1
          NND(d2,d4) = NND(d2,d4) + 1
!           Find unique pairs of super-diagonal neighbors
!           and increment counters
          NNS(d0,d7) = NNS(d0,d7) + 1
          NNS(d1,d6) = NNS(d1,d6) + 1
          NNS(d2,d5) = NNS(d2,d5) + 1
          NNS(d3,d4) = NNS(d3,d4) + 1
       enddo
     enddo
   enddo
!    Make matrices symmetric
   NNL = NNL + transpose(NNL)
   NND = NND + transpose(NND)
   NNS = NNS + transpose(NNS) 
!    Normalize by total number of pairs
   PNNL = dfloat(NNL)/dfloat( 6*ngrid)
   PNND = dfloat(NND)/dfloat(12*ngrid)
   PNNS = dfloat(NNS)/dfloat( 8*ngrid)
   dvoxel = (L(1)+L(2)+L(3))/3.
!    Output results to file
   write(9,*)
   write(9,*)'Pairwise probabilities (assumes cubic voxels)'
   write(9,*)'----------------------------------------------------'
   write(9,*)'             d=0         d=L   d=sqrt2*L   d=sqrt3*L'
   write(9,*)'Pair        Self     Lateral    Diagonal  Super-diag'
   write(9,fmt='(A5,4F12.7)')'00',frac(0),PNNL(0,0), &
                                PNND(0,0),PNNS(0,0)
   write(9,fmt='(A5,4F12.7)')'01',0.,2.*PNNL(0,1), &
                        2.*PNND(0,1),2.*PNNS(0,1)
   write(9,fmt='(A5,4F12.7)')'02',0.,2.*PNNL(0,2), &
                        2.*PNND(0,2),2.*PNNS(0,2)
   write(9,fmt='(A5,4F12.7)')'11',frac(1),PNNL(1,1), &
                                PNND(1,1),PNNS(1,1)
   write(9,fmt='(A5,4F12.7)')'12',0.,2.*PNNL(1,2), &
                        2.*PNND(1,2),2.*PNNS(1,2)
   write(9,fmt='(A5,4F12.7)')'22',frac(2),PNNL(2,2), &
                                PNND(2,2),PNNS(2,2)
   write(9,*)
   write(9,*)'Radial distribution functions (assumes cubic voxels)'
   write(9,*)'----------------------------------------------------'
   write(9,*)'             Pair'
   write(9,*)'Dist (um)       00        01        02        11',&
                                          '        12        22'
   write(9,fmt='(F9.5,6F10.4)') 0.,1./frac(0), &
                           0., 0., 1./frac(1), 0., 1./frac(2)
   write(9,fmt='(F9.5,6F10.4)') dvoxel, &
                         PNNL(0,0)/frac(0)/frac(0), &
                         PNNL(0,1)/frac(0)/frac(1), &
                         PNNL(0,2)/frac(0)/frac(2), &
                         PNNL(1,1)/frac(1)/frac(1), &
                         PNNL(1,2)/frac(1)/frac(2), &
                         PNNL(2,2)/frac(2)/frac(2)
   write(9,fmt='(F9.5,6F10.4)') dvoxel*sqrt(2.), &
                         PNND(0,0)/frac(0)/frac(0), &
                         PNND(0,1)/frac(0)/frac(1), &
                         PNND(0,2)/frac(0)/frac(2), &
                         PNND(1,1)/frac(1)/frac(1), &
                         PNND(1,2)/frac(1)/frac(2), &
                         PNND(2,2)/frac(2)/frac(2)
   write(9,fmt='(F9.5,6F10.4)') dvoxel*sqrt(3.), &
                         PNNS(0,0)/frac(0)/frac(0), &
                         PNNS(0,1)/frac(0)/frac(1), &
                         PNNS(0,2)/frac(0)/frac(2), &
                         PNNS(1,1)/frac(1)/frac(1), &
                         PNNS(1,2)/frac(1)/frac(2), &
                         PNNS(2,2)/frac(2)/frac(2)
   write(9,*)

!**************************************************************
!      output results of conductivity calculations

   write(9,*)
   write(9,*)
   write(9,*) 'BC = ', bctype(BC), ' maxpass =',maxpass
   write(9,*) '***********************************************'
   write(6,*) '***********************************************'

!         Calculate mean conductivity for 3 directions
!         and standard deviation of mean. Use conductivites
!         and boundary conditions previously specified.

!         Then output results.

   do carrier = 1, 2 ! loop over ionic and electronic types

!          ***the heart of the calculation***
    call k_FVM(carrier, steps, keff, kerr)	

    do pass = CG, maxpass
!           detemine isotropic average and uncertainty
      k1avg(pass) = (keff(1,pass) + keff(2,pass) + keff(3,pass))/3.
      k2avg(pass) = (keff(1,pass)**2 + keff(2,pass)**2 &
                                     + keff(3,pass)**2)/3.
      k2avg(pass) = sqrt(abs(k2avg(pass) &
                  - k1avg(pass)*k1avg(pass))/3.)

!          compute number of coarse-grain voxels in each direction
      skip2 = 2**pass
      cnode = (nnode + skip2 - 1)/skip2


!          output mean and standard dev of mean (standard error)
      write(6,*)
      write(6,*) 'Summary for ',ctype(carrier), ' PASS =',pass
      write(6,*) '    BC = ', bctype(BC)
      write(6,*) '    k_avg =', k1avg(pass),' +-', k2avg(pass)
                      
      write(6,*)
      write(6,*) '***************************************'

      write(9,*)
      write(9,*) ctype(carrier), '  PASS =',pass
      write(9,fmt='(2(A,F8.4),A,I7)') '   k_x   = ',keff(1,pass), &
                      ' +-',kerr(1,pass),'     steps =',steps(1,pass)
      write(9,fmt='(2(A,F8.4),A,I7)') '   k_y   = ',keff(2,pass), &
                      ' +-',kerr(2,pass),'     steps =',steps(2,pass)
      write(9,fmt='(2(A,F8.4),A,I7)') '   k_z   = ',keff(3,pass), &
                      ' +-',kerr(3,pass),'     steps =',steps(3,pass)
      write(9,fmt='(2(A,F8.4))') '   k_avg = ',k1avg(pass), &
                                         ' +-',k2avg(pass)

    enddo !end loop over passes


    write(9,*)
    write(9,fmt='(2(A,F8.4),A)') '   k_min = ',kmin(carrier), &
                                 '   (theoretical)'
    write(9,fmt='(2(A,F8.4),A)') '   k_med = ',kmed(carrier), &
                                 '   (theoretical)'
    write(9,fmt='(2(A,F8.4),A)') '   k_max = ',kmax(carrier), &
                                 '   (theoretical)'
    write(9,*)
   enddo
   close(9)
   stop
 end program


!********************************************************************
!**                                                                **
!**    Subroutine to coarse-grain a grid of spatially varying      **
!**    conductivity into a single effective conductivity for       **
!**    flow of current in a particular direction. Procedure is     **
!**    based on the renormalization scheme of M.R. Karim and       **
!**    K. Krabbenhoft, Transp. Porous Med. 85, 677?90 (2010).     **
!**                                                                **
!**    The original KK method is altered/generalized to compute    **
!**    an overall conductivity for an arbitrary-size rectangular   **
!**    3D grid (i.e. any number of nodes in each direction is      **
!**    allowed). Non-cubic voxels are likewise allowed, but it is  **
!**    assumed that each voxel of the input grid has the same      **
!**    size/shape. The method divides the grid into small blocks   **
!**    of nodes/voxels (each size 2x2x2 or smaller). The           **
!**    effective conductance of each block is then computed.       **
!**    These blocks are recursively combined into larger blocks    **
!**    until the original grid is reduced to one large voxel.      **
!**                                                                **
!**    Note that conductance (which incorporates geometry of       **
!**    the grid) rather than conductivity is used for inter-       **
!**    mediate calculations, and is returned by this subroutine.   ** 
!**                                                                **
!********************************************************************
 subroutine k_coarse(pass,L,NX,NY,NZ,kgrid)
   implicit none
   integer, intent(in) :: pass      !number of coarse-grain passes
   integer, intent(in) :: NX,NY,NZ  !original dimensions of grid
   real,    intent(in) :: L(3)      !dimensions of single voxel
   real    :: kgrid(NX,NY,NZ,3)     !input/output conductivities

!        local variables
   real    :: eps = 1.e-16          !small param for num stability
   real    :: k8(2,2,2)             !initial block conductances
   real    :: kL,kU1,kU2,kU         !intermed block conductances
   integer :: blockx,blocky,blockz  !dimensions of block
   integer :: ix,iy,iz,jx,jy,jz     !counters and indices
   integer :: skip,d,it,kx,ky,kz    !counters and indices

!******************************************************************
!        initialize number of nodes in each direction and number
!        of coarse-grain passes 

!        Perform multiple coarse-grain passes
   do it = 1, pass
      skip = 2**(it-1)

!        loop over x,y,z directions, skipping over
!        nodes/voxels not needed in this pass
!        and computing size and location of each block
      do iz = 1, NZ, 2*skip
        jz = iz + skip
        kz = min(iz + 2*skip - 1, NZ)
        blockz = 2
        if (jz > NZ) then
           jz = iz
           blockz = 1
        endif
        do iy = 1, NY, 2*skip
          jy = iy + skip
          ky = min(iy + 2*skip - 1, NY)
          blocky = 2
          if (jy > NY) then
             jy = iy
             blocky = 1
          endif
          do ix = 1, NX, 2*skip
            jx = ix + skip
            kx = min(ix + 2*skip - 1, NX)
            blockx = 2
            if (jx > NX) then
               jx = ix
               blockx = 1
            endif
!              Now get conductance of block of voxels
!              in the designated direction, accounting
!              for fact that block may not be of full 2x2x2 size.
!              Results are stored in original array of node
!              conductances and therefore replace original values.
!              select appropriate formulas based on direction,
!              computing lower and upper conductance bounds
!              by series and parallel resistance formulas.
!              Then store single block conductance that is
!              geometric average of lower and upper bounds.
!         x direction
            k8(:,:,:) = 0.
            k8(1:blockx,1:blocky,1:blockz) = &
                       kgrid(ix:jx:skip,iy:jy:skip,iz:jz:skip,1)
            if (blockx == 1) k8(2,1:2,1:2) = 1./eps
            kL = k8(1,1,1)*k8(2,1,1)/(k8(1,1,1)+k8(2,1,1)+eps) &
               + k8(1,2,1)*k8(2,2,1)/(k8(1,2,1)+k8(2,2,1)+eps) &
               + k8(1,1,2)*k8(2,1,2)/(k8(1,1,2)+k8(2,1,2)+eps) &
               + k8(1,2,2)*k8(2,2,2)/(k8(1,2,2)+k8(2,2,2)+eps)
            kU1 = k8(1,1,1) + k8(1,2,1) + k8(1,1,2) + k8(1,2,2)
            kU2 = k8(2,1,1) + k8(2,2,1) + k8(2,1,2) + k8(2,2,2)
            kU  = kU1*kU2/(kU1+kU2+eps)
            kgrid(ix:kx,iy:ky,iz:kz,1) = sqrt(kL*kU)
!         y direction
            k8(:,:,:) = 0.
            k8(1:blockx,1:blocky,1:blockz) = &
                       kgrid(ix:jx:skip,iy:jy:skip,iz:jz:skip,2)
            if (blocky == 1) k8(1:2,2,1:2) = 1./eps
            kL = k8(1,1,1)*k8(1,2,1)/(k8(1,1,1)+k8(1,2,1)+eps) &
               + k8(2,1,1)*k8(2,2,1)/(k8(2,1,1)+k8(2,2,1)+eps) &
               + k8(1,1,2)*k8(1,2,2)/(k8(1,1,2)+k8(1,2,2)+eps) &
               + k8(2,1,2)*k8(2,2,2)/(k8(2,1,2)+k8(2,2,2)+eps)
            kU1 = k8(1,1,1) + k8(2,1,1) + k8(1,1,2) + k8(2,1,2)
            kU2 = k8(1,2,1) + k8(2,2,1) + k8(1,2,2) + k8(2,2,2)
            kU  = kU1*kU2/(kU1+kU2+eps)
            kgrid(ix:kx,iy:ky,iz:kz,2) = sqrt(kL*kU)
!         z direction
            k8(:,:,:) = 0.
            k8(1:blockx,1:blocky,1:blockz) = &
                       kgrid(ix:jx:skip,iy:jy:skip,iz:jz:skip,3)
            if (blockz == 1) k8(1:2,1:2,2) = 1./eps
            kL = k8(1,1,1)*k8(1,1,2)/(k8(1,1,1)+k8(1,1,2)+eps) &
               + k8(2,1,1)*k8(2,1,2)/(k8(2,1,1)+k8(2,1,2)+eps) &
               + k8(1,2,1)*k8(1,2,2)/(k8(1,2,1)+k8(1,2,2)+eps) &
               + k8(2,2,1)*k8(2,2,2)/(k8(2,2,1)+k8(2,2,2)+eps)
            kU1 = k8(1,1,1) + k8(2,1,1) + k8(1,2,1) + k8(2,2,1)
            kU2 = k8(1,1,2) + k8(2,1,2) + k8(1,2,2) + k8(2,2,2)
            kU  = kU1*kU2/(kU1+kU2+eps)
            kgrid(ix:kx,iy:ky,iz:kz,3) = sqrt(kL*kU)
         enddo
        enddo
      enddo
   enddo
!     At this point the original grid has been transformed to a 
!     reduced number of independent nodes. The original conduc-
!     tivities have been replaced with conductances, meaning 
!     geometry is incorporated and the last nodes in the x,y,z 
!     directions may be smaller than the remaining nodes. Now
!     make final corrections based on shape of original voxels 
!     and add small parameter for numerical stability.
   do d = 1, 3
      kgrid(:,:,:,d) = kgrid(:,:,:,d)*L(1)*L(2)*L(3)/L(d)**2 + eps
   enddo
 end subroutine k_coarse


!********************************************************************
!**                                                                **
!**  Subroutine to loop over all possible levels of coarse-        **
!**  graining to make representative conductivity color maps       **
!**                                                                **
!********************************************************************
 subroutine k_maps(carrier,kdmax,kdmed)

   use nodes
   implicit none

!        subroutine input/output data
   integer, intent(in)  :: carrier  !type of conduction
   real,    intent(in)  :: kdmax    !max domain cond
   real,    intent(in)  :: kdmed    !median domain cond
   
!        local variables
   real     :: kexpon             !exponent for scaling cond
   real     :: kgrid(NX,NY,NZ,3)  !grid conductivities
   real     :: color(NX,NY)       !grid cond converted to color
   real     :: delx,dely,delz     !dimensions of coarse-grained node
   integer  :: i,j,k,iZ           !pixel indices
   integer  :: pass,skip2         !coarse-grain indices

! arbitrary z value to make all maps with; must be (1...NZ)
   iZ = 1

! Calculate exponent (0...1) for nonlinear conductivity-to-color 
! encoding. This is equivalent to "gamma correction" used to  
! tune display of luminance in computer images.
   kexpon = log(0.5) / log(kdmed/kdmax)
   write(6,fmt='(2(A,F8.4),A)') '   expon = ',kexpon
   write(9,fmt='(2(A,F8.4),A)') '   expon = ',kexpon


!  loop over all coarse-graining level
   do pass = 0, maxpass
!        initialize isotropic kgrid then coarse-grain it, 
!        making it anisotropic
     do k = 1, NZ
       do j = 1, NY
         do i = 1, NX
            kgrid(i,j,k,1) = kdomain(domain(i,j,k),carrier)
            kgrid(i,j,k,2) = kdomain(domain(i,j,k),carrier)
            kgrid(i,j,k,3) = kdomain(domain(i,j,k),carrier)
         enddo
       enddo
     enddo
     call k_coarse(pass,L,NX,NY,NZ,kgrid)

!    rescale conductances to turn back into conductivities
     skip2 = 2**pass
     delz = L(3)*min(skip2,NZ+1-iZ+mod(iZ-1,skip2))
     do j = 1, NY
        dely = L(2)*min(skip2,NY+1-j+mod(j-1,skip2))
        do i = 1, NX
!          convert node conductance into conductivity;
!          average node conductivity in 3 directions;
!          scale by max domain conductivity; then map
!          to color scale by exponent
           delx = L(1)*min(skip2,NX+1-i+mod(i-1,skip2))
           color(i,j) = ( (kgrid(i,j,iZ,1)*delx**2 &
                          +kgrid(i,j,iZ,2)*dely**2 &
                          +kgrid(i,j,iZ,3)*delz**2)&
                         /(3*kdmax*delx*dely*delz) )**kexpon
        enddo
     enddo
!    create BMP graphic file of conductivity map for one z value
     call makebmp(color,NX,NY,"kmap",carrier,pass)
   enddo
   end subroutine k_maps


!********************************************************************
!**                                                                **
!**   Subroutine to make Windows 3.x type BMP image with custom    **
!**   8-bit color palette and run-length encoding (RLE) for        **
!**   compression. The input image array light(NX,NY) contains     **
!**   real intensity values from 0 to 1. These are scaled and      **
!**   discretized to integer 0...255, which corresponds to the     **
!**   RGB color map value. RLE is then done. Header and pixel      **
!**   data are finally written to the .bmp file.                   **
!**   By changing 'dilate' each pixel can be repeated in horiz     **
!**   and vert directions in order to enlarge native size of       **
!**   image.                                                       **
!**             written by Dean Wheeler, 2013-2015                 **
!**                                                                **
!********************************************************************
      subroutine makebmp(light,NX,NY,fname,id1,id2)
      implicit none

!* interface arguments
      real, intent(in) :: light(NX,NY) !scaled cond array (0...1)
      integer, intent(in) :: NX,NY     !dimensions of grid
      character(len=4),intent(in):: fname !image base filename
      integer, intent(in) :: id1,id2   !image 1- and 2-digit indices

!* local variables
      integer    :: i,j,h,k,ic,run,maxp    !counters and indices
      integer    :: xp,yp,dilate,pad       !numbers of pixels
      integer    :: color,colornxt         !color storage  
      character(len=12)      :: fnameout   !BMP file name
      character(len=1078)    :: header     !BMP file header
      character(len=9*NX*NY) :: RLEpixels  !compressed pixel data
      character(len=4)       :: byt4       !4-byte integer
  
!* stretch ratio for enlarging image (user can change)
      dilate = 3

!* pixels in each direction
     xp = NX*dilate              !image width
     pad = 4*int((xp+3)/4) - xp  !padding so line is multpl of 4
     yp = NY*dilate              !image height

!* convert input image array to 1-byte color value from color table;
!* and do run-length encoding on each line
   k = 1
   do j = 1, NY ! loop over lines
     do h = 1, dilate !loop over repeats of each line
       run = 0
       color = int(255.999*light(1,j)) !color of leftmost pixel
       do i = 1, NX  !loop over pixels in one line
         colornxt = int(255.999*light(i,j)) !store next pixel
!          check for end of run: either pixel color
!          change or too many pixels in one run
         if ((colornxt /= color) .or. &
                   (run + dilate > 255) ) then
            RLEpixels(k:k+1) = char(run)//char(color)
            k = k + 2  !increment byte counter
            run = 0
            color = colornxt
         endif
!          increment run-length counter
         run = run + dilate
       enddo
  !      finish the line by ending run and using end-of-line marker
       RLEpixels(k:k+1)   = char(run)//char(color)
       RLEpixels(k+2:k+3) = char(0)//char(0)
       k = k + 4  !increment byte counter
     enddo
   enddo
!    finish last line
   RLEpixels(k-1:k-1) = char(1)  !end-of-bitmap marker 
   maxp = k-1              !total number of bytes in pixel data
   if (maxp > 9*NX*NY) write(6,*) &
                     'must increase length of RLEpixels variable'

!* BMP header (bytes 1-54)
      header( 1: 2) = 'BM'            !declare this is BMP file
      header( 3: 6) = byt4(1078+maxp) !tot file size = header+data
      header( 7:10) = byt4(0)         !4 reserved bytes     
      header(11:14) = byt4(1078)      !header size, incl color table
      header(15:18) = byt4(40)        !bit-map headr size must = 40
      header(19:22) = byt4(xp)        !image width in pixels
      header(23:26) = byt4(yp)        !image height in pixels
      header(27:30) = byt4(8*256**2+1)!set 8 color bits per pixel
      header(31:34) = byt4(1)         !use RLE-8 compression
      header(35:38) = byt4(maxp)      !size of image pixel data
      header(39:42) = byt4(3937)      !horiz resolution (pixels/m)
      header(43:46) = byt4(3937)      !vert resolution (pixels/m)
      header(47:50) = byt4(256)       !color table has 256 colors
      header(51:54) = byt4(0)         !num of important colors
!* BMP color table for 256 colors, comprising bytes 55-1078
!* For each color give blue,green,red,alpha byte values (0-255)
      do ic = 0, 127
!          linearly interpolate between black and red for 
!          lower intensity values
         header(55+ic*4:58+ic*4) = &
                   char(0)//char(0)//char(ic*2)//char(0)
!          linearly interpolate between red and white for 
!          higher intensity values
         header(567+ic*4:570+ic*4) = &
                   char(ic*2)//char(ic*2)//char(255)//char(0)
      enddo
!* generate full filename
      write(fnameout,'(A,i1.1,''-''i2.2,''.bmp'')') fname,id1,id2
!* write header and pixel data to bmp file
      open(unit=12,file=fnameout,status='unknown')
      write(12,'(A)',advance='no') header//RLEpixels(1:maxp)
      close(12)
      return
      end subroutine makebmp

!************************************************************
!**  function to convert integer value 
!**  to 4 bytes in little-endian format

       character*4 function byt4(i)
       integer i
       byt4 = char(mod(i,256))//char(mod(i/256,256)) &
            //char(mod(i/256**2,256))//char(i/256**3)
       return
       end



!********************************************************************
!**                                                                **
!**  Subroutine to calculate the effective conductivity on a grid  **
!**  for a particular charge carrier (1,2) in a given cartesian    **
!**  direction (1,2,3) and for given BCs (1,2).                    **
!**                                                                **
!**  The face-centered finite volume method is used. An iterative  **
!**  scheme is used to solve the resulting matrix equation. The    **
!**  user can adjust how rapidly/aggressively to attempt the       **
!**  convergence.                                                  **
!**                                                                **
!**  The computed potential at each node is taken as relative to   **
!**  the average linear potential profile induced by an external   **
!**  field. For example, if the field is in the x direction:       **
!**       potl_tot = potl + field*x                                **
!**                                                                **
!********************************************************************
 subroutine k_FVM(carrier,steps,keff,kerr)

   use nodes
   implicit none

!        subroutine input/output data
   integer, intent(in)  :: carrier          !ionic=1, electronic=2
   integer, intent(out) :: steps(3,0:maxpass) !num of matrx iter
   real,    intent(out) :: keff(3,0:maxpass)  !effective conductivity 
   real,    intent(out) :: kerr(3,0:maxpass)  !std dev of mean
   
!        local variables
   real     :: kgrid(NX,NY,NZ,3)  !grid conductivities
   real     :: pgrid(NX,NY,NZ,3)  !grid potential
   real     :: color(NX,NY)       !grid cond converted to color
   real     :: delx,dely,delz     !dimensions of coarse-grained node
   real*8   :: keff1,keff2,keff3,shape,kp(6),ksite(3),knabor,ktotal
   real     :: field(6),potlsat,potlold,threshold,TS,TSi,TSf,rate,pn
   integer  :: dir,site,site1,j1,stepsmax,nsoln,unbal
   integer  :: post1,post2,post3
   integer  :: i,j,k
   logical  :: stopflag
   integer  :: in,jn,kn,d,lo,hi,ngridc,aerr,pass,skip2
   real     :: xc,yc,zc
   integer  :: x0,x1,y0,y1,z0,z1

   integer,allocatable :: inode(:),nabor(:,:)
   real,allocatable    :: potl(:,:,:),ktot(:,:,:),kpair(:,:,:,:)
   real,allocatable    :: boundary(:,:,:)

   real    :: eps = 1.e-16          !small param for num stability
   integer :: cnode(3)              !max node after coarse-graining
   integer :: maxcnode              !largest value of cnode

!*******************************************************************
!     CHANGE THESE VALUES TO ADJUST CONVERGENCE OF SOLUTION
!     *****************************************************

   stepsmax = 10       ! max number of matrix iterations per element
   threshold = 1.E-4   ! error threshold to stop interations
   TSi = 0.1           ! initial (minimum) timestep parameter
   TSf = 1.5           ! final (maximum) timestep parameter
   rate = 0.8          ! how fast to adjust timestep value

!      TSi and TSf (timestep) values must be greater than 0, and 
!      control how much change is made to node potentials during 
!      a given iteration. Lower values provide more convergence  
!      stability and larger values provide convergence 
!      acceleration. Adjust TSi, TSf, and 'rate' to balance 
!      computational time vs. stability. Usually you want to 
!      start off with a small timestep, then speed up to a larger
!      value. Also change 'threshold' value to balance 
!      computational time vs. accuracy of result. If the 
!      conductivity result is too sensitive to the exact 
!      convergence parameters chosen, consider using a more
!      coarse-grained grid.
!
!*******************************************************************
!      initialize potential of elements on full-size grid
   pgrid = 0. 
   do i = 1, nnode(1)
      pgrid(i,:,:,1) = i - 0.5
   enddo
   do j = 1, nnode(2)
      pgrid(:,j,:,2) = j - 0.5
   enddo
   do k = 1, nnode(3)
      pgrid(:,:,k,3) = k - 0.5
   enddo

!      major loop over level of coarse-graining

   CGloop: do pass = maxpass, CG, -1

   skip2 = 2**pass
   cnode = (nnode + skip2 - 1)/skip2
   maxcnode = maxval(cnode)
   ngridc = cnode(1)*cnode(2)*cnode(3)

!       allocate variable-size arrays for storing node
!       information on coarse-grain grid

   allocate (potl(0:cnode(1),0:cnode(2),0:cnode(3)),STAT=aerr)
   allocate (ktot(0:cnode(1),0:cnode(2),0:cnode(3)),STAT=aerr)
   allocate (kpair(6,0:cnode(1),0:cnode(2),0:cnode(3)),STAT=aerr)
   allocate (boundary(0:cnode(1),0:cnode(2),0:cnode(3)),STAT=aerr)
   allocate (inode(0:maxcnode),STAT=aerr)
   allocate (nabor(6,maxcnode),STAT=aerr)
   if (aerr /= 0) then
      write(6,*) 'not enough memory'
      stop
   endif

!*********************************************************************
!        initialize isotropic kgrid then coarse-grain it, 
!        making it anisotropic

   do k = 1, NZ
     do j = 1, NY
       do i = 1, NX
          kgrid(i,j,k,1) = kdomain(domain(i,j,k),carrier)
          kgrid(i,j,k,2) = kdomain(domain(i,j,k),carrier)
          kgrid(i,j,k,3) = kdomain(domain(i,j,k),carrier)
       enddo
     enddo
   enddo

   call k_coarse(pass,L,NX,NY,NZ,kgrid)

!    map position from coarse-grain grid to original grid 

   inode(0) = 0
   do i = 1, maxcnode
      inode(i) = 1 + (i-1)*skip2
   enddo

!    Boundary Conditions are integrated with a neighbor
!    list in each of six directions.
!    If BC=1 set up insulating boundaries 
!    on lateral sides of cuboid and constant-potential
!    planes on end of cuboid that are normal to field.
!    If BC=2 set up periodic boundary conditions on
!    outside of cuboid.

   do d = 1, 3 !loop over x,y,z axes
      do i = 1, cnode(d) !independent nodes on this axis

         if (BC == 1) then 
!             implement insulating BCs in pos and neg directions
            nabor(2*d-1,i) = min(i+1,cnode(d))  
            nabor(2*d,i)   = max(i-1,1)
         else 
!             implement periodic BCs in pos and neg directions
            nabor(2*d-1,i) = mod(i,cnode(d))+1
            nabor(2*d,i)   = mod(i-2+cnode(d),cnode(d))+1
         endif
      enddo
   enddo

!********************************************************************* 
!  major loop over 3 coordinate directions to impose electric field
   directionloop: do dir = 1, 3

!       write progress report to screen
   write(6,*)
   write(6,*) '   *calculating ',trim(ctype(carrier)), &
                             '  dir=',dir,' pass=',pass


!    If only one node in field direction, there is a 
!    better way to compute effective conductivity by
!    summing up all parallel conductances and correcting
!    for geometry where 
!           keff = shape*(sum of conductances)
!          shape = L/A = L^2/Vol
!    After this calculation, terminate subroutine

   if (cnode(dir) == 1) then
     steps(dir,pass) = 0
     keff(dir,pass) = 0.
     kerr(dir,pass) = 0.
     shape = (nnode(dir)*L(dir))**2/(ngrid*L(1)*L(2)*L(3))
     do k = 1, cnode(3)
       kn = inode(k)
       do j = 1, cnode(2)
         jn = inode(j)
         do i = 1, cnode(1)
            in = inode(i)
            keff(dir,pass) = keff(dir,pass) &
                           + shape*kgrid(in,jn,kn,dir)
         enddo ! end loop over nodes
       enddo
     enddo

!      output results to screen
     write(6,fmt='(2(A,F8.4))') '     k = ',keff(dir,pass), &
                                      ' +-',kerr(dir,pass)

!      end this cycle over coordinate directions
     cycle
   endif

! else more than one node in field direction ...
!
!  ********************************************
!     interpolate INITIAL potential of nodes

   potl = 0. 
   do k = 1, cnode(3)
      zc = 0.5*(inode(k) + min(k*skip2, nnode(3)))
      z0 = max(1,floor(zc))
      z1 = min(nnode(3),ceiling(zc))
      do j = 1, cnode(2)
         yc = 0.5*(inode(j) + min(j*skip2, nnode(2)))
         y0 = max(1,floor(yc))
         y1 = min(nnode(2),ceiling(yc))
         do i = 1, cnode(1)
            xc = 0.5*(inode(i) + min(i*skip2, nnode(1)))
            x0 = max(1,floor(xc))
            x1 = min(nnode(1),ceiling(xc))
            potl(i,j,k) = (pgrid(x0,y0,z0,dir)+pgrid(x1,y0,z0,dir)&
                    + pgrid(x0,y1,z0,dir) + pgrid(x1,y1,z0,dir) &
                    + pgrid(x0,y0,z1,dir) + pgrid(x1,y0,z1,dir) &
                    + pgrid(x0,y1,z1,dir) + pgrid(x1,y1,z1,dir))/8.
         enddo
      enddo
   enddo

!  *********************************************
!     pre-calculate neighbor effective pair
!     conductances as harmonic mean of respective 
!     node conductances then normalize by total 
!     for all neighbors
   boundary = 0.
   do k = 1, cnode(3)
      kn = inode(k)
      do j = 1, cnode(2)
         jn = inode(j)
         do i = 1, cnode(1)
            in = inode(i)
            ksite(1) = kgrid(in,jn,kn,1) !central node conduct
            ksite(2) = kgrid(in,jn,kn,2)
            ksite(3) = kgrid(in,jn,kn,3)
            knabor = kgrid(inode(nabor(1,i)),jn,kn,1)
            kp(1)  = 2.d0/(1.d0/ksite(1) + 1.d0/knabor)
            knabor = kgrid(inode(nabor(2,i)),jn,kn,1)
            kp(2)  = 2.d0/(1.d0/ksite(1) + 1.d0/knabor)
            knabor = kgrid(in,inode(nabor(3,j)),kn,2)
            kp(3)  = 2.d0/(1.d0/ksite(2) + 1.d0/knabor)
            knabor = kgrid(in,inode(nabor(4,j)),kn,2)
            kp(4)  = 2.d0/(1.d0/ksite(2) + 1.d0/knabor)
            knabor = kgrid(in,jn,inode(nabor(5,k)),3)
            kp(5)  = 2.d0/(1.d0/ksite(3) + 1.d0/knabor)
            knabor = kgrid(in,jn,inode(nabor(6,k)),3)
            kp(6)  = 2.d0/(1.d0/ksite(3) + 1.d0/knabor)

!             If BC=1 and if node is on leading or trailing 
!             surface in direction of field, then one neighbor
!             is constant-potential plane. Therefore need to
!             to reassign this neighbor to node zero and make
!             corrections to pair properties. Likewise if BC=2
!             then need to introduce a potential offset.
!             In any case, the imposed voltage difference across
!             the grid is equal to the number of nodes in the
!             field direction.
!             These changes must be at this point in the code.
!                 post1 = position in direction of field
!                 post2,post3 = position in other directions
                   
            if (dir == 1) then
               post1 = i
               post2 = j
               post3 = k
            elseif (dir == 2) then
               post1 = j
               post2 = k
               post3 = i
            else
               post1 = k
               post2 = i
               post3 = j
            endif
            if (BC == 1) then !mixed boundaries
!              Check for active boundary and if found set
!              to conducting boundary of fixed potential
!              equal to zero or +nnode.
               if (post1 == 1) then
                  nabor(2*dir,1) = 0   !reassign neighbor
                  kp(2*dir) = 2.d0*ksite(dir)
               elseif (post1 == cnode(dir)) then
                  nabor(2*dir-1,post1) = 0   !reassign neighbor
                  kp(2*dir-1) = 2.d0*ksite(dir)
                  boundary(i,j,k) = 2.d0*ksite(dir)*nnode(dir)
               endif
          
!                  Check for insulating boundaries.
!                  Above implementation of insulating BCs does
!                  work, but fix below reduces number of steps
!                  required to converge to solution because it
!                  enforces BC each step.
               if (post2 == 1) then
                  kp(2*mod(dir,3)+2) = 0.d0
               elseif (post2 == cnode(mod(dir,3)+1)) then
                  kp(2*mod(dir,3)+1) = 0.d0
               endif
               if (post3 == 1) then !insulating boundary
                  kp(2*mod(dir+1,3)+2) = 0.d0
               elseif (post3 == cnode(mod(dir+1,3)+1)) then
                  kp(2*mod(dir+1,3)+1) = 0.d0
               endif
            else !periodic boundaries
!              Check for active boundary and if found set 
!              potential jump across boundary equal to +-nnode
               if (post1 == 1) then !on active boundary
                  boundary(i,j,k) = -kp(2*dir)*nnode(dir)
               elseif (post1 == cnode(dir)) then !active boundary
                  boundary(i,j,k) = kp(2*dir-1)*nnode(dir)
               endif
            endif  

!              Normalize and store kp and boundary values;
!              store ktot
            ktotal = kp(1) +kp(2) +kp(3) +kp(4) +kp(5) +kp(6)
            do j1 = 1, 6
               kpair(j1,i,j,k) = kp(j1)/ktotal
            enddo
            boundary(i,j,k) = boundary(i,j,k)/ktotal
            ktot(i,j,k) = ktotal
         enddo
      enddo
   enddo

!     **********************************************
!     perform iterations to solve for relative node
!     potentials using finite difference equations

   stopflag = .false.
   steps(dir,pass) = 0
   TS = TSi
   do while (.not. stopflag)
     steps(dir,pass) = steps(dir,pass) + 1
     unbal = 0
!         update timestep parameter
     TS = TSf +(TS-TSf)*(1.-rate*float(999+ngridc)**(-0.5))

!       loop over nodes

     do k = 1, cnode(3)
       do j = 1, cnode(2)
         do i = 1, cnode(1)

!            Potential of site that would exactly satisfy the
!            node equation; it is a weighted average of neighbor
!            potentials plus a term to account for BCs
           potlsat = kpair(1,i,j,k)*potl(nabor(1,i),j,k) &
                   + kpair(2,i,j,k)*potl(nabor(2,i),j,k) &
                   + kpair(3,i,j,k)*potl(i,nabor(3,j),k) &
                   + kpair(4,i,j,k)*potl(i,nabor(4,j),k) &
                   + kpair(5,i,j,k)*potl(i,j,nabor(5,k)) &
                   + kpair(6,i,j,k)*potl(i,j,nabor(6,k)) &
                   + boundary(i,j,k)

!            update site potential with use of timestep parameter
           potlold = potl(i,j,k)
           potl(i,j,k) = TS*potlsat + (1.-TS)*potlold

!            count number of nodes with potential errors
!            above threshold
           if (abs(potlsat - potlold) >= threshold) &
                                       unbal = unbal + 1
         enddo ! end loop over nodes
       enddo ! end loop over nodes
     enddo ! end loop over nodes

!       continue iterations if still have node errors
!       and maximum iterations not yet reached
     stopflag = .true.
     if ((unbal > 0) .and. (steps(dir,pass) < stepsmax*ngridc)) &
                                      stopflag =.false.

!       keep track of convergence progress
     if ((mod(steps(dir,pass),300) == 0) .or. stopflag) then
        write(6,'(A,I3,A,F5.2)')'     balanced=', &
             nint(100.-100.*unbal/ngridc),'%   TS=',TS
     endif 
   enddo !end loop over iterations on potential  
   write(6,*)'    iterations stopped at ',steps(dir,pass)

!     ***********************************************
!     calculate effective conductivity for three 
!     sampling planes normal to the field direction 
!     where keff = shape*I
!          shape = L/(A*deltaV) = L^2/(Vol*deltaV)

   keff1 = 0.d0
   keff2 = 0.d0
   keff3 = 0.d0
   j1 = 2*dir-1  !neighbor across/opposite the plane
   shape = nnode(dir)*L(dir)**2 /(ngrid*L(1)*L(2)*L(3))
   do k = 1, cnode(3)
    do j = 1, cnode(2)
     do i = 1, cnode(1)
      if (dir == 1) then
         post1 = i
         pn = potl(nabor(j1,i),j,k)
       elseif (dir == 2) then
         post1 = j
         pn = potl(i,nabor(j1,j),k)
      else
         post1 = k
         pn = potl(i,j,nabor(j1,k))
      endif
!        first plane for calculating current
      if (post1 == max(1,nint(cnode(dir)/6.d0))) then
        keff1 = keff1 + ktot(i,j,k)*shape*( kpair(j1,i,j,k) &
              *(pn - potl(i,j,k)) + max(0.,boundary(i,j,k)) )
      endif
!        second plane for calculating current
      if (post1 == max(1,nint(cnode(dir)/2.d0))) then
        keff2 = keff2 + ktot(i,j,k)*shape*( kpair(j1,i,j,k) &
              *(pn - potl(i,j,k)) + max(0.,boundary(i,j,k)) )
      endif
!        third plane for calculating current
      if (post1 == nint(cnode(dir)*5.d0/6.d0)) then
        keff3 = keff3 + ktot(i,j,k)*shape*( kpair(j1,i,j,k) &
              *(pn - potl(i,j,k)) + max(0.,boundary(i,j,k)) )
      endif
     enddo ! end loop over nodes
    enddo
   enddo

!     calculate mean and std dev of mean for 3 sampling planes
!   
   keff(dir,pass) = (keff1 + keff2 + keff3)/3.d0
   kerr(dir,pass) = sqrt(abs((keff1**2+keff2**2+keff3**2)/3.d0 &
                    -keff(dir,pass)**2)/3.d0)

!       warn if keff not fully converged/consistent

   if (kerr(dir,pass) > (0.01*keff(dir,pass) + 10.*threshold)) then
      write(6,*)'    conductivity not well-converged:'
      write(6,fmt='(3(A,F9.4))')'     k1=',keff1, &
                          ' k2=',keff2,' k3=',keff3
   endif

!  output to screen results
   write(6,fmt='(2(A,F8.4))') '     k = ',keff(dir,pass), &
                                    ' +-',kerr(dir,pass)

!  *******************************************************
!    Transfer potentials of nodes back onto original grid

!      (this needs to be improved with some kind of interpolation
!       algorithm to make more smooth representation of potentials
!       on original (finer) grid to speed up finite volume calc)

   do k = 1, cnode(3)
      z0 = inode(k)
      z1 = min(k*skip2,nnode(3))
      do j = 1, cnode(2)
         y0 = inode(j)
         y1 = min(j*skip2,nnode(2))
         do i = 1, cnode(1)
            x0 = inode(i)
            x1 = min(i*skip2,nnode(1))
            pgrid(x0:x1,y0:y1,z0:z1,dir) = potl(i,j,k)
         enddo
      enddo
   enddo
!  ***************************************************

!  end major loop over coordinate directions
   enddo directionloop

!  deallocate temporary arrays for node properties
   deallocate (potl)
   deallocate (ktot)
   deallocate (kpair)
   deallocate (boundary)
   deallocate (inode)
   deallocate (nabor)

!  end major loop over degree of coarse-graining
   enddo CGloop

   return
 end subroutine k_FVM
