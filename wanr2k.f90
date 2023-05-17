      Program Wannier_band_structure
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=200,nblocks=10,nE = 1000
      real,parameter::Emin = -2.0d0,Emax = 2.0d0,eta = 1e-2
      complex, parameter:: icomp=(0,1)
!------------------------------------------------------
      integer*4 ik,ikmax,iz,row_offset,col_offset,iblock,mb,ib,jb,iE,iside
      real*8 kz,delE,E
      character(len=30)::klabel(nkpath)
      character(len=80) hamil_file,nnkp,line
      integer*4,parameter::nk=(nkpath-1)*np+1
      integer*4 i,j,k,nr,i1,i2,nb,lwork,info
      real*8,parameter::third=1d0/3d0!,kz=0d0
      real*8 phase,pi2,jk,a,b,ef
      real*8 klist(3,1:nk),xk(nk),kpath(3,np),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),Aspec(nE,nk,2)
      real*8,allocatable:: rvec(:,:),ene(:,:),rwork(:),w(:,:),layer(:,:,:)
      integer*4,allocatable:: ndeg(:)
      complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:)
      complex*16 temp1,temp2
!------------------------------------------------------
      write(hamil_file,'(a,a)')trim(adjustl(prefix)),"_hr.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

      pi2=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
      read(98,*)bvec
!---------------kpath
      data kpath(:,1) /     0.2d0,      0.0d0,    0.0d0/  !M
      data kpath(:,2) /     0.0d0,      0.0d0,    0.0d0/  !G
      data kpath(:,3) /     0.1d0,      0.1d0,    0.0d0/  !K
      ! data kpath(:,4) /     0.0d0,      0.0d0,    0.0d0/  !G
      ! data kpath(:,5) /     0.0d0,      0.0d0,    0.5d0/  !A
      ! data kpath(:,6) /     0.5d0,      0.0d0,    0.5d0/  !L
      ! data kpath(:,7) /     third,      third,    0.5d0/  !H
      ! data kpath(:,8) /     0.0d0,      0.0d0,    0.5d0/  !A

      data klabel     /'M','Gamma','K'/

      ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

      xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
      xkl(1)=xk(1)
      

      k=0
      ktemp1=0d0
      do i=1,nkpath-1
       do j=1,np
        k=k+1
        jk=dfloat(j-1)/dfloat(np)
        klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
        ktemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
        if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
        if(j.eq.1) xkl(i)=xk(k)
        ktemp1=ktemp2
       enddo
      enddo
      klist(:,nk)=kpath(:,nkpath)
      ktemp2=klist(1,nk)*bvec(:,1)+klist(2,nk)*bvec(:,2)+klist(3,nk)*bvec(:,3)
      xk(nk)=xk(nk-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
      xkl(nkpath)=xk(nk)
!      write(*,*)klist
      klist=klist*pi2

!------read H(R)
      open(99,file=trim(adjustl(hamil_file)),err=444)
      open(100,file='band.dat')
      read(99,*)
      read(99,*)nb,nr
      mb = nb*nblocks
      allocate(rvec(3,nr),Hk(mb,mb),Hamr(nb,nb,nr),ndeg(nr),ene(mb,nk))
      read(99,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)rvec(1,k),rvec(2,k),rvec(3,k),i1,i2,a,b
               hamr(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
     lwork=max(1,2*mb-1)
     allocate(work(max(1,lwork)),rwork(max(1,3*mb-2)))
     allocate(w(mb,nk),layer(mb,nk,2))

!---- Fourrier transform H(R) to H(k)
      ene=0d0
      do k=1,nk
         HK=(0d0,0d0)
         do j=1,nr
            iz = rvec(3,j)
            phase=0.0d0
            do i=1,2
               phase=phase+klist(i,k)*rvec(i,j)
            enddo

            row_offset=((iz + abs(iz))/2)*nb
            col_offset=(abs(iz - abs(iz))/2)*nb

            do ib=1,nb
               do jb=1,nb
                  do iblock=0,nblocks-abs(iz)-1
                        i1 = ib + row_offset + iblock*nb
                        i2 = jb + col_offset + iblock*nb
                        Hk(i1,i2)=Hk(i1,i2)+Hamr(ib,jb,j)* &
                        dcmplx(cos(phase),-sin(phase))/float(ndeg(j))

                  enddo
               enddo
            enddo

         enddo

         call zheev('V','U',mb,Hk,mb,ene(:,k),work,lwork,rwork,info)
         !ef = (minval(ene(13*nblocks,:)) + maxval(ene(12*nblocks,:)))/2
         ef = 6.0d0
         ! extract required spectral function
         delE = (Emax - Emin)/(nE - 1)
         do i=1,mb
            w(:,k) = dconjg(Hk(:,i)) * Hk(:,i) 
            layer(i, k, 1) = sum(w(1:3*nb, k))  !changed stuff here!
            layer(i, k, 2) = sum(w(mb-3*(nb-1):mb, k))
         enddo
         do i=1,mb
            do iE=1,nE
                  E = Emin + iE * delE
                  do iside=1,2
                        Aspec(iE,k,iside) = Aspec(iE,k,iside) - aimag((layer(i,k,iside))/((ene(i,k)-ef) - E + icomp*eta))
                  enddo
            enddo
         enddo

         
      enddo
      do iE=1,nE
         E = Emin + iE * delE
         do k=1,nk
           write(100,'(4(x,f12.6))') xk(k),E,Aspec(iE,k,1),Aspec(iE,k,2)
         enddo
           !write(100,*)
           !write(100,*)
      enddo
      call write_plt(nkpath,xkl,klabel,ef)
      stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file)),' not found'
      stop

      end

     subroutine write_plt(nkp,xkl,kl,ef)
     implicit none
     integer nkp,i
     real*8 xkl(nkp), ef
     character(len=30)kl(nkp), klabel
     
     open(99,file='spectral.plt')
     write(99,'(a,f12.8)')'ef=',ef
     write(99,'(a)') 'set xtics ( \'
     do i=1,nkp
        klabel=adjustl(kl(i))
        if(klabel(1:1).eq.'g'.or.klabel(1:1).eq.'G')kl(i)="{/Symbol \107}"
        if(i.ne.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i),", \"
        if(i.eq.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i)," )"
     enddo
     write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [',xkl(1),':',xkl(nkp),']'
     write(99,'(a)') &
          'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 7.00in, 7.50in'
     write(99,'(a,f4.2,a)')'set output "spectral.pdf"'
     write(99,'(10(a,/),a)') &
          'set encoding iso_8859_1',&
          'set size ratio 0 1.0,1.0',&
          'unset colorbox',&
          'set ylabel "E-E_{f} (eV)"',&
          'set yrange [ -1 : 1.0 ]',&
          'unset key',&
          'set ytics 0.5 scale 1 nomirror out',&
          'set mytics 2',&
          'set parametric',&
          'set trange [-10:10]',&
          'plot "band.dat" u 1:2:4 with image,\'
!     do i=2,nkp-1
!       write(99,'(f12.6,a)') xkl(i),',t with l lt 2  lc -1,\'
!     enddo
!     write(99,'(a)') 't,0 with l lt 2  lc -1'
    end subroutine write_plt
