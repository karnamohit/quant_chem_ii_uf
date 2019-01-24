program scf_KR
use glovar
implicit none

real*8,parameter :: econv=1.0E-7,densconv=1.0E-5,maxstep=50
real*8,parameter :: nsing=1,ntrip=1

print *,'convergence parameters for: energy=',econv,'density residue=',densconv

call aces_init_rte
call aces_ja_init

!*****Reading JOBARC (after running ACES II)*****
        call getrec(10, "JOBARC", "NATOMS", 1, natom)
        call getrec(20, "JOBARC", "NBASTOT", 1, nbasorb)
        call getrec(30, "JOBARC", "NMPROTON", 1, nmolocc)
        call getrec(40, "JOBARC", "NUCREP", 1, vnn)

!        nmolocc=2

        print *,'number of basis functions used =',nbasorb
        print *,'number of atoms in the molecule =',natom
        print *,'number of electrons =',nmolocc
        print *,'nuclear repulsion, Vnn =',vnn,'Ha'

        print *,'primary test concluded.'

        n=nbasorb*(nbasorb+1)/2
        nmolhlfocc=nmolocc/2

        allocate(core1d(n),ovrlp1d(n))
        allocate(core2d(nbasorb,nbasorb),ovrlp2d(nbasorb,nbasorb))
        allocate(veff4d(nbasorb,nbasorb,nbasorb,nbasorb))
        allocate(dummy1(600),dummy2(600),dummy3(n,n))
        allocate(transform(nbasorb,nbasorb))
        allocate(fock(nbasorb,nbasorb),fockp(nbasorb,nbasorb))
        allocate(constp(nbasorb,nbasorb),const(nbasorb,nbasorb))
        allocate(dens2d(nbasorb,nbasorb),ovrlpinv2d(nbasorb,nbasorb))
        allocate(densold(nbasorb,nbasorb),densnew(nbasorb,nbasorb))
        allocate(densresm(nbasorb,nbasorb))
        allocate(onee(nbasorb,nbasorb),twoe(nbasorb,nbasorb))
        allocate(eigen(nbasorb,nbasorb))

        core1d=0
        ovrlp1d=0
        core2d=0
        ovrlp2d=0
        veff4d=0
        dummy1=0
        dummy2=0
        dummy3=0
        transform=0
        fock=0
        fockp=0
        const=0
        constp=0
        dens2d=0
        ovrlpinv2d=0
        densold=0
        densnew=0
        densresm=0

!*****Getting overlap, core hamiltonian matrices*****
        call Get1EInt(core1d,ovrlp1d,dummy1,dummy2,n)
        call EXPND2(core1d,core2d,nbasorb)
        call EXPND2(ovrlp1d,ovrlp2d,nbasorb)
        call Get2EInt(veff4d,dummy3,dummy1,dummy2,nbasorb,n)

        print *,'subroutine call check 1'

        deallocate(core1d,ovrlp1d)
        deallocate(dummy1,dummy2,dummy3)

        transform=0.0
        eigen=ovrlp2d

!*****Diagonalizing the overlap matrix, S, and getting the transformation matrix, U*****
        call eig(eigen,transform,24,nbasorb,1)

        print *,'subroutine call check 2'

        allocate(dummy4(nbasorb,nbasorb))

        dummy4=0

!*****Forming S^(-1/2) matrix*****
        do i=1,nbasorb

         dummy4(i,i)=eigen(i,i)**(-0.5)

        end do

        ovrlpinv2d=matmul(matmul(transform,dummy4),transpose(transform))

        deallocate(dummy4)

        print *,'S^(-1/2) matrix formed'

!*****Building an initial Fock matrix*****
        do mm=1,nbasorb
         do nn=1,nbasorb

          fock(mm,nn)=core2d(mm,nn)

         end do
        end do

        print *,'intial Fock matrix formed'

!*****Forming F'=S^(-1/2)*F, from F'C'=C'E *****
        fockp=matmul(ovrlpinv2d,matmul(fock,ovrlpinv2d))

        print *,'transformed Fock matrix formed'

!*****Diagonalizing F' and getting eigenvalues and the coeff. matrix*****
        call eig(fockp,constp,23,nbasorb,0)

        print *,'subroutine call check 3'

!*****Forming C=S^(-1/2)*C' *****
        const=matmul(ovrlpinv2d,constp)

        allocate(dummy4(nbasorb,nmolhlfocc),dummy3(nmolhlfocc,nmolhlfocc))

        do i=1,nbasorb
         do j=1,nmolhlfocc

          dummy4(i,j)=const(i,j)

         end do
        end do

        dummy3=matmul(transpose(dummy4),matmul(ovrlp2d,dummy4))

        deallocate(dummy3,dummy4)

!*****Building the charge density matrix*****
        do i=1,nbasorb
         do j=1,nbasorb
          do k=1,nmolhlfocc

           dens2d(i,j)=dens2d(i,j)+2*const(i,k)*const(j,k)

          end do
         end do
        end do

!*****Calculating the initial HF energy*****
        ehf=vnn

        do mm=1,nbasorb
         do nn=1,nbasorb

          ehf=ehf+0.5*dens2d(nn,mm)*core2d(mm,nn)

         end do
        end do

        do i=1,nmolhlfocc

         ehf=ehf+fockp(i,i)

        end do

        print *,'initial energy is',ehf

!*****HF energy and density matrix update for the SCF procedure*****
        ehfold=ehf
        densold=dens2d
        de=abs(ehfold)

        step=0

        print *,'Starting the SCF procedure...'

!*****Self-Consistent Field (SCF) procedure*****
        do while (de.gt.econv.OR.densres.gt.densconv.AND.step.lt.maxstep)

                step=step+1

                print *,' '

                print *,'SCF step',step

                onee=0.0
                twoe=0.0

        !*****Building a Fock matrix*****
                do mm=1,nbasorb
                 do nn=1,nbasorb

                  onee(mm,nn)=core2d(mm,nn)

                 end do
                end do

                do mm=1,nbasorb
                 do nn=1,nbasorb
                  do ll=1,nbasorb
                   do dd=1,nbasorb

                    twoe(mm,nn)=twoe(mm,nn)+densold(dd,ll)*(veff4d(mm,nn,ll,dd)-0.5*veff4d(mm,dd,ll,nn))

                   end do
                  end do
                 end do
                end do

                fock=onee+twoe

        !*****Forming F' *****
                fockp=matmul(ovrlpinv2d,matmul(fock,ovrlpinv2d))

        !*****Diagonalizing F' and getting eigenvalues and the coeff. matrix*****
                i=23+step

                call eig(fockp,constp,i,nbasorb,0)

                print *,'******************************************'
                print *,'                Eigenvalues               '
                print *,'******************************************'
                do i=1,nmolhlfocc
                 print *,i,' ',fockp(i,i)
                end do
                print *,'------------------------------------------'
                do i=(nmolhlfocc+1),nbasorb
                 print *,i,' ',fockp(i,i)
                end do
                print *,'******************************************'

        !*****Forming C=S^(-1/2)*C' *****
                const=matmul(ovrlpinv2d,constp)

        !*****Building the new charge density matrix*****
                print *,'*************************'
                print *,'     Density Matrix      '
                print *,'*************************'

                do mm=1,nbasorb
                 do nn=1,nbasorb

                  densnew(mm,nn)=0

                  do i=1,nmolhlfocc

                   densnew(mm,nn)=densnew(mm,nn)+(2.0*const(mm,i)*const(nn,i))

                  end do

                  if (nn.eq.nbasorb) then

                   print *,densnew(mm,nn)

                  else

                   write(*,"(f8.5,1x)",advance="no")densnew(mm,nn)

                  end if

                 end do
                end do

                print *,'*************************'

        !*****Building density residue matrix*****
                do mm=1,nbasorb
                 do nn=1,nmolhlfocc

                  densresm(mm,nn)=abs(densold(mm,nn)-densnew(mm,nn))

                 end do
                end do

        !*****Calculating the new 1-e, 2-e contributions and the HF energy*****
                ehfnew=vnn

                do mm=1,nbasorb
                 do nn=1,nbasorb

                  ehfnew=ehfnew+densold(nn,mm)*(onee(mm,nn)+0.5*twoe(mm,nn))

                 end do
                end do

        !*****Convergence check and energy, density updates*****
                de=abs(ehfnew-ehfold)
                densres=maxval(densresm)

                densold=densnew
                ehfold=ehfnew

                print *,'delta(E(HF)) =',de
                print *,'max(P(res)) =',densres 
                print *,'E(HF) for step',step,'is =',ehfnew

        end do

        if (de.lt.econv.and.densres.lt.densconv) then

         print *,'SCF converged in',step,'steps!'

!*****Converged HF density matrix and energy update*****
         dens2d=densold
         ehf=ehfold

         print *,'At convergence:'
         print *,'E(tot) =',ehf,'Ha'
         print *,'E(HF) =',(ehf-vnn),'Ha'
         print *,'largest density residue matrix element, max(P(res)) =',densres

!*****Partial transformation of the 2-e integrals from AO to MO basis*****
        xx=nmolhlfocc
        yy=nbasorb
        zz=yy-xx

        allocate(ijabD(xx,(xx+1):yy,xx,(xx+1):yy))
        allocate(ilndD(xx,yy,yy,yy))
        allocate(iladD(xx,(xx+1):yy,yy,yy))
        allocate(ijadD(xx,(xx+1):yy,xx,yy))

        ijabD=0.0
        ilndD=0.0
        iladD=0.0
        ijadD=0.0

!        do i=1,xx
!         do a=(xx+1),yy
!          do j=1,xx
!           do b=(xx+1),yy
!
!            ii=0.0
!
!            do mm=1,nbasorb
!             do nn=1,nbasorb
!              do ll=1,nbasorb
!               do dd=1,nbasorb
!
!                ii=ii+const(dd,b)*const(ll,j)*const(nn,a)*const(mm,i)*veff4d(mm,nn,ll,dd)
!
!               end do
!              end do
!             end do
!            end do
!           end do
!          end do
!         end do
!        end do

        do i=1,xx
         do nn=1,yy
          do ll=1,yy
           do dd=1,yy

            ii=0.0

            do mm=1,yy

             ii=ii+const(mm,i)*veff4d(mm,nn,ll,dd)

            end do

            ilndD(i,nn,ll,dd)=ii

           end do
          end do
         end do
        end do

        do a=(xx+1),yy
         do i=1,xx
          do ll=1,yy
           do dd=1,yy

            aa=0.0

            do nn=1,yy

             aa=aa+const(nn,a)*ilndD(i,nn,ll,dd)

            end do

            iladD(i,a,ll,dd)=aa

           end do
          end do
         end do
        end do

        do j=1,xx
         do a=(xx+1),yy
          do i=1,xx
           do dd=1,yy

            jj=0.0

            do ll=1,yy

             jj=jj+const(ll,j)*iladD(i,a,ll,dd)

            end do

            ijadD(i,a,j,dd)=jj

           end do
          end do
         end do
        end do

        do b=(xx+1),yy
         do j=1,xx
          do a=(xx+1),yy
           do i=1,xx

            bb=0.0

            do dd=1,yy

             bb=bb+const(dd,b)*ijadD(i,a,j,dd)

            end do

            ijabD(i,a,j,b)=bb

           end do
          end do
         end do
        end do

!*****MBPT II order energy correction*****
        ehfnew=0.0
        ii=0.0

        do i=1,xx
         do a=(xx+1),yy
          do j=1,xx
           do b=(xx+1),yy

            de=fockp(i,i)+fockp(j,j)-fockp(a,a)-fockp(b,b)
            ii=2.0*(ijabD(i,a,j,b)**(2.0))

            ehfnew=ehfnew+ii/de

           end do
          end do
         end do
        end do

        do i=1,xx
         do a=(xx+1),yy
          do j=1,xx
           do b=(xx+1),yy

            de=fockp(i,i)+fockp(j,j)-fockp(a,a)-fockp(b,b)
            ii=ijabD(i,a,j,b)*ijabD(i,b,j,a)

            ehfnew=ehfnew-ii/de

           end do
          end do
         end do
        end do

        print *,'E(MP2) =',ehfnew
        print *,'E(tot) =',(ehfnew+ehf)

        deallocate(ilndD,iladD,ijadD)

!*****Some more partial transformation of 2-e integrals from  AO to MO*****
        allocate(ajibD((xx+1):yy,xx,xx,(xx+1):yy))
        allocate(ajbiD((xx+1):yy,(xx+1):yy,xx,xx))
        allocate(alndD((xx+1):yy,yy,yy,yy))
        allocate(alidD((xx+1):yy,yy,xx,yy))
        allocate(ajidD((xx+1):yy,xx,xx,yy))

        do a=(xx+1),yy
         do nn=1,yy
          do ll=1,yy
           do dd=1,yy

            aa=0.0

            do mm=1,yy

             aa=aa+const(mm,a)*veff4d(mm,nn,ll,dd)

            end do

            alndD(a,nn,ll,dd)=aa

           end do
          end do
         end do
        end do

        do i=1,xx
         do a=(xx+1),yy
          do ll=1,yy
           do dd=1,yy

            ii=0.0

            do nn=1,yy

             ii=ii+const(nn,i)*alndD(a,nn,ll,dd)

            end do

            alidD(a,i,ll,dd)=ii

           end do
          end do
         end do
        end do

        do j=1,xx
         do i=1,xx
          do a=(xx+1),yy
           do dd=1,yy

            jj=0.0

            do ll=1,yy

             jj=jj+const(ll,j)*alidD(a,i,ll,dd)

            end do

            ajidD(a,i,j,dd)=jj

           end do
          end do
         end do
        end do

        do b=(xx+1),yy
         do j=1,xx
          do i=1,xx
           do a=(xx+1),yy

            bb=0.0

            do dd=1,yy

             bb=bb+const(dd,b)*ajidD(a,i,j,dd)

            end do

            ajibD(a,i,j,b)=bb

           end do
          end do
         end do
        end do

        deallocate(ajidD,alidD)

        allocate(albdD((xx+1):yy,xx,yy,yy))
        allocate(ajbdD((xx+1):yy,(xx+1):yy,xx,yy))
        
        do b=(xx+1),yy
         do a=(xx+1),yy
          do ll=1,yy
           do dd=1,yy

            bb=0.0

            do nn=1,yy

             bb=bb+const(nn,b)*alndD(a,nn,ll,dd)

            end do

            albdD(a,b,ll,dd)=bb

           end do
          end do
         end do
        end do

        do j=1,xx
         do b=(xx+1),yy
          do a=(xx+1),yy
           do dd=1,yy

            jj=0.0

            do ll=1,yy

             jj=jj+const(ll,j)*albD(a,b,ll,dd)

            end do

            ajbdD(a,b,j,dd)=jj

           end do
          end do
         end do
        end do


        do i=1,xx
         do j=1,xx
          do b=(xx+1),yy
           do a=(xx+1),yy

            ii=0.0

            do dd=1,yy

             ii=ii+const(dd,i)*ajbdD(a,b,j,dd)

            end do

            ajbiD(a,b,j,i)=ii

           end do
          end do
         end do
        end do

        deallocate(ajbdD,albdD,alndD)

!*****CI Singles Hamiltonian*****
        allocate(AmatSing(xx*zz,xx*zz))
        allocate(AmatTrip(xx*zz,xx*zz))
        allocate(cisSing(xx*zz,xx*zz))
        allocate(cisTrip(xx*zz,xx*zz))

        ia=0

        do i=1,xx
         do a=xx+1,zz

          ia=ia+1
          jb=0

          do j=1,xx
           do b=xx+1,zz

            jb=jb+1

            if (i.eq.j.AND.a.eq.b) then

             de=(fockp(a,a)-fockp(i,i))

            else

             de=0

            end if

            AmatSing(ia,jb)=de+2*ajibD(a,i,j,b)-ajbiD(a,b,j,i)
            AmatTrip(ia,jb)=de-ajbiD(a,b,j,i)

           end do
          end do

         end do
        end do

!*****Energies for spin-adapted CIS excitations, triplets appear thrice*****
        call eig(AmatSing,cisSing,66,xx*zz,0)
        call eig(AmatTrip,cisTrip,67,xx*zz,0)

        print *,'first',nsing,'singlet excitations'

        do i=1,nsing

         print *,'level #',i,'E(CIS)=',AmatSing(i,i)

        end do

        print *,'first',ntrip,'triplet excitations'

        do i=1,ntrip

         print *,'level #',i,'E(CIS)=',AmatTrip(i,i)

        end do

!*****SCF fails to converge*****
        else

         print *,'SCF did not converge!'

        end if

!*****Deallocation of memory*****
        deallocate(core2d,ovrlp2d)
        deallocate(veff4d)
        deallocate(ijabD,ajibD,ajbiD)
        deallocate(AmatSing,AmatTrip)
        deallocate(cisSing,cisTrip)
        deallocate(transform)
        deallocate(fock,fockp)
        deallocate(constp,const)
        deallocate(dens2d,ovrlpinv2d)
        deallocate(densold,densnew)
        deallocate(densresm)
        deallocate(onee,twoe)
        deallocate(eigen)

        print *,"final deallocation check"

call aces_ja_fin
call aces_fin

end program
