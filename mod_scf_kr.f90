module glovar
integer :: dd,ll,mm,nn,i,j,k,l,a,b,c,d,n,step
real*8 :: xx,yy,zz
integer :: nbasorb,nmolhlfocc,nmolhlfunocc,nmolocc,nmolunocc,natom
real*8,dimension(:,:),allocatable ::core2d,ovrlp2d,dens2d,fock,eigen,fockp,const
real*8,dimension(:,:),allocatable :: constp,ovrlpinv2d,dummy3,dummy4,transform
real*8,dimension(:,:),allocatable :: densnew,densold,densresm,onee,twoe
real*8,dimension(:,:,:,:),allocatable :: veff4d
real*8,dimension(:),allocatable :: core1d,ovrlp1d
real*8,dimension(:),allocatable :: dummy1
integer,dimension(:),allocatable :: dummy2
real*8 :: vnn,ehf,de,ehfnew,ehfold,densres
end module
