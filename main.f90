program main

    use omp_lib

    implicit none

    integer:: N1,N2,N3
    integer:: i,ii,j,j1,j2,j3
    character(len=60):: file_in,file_out,cwd
    real(8),parameter:: pi=4.*atan(1.)
    real(8),dimension(3):: xmin,xmax,node
    real(8):: dx
    integer:: Ngrain,istat
    real(8),dimension(:),allocatable:: t1
    integer,dimension(:,:,:),allocatable:: cellid
    real(8),dimension(:,:),allocatable:: X
    integer,dimension(:),allocatable:: voxelid
    integer,dimension(:,:),allocatable:: badgrain
    integer,dimension(:),allocatable:: gc,neighb
    integer,dimension(:,:),allocatable:: gneighbor
    integer,dimension(:,:),allocatable:: grain_remain
    real(8),dimension(:),allocatable:: SG_vol
    integer:: pt_gc,pt_voxelid,pt_gc0
    integer:: cnt2
    logical:: icon,ineigh
    integer:: c,NgrainSG
    real(8),parameter:: big=1000000.
    integer:: N1_start,N1_end,dN1
    integer:: Ngrain_start,Ngrain_end,dNgrain
    integer:: nthreads,tid

    call getarg(1,file_in)
    file_in=trim(file_in)
    call getarg(2,file_out)
    file_out=trim(file_out)

    ! Start reading input parameters
    !file_in='Small_IN100_20140317_XSG_BF.csv'
    call getcwd(cwd)
    open(unit=10,file=file_in,status='old',action='read',iostat=istat)
    !open(unit=10,file=trim(cwd)//'/Data/'//trim(file_in),status='old',action='read',iostat=istat)
    if (istat/=0) stop "*** The input file could not be opened ***"
    read(10,*) Ngrain
    read(10,*) dx
    read(10,*) (xmin(i),i=1,3)
    read(10,*) (xmax(i),i=1,3)

    ! Memory allocation
    allocate(X(Ngrain,4),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for X ***"

    ! Finish reading input paramters
    do i=1,Ngrain
        read(10,*) (X(i,j),j=1,4)
    enddo
    close(10)

    ! Calculate the number of voxels in each direction
    N1=nint((xmax(1)-xmin(1))/dx)
    N2=nint((xmax(2)-xmin(2))/dx)
    N3=nint((xmax(3)-xmin(3))/dx)

    ! Memory allocation
    allocate(cellid(N1,N2,N3),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for cellid ***"
    allocate(t1(Ngrain),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for t1 ***"

    ! Initialization
    cellid(:,:,:)=0
    t1(:)=0.

    !$omp parallel &
    !$omp shared (cellid) &
    !$omp private (j1,j2,j3,N1_start,N1_end,node,t1)

    nthreads=omp_get_num_threads()

    dN1=N1/nthreads

    tid=omp_get_thread_num()

    N1_start=tid*dN1+1

    if (tid==(nthreads-1)) then
        N1_end=N1
    else
        N1_end=N1_start+dN1
    end if

    ! For each voxel, find the grain it is reached by first
    do j1=N1_start,N1_end
        do j2=1,N2
            do j3=1,N3
                node=(/ (j1-.5)*dx+xmin(1),(j2-.5)*dx+xmin(2),(j3-.5)*dx+xmin(3) /)
                t1=(node(1)-X(:,1))**2+(node(2)-X(:,2))**2+(node(3)-X(:,3))**2
                t1=t1/X(:,4)**2
                cellid(j1,j2,j3)=minloc(t1,1)
            enddo
        enddo
    enddo

    !$omp end parallel

    ! Memory allocation
    allocate(voxelid(N1*N2*N3),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for voxelid ***"
    allocate(gc(N1*N2*N3),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for gc ***"
    allocate(neighb(N1*N2*N3),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for neighb ***"
    allocate(badgrain(N1*N2*N3,15),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for badgrain ***"

    ! Initialization
    badgrain(:,:)=0
    gc(:)=0

    c=0
    pt_voxelid=N1*N2*N3
    do while (pt_voxelid/=0)
        c=c+1
        write(*,*)'c = ',c,', pt_voxelid = ',pt_voxelid
        voxelid(:)=0
        pt_voxelid=0

        !$omp parallel &
        !$omp shared (pt_voxelid,voxelid,badgrain) &
        !$omp private (tid,Ngrain_start,Ngrain_end,j,pt_gc,pt_gc0,gc,neighb,j1,j2,j3,cnt2,i,ii,icon)

        nthreads=omp_get_num_threads()

        dNgrain=Ngrain/nthreads

        tid=omp_get_thread_num()

        Ngrain_start=tid*dNgrain+1

        if (tid==(nthreads-1)) then
            Ngrain_end=Ngrain
        else
            Ngrain_end=Ngrain_start+dNgrain
        end if


        ! Loop over all the grains
        do j=Ngrain_start,Ngrain_end
            pt_gc=0
            pt_gc0=0
            gc(:)=0
            neighb(:)=0
            j1=ceiling((X(j,1)-xmin(1))/dx)
            j2=ceiling((X(j,2)-xmin(2))/dx)
            j3=ceiling((X(j,3)-xmin(3))/dx)
            neighb(N1*N2*(j3-1)+N1*(j2-1)+j1)=1
            ! Collect all the voxels associated to grain j and store them in gc
            do j1=1,N1
                do j2=1,N2
                    do j3=1,N3
                        if (cellid(j1,j2,j3)==j.and.neighb(N1*N2*(j3-1)+N1*(j2-1)+j1)/=1) then
                            pt_gc=pt_gc+1
                            gc(pt_gc)=N1*N2*(j3-1)+N1*(j2-1)+j1
                        endif
                    enddo !j3=1,N3
                enddo !j2=1,N2
            enddo !j1=1,N1
            cnt2=1
            ! Loop over the voxels in gc(:,:) until each voxel is found to be either:
            !   - simply connected with the grain and thus stored in neighb(:,:) or,
            !   - disconnected from the remaining of the grain and thus stored in voxelid(:,:)
            do while (cnt2/=0)
                cnt2=0
                pt_gc0=pt_gc
                ! Loop over the voxels in gc(:)
                do i=1,pt_gc0
                    ii=i+pt_gc-pt_gc0
                    if (ii>pt_gc) then
                        exit
                    else
                        icon=.false.
                        ! Verify if the voxel gc(i) is connected to the growing set neighb(:) which currently
                        ! contains some of the voxels that are simply connected to the center of the grain
                        if ((gc(ii)<N1*N2*N3).and.(neighb(gc(ii)+1)==1)) then
                            icon=.true.
                        elseif ((gc(ii)>1).and.(neighb(gc(ii)-1)==1)) then
                            icon=.true.
                        elseif ((gc(ii)<(N1*(N2*N3-1)+1)).and.(neighb(gc(ii)+N1)==1)) then
                            icon=.true.
                        elseif ((gc(ii)>N1).and.(neighb(gc(ii)-N1)==1)) then
                            icon=.true.
                        elseif ((gc(ii)<(N1*N2*(N3-1)+1)).and.(neighb(gc(ii)+N1*N2)==1)) then
                            icon=.true.
                        elseif ((gc(ii)>N1*N2).and.(neighb(gc(ii)-N1*N2)==1)) then
                            icon=.true.
                        endif
                        ! If the voxel gc(i) is connected the current set of identified connected voxels, then add it to neighb(:)
                        if (icon) then
                            cnt2=cnt2+1
                            neighb(gc(ii))=1
                            pt_gc=pt_gc-1
                            gc(ii:pt_gc)=gc((ii+1):(pt_gc+1))
                        endif
                    endif
                enddo !i=1,pt_gc0
            enddo !while (cnt2/=0)
            do i=1,pt_gc
                !$omp critical
                pt_voxelid=pt_voxelid+1
                voxelid(pt_voxelid)=gc(i)
                !$omp end critical
                badgrain(gc(i),c)=j
            enddo !i=1,pt_gc
        enddo !j=1,Ngrain

        !$omp end parallel

        ! Loop over the pt_voxelid voxels which are not simply connected to the remaining of the
        ! grain they were associated to
        do j=1,pt_voxelid
            j3=int((voxelid(j)-1)/float(N1*N2))+1
            j2=int((voxelid(j)-1-N1*N2*(j3-1))/float(N1))+1
            j1=voxelid(j)-N1*N2*(j3-1)-N1*(j2-1)
            node=(/ (j1-.5)*dx+xmin(1),(j2-.5)*dx+xmin(2),(j3-.5)*dx+xmin(3) /)
            t1=(node(1)-X(:,1))**2+(node(2)-X(:,2))**2+(node(3)-X(:,3))**2
            t1=t1/X(:,4)**2
            do i=1,c
                if (badgrain(voxelid(j),i)/=0) then
                    t1(badgrain(voxelid(j),i))=big
                end if
            enddo !i=1,c
            cellid(j1,j2,j3)=minloc(t1,1)
        enddo !j=1,pt_voxelid
    enddo !while (pt_voxelid/=0)

    ! Memory deallocation
    deallocate(t1,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for t1 ***"
    deallocate(voxelid,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for voxelid ***"
    deallocate(gc,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for gc ***"
    deallocate(neighb,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for neighb ***"
    deallocate(badgrain,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for badgrain ***"

    ! Memory allocation
    allocate(gneighbor(Ngrain,Ngrain),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for gneighbor ***"

    ! Initialization
    gneighbor(:,:)=0

    ! Identify the grains that are connected to one grain only
    do j1=1,(N1-1)
        do j2=1,(N2-1)
            do j3=1,(N3-1)
                if (cellid(j1,j2,j3)/=cellid(j1+1,j2,j3)) then
                    ineigh=.false.
                    do i=2,(1+gneighbor(cellid(j1,j2,j3),1))
                        if (gneighbor(cellid(j1,j2,j3),i)==cellid(j1+1,j2,j3)) then
                            ineigh=.true.
                            exit
                        endif
                    enddo
                    if (.not.ineigh) then
                        gneighbor(cellid(j1,j2,j3),1)=gneighbor(cellid(j1,j2,j3),1)+1
                        gneighbor(cellid(j1,j2,j3),gneighbor(cellid(j1,j2,j3),1)+1)=cellid(j1+1,j2,j3)
                    endif
                    ineigh=.false.
                    do i=2,(1+gneighbor(cellid(j1+1,j2,j3),1))
                        if (gneighbor(cellid(j1+1,j2,j3),i)==cellid(j1,j2,j3)) then
                            ineigh=.true.
                            exit
                        endif
                    enddo
                    if (.not.ineigh) then
                        gneighbor(cellid(j1+1,j2,j3),1)=gneighbor(cellid(j1+1,j2,j3),1)+1
                        gneighbor(cellid(j1+1,j2,j3),gneighbor(cellid(j1+1,j2,j3),1)+1)=cellid(j1,j2,j3)
                    endif
                endif
                if (cellid(j1,j2,j3)/=cellid(j1,j2+1,j3)) then
                    ineigh=.false.
                    do i=2,(1+gneighbor(cellid(j1,j2,j3),1))
                        if (gneighbor(cellid(j1,j2,j3),i)==cellid(j1,j2+1,j3)) then
                            ineigh=.true.
                            exit
                        endif
                    enddo
                    if (.not.ineigh) then
                        gneighbor(cellid(j1,j2,j3),1)=gneighbor(cellid(j1,j2,j3),1)+1
                        gneighbor(cellid(j1,j2,j3),gneighbor(cellid(j1,j2,j3),1)+1)=cellid(j1,j2+1,j3)
                    endif
                    ineigh=.false.
                    do i=2,(1+gneighbor(cellid(j1,j2+1,j3),1))
                        if (gneighbor(cellid(j1,j2+1,j3),i)==cellid(j1,j2,j3)) then
                            ineigh=.true.
                            exit
                        endif
                    enddo
                    if (.not.ineigh) then
                        gneighbor(cellid(j1,j2+1,j3),1)=gneighbor(cellid(j1,j2+1,j3),1)+1
                        gneighbor(cellid(j1,j2+1,j3),gneighbor(cellid(j1,j2+1,j3),1)+1)=cellid(j1,j2,j3)
                    endif
                endif
                if (cellid(j1,j2,j3)/=cellid(j1,j2,j3+1)) then
                    ineigh=.false.
                    do i=2,(1+gneighbor(cellid(j1,j2,j3),1))
                        if (gneighbor(cellid(j1,j2,j3),i)==cellid(j1,j2,j3+1)) then
                            ineigh=.true.
                            exit
                        endif
                    enddo
                    if (.not.ineigh) then
                        gneighbor(cellid(j1,j2,j3),1)=gneighbor(cellid(j1,j2,j3),1)+1
                        gneighbor(cellid(j1,j2,j3),gneighbor(cellid(j1,j2,j3),1)+1)=cellid(j1,j2,j3+1)
                    endif
                    ineigh=.false.
                    do i=2,(1+gneighbor(cellid(j1,j2,j3+1),1))
                        if (gneighbor(cellid(j1,j2,j3+1),i)==cellid(j1,j2,j3)) then
                            ineigh=.true.
                            exit
                        endif
                    enddo
                    if (.not.ineigh) then
                        gneighbor(cellid(j1,j2,j3+1),1)=gneighbor(cellid(j1,j2,j3+1),1)+1
                        gneighbor(cellid(j1,j2,j3+1),gneighbor(cellid(j1,j2,j3+1),1)+1)=cellid(j1,j2,j3)
                    endif
                endif
            enddo
        enddo
    enddo

    ! Memory allocation
    allocate(grain_remain(Ngrain,2),stat=istat)
    if (istat/=0) stop "*** Could not allocate memory for grain_remain ***"

    ! Initialization
    grain_remain(:,:)=0

    NgrainSG=0
    do j=1,Ngrain
        if (gneighbor(j,1)>1) then
            NgrainSG=NgrainSG+1
            grain_remain(NgrainSG,1)=j
            grain_remain(j,2)=NgrainSG
        else
            do j1=1,N1
                do j2=1,N2
                    do j3=1,N3
                        if (cellid(j1,j2,j3)==j) then
                            cellid(j1,j2,j3)=gneighbor(j,2)
                        endif
                    enddo
                enddo
            enddo
        endif
    enddo

    ! Memory allocation
    allocate(SG_vol(NgrainSG),stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for SG_vol ***"

    ! Initialize
    SG_vol(:)=0.

    do j1=1,N1
        do j2=1,N2
            do j3=1,N3
                SG_vol(grain_remain(cellid(j1,j2,j3),2))=SG_vol(grain_remain(cellid(j1,j2,j3),2))+1
            enddo
        enddo
    enddo
    SG_vol=SG_vol*dx**3

    ! Memory deallocation
    deallocate(cellid,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for cellid ***"

    !file_out='Small_IN100_20140317_XSG_BF.out'
    open(unit=20,file=file_out,status='replace',action='write',iostat=istat)
    !open(unit=20,file=trim(cwd)//'/Data/'//trim(file_out),status='replace',action='write',iostat=istat)
    if (istat/=0) stop "*** The output file could not be created ***"
    write(20,*)NgrainSG
    do i=1,NgrainSG
        write(20,*)SG_vol(i),grain_remain(i,1),gneighbor(grain_remain(i,1),1),&
        (gneighbor(grain_remain(i,1),j),j=2,gneighbor(grain_remain(i,1),1)+1)
    enddo
    close(20)

    ! Memory deallocation
    deallocate(gneighbor,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for neighbor ***"
    deallocate(grain_remain,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for grain_remain ***"

    ! Memory deallocation
    deallocate(SG_vol,stat=istat)
    if (istat/=0) stop "*** Could not deallocate memory for SG_vol ***"

end program main
