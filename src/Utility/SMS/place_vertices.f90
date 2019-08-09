! Prototype to read in SMS 11.xx map file and manipulate vertices 
! Input: old.map
! Output: new.map

! pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o place_vertices place_vertices.f90

  integer, parameter :: mx_vert=10000 !max # of vertices along any arc
  integer, parameter :: mx_nodes=10000 !max # of nodes
  integer, parameter :: mx_arcs=10000 !max # of arcs

  character(len=100) :: line_str,str_tmp,str_tmp2
  character(len=2) :: char2
  character(len=5) :: char5
  integer :: line_node(2),line_arc(2),iarc_nodes(2,mx_arcs),nvert(mx_arcs),node_id(mx_nodes)
  real :: xy_nodes(2,mx_nodes),xy_vert(2,mx_vert,mx_arcs)

  !Test memory
  xy_vert(2,mx_vert,mx_arcs)=0
 
  ! Scan 
  open(16,file='new.map',status='replace')
  open(15,file='old.map',status='old')

  nline=0
  line_node(:)=0 !start and end line # for reading node section
  line_arc(:)=0
  ifl_node=0 !flag
  ifl_arc=0 !flag
  do
    read(15,'(a)',end=99)line_str
    nline=nline+1
    line_str=adjustl(line_str) !place blanks at end
    len_str=len_trim(line_str)

    if(len_str>=4.and.line_str(1:4).eq."NODE".and.ifl_node==0) then
      ifl_node=1 !only read the 1st
      line_node(1)=nline
    endif
    
    if(len_str>=3.and.line_str(1:3).eq."ARC".and.ifl_arc==0) then
      ifl_arc=1 !only read the 1st
      line_arc(1)=nline
      line_node(2)=nline-1 !assuming 2 sections r rite after each other
    endif

    if(len_str>=6.and.line_str(1:6).eq."ENDCOV") line_arc(2)=nline-1
  enddo

99   close(15)
  print*, nline,' lines read'
  print*, '1st n last lines of node list are ',line_node
  print*, '1st n last lines of arc list are ',line_arc

  if(line_node(1)*line_node(2)*line_arc(1)*line_arc(2)==0) then
    print*, 'Not found:',line_arc,line_node
    stop
  endif

  !Parse node list
  open(15,file='old.map',status='old')
  rewind(15)

  if(mod(line_node(2)-line_node(1)+1,4)/=0) stop 'Wrong node list'
  nnodes=(line_node(2)-line_node(1)+1)/4
  if(nnodes>mx_nodes) then
    print*, 'Too many nodes(1):',i
    stop 
  endif

  do i=1,line_node(1)-1
    read(15,*)
  enddo !i
    
  do i=1,nnodes
    read(15,*) !NODE
    read(15,*)char2,xy_nodes(1:2,i)
!    print*, char2,xy_nodes(1:2,nnodes)
    !Note that SMS does not number nodes consecutively
    read(15,*)char2,node_id(i) !ID
    read(15,*) !END
  enddo !i

  !Parse arc list
  rewind(15)
  do i=1,line_arc(1)-1
    read(15,*)
  enddo !i

  !# of arcs
  narcs=0  
  do i=line_arc(1),line_arc(2)
    read(15,*)line_str
    line_str=adjustl(line_str) !place blanks at end
    len_str=len_trim(line_str)
    if(len_str>=12.and.line_str(1:len_str).eq."ARCELEVATION") then
      narcs=narcs+1
      if(narcs>mx_arcs) then
        print*, 'Too many arcs:',i
        stop
      endif
    endif
  enddo !i
  print*, narcs,' arcs found'

  !Redo
  rewind(15)
  do i=1,line_arc(1)-1
    read(15,*)
  enddo !i

  nvert=0
  do i=1,narcs
    read(15,*) !ARC
    read(15,*) !ID
    read(15,*) !ARCELEVATION
    read(15,*)char5,iarc_nodes(1:2,i) !NODES
    print*, char5,iarc_nodes(1:2,i)
    read(15,'(a)')line_str 
    line_str=adjustl(line_str) !place blanks at end
    len_str=len_trim(line_str)
    print*, line_str !,line_str(1:len_str),len_str
    if(len_str>=3.and.line_str(1:3).eq."END") then
    else !ARCVERTICES
      loc=index(line_str,' ')
      read(line_str(loc+1:len_str),*)nvert(i)
      if(nvert(i)>mx_vert) then
        print*, 'Too many vertices:',i
        stop
      endif
      do k=1,nvert(i)
        read(15,*)xy_vert(1:2,k,i)
        print*, xy_vert(1:2,k,i)
      enddo !k
      read(15,*) !END
    endif
  enddo !i

  !Check
  rewind(15)
  do i=1,line_node(1)-1
    read(15,'(a)')line_str
    write(16,*)line_str
  enddo

  do i=1,nnodes
    write(16,*)'NODE'
    write(16,*)'XY',xy_nodes(1:2,i),0.
    write(16,*)'ID',node_id(i)
    write(16,*)'END'
  enddo !i

  do i=1,narcs
    write(16,*)'ARC'
    write(16,*)'ID',i
    write(16,*)'ARCELEVATION',0.
    write(16,*)'NODES',iarc_nodes(1:2,i)
    if(nvert(i)>0) then
      write(16,*)'ARCVERTICES',nvert(i)
      do k=1,nvert(i)
        write(16,*)xy_vert(1:2,k,i),0.
      enddo !k
    endif
    write(16,*)'END'
  enddo !i

  do i=line_node(1),line_arc(2)
    read(15,*)
  enddo !i
  do i=line_arc(2)+1,nline
    read(15,'(a)')line_str
    write(16,*)line_str
  enddo !i

  stop
  end
