program test_fusion_table_digest
 use iso_c_binding
 use fusion_birth_table_fortran
 use fusion_beam_birth_table_fortran
 implicit none
 integer(c_int8_t),allocatable,target::bytes(:)
 integer(c_int8_t)::digest(32)
 integer(c_size_t)::n,packed
 integer(c_int)::status
 type(c_ptr)::table
 integer::u,j,k
 character(4096)::path
 character(64)::hex
 do k=1,2
  call get_command_argument(k,path)
  open(newunit=u,file=trim(path),access='stream',form='unformatted',status='old')
  inquire(unit=u,size=n);allocate(bytes(n));read(u)bytes;close(u)
  if(k==1)then
   call fusion_birth_table_unpack(bytes,n,table,status)
  else
   call fusion_beam_birth_table_unpack(bytes,n,table,status)
  endif
  if(status/=0)error stop 'unpack'
  if(k==1)then
   call fusion_birth_table_content_digest(table,digest,packed,status)
   call fusion_birth_table_destroy(table)
  else
   call fusion_beam_birth_table_content_digest(table,digest,packed,status)
   call fusion_beam_birth_table_destroy(table)
  endif
  if(status/=0.or.packed/=n)error stop 'digest'
  do j=1,32
   write(hex(2*j-1:2*j),'(z2.2)')iand(int(digest(j)),255)
  enddo
  print '(i0,1x,i0,1x,a)',k,packed,hex
  digest=9;packed=9
  if(k==1)then
   call fusion_birth_table_content_digest(c_null_ptr,digest,packed,status)
  else
   call fusion_beam_birth_table_content_digest(c_null_ptr,digest,packed,status)
  endif
  if(status==0.or.packed/=0.or.any(digest/=0))error stop 'failure outputs'
  deallocate(bytes)
 enddo
 print *,'FORTRAN_TABLE_DIGEST_PASS'
end program
