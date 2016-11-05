! io_pnetcdf.F

! parallel NetCDF I/O

!_______________________________________________________________________
      subroutine def_var_pnetcdf(ncid,name,nvdims,vdims,varid,
     $                           long_name,units,coords,lcoords)
# include "mpif.h"
# include "pnetcdf.inc"
      integer vdims(nvdims)
      integer ncid,nvdims,varid
      logical lcoords
      character*(*) name,long_name,units,coords
      integer status
      integer(mpi_offset_kind) length
! define variable
      status=nfmpi_def_var(ncid,name,nf_double,nvdims,vdims,varid)
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
! define attributes
      length=len(trim(long_name))
      status=nfmpi_put_att_text(ncid,varid,'long_name',
     $                                           length,trim(long_name))
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

      length=len(trim(units))
      status=nfmpi_put_att_text(ncid,varid,'units',length,trim(units))
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
! add coordinates attribute, if necessary
      if(lcoords) then
        length=len(trim(coords))
        status=nfmpi_put_att_text(ncid,varid,'coordinates',
     $                                              length,trim(coords))
        call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
      end if
      return
      end

!_______________________________________________________________________
      subroutine handle_error_pnetcdf(routine,status,nf_noerr)
      implicit none
      include 'pom.h'
      integer status,nf_noerr
      character*(*) routine
      if(status.ne.nf_noerr) then
        error_status=1
        if(my_task.eq.master_task) write(*,'(/i4,a,a)')
     $                         status,'] Error: NetCDF routine ',routine
      end if
      return
      end

!_______________________________________________________________________
      subroutine write_output_pnetcdf
! write output data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character*120 netcdf_out_file,str_tmp
      integer nprint
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid,dum_varid,dvm_varid
      integer time_varid,uab_varid,vab_varid,elb_varid,u_varid,v_varid,
     $        w_varid,t_varid,s_varid,rho_varid
      integer ncid,status
      integer vdims(4)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(4),edge(4)

! create netcdf file
      nprint=(iint+int(time0*86400./dti))/iprint
      write(netcdf_out_file,'(a,''out/'',a,''.'',i4.4,''.nc'')')
     $                            trim(wrk_pth),trim(netcdf_file),nprint
      if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $                                             trim(netcdf_out_file)
      status=nfmpi_create(pom_comm,netcdf_out_file,
     $                    nf_clobber+nf_64bit_offset,mpi_info_null,ncid)
      call handle_error_pnetcdf('nf_create: '//netcdf_out_file,
     $                                                  status,nf_noerr)

! define global attributes
      length=len(trim(title))
      status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                                                      trim(title))
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

      str_tmp='output file'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                                                    trim(str_tmp))
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

! define dimensions
      length=1
      status=nfmpi_def_dim(ncid,'time',length,time_dimid)
      call handle_error_pnetcdf('nf_def_dim: time',status,nf_noerr)
      length=kb
      status=nfmpi_def_dim(ncid,'z',length,z_dimid)
      call handle_error_pnetcdf('nf_def_dim: z',status,nf_noerr)
      length=jm_global
      status=nfmpi_def_dim(ncid,'y',length,y_dimid)
      call handle_error_pnetcdf('nf_def_dim: y',status,nf_noerr)
      length=im_global
      status=nfmpi_def_dim(ncid,'x',length,x_dimid)
      call handle_error_pnetcdf('nf_def_dim: x',status,nf_noerr)

! define variables and their attributes
      vdims(1)=time_dimid
      str_tmp='days since '//time_start
      call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,
     $                                       'time',str_tmp,' ',.false.)

      vdims(1)=z_dimid
      call def_var_pnetcdf(ncid,'z',1,vdims,z_varid,
     $                   'sigma of cell face','sigma_level',' ',.false.)
      length=22
      status=nfmpi_put_att_text(ncid,z_varid,'standard_name',
     $                                  length,'ocean_sigma_coordinate')
      length=26
      status=nfmpi_put_att_text(ncid,z_varid,'formula_terms',
     $                              length,'sigma: z eta: elb depth: h')
      call handle_error_pnetcdf('nf_put_att_text: z',status,nf_noerr)

      call def_var_pnetcdf(ncid,'zz',1,vdims,zz_varid,
     $                'sigma of cell centre','sigma_level', ' ',.false.)
      length=22
      status=nfmpi_put_att_text(ncid,zz_varid,'standard_name',
     $                                  length,'ocean_sigma_coordinate')
      length=27
      status=nfmpi_put_att_text(ncid,zz_varid,'formula_terms',
     $                             length,'sigma: zz eta: elb depth: h')
      call handle_error_pnetcdf('nf_put_att_text: zz',status,nf_noerr)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      call def_var_pnetcdf(ncid,'dx',2,vdims,dx_varid,
     $            'grid increment in x','metre','east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'dy',2,vdims,dy_varid,
     $            'grid increment in y','metre','east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'east_u',2,vdims,east_u_varid,
     $            'easting of u-points','metre','east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'east_v',2,vdims,east_v_varid,
     $            'easting of v-points','metre','east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'east_e',2,vdims,east_e_varid,
     $    'easting of elevation points','metre','east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'east_c',2,vdims,east_c_varid,
     $        'easting of cell corners','metre','east_c north_c',.true.)
      call def_var_pnetcdf(ncid,'north_u',2,vdims,north_u_varid,
     $           'northing of u-points','metre','east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'north_v',2,vdims,north_v_varid,
     $           'northing of v-points','metre','east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'north_e',2,vdims,north_e_varid,
     $   'northing of elevation points','metre','east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'north_c',2,vdims,north_c_varid,
     $       'northing of cell corners','metre','east_c north_c',.true.)
      call def_var_pnetcdf(ncid,'rot',2,vdims,rot_varid,
     $                    'Rotation angle of x-axis wrt. east','degree',
     $                                          'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'h',2,vdims,h_varid,
     $        'undisturbed water depth','metre','east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'fsm',2,vdims,fsm_varid,
     $      'free surface mask','dimensionless','east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'dum',2,vdims,dum_varid,
     $        'u-velocity mask','dimensionless','east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'dvm',2,vdims,dvm_varid,
     $        'v-velocity mask','dimensionless','east_v north_v',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=time_dimid
      call def_var_pnetcdf(ncid,'uab',3,vdims,uab_varid,
     $           'depth-averaged u','metre/sec','east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'vab',3,vdims,vab_varid,
     $           'depth-averaged v','metre/sec','east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'elb',3,vdims,elb_varid,
     $              'surface elevation','metre','east_e north_e',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=z_dimid
      vdims(4)=time_dimid
      call def_var_pnetcdf(ncid,'u',4,vdims,u_varid,
     $              'x-velocity','metre/sec','east_u north_u zz',.true.)
      call def_var_pnetcdf(ncid,'v',4,vdims,v_varid,
     $             'y-velocity','metre/sec', 'east_v north_v zz',.true.)
      call def_var_pnetcdf(ncid,'w',4,vdims,w_varid,
     $               'z-velocity','metre/sec','east_e north_e z',.true.)
      call def_var_pnetcdf(ncid,'t',4,vdims,t_varid,
     $           'potential temperature','K','east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'s',4,vdims,s_varid,
     $       'salinity x rho / rhoref','PSS','east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'rho',4,vdims,rho_varid,
     $                          '(density-1000)/rhoref','dimensionless',
     $                                       'east_e north_e zz',.true.)

! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf('nf_enddef: output',status,nf_noerr)

! write data
      start(1)=1
      edge(1)=1
      status=nfmpi_put_vara_double_all(ncid,time_varid,start,edge,time)
      call handle_error_pnetcdf('nf_put_vara_double:time',
     $                                                  status,nf_noerr)

      start(1)=1
      edge(1)=kb
      status=nfmpi_put_vara_double_all(ncid,z_varid,start,edge,z)
      call handle_error_pnetcdf('nf_put_var_double: z',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,zz_varid,start,edge,zz)
      call handle_error_pnetcdf('nf_put_var_double: zz',status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_put_vara_double_all(ncid,dx_varid,start,edge,dx)
      call handle_error_pnetcdf('nf_put_var_double: dx',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,dy_varid,start,edge,dy)
      call handle_error_pnetcdf('nf_put_var_double: dy',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,east_u_varid,start,edge,
     $                                                           east_u)
      call handle_error_pnetcdf('nf_put_var_double: east_u',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,east_v_varid,start,edge,
     $                                                           east_v)
      call handle_error_pnetcdf('nf_put_var_double: east_v',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,east_e_varid,start,edge,
     $                                                           east_e)
      call handle_error_pnetcdf('nf_put_var_double: east_e',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,east_c_varid,start,edge,
     $                                                           east_c)
      call handle_error_pnetcdf('nf_put_var_double: east_c',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,north_u_varid,start,edge,
     $                                                          north_u)
      call handle_error_pnetcdf('nf_put_var_double: north_u',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,north_v_varid,start,edge,
     $                                                          north_v)
      call handle_error_pnetcdf('nf_put_var_double: north_v',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,north_e_varid,start,edge,
     $                                                          north_e)
      call handle_error_pnetcdf('nf_put_var_double: north_e',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,north_c_varid,start,edge,
     $                                                          north_c)
      call handle_error_pnetcdf('nf_put_var_double: north_c',status,
     $                                                         nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,rot_varid,start,edge,rot)
      call handle_error_pnetcdf('nf_put_var_double: rot',
     $                                                  status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,h_varid,start,edge,h)
      call handle_error_pnetcdf('nf_put_var_double: h',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,fsm_varid,start,edge,fsm)
      call handle_error_pnetcdf('nf_put_var_double: fsm',
     $                                                  status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,dum_varid,start,edge,dum)
      call handle_error_pnetcdf('nf_put_var_double: dum',
     $                                                  status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,dvm_varid,start,edge,dvm)
      call handle_error_pnetcdf('nf_put_var_double: dvm',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=1
      status=nfmpi_put_vara_double_all(ncid,uab_varid,start,edge,uab)
      call handle_error_pnetcdf('nf_put_vara_double: uab',
     $                                                  status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,vab_varid,start,edge,vab)
      call handle_error_pnetcdf('nf_put_vara_double: vab',
     $                                                  status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,elb_varid,start,edge,elb)
      call handle_error_pnetcdf('nf_put_vara_double: elb',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      start(4)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      edge(4)=1
      status=nfmpi_put_vara_double_all(ncid,u_varid,start,edge,u)
      call handle_error_pnetcdf('nf_put_vara_double: u',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,v_varid,start,edge,v)
      call handle_error_pnetcdf('nf_put_vara_double: v',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,w_varid,start,edge,w)
      call handle_error_pnetcdf('nf_put_vara_double: w',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,t_varid,start,edge,t)
      call handle_error_pnetcdf('nf_put_vara_double: t',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,s_varid,start,edge,s)
      call handle_error_pnetcdf('nf_put_vara_double: s',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,rho_varid,start,edge,rho)
      call handle_error_pnetcdf('nf_put_vara_double: rho',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: output',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine write_restart_pnetcdf
! write data for a seamless restart
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character*120 netcdf_out_file,str_tmp
      integer nprint
      integer ncid,status
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer time_varid,wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)

! create netcdf restart file
      nprint=(iint+int(time0*86400./dti))/irestart
      write(netcdf_out_file,'(a,''out/'',a,''.'',i4.4,''.nc'')')
     $                         trim(wrk_pth),trim(write_rst_file),nprint
      if(my_task.eq.0) write(*,'(/''writing file '',a)')
     $                                             trim(netcdf_out_file)
      status=nfmpi_create(pom_comm,netcdf_out_file,
     $                    nf_clobber+nf_64bit_offset,mpi_info_null,ncid)
      call handle_error_pnetcdf('nf_create',status,nf_noerr)

! define global attributes
      length=len(trim(title))
      status=nfmpi_put_att_text(ncid,nf_global,'title',length,
     $                                                      trim(title))
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)

      str_tmp='restart file'
      length=len(trim(str_tmp))
      status=nfmpi_put_att_text(ncid,nf_global,'description',length,
     $                                                    trim(str_tmp))
      call handle_error_pnetcdf('nf_put_att_text',status,nf_noerr)
! define dimensions
      length=1
      status=nfmpi_def_dim(ncid,'time',length,time_dimid)
      call handle_error_pnetcdf('nf_def_dim',status,nf_noerr)

      length=kb
      status=nfmpi_def_dim(ncid,'z',length,z_dimid)
      call handle_error_pnetcdf('nf_def_dim',status,nf_noerr)

      length=jm_global
      status=nfmpi_def_dim(ncid,'y',length,y_dimid)
      call handle_error_pnetcdf('nf_def_dim',status,nf_noerr)

      length=im_global
      status=nfmpi_def_dim(ncid,'x',length,x_dimid)
      call handle_error_pnetcdf('nf_def_dim',status,nf_noerr)

! define variables and their attributes
      vdims(1)=time_dimid
      str_tmp='days since '//time_start
      call def_var_pnetcdf(ncid,'time',1,vdims,time_varid,
     $                                       'time',str_tmp,' ',.false.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      call def_var_pnetcdf(ncid,'wubot',2,vdims,wubot_varid,
     $                     'x-momentum flux at the bottom',
     $                     'metre^2/sec^2',
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'wvbot',2,vdims,wvbot_varid,
     $                     'y-momentum flux at the bottom',
     $                     'metre^2/sec^2',
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'aam2d',2,vdims,aam2d_varid,
     $                     'vertical average of aam',
     $                     'metre^2/sec',
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'ua',2,vdims,ua_varid,
     $                     'vertical mean of u',
     $                     'metre/sec',
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'uab',2,vdims,uab_varid,
     $                     'vertical mean of u at time -dt',
     $                     'metre/sec',
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'va',2,vdims,va_varid,
     $                     'vertical mean of v',
     $                     'metre/sec',
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'vab',2,vdims,vab_varid,
     $                     'vertical mean of v at time -dt',
     $                     'metre/sec',
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'el',2,vdims,el_varid,
     $                     'surface elevation in external mode',
     $                     'metre',
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'elb',2,vdims,elb_varid,
     $                     'surface elevation in external mode at -dt',
     $                     'metre',
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'et',2,vdims,et_varid,
     $                     'surface elevation in internal mode',
     $                     'metre',
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'etb',2,vdims,etb_varid,
     $                     'surface elevation in internal mode at -dt',
     $                     'metre',
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'egb',2,vdims,egb_varid,
     $                     'surface elevation for pres. grad. at -dt',
     $                     'metre',
     $                     'east_e north_e',.true.)
      call def_var_pnetcdf(ncid,'utb',2,vdims,utb_varid,
     $                     'ua time averaged over dti',
     $                     'metre/sec',
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'vtb',2,vdims,vtb_varid,
     $                     'va time averaged over dti',
     $                     'metre/sec',
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'adx2d',2,vdims,adx2d_varid,
     $                     'vertical integral of advx',
     $                     '-',
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'ady2d',2,vdims,ady2d_varid,
     $                     'vertical integral of advy',
     $                     '-',
     $                     'east_v north_v',.true.)
      call def_var_pnetcdf(ncid,'advua',2,vdims,advua_varid,
     $                     'sum of 2nd, 3rd and 4th terms in eq (18)',
     $                     '-',
     $                     'east_u north_u',.true.)
      call def_var_pnetcdf(ncid,'advva',2,vdims,advva_varid,
     $                     'sum of 2nd, 3rd and 4th terms in eq (19)',
     $                     '-',
     $                     'east_v north_v',.true.)

      vdims(1)=x_dimid
      vdims(2)=y_dimid
      vdims(3)=z_dimid
      call def_var_pnetcdf(ncid,'u',3,vdims,u_varid,
     $                     'x-velocity',
     $                     'metre/sec',
     $                     'east_u north_u zz',.true.)
      call def_var_pnetcdf(ncid,'ub',3,vdims,ub_varid,
     $                     'x-velocity at time -dt',
     $                     'metre/sec',
     $                     'east_u north_u zz',.true.)
      call def_var_pnetcdf(ncid,'v',3,vdims,v_varid,
     $                     'y-velocity',
     $                     'metre/sec',
     $                     'east_v north_v zz',.true.)
      call def_var_pnetcdf(ncid,'vb',3,vdims,vb_varid,
     $                     'y-velocity at time -dt',
     $                     'metre/sec',
     $                     'east_v north_v zz',.true.)
      call def_var_pnetcdf(ncid,'w',3,vdims,w_varid,
     $                     'sigma-velocity',
     $                     'metre/sec',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'t',3,vdims,t_varid,
     $                     'potential temperature',
     $                     'K',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'tb',3,vdims,tb_varid,
     $                     'potential temperature at time -dt',
     $                     'K',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'s',3,vdims,s_varid,
     $                     'salinity x rho / rhoref',
     $                     'PSS',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'sb',3,vdims,sb_varid,
     $                     'salinity x rho / rhoref at time -dt',
     $                     'PSS',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'rho',3,vdims,rho_varid,
     $                     '(density-1000)/rhoref',
     $                     'dimensionless',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'km',3,vdims,km_varid,
     $                     'vertical kinematic viscosity',
     $                     'metre^2/sec',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'kh',3,vdims,kh_varid,
     $                     'vertical diffusivity',
     $                     'metre^2/sec',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'kq',3,vdims,kq_varid,
     $                     'kq',
     $                     'metre^2/sec',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'l',3,vdims,l_varid,
     $                     'turbulence length scale',
     $                     '-',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2',3,vdims,q2_varid,
     $                     'twice the turbulent kinetic energy',
     $                     'metre^2/sec^2',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2b',3,vdims,q2b_varid,
     $                     'twice the turbulent kinetic energy at -dt',
     $                     'metre^2/sec^2',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'aam',3,vdims,aam_varid,
     $                     'horizontal kinematic viscosity',
     $                     'metre^2/sec',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2l',3,vdims,q2l_varid,
     $                     'q2 x l',
     $                     'metre^3/sec^2',
     $                     'east_e north_e zz',.true.)
      call def_var_pnetcdf(ncid,'q2lb',3,vdims,q2lb_varid,
     $                     'q2 x l at time -dt',
     $                     'metre^3/sec^2',
     $                     'east_e north_e zz',.true.)

! end definitions
      status=nfmpi_enddef(ncid)
      call handle_error_pnetcdf('nf_enddef',status,nf_noerr)

! write data
      start(1)=1
      edge(1)=1
      status=nfmpi_put_vara_double_all(ncid,time_varid,start,edge,time)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_put_vara_double_all(ncid,wubot_varid,
     $                                                 start,edge,wubot)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,wvbot_varid,
     $                                                 start,edge,wvbot)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,aam2d_varid,
     $                                                 start,edge,aam2d)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,ua_varid,start,edge,ua)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,uab_varid,start,edge,uab)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,va_varid,start,edge,va)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,vab_varid,start,edge,vab)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,el_varid,start,edge,el)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,elb_varid,start,edge,elb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,et_varid,start,edge,et)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,etb_varid,start,edge,etb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,egb_varid,start,edge,egb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,utb_varid,start,edge,utb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,vtb_varid,start,edge,vtb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,adx2d_varid,
     $                                                 start,edge,adx2d)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,ady2d_varid,
     $                                                 start,edge,ady2d)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,advua_varid,
     $                                                 start,edge,advua)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,advva_varid,
     $                                                 start,edge,advva)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      status=nfmpi_put_vara_double_all(ncid,u_varid,start,edge,u)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,ub_varid,start,edge,ub)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,v_varid,start,edge,v)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,vb_varid,start,edge,vb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,w_varid,start,edge,w)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,t_varid,start,edge,t)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,tb_varid,start,edge,tb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,s_varid,start,edge,s)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,sb_varid,start,edge,sb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,rho_varid,start,edge,rho)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,km_varid,start,edge,km)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,kh_varid,start,edge,kh)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,kq_varid,start,edge,kq)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,l_varid,start,edge,l)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,q2_varid,start,edge,q2)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,q2b_varid,start,edge,q2b)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,aam_varid,start,edge,aam)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,q2l_varid,start,edge,q2l)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)
      status=nfmpi_put_vara_double_all(ncid,q2lb_varid,start,edge,q2lb)
      call handle_error_pnetcdf('nf_put_vara_double',status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_grid_pnetcdf
! read grid data
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character*120 netcdf_grid_file
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid,dum_varid,dvm_varid
      integer ncid,status
      integer vdims(2)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf file
      write(netcdf_grid_file,'(a,''in/'',a,''.grid.nc'')')
     $                                   trim(wrk_pth),trim(netcdf_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $  trim(netcdf_grid_file)
      status=nfmpi_open(pom_comm,netcdf_grid_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_grid_file,
     $                          status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'zz',zz_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: zz'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'dx',dx_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: dx'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'dy',dy_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: dy'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'east_u',east_u_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_u'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'east_v',east_v_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_v'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'east_e',east_e_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_e'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'east_c',east_c_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: east_c'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'north_u',north_u_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_u'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'north_v',north_v_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_v'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'north_e',north_e_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_e'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'north_c',north_c_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: north_c'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'rot',rot_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: rot'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'h',h_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: h'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'fsm',fsm_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: fsm'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'dum',dum_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: dum'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'dvm',dvm_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: dvm'
     $                          ,status,nf_noerr)
! get data
      start(1)=1
      edge(1)=kb
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,z)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,zz_varid,start,edge,zz)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_double_all(ncid,dx_varid,start,edge,dx)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,dy_varid,start,edge,dy)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,east_u_varid,start,edge,
     $                               east_u)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,east_v_varid,start,edge,
     $                               east_v)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,east_e_varid,start,edge,
     $                               east_e)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,east_c_varid,start,edge,
     $                               east_c)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,north_u_varid,start,edge,
     $                               north_u)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,north_v_varid,start,edge,
     $                               north_v)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,north_e_varid,start,edge,
     $                               north_e)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,north_c_varid,start,edge,
     $                               north_c)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,rot_varid,start,edge,rot)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,h_varid,start,edge,h)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,fsm_varid,start,edge,fsm)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,dum_varid,start,edge,dum)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,dvm_varid,start,edge,dvm)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

! close file:
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close: grid',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_restart_pnetcdf
! read data for a seamless restart
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character*120 netcdf_in_file
      integer ncid,status
      integer time_varid,wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)
      integer i,j

! open netcdf restart file
      write(netcdf_in_file,'(a,''in/'',a)')
     $                                 trim(wrk_pth),trim(read_rst_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $ trim(netcdf_in_file)
      status=nfmpi_open(pom_comm,netcdf_in_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//trim(read_rst_file),
     $                          status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'time',time_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: time'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'wubot',wubot_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: wubot'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'wvbot',wvbot_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: wvbot'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'aam2d',aam2d_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: aam2d'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'ua',ua_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ua'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'uab',uab_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: uab'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'va',va_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: va'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'vab',vab_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vab'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'el',el_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: el'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'elb',elb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: elb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'et',et_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: et'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'etb',etb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: etb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'egb',egb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: egb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'utb',utb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: utb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'vtb',vtb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vtb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'u',u_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: u'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'ub',ub_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ub'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'v',v_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: v'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'vb',vb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'w',w_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: w'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'t',t_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: t'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'tb',tb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'s',s_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: s'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'sb',sb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: sb'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'rho',rho_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: rho'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'adx2d',adx2d_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: adx2d'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'ady2d',ady2d_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ady2d'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'advua',advua_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: advua'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'advva',advva_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: advva'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'km',km_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: km'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'kh',kh_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: kh'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'kq',kq_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: kq'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'l',l_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: l'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'q2',q2_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: q2'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'q2b',q2b_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: q2b'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'aam',aam_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: aam'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'q2l',q2l_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: q2l'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'q2lb',q2lb_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: q2lb'
     $                          ,status,nf_noerr)

! get data
      start(1)=1
      edge(1)=1
      status=nfmpi_get_vara_double_all(ncid,time_varid,start,edge,time0)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_double_all(ncid,wubot_varid,
     $                                                 start,edge,wubot)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,wvbot_varid,
     $                                                 start,edge,wvbot)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,aam2d_varid,
     $                                                 start,edge,aam2d)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,ua_varid,start,edge,ua)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,uab_varid,start,edge,uab)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,va_varid,start,edge,va)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,vab_varid,start,edge,vab)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,el_varid,start,edge,el)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,elb_varid,start,edge,elb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,et_varid,start,edge,et)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,etb_varid,start,edge,etb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,egb_varid,start,edge,egb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,utb_varid,start,edge,utb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,vtb_varid,start,edge,vtb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,adx2d_varid,
     $                                                 start,edge,adx2d)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,ady2d_varid,
     $                                                 start,edge,ady2d)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,advua_varid,
     $                                                 start,edge,advua)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,advva_varid,
     $                                                 start,edge,advva)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=kb
      status=nfmpi_get_vara_double_all(ncid,u_varid,start,edge,u)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,ub_varid,start,edge,ub)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,v_varid,start,edge,v)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,vb_varid,start,edge,vb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,w_varid,start,edge,w)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,t_varid,start,edge,t)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,tb_varid,start,edge,tb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,s_varid,start,edge,s)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,sb_varid,start,edge,sb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,rho_varid,start,edge,rho)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,km_varid,start,edge,km)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,kh_varid,start,edge,kh)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,kq_varid,start,edge,kq)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,l_varid,start,edge,l)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,q2_varid,start,edge,q2)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,q2b_varid,start,edge,q2b)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,aam_varid,start,edge,aam)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,q2l_varid,start,edge,q2l)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,q2lb_varid,start,edge,q2lb)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

! initialize elevation
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do

! initialize time
      time=time0

      return
      end

!_______________________________________________________________________
      subroutine read_initial_ts_pnetcdf(k,zlev,temp,salt)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer k
      double precision zlev(k),temp(im,jm,k),salt(im,jm,k)
      character*120 netcdf_ic_file
      integer ncid,status
      integer z_varid,t_varid,s_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)

! open netcdf ic file
      write(netcdf_ic_file,'(a,''in/'',a,''.ts_initial.nc'')')
     $                                   trim(wrk_pth),trim(netcdf_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                              trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_ic_file,
     $                                                  status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'t',t_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: t',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'s',s_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: s',status,nf_noerr)

! get data
      start(1)=1
      edge(1)=k
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,zlev)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=k
      status=nfmpi_get_vara_double_all(ncid,t_varid,start,edge,temp)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,s_varid,start,edge,salt)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_clim_ts_pnetcdf(k,zlev,temp,salt)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer k
      double precision zlev(k),temp(im,jm,k),salt(im,jm,k)
      character*120 netcdf_ic_file
      integer ncid,status
      integer z_varid,t_varid,s_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)

! open netcdf file
      write(netcdf_ic_file,'(a,''in/'',a,''.ts_clim.nc'')')
     $                                   trim(wrk_pth),trim(netcdf_file)
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                              trim(netcdf_ic_file)
      status=nfmpi_open(pom_comm,netcdf_ic_file,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_ic_file,
     $                                                  status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'t',t_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: t',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'s',s_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: s',status,nf_noerr)

! get data
      start(1)=1
      edge(1)=k
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,zlev)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=k
      status=nfmpi_get_vara_double_all(ncid,t_varid,start,edge,temp)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,s_varid,start,edge,salt)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_wind_pnetcdf(n,wu,wv)
! read wind stress
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer n
      double precision wu(im,jm),wv(im,jm)
      integer i,j
      character*120 netcdf_file2
      integer ncid,status
      integer wu_varid,wv_varid
      integer vdims(2)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf file
      write(netcdf_file2,'(a,''in/'',a,''.wind.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                                trim(netcdf_file2)
      status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_file2,
     $                                                 status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'wu',wu_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: wu',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'wv',wv_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: wv',status,nf_noerr)

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_double_all(ncid,wu_varid,start,edge,wu)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,wv_varid,start,edge,wv)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_restore_tau_interior_pnetcdf(n,k,zlev,tau)
! read restore temperature interior
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer n,k
      double precision zlev(k),tau(im,jm,k)
      character*120 netcdf_file2
      integer ncid,status
      integer z_varid,t_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)

! open netcdf file
      write(netcdf_file2,
     $          '(a,''in/'',a,''.restore_tau_interior.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                                trim(netcdf_file2)
      status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_file2,
     $                                                  status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'tau',t_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tau',status,nf_noerr)

! get data
      start(1)=1
      edge(1)=k
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,zlev)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=k
      status=nfmpi_get_vara_double_all(ncid,t_varid,start,edge,tau)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_heat_pnetcdf(n,shf,swr)
! read heat flux
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      double precision shf(im,jm),swr(im,jm)
      integer i,j,n
      character*120 netcdf_heat_file
      integer ncid,status
      integer shf_varid,swr_varid
      integer vdims(2)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf heat file
      write(netcdf_heat_file,'(a,''in/'',a,''.heat.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $                                            trim(netcdf_heat_file)
      status=nfmpi_open(pom_comm,netcdf_heat_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_heat_file,
     $                          status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'shf',shf_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: shf'
     $                          ,status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'swr',swr_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: swr'
     $                          ,status,nf_noerr)

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_double_all(ncid,shf_varid,start,edge,shf)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,swr_varid,start,edge,swr)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_water_pnetcdf(n,wat)
! read water flux
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      double precision wat(im,jm)
      integer i,j,n
      character*120 netcdf_water_file
      integer ncid,status
      integer wat_varid
      integer vdims(2)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf heat file
      write(netcdf_water_file,'(a,''in/'',a,''.water.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $                                           trim(netcdf_water_file)
      status=nfmpi_open(pom_comm,netcdf_water_file,nf_nowrite,
     $                  mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_water_file,
     $                          status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'w',wat_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: water'
     $                          ,status,nf_noerr)

! get data
      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      status=nfmpi_get_vara_double_all(ncid,wat_varid,start,edge,wat)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                          status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_restore_t_interior_pnetcdf(n,k,zlev,temp)
! read restore temperature interior
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer n,k
      double precision zlev(k),temp(im,jm,k)
      character*120 netcdf_file2
      integer ncid,status
      integer z_varid,t_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)

! open netcdf ic file
      write(netcdf_file2,
     $            '(a,''in/'',a,''.restore_t_interior.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                                trim(netcdf_file2)
      status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_file2,
     $                                                  status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'t',t_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: t',status,nf_noerr)

! get data
      start(1)=1
      edge(1)=k
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,zlev)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=k
      status=nfmpi_get_vara_double_all(ncid,t_varid,start,edge,temp)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_restore_s_interior_pnetcdf(n,k,zlev,salt)
! read restore salinity interior
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      integer n,k
      double precision zlev(k),salt(im,jm,k)
      character*120 netcdf_file2
      integer ncid,status
      integer z_varid,s_varid
      integer vdims(3)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(3),edge(3)

! open netcdf ic file
      write(netcdf_file2,
     $            '(a,''in/'',a,''.restore_s_interior.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                                trim(netcdf_file2)
      status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_file2,
     $                                                  status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'s',s_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: s',status,nf_noerr)

! get data
      start(1)=1
      edge(1)=k
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,zlev)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=j_global(1)
      start(3)=1
      edge(1)=im
      edge(2)=jm
      edge(3)=k
      status=nfmpi_get_vara_double_all(ncid,s_varid,start,edge,salt)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_boundary_conditions_pnetcdf(n,k,zlev,t_w,s_w,u_w,
     $                             t_e,s_e,u_e,t_n,s_n,v_n,t_s,s_s,v_s)
! read boundary conditions
      implicit none
# include "mpif.h"
# include "pnetcdf.inc"
      include 'pom.h'
      character*120 netcdf_file2
      integer n,k
      double precision zlev(k)
      double precision t_w(jm,k),s_w(jm,k),u_w(jm),t_e(jm,k),s_e(jm,k),
     $   u_e(jm),t_n(im,k),s_n(im,k),v_n(im),t_s(im,k),s_s(im,k),v_s(im)
      integer ncid,status
      integer z_varid
      integer tw_varid,sw_varid,uw_varid
      integer te_varid,se_varid,ue_varid
      integer tn_varid,sn_varid,vn_varid
      integer ts_varid,ss_varid,vs_varid
      integer vdims(2)
      integer(mpi_offset_kind) length
      integer(mpi_offset_kind) start(2),edge(2)

! open netcdf ic file
      write(netcdf_file2,'(a,''in/'',a,''.bc.'',i4.4,''.nc'')')
     $                                 trim(wrk_pth),trim(netcdf_file),n
      if(my_task.eq.0) write(*,'(/''reading file '',a)')
     $                                                trim(netcdf_file2)
      status=nfmpi_open(pom_comm,netcdf_file2,nf_nowrite,
     $                                               mpi_info_null,ncid)
      call handle_error_pnetcdf('nfmpi_open: '//netcdf_file2,
     $                                                  status,nf_noerr)

! get variables
      status=nfmpi_inq_varid(ncid,'z',z_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: z',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'tw',tw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tw',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'sw',sw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: sw',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'uw',uw_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: uw',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'te',te_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: te',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'se',se_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: se',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'ue',ue_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ue',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'tn',tn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: tn',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'sn',sn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: sn',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'vn',vn_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vn',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'ts',ts_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ts',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'ss',ss_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: ss',status,nf_noerr)
      status=nfmpi_inq_varid(ncid,'vs',vs_varid)
      call handle_error_pnetcdf('nfmpi_inq_varid: vs',status,nf_noerr)

! get data
      start(1)=1
      edge(1)=k
      status=nfmpi_get_vara_double_all(ncid,z_varid,start,edge,zlev)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all',
     $                                                  status,nf_noerr)

      start(1)=j_global(1)
      edge(1)=jm
      status=nfmpi_get_vara_double_all(ncid,uw_varid,start,edge,u_w)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: u_w',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,ue_varid,start,edge,u_e)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: u_e',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      edge(1)=im
      status=nfmpi_get_vara_double_all(ncid,vn_varid,start,edge,v_n)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: v_n',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,vs_varid,start,edge,v_s)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: v_s',
     $                                                  status,nf_noerr)

      start(1)=j_global(1)
      start(2)=1
      edge(1)=jm
      edge(2)=k
      status=nfmpi_get_vara_double_all(ncid,tw_varid,start,edge,t_w)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: t_w',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,sw_varid,start,edge,s_w)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: s_w',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,te_varid,start,edge,t_e)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: t_e',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,se_varid,start,edge,s_e)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: s_e',
     $                                                  status,nf_noerr)

      start(1)=i_global(1)
      start(2)=1
      edge(1)=im
      edge(2)=k
      status=nfmpi_get_vara_double_all(ncid,tn_varid,start,edge,t_n)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: t_n',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,sn_varid,start,edge,s_n)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: s_n',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,ts_varid,start,edge,t_s)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: t_s',
     $                                                  status,nf_noerr)
      status=nfmpi_get_vara_double_all(ncid,ss_varid,start,edge,s_s)
      call handle_error_pnetcdf('nfmpi_get_vara_double_all: s_s',
     $                                                  status,nf_noerr)

! close file
      status=nfmpi_close(ncid)
      call handle_error_pnetcdf('nf_close',status,nf_noerr)

      return
      end

