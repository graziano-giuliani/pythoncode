
#include <stdio.h>
#include <netcdf.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

int main(int argc,char *argv[])
{
  float *sfgrd, *emrf;
  long nsf;
  char *fname, *p, *e, foutname[PATH_MAX], line[256];
  int i, ncres, fres, ncid, idim, ivar[2];
  FILE *fp;
  char units[2][16] = { "cm-1", "1" };

  fprintf(stdout, "Convert Emissivity ASCII file to NetCDF file\n");

  if (argc < 2)
  {
    fprintf(stderr,"Not enough arguments. Need name of input file.\n");
    fprintf(stderr,"Usage :\n");
    fprintf(stderr,"%s filename\n",argv[0]);
    return -1;
  }
  fname = strdup(argv[1]);
  fp = fopen(fname,"r");
  if ( fp == NULL )
  {
    fprintf(stderr,"Cannot open file %s\n",fname);
    return -2;
  }
  p = strrchr(fname, '/');
  if ( p == NULL ) p = fname;
  else p = p+1;
  e = strrchr(fname, '.');
  if ( e != NULL ) *e = 0;
  sprintf(foutname,"%s.nc", p);
  ncres = nc_create(foutname,NC_CLOBBER,&ncid);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Cannot open output file %s\n",foutname);
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -3;
  }
  fres = fscanf(fp, "%ld\n", &nsf);
  if ( fres != 1 )
  {
    fprintf(stderr,"Error reading number of hinge points from input file\n");
    return -4;
  }
  if ( nsf < 0 || nsf > 1000000000 )
  {
    fprintf(stderr,"Error reading nsf from input file: %ld\n",nsf);
    return -4;
  }
  sfgrd = (float *) malloc(nsf*sizeof(float));
  emrf = (float *) malloc(nsf*sizeof(float));
  if ( sfgrd == NULL || emrf == NULL )
  {
    fprintf(stderr,"Error malloc data space: %ld",nsf*2*sizeof(float));
    return -5;
  }
  ncres = nc_def_dim(ncid,"nspectral",nsf,&idim);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error creating dim nspectral\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_def_var(ncid,"SfGrd",NC_FLOAT,1,&idim,ivar);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error creating SfGrd var\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_def_var(ncid,"EmRf",NC_FLOAT,1,&idim,ivar+1);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error creating EmRf var\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_put_att_text(ncid,ivar[0],"units",strlen(units[0])+1,units[0]);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error adding units\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_put_att_text(ncid,ivar[1],"units",strlen(units[1])+1,units[1]);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error adding units\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_enddef(ncid);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error finalizing outout file\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -7;
  }
  for ( i = 0; i < nsf; i ++ )
  {
    fres = fscanf(fp,"%f %f\n",sfgrd+i,emrf+i);
    if ( fres != 2 )
    {
      fprintf(stderr,"Error reading input data at record %d\n",i);
      return -8;
    }
  }
  ncres = nc_put_var_float(ncid,ivar[0],sfgrd);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error write var SfGrd\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -9;
  }
  ncres = nc_put_var_float(ncid,ivar[1],emrf);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error write var EmRf\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -9;
  }
  (void) fclose(fp);
  ncres = nc_close(ncid);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error closing output file!\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -10;
  }
  return 0;
}
