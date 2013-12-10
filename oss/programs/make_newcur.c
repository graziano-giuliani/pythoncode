
#include <stdio.h>
#include <netcdf.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

int main(int argc,char *argv[])
{
  float *wvl, *irr;
  long nspi;
  char *fname, *p, *e, foutname[PATH_MAX], line[256];
  int i, ncres, fres, ncid, idim, ivar[2];
  FILE *fp;
  char units[2][16] = { "cm-1", "mW m-2 / cm-1" };
  const float solcon = 1368.0;

  fprintf(stdout, "Convert ASCII Solar Irradiance to NetCDF file\n");

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
  ncres = nc_put_att_float(ncid,NC_GLOBAL,"solar_constant",NC_FLOAT,1,&solcon);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error adding global attribute\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -3;
  }
  fres = fscanf(fp, "%ld\n", &nspi);
  if ( fres != 1 )
  {
    fprintf(stderr,"Error reading nspi from input file\n");
    return -4;
  }
  if ( nspi < 0 || nspi > 1000000000 )
  {
    fprintf(stderr,"Error reading nspi from input file: %ld\n",nspi);
    return -4;
  }
  wvl = (float *) malloc(nspi*sizeof(float));
  irr = (float *) malloc(nspi*sizeof(float));
  if ( wvl == NULL || irr == NULL )
  {
    fprintf(stderr,"Error malloc data space: %ld",nspi*2*sizeof(float));
    return -5;
  }
  ncres = nc_def_dim(ncid,"nspectral",nspi,&idim);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error creating dim nspectral\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_def_var(ncid,"FREQ",NC_FLOAT,1,&idim,ivar);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error creating FREQ var\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -6;
  }
  ncres = nc_def_var(ncid,"IRRADIANCE",NC_FLOAT,1,&idim,ivar+1);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error creating IRRADIANCE var\n");
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
  /* Skip 2 lines header */
  (void) fgets(line,256,fp);
  (void) fgets(line,256,fp);
  for ( i = 0; i < nspi; i ++ )
  {
    fres = fscanf(fp,"%f %f\n",wvl+i,irr+i);
    if ( fres != 2 )
    {
      fprintf(stderr,"Error reading input data at record %d\n",i);
      return -8;
    }
    /* Convert units of measure */
    irr[i] *= 1.0e+7;
  }
  ncres = nc_put_var_float(ncid,ivar[0],wvl);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error write var FREQ\n");
    fprintf(stderr,"%s\n",nc_strerror(ncres));
    return -9;
  }
  ncres = nc_put_var_float(ncid,ivar[1],irr);
  if ( ncres != NC_NOERR )
  {
    fprintf(stderr,"Error write var IRRADIANCE\n");
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
