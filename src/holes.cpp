// Authors:	Sven Wiese
//-----------------------------------------------------

#include <cstdio>
#include <cstdlib>

#include "holes.hpp"
#include "utils.hpp"

int
readholes (const char *filename,
           int        *nvarswithholes_p,
           double     **lbvarswithholes_p,
           double     **ubvarswithholes_p,
           int        **nholesvarswithholes_p,
           char       ***namevarswithholes_p,
           double     ***lbholes_p,
           double     ***ubholes_p)
{
    int status=0;
    
    FILE *fin  = NULL;
    
    int nvarswithholes=0;
    double *lbvarswithholes=NULL;
    double *ubvarswithholes=NULL;
    int *nholesvarswithholes=NULL;
    char **namevarswithholes=NULL;
    double **lbholes=NULL;
    double **ubholes=NULL;
    
    int i, j;
    
    fin = fopen (filename, "r");
    if ( fin == NULL ) {
        printf ("ERROR1 %s\n", filename);
        status = ERR_BADARGUMENT;
        goto TERMINATE;
    }
    
    /*read number of variables with holes*/
    if ( fscanf (fin, "%d", &nvarswithholes) != 1 ) {
        status = ERR_BADFILEFORMAT;
        goto TERMINATE;
    }
    
    /*Allocate memory*/
    lbvarswithholes=(double*)malloc(nvarswithholes*sizeof(double));
    ubvarswithholes=(double*)malloc(nvarswithholes*sizeof(double));
    nholesvarswithholes=(int*)malloc(nvarswithholes*sizeof(int));
    namevarswithholes=(char**)malloc(nvarswithholes*sizeof(char*));
    lbholes=(double**)malloc(nvarswithholes*sizeof(double*));
    ubholes=(double**)malloc(nvarswithholes*sizeof(double*));
    
    if ( lbvarswithholes == NULL ||
         ubvarswithholes == NULL ||
         nholesvarswithholes == NULL ||
         namevarswithholes == NULL ||
         lbholes == NULL ||
         ubholes == NULL ) {
        status = ERR_NOMEMORY;
        goto TERMINATE;
    }
    
    for(i=0;i<nvarswithholes;i++){
        namevarswithholes[i]=(char*)malloc(BUFFERSIZE*sizeof(char));
        if ( namevarswithholes[i] == NULL) {
            status = ERR_NOMEMORY;
            goto TERMINATE;
        }
    }
    
    /*read hole information row by row*/
    for(i=0;i<nvarswithholes;i++){
        int coef1,coef2;
        if(fscanf(fin, "%s",namevarswithholes[i]) != 1){
            status = ERR_BADFILEFORMAT;
            goto TERMINATE;
        }
        if(fscanf(fin, "%d",&coef1) != 1){
            status = ERR_BADFILEFORMAT;
            goto TERMINATE;
        }
        if(fscanf(fin, "%d",&coef2) != 1){
            status = ERR_BADFILEFORMAT;
            goto TERMINATE;
        }
        if(fscanf(fin, "%d",&nholesvarswithholes[i]) != 1){
            status = ERR_BADFILEFORMAT;
            goto TERMINATE;
        }
        lbvarswithholes[i]=(double)coef1;
        ubvarswithholes[i]=(double)coef2;
        lbholes[i]=(double*)malloc(nholesvarswithholes[i]*sizeof(double));
        ubholes[i]=(double*)malloc(nholesvarswithholes[i]*sizeof(double));
        if ( lbholes[i] == NULL ||
             ubholes[i] == NULL) {
            status = ERR_NOMEMORY;
            goto TERMINATE;
        }
        for(j=0;j<nholesvarswithholes[i];j++){
            if(fscanf(fin, "%d",&coef1) != 1){
                status = ERR_BADFILEFORMAT;
                goto TERMINATE;
            }
            if(fscanf(fin, "%d",&coef2) != 1){
                status = ERR_BADFILEFORMAT;
                goto TERMINATE;
            }
            lbholes[i][j]=(double)coef1;
            ubholes[i][j]=(double)coef2;
        }
    }
    
#ifdef OUTPL
	printf("----- holes read -----\n");
	for (int i=0; i<nvarswithholes; i++){
		printf("%s %d %d %d",namevarswithholes[i],nholesvarswithholes[i],(int)lbvarswithholes[i],(int)ubvarswithholes[i]);
		for (int j=0; j<nholesvarswithholes[i]; j++){
			printf(" %d %d",(int)lbholes[i][j],(int)ubholes[i][j]);
		}
		printf("\n");
	}

#endif

    /* Install the hole information */
    *nvarswithholes_p=nvarswithholes;
    *lbvarswithholes_p=lbvarswithholes;
    *ubvarswithholes_p=ubvarswithholes;
    *namevarswithholes_p=namevarswithholes;
    *nholesvarswithholes_p=nholesvarswithholes;
    *lbholes_p=lbholes;
    *ubholes_p=ubholes;
    
    lbvarswithholes=NULL;
    ubvarswithholes=NULL;
    namevarswithholes=NULL;
    nholesvarswithholes=NULL;
    lbholes=NULL;
    ubholes=NULL;
    
TERMINATE:
    
    if ( fin != NULL )
        fclose (fin);

   FREEN (&lbvarswithholes);
   FREEN (&ubvarswithholes);
   FREEN_mat (&lbholes,nvarswithholes);
   FREEN_mat (&ubholes,nvarswithholes);
   FREEN (&nholesvarswithholes);
   FREEN_mat (&namevarswithholes,nvarswithholes);
    
    return status;
}/*END readholes*/
