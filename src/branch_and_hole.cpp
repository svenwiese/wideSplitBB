// Authors:	Sven Wiese
//-----------------------------------------------------

#include <cassert>
#include <iostream>
#include <string>
#include <map>
#include <cfloat>
#include <cstring>
#include <algorithm>

#include "cplex.h"

#include "holes.hpp"
#include "utils.hpp"

struct Options {

  typedef std::map<std::string, std::pair<bool, std::string> > b_opt_t;
  typedef std::map<std::string, std::pair<int, std::string> > i_opt_t;
  typedef std::map<std::string, std::pair<double, std::string > > d_opt_t;
  typedef std::map<std::string, std::pair<std::string, std::string> > s_opt_t;
   
  void usage(){
    std::cerr <<"  usage:" <<std::endl;
    std::cerr <<"    "<<"branch_and_hole " <<"[options]   instance{.mps , .lp} " <<std::endl;
    std::cerr <<"   OPTIONS:"<<std::endl;
    for(b_opt_t::iterator k = b_opt.begin() ; k != b_opt.end() ; k++){
      std::cerr<<"    -"<<k->first<<" [Boolean option]: "<<k->second.second<<" (default "<<k->second.first<<")"<<std::endl;
    }
    for(i_opt_t::iterator k = i_opt.begin() ; k != i_opt.end() ; k++){
      std::cerr<<"    -"<<k->first<<" [Integer option]: "<<k->second.second<<" (default "<<k->second.first<<")"<<std::endl;
    }
    for(d_opt_t::iterator k = d_opt.begin() ; k != d_opt.end() ; k++){
      std::cerr<<"    -"<<k->first<<" [Double option]: "<<k->second.second<<" (default "<<k->second.first<<")"<<std::endl;
    }
    for(s_opt_t::iterator k = s_opt.begin() ; k != s_opt.end() ; k++){
      std::cerr<<"    -"<<k->first<<" [String option]: "<<k->second.second<<" (default "<<k->second.first<<")"<<std::endl;
    }
  }
  enum KeyTypes{
    BoolKey,
    IntKey,
    DoubleKey,
    StringKey};

  b_opt_t b_opt;
  i_opt_t i_opt;
  d_opt_t d_opt;
  s_opt_t s_opt;

  std::map<std::string, KeyTypes> keys;

  std::map<std::string, bool> stringset;

  bool isOpt(const char * a){
   return a[0]=='-';}

  Options(){
    
    keys["inc"] = BoolKey;
    b_opt["inc"] = std::make_pair(false, "Incumbent callback switch");
    keys["br"] = BoolKey;
    b_opt["br"] = std::make_pair(false, "Branch callback switch");
    keys["c"] = BoolKey;
    b_opt["c"] = std::make_pair(false, "Turn off cplex cuts");
    keys ["filter"] = BoolKey;
    b_opt["filter"] = std::make_pair(false, "Filter cuts (CPX_PURGE otherwise)");
    
    keys["n"] = IntKey;
    i_opt["n"] = std::make_pair(0, "Number of rounds of wide split cuts at root node");
    keys["log"] = IntKey;
    i_opt["log"] = std::make_pair(1, "log level");
    
    keys["bab_time"] = DoubleKey;
    d_opt["bab_time"] = std::make_pair(DBL_MAX, "time limit in branch-and-bound");

    keys["hfile"] = StringKey;
    stringset["hfile"] = false;
    s_opt["hfile"] = std::make_pair("instance_holes.txt", "holefile");

    keys["logfile"] = StringKey;
    stringset["logfile"] = false;
    s_opt["logfile"] = std::make_pair("instance.bblog", "logfile");
  }

  int process(int argc, const char ** argv){
    int num = 1;
    while(argv[num] != NULL){
      if(!isOpt(argv[num])){
         return num;
      }
      const char * k = &argv[num][1];
      if(keys.find(k) == keys.end()){
         printf("Unknowns option %s\n", argv[num]);
         usage();
         exit(0);
      }
      switch(keys[k]){
        case BoolKey:
         std::cout<<"Boolean key "<<k<<std::endl;
         b_opt[k].first = true;
         break;
       case IntKey:
         num++;
         assert(argv[num] != NULL);
         i_opt[k].first = atoi(argv[num]);
         break;
       case DoubleKey:
         num++;
         assert(argv[num] != NULL);
         d_opt[k].first = atof(argv[num]);
         break;
       case StringKey:
         num++;
         assert(argv[num] != NULL);
         s_opt[k].first = argv[num];
         stringset[k] = true;
         break;
       default:
         printf("Unknowns option %s\n", argv[num]);
         usage();
         exit(0);
       }
      num++;
    }
    return num;
  }
};

static int inccalled = 0;
static int increjected = 0;
static int brstr = 0;
static int brchange = 0;
static int npass = 0;
static int cuts_tot = 0;

static bool inc_rejected = false;
static int varindex_found = 0;
static int cpxvarindex_found = 0;
static int holeindex_found = 0;

struct incumbentdata {

   incumbentdata() : 	nvarswithholes(0),
			lbvarswithholes(NULL),
			ubvarswithholes(NULL),
			nholesvarswithholes(NULL),
			namevarswithholes(NULL),
			lbholes(NULL),
			ubholes(NULL),
			ncols(0),
			loglevel(0){}

   int nvarswithholes;
   double *lbvarswithholes;
   double *ubvarswithholes;
   int *nholesvarswithholes;
   char **namevarswithholes;
   double **lbholes;
   double **ubholes;
   std::map<std::string, int> name_to_index;
   int ncols;
   int loglevel;
};

struct branchdata {

	branchdata() : 	nvarswithholes(0),
			lbvarswithholes(NULL),
			ubvarswithholes(NULL),
			nholesvarswithholes(NULL),
			namevarswithholes(NULL),
			lbholes(NULL),
			ubholes(NULL),
			ncols(0),
			loglevel(0),
			has_hole(NULL){}

   int nvarswithholes;
   double *lbvarswithholes;
   double *ubvarswithholes;
   int *nholesvarswithholes;
   char **namevarswithholes;
   double **lbholes;
   double **ubholes;
   std::map<std::string, int> name_to_index;
   std::map<int, std::string> index_to_name;
   std::map<int, int> index_to_index;
   int ncols;
   int loglevel;
   bool *has_hole;
};

struct cutdata {

cutdata() : 	nvarswithholes(0),
		lbvarswithholes(NULL),
		ubvarswithholes(NULL),
		nholesvarswithholes(NULL),
		namevarswithholes(NULL),
		lbholes(NULL),
		ubholes(NULL),
		ncols(0),
		maxpass(0),
		ctype(NULL),
		filter(0),
		fout(NULL),
		loglevel(0){}

   int nvarswithholes;
   double *lbvarswithholes;
   double *ubvarswithholes;
   int *nholesvarswithholes;
   char **namevarswithholes;
   double **lbholes;
   double **ubholes;
   std::map<std::string, int> name_to_index;
   int ncols;
   int maxpass;
   char *ctype;
   int filter;
   FILE *fout;
   int loglevel;
};

int CPXPUBLIC
 hole_incumbentcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           double objval,
           double *x,
           int *isfeas_p,
           int *useraction_p){

	incumbentdata *inc = (incumbentdata*) cbhandle;

	if (inc->loglevel>=3) printf("incumbentcallback called.\n");

	inccalled++;

	int i = 0;
	int j = 0;
	int index = 0;
	bool found = false;
	// check whether any x[j] is in a hole
	for(i=0; i<inc->nvarswithholes; i++){
		index = inc->name_to_index[inc->namevarswithholes[i]];
		for(j=0; j<inc->nholesvarswithholes[i]; j++){
			if (x[index] > inc->lbholes[i][j] - 1 + EPSVIOL && x[index] < inc->ubholes[i][j] + 1 -EPSVIOL){
				found = true;
				break;
			}
		}
		if (found) break;
	}

	if (found) {
		*isfeas_p = 0;
		increjected++;
		if (wherefrom == CPX_CALLBACK_MIP_INCUMBENT_NODESOLN){
			inc_rejected = true;
			varindex_found = i;
			cpxvarindex_found = index;
			holeindex_found = j;
		}
		if (inc->loglevel>=1) printf("integer solution rejected.\n");
	}

	return 0;
}

int CPXPUBLIC
 hole_branchcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           int type,
           int sos,
           int nodecnt,
           int bdcnt,
           const int *nodebeg,
           const int *indices,
           const char *lu,
           const double *bd,
           const double *nodeest,
           int *useraction_p){

	*useraction_p = CPX_CALLBACK_DEFAULT;

	branchdata *branch = (branchdata*) cbhandle;

	// get a local copy of inc_rejected and reset it
	bool rejected = inc_rejected;
	inc_rejected = false;

	if (branch->loglevel>=3) printf("branchcallback called, branchtype: %c. inc_rejected: %d.\n",type,rejected);

	int status = 0;

	int *mynodebeg = NULL;

	bool enforce = false;
	int cnt = 0;

	char *varlu = NULL;
	double *varbd = NULL;

	double *x = NULL;

	if (rejected){ // branch on the variable detected in the inccb
		varlu = (char*) malloc(2*sizeof(char));
		varbd = (double*) malloc(2*sizeof(double));

		varlu[0] = 'U';
		varlu[1] = 'L';
		varbd[0] = branch->lbholes[varindex_found][holeindex_found]-1;
		varbd[1] = branch->ubholes[varindex_found][holeindex_found]+1;
		double est = 0; //adjust this
		status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, &est);
		if (status) goto TERMINATE;
		for (int c=0; c<2; c++){
			int seqnr = 0; //adjust this
			status = CPXbranchcallbackbranchbds( env, cbdata, wherefrom, 1, &cpxvarindex_found, varlu+c, varbd+c, est, NULL, &seqnr);
			if (status) goto TERMINATE;
		}
		brchange++;
		if (branch->loglevel>=2){
					printf("branching on variable with hole, setting new bounds to %.2lf / %.2lf\n",
					varbd[0],varbd[1]);
		}
	} else { // see whether the cplex branching can be strengthened
		if (type != CPX_TYPE_VAR || nodecnt == 0) goto TERMINATE;

		/* copy the branching data */
		varlu = (char*) malloc(bdcnt*sizeof(char));
		varbd = (double*) malloc(bdcnt*sizeof(double));
		for (int k=0; k<bdcnt; k++){
			varbd[k] = bd[k];
			varlu[k] = lu[k];
		}

		mynodebeg = (int *) malloc ((nodecnt+1) * sizeof (int));
		if (mynodebeg == NULL) goto TERMINATE;
		for (int c=0; c<nodecnt; c++){
			mynodebeg[c] = nodebeg[c];
		}
		mynodebeg[nodecnt] = bdcnt;

		/* now check whether any of the new bounds falls in a hole
	 	 * and can therefore be enforced */
		for (int c=0; c<nodecnt; c++){
			for (int k=mynodebeg[c]; k<mynodebeg[c+1]; k++){
				if (lu[k] == 'B') continue;
				int varind = indices[k];
				if (branch->has_hole[varind] == TRUE) {
					int i = branch->index_to_index[varind];
					for(int j=0; j<branch->nholesvarswithholes[i]; j++){
						if ( bd[k] > branch->lbholes[i][j]-1 && bd[k] < branch->ubholes[i][j]+1){
							varbd[k] = (lu[k] == 'L') ? branch->ubholes[i][j]+1 : branch->lbholes[i][j]-1;
							varlu[k] = lu[k];
							enforce = true;
							brstr++;
							if (branch->loglevel>=2){
								printf("cplex branching on variable with hole: %s, suggested new %c to %.2lf,\n",
								branch->namevarswithholes[i],lu[k],bd[k]);
								printf("\tthat's in a hole. setting new bound to %.2lf\n",varbd[k]);
							}
							break;
						}
					}
				}
			}
		}

		if (enforce == true){
			cnt = 0;
			for (int c=0; c<nodecnt; c++){
				int seqnr = 0;
				status = CPXbranchcallbackbranchbds( env, cbdata, wherefrom, mynodebeg[c+1]-mynodebeg[c],
							 indices+cnt, varlu+cnt, varbd+cnt, nodeest[c], NULL, &seqnr);
				if (status) goto TERMINATE;
				cnt += mynodebeg[c+1]-mynodebeg[c];
			}
		} else {
			goto TERMINATE;
		}
	}

	*useraction_p = CPX_CALLBACK_SET;

TERMINATE:

	FREEN(&varlu);
	FREEN(&varbd);
	FREEN(&mynodebeg);
	FREEN(&x);

	return status;
}

int CPXPUBLIC
 hole_cutcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           int *useraction_p){

	int status = 0;
     
    /* insert your cut separation here */

TERMINATE:

	return status;
}

int CPXPUBLIC
 empty_incumbentcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           double objval,
           double *x,
           int *isfeas_p,
           int *useraction_p){

	return 0;
}

int CPXPUBLIC
 empty_branchcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           int type,
           int sos,
           int nodecnt,
           int bdcnt,
           const int *nodebeg,
           const int *indices,
           const char *lu,
           const double *bd,
           const double *nodeest,
           int *useraction_p){

	return 0;
}

int CPXPUBLIC
 empty_cutcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           int *useraction_p){

	return 0;
}
  
int main(int argc, const char *argv[])
{
  Options opt;
  int num = opt.process(argc, argv);

  if ( num != argc - 1 ) {
    std::cerr <<"Incorrect number of command line parameters." <<std::endl;
    opt.usage();
    exit(1);
  }
    
  if (opt.i_opt["n"].first > 0) {
      std::cerr <<"WARNING: Chose to separate user cuts, but no separation is currently implemented." <<std::endl;
  }
    
  int status = 0;
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;

  int norigrows = 0;
  char *ctype = NULL;
  double objval;
  double cutoff;
  std::string opt_stat;

  int filter = opt.b_opt["filter"].first ? CPX_USECUT_FILTER : CPX_USECUT_PURGE;

  std::string FileName = argv[num];
  std::string logFileName = (FileName);
  std::string hFileName (FileName);
  if (opt.stringset["hfile"]){
	hFileName = opt.s_opt["hfile"].first;
  } else {
	size_t pos = hFileName.find(".mps");
	if (pos != std::string::npos){
		hFileName.replace(pos,4,"_holes.txt");
	} else {
		pos = hFileName.find(".lp");
		hFileName.replace(pos,3,"_holes.txt");
	}
  }
  if (opt.stringset["logfile"]){
	logFileName = opt.s_opt["logfile"].first;
  } else {
	size_t pos = logFileName.find(".mps");
	if (pos != std::string::npos){
		logFileName.replace(pos,4,".bblog");
	} else {
		pos = logFileName.find(".lp");
		logFileName.replace(pos,3,".bblog");
	}
  }


  std::cout<<FileName<<", "<<hFileName<<", "<<logFileName<<std::endl;

  FILE *fout = NULL;

   int nvarswithholes=0;
   double *lbvarswithholes=NULL;
   double *ubvarswithholes=NULL;
   int *nholesvarswithholes=NULL;
   char **namevarswithholes=NULL;
   double **lbholes=NULL;
   double **ubholes=NULL;
   int ncols = 0;
   std::map<std::string, int> name_to_index;
   std::map<int, std::string> index_to_name;
   std::map<int, int> index_to_index;
   bool *has_hole = NULL;

  incumbentdata inc;
  branchdata branch;
  cutdata cut;

  env = CPXopenCPLEX(&status); 
  if (status) goto TERMINATE;

  fout = fopen (logFileName.c_str(),"a");
  if (fout==NULL){
	status = ERR_OPENFILE;
	goto TERMINATE;
  }

    // Read file describing problem
    lp = CPXcreateprob(env, &status, FileName.c_str());
    status = CPXreadcopyprob(env, lp, FileName.c_str(),NULL);
    if(status != 0){
      printf("Could not read file %s error %i\n", FileName.c_str(), status);
      exit(1);
    }

  CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
  // turn off nonlinear and dual reductions (should be implied by presence of incumbentcallback)
  CPXsetintparam(env, CPX_PARAM_PRELINEAR, CPX_OFF);
  CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
  // also turn off primal reductions because we use the function CPXgetcallbacknodelp
  CPXsetintparam(env, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
  // let callback ork on original problems
  CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);

  CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-08);
  CPXsetdblparam(env, CPX_PARAM_EPAGAP, 1e-08);
  CPXsetintparam(env, CPX_PARAM_THREADS, 1);
  CPXsetdblparam(env, CPX_PARAM_TILIM, opt.d_opt["bab_time"].first);

  if (opt.b_opt["c"].first){
	CPXsetintparam(env, CPXPARAM_MIP_Cuts_LiftProj, -1);
	CPXsetintparam(env, CPXPARAM_MIP_Cuts_Gomory, -1);
	CPXsetintparam(env, CPXPARAM_MIP_Limits_EachCutLimit, 0);
  }

  ncols = CPXgetnumcols(env, lp);
  norigrows = CPXgetnumrows(env, lp);

  // read hole information
  has_hole = (bool*)malloc(ncols*sizeof(bool));
  if ( has_hole == NULL ) {
      status = ERR_NOMEMORY;
      goto TERMINATE;
   }
    
  if (opt.b_opt["inc"].first || opt.b_opt["bra"].first || opt.i_opt["n"].first > 0) {
	status = readholes (hFileName.c_str(),&nvarswithholes,&lbvarswithholes,
                           &ubvarswithholes,&nholesvarswithholes,
                           &namevarswithholes,&lbholes,&ubholes);
  if ( status ) goto TERMINATE;
      
    // map variable names to column indices and vice versa
	char* name[1] = {NULL};
	char* namestore = NULL;
	namestore = (char*)malloc(BUFFERSIZE*sizeof(char));
	if ( namestore == NULL) {
            status = ERR_NOMEMORY;
            goto TERMINATE;
        }
	int surplus = 0;
	for (int j=0; j<ncols; j++) has_hole[j] = FALSE;
	for (int i=0; i<nvarswithholes; i++){
		for (int j=0; j<ncols; j++){
			status = CPXgetcolname (env, lp, name, namestore, BUFFERSIZE, &surplus, j, j);
			if ( status ) goto TERMINATE;
			if (!strcmp(namestore,namevarswithholes[i])){
				name_to_index[namestore] = j;
				index_to_name[j] = namestore;
				index_to_index[j] = i;
				has_hole[j] = TRUE;
				break;
			}
		}
	}
	FREEN(&namestore);
  }

  // get ctype information
   ctype = (char *) malloc (ncols * sizeof (char));
   if ( ctype == NULL ) {
      status = ERR_NOMEMORY;
      goto TERMINATE;
   }
   status = CPXgetctype (env, lp, ctype, 0, ncols-1);
   if ( status ) goto TERMINATE;

   inc.nvarswithholes=nvarswithholes;
   inc.lbvarswithholes=lbvarswithholes;
   inc.ubvarswithholes=ubvarswithholes;
   inc.nholesvarswithholes=nholesvarswithholes;
   inc.namevarswithholes=namevarswithholes;
   inc.lbholes=lbholes;
   inc.ubholes=ubholes;
   inc.name_to_index=name_to_index;
   inc.ncols=ncols;
   inc.loglevel=opt.i_opt["log"].first;

   branch.nvarswithholes=nvarswithholes;
   branch.lbvarswithholes=lbvarswithholes;
   branch.ubvarswithholes=ubvarswithholes;
   branch.nholesvarswithholes=nholesvarswithholes;
   branch.namevarswithholes=namevarswithholes;
   branch.lbholes=lbholes;
   branch.ubholes=ubholes;
   branch.name_to_index=name_to_index;
   branch.index_to_name=index_to_name;
   branch.index_to_index=index_to_index;
   branch.ncols=ncols;
   branch.loglevel=opt.i_opt["log"].first;
   branch.has_hole=has_hole;

   cut.nvarswithholes=nvarswithholes;
   cut.lbvarswithholes=lbvarswithholes;
   cut.ubvarswithholes=ubvarswithholes;
   cut.nholesvarswithholes=nholesvarswithholes;
   cut.namevarswithholes=namevarswithholes;
   cut.lbholes=lbholes;
   cut.ubholes=ubholes;
   cut.name_to_index=name_to_index;
   cut.ncols=ncols;
   cut.maxpass=opt.i_opt["n"].first;
   cut.ctype=ctype;
   cut.filter=filter;
   cut.fout=fout;
   cut.loglevel=opt.i_opt["log"].first;

  if(opt.b_opt["inc"].first)
    CPXsetincumbentcallbackfunc(env, hole_incumbentcallback, &inc);
  else
    CPXsetincumbentcallbackfunc(env, empty_incumbentcallback, &inc);
  if(opt.b_opt["br"].first)
    CPXsetbranchcallbackfunc(env, hole_branchcallback, &branch);
  else
    CPXsetbranchcallbackfunc(env, empty_branchcallback, &branch);
  if(opt.i_opt["n"].first > 0)
    CPXsetusercutcallbackfunc(env, hole_cutcallback, &cut);
  else
    CPXsetusercutcallbackfunc(env, empty_cutcallback, &cut);

  // set CPX_WORKDIR on cluster!!!
  status = CPXsetstrparam (env, CPX_PARAM_WORKDIR, "/mnt/cluster-tmp/sven/");
  if (status) goto TERMINATE;
  // print log line
  fprintf(fout,"\n%s, inccb %d, brcb %d, maxrounds %4d, cpxcuts off %d, filter %d, time %.2lf\n-------------------------------------------------------------------------------\n",
	  FileName.c_str(),opt.b_opt["inc"].first,opt.b_opt["br"].first,opt.i_opt["n"].first,
	  opt.b_opt["c"].first,opt.b_opt["filter"].first,opt.d_opt["bab_time"].first);
  CPXmipopt(env, lp); 
  // print log lines
  if (opt.i_opt["n"].first>0) fprintf(fout,"... total: %d\n",cuts_tot);
  if (opt.b_opt["inc"].first) fprintf(fout,"-------------------------------\nincumbent rejected/called: %d/%d\n",increjected,inccalled);
  if (opt.b_opt["br"].first) fprintf(fout,"-------------------------------\nbranch strengthenings/changes: %d/%d\n",brstr,brchange);
  CPXgetbestobjval(env, lp, &cutoff);
  CPXgetobjval(env, lp, &objval);

  status = CPXgetstat(env, lp);
  switch (status){
    case CPXMIP_OPTIMAL:
    case CPXMIP_OPTIMAL_TOL:
      opt_stat = "OPTIMAL";
      break;
    case CPXMIP_TIME_LIM_FEAS:
    case CPXMIP_TIME_LIM_INFEAS:
      opt_stat = "TIME_LIMIT";
      break;
    default:
      opt_stat = "OTHER_EXIT";
  }
  printf("\nBRANCH-AND-HOLE: nodes %i bound %.10f sol %.10f %s cuts %d incs %d/%d branch %d/%d\n", CPXgetnodecnt(env, lp), cutoff, objval, opt_stat.c_str(), cuts_tot, increjected, inccalled, brstr, brchange);
  fprintf(fout,"-------------------------------\nFINAL: nodes %i bound %.10f sol %.10f %s\n", CPXgetnodecnt(env, lp), cutoff, objval, opt_stat.c_str());

TERMINATE:

   FREEN (&ctype);

   FREEN (&lbvarswithholes);
   FREEN (&ubvarswithholes);
   FREEN_mat (&lbholes,nvarswithholes);
   FREEN_mat (&ubholes,nvarswithholes);
   FREEN (&nholesvarswithholes);
   FREEN_mat (&namevarswithholes,nvarswithholes);
   FREEN (&has_hole);

   /* Close files */
   if ( fout != NULL ) 
      fclose (fout);

   /* Free LP, if necessary */
   if ( lp != NULL ) {
      int xstatus = CPXfreeprob (env, &lp);
      if ( xstatus )
         printf ("Error in CPXfreeprob, status=%d\n", xstatus);
      if ( !status ) status = xstatus;
   }

   /* Free CPLEX environment, if necessary */
   if ( env != NULL ) {
      int xstatus = CPXcloseCPLEX (&env);
      if ( xstatus )
         printf ("Error in CPXcloseCPLEX, status=%d\n", xstatus);
      if ( !status ) status = xstatus;
   }

   if ( status ) 
      printf ("FINAL STATUS: %d\n", status);

   return status;
}


