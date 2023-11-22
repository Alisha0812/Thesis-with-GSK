/***********************************************
* CMI Data Analysis and Strategy Working group *
*                                              *
* Inferential model for CMI data analsyis      *
* NOV2022                                      *
************************************************/

/* Import the data */

proc import file="C:/Users/fs696631/OneDrive - GSK/PROJECTS/ProjectsCollaborations/CMIworkinggroup/scripts/data/dat_test.csv"
    out=work.dat
    dbms=csv;
run;

/* Sort by group and visit */

proc sort data=dat ;
by gr visit;
run;


/* Compute mean on pre and bkg */
proc means data=dat; 
var y0 x;
output out=desc1 N=n MEAN(y0)=m_pre MEAN(x)=m_bkg;
run;
data desc1;
 set desc1;
 gm_pre= round(10**m_pre,0.0001);
 gm_bkg= round(10**m_bkg, 0.0001);
 key=1;
run;


/* Run the model */

PROC MIXED DATA=dat/*order=data*/; 
   CLASS gr subjid ;
   MODEL y = gr y0 x / DDFM=KR cl outp=pred ;
   /*REPEATED visit / SUBJECT=usubjid TYPE=UNR GROUP=Adj;*/ /* To run when multiple visits are to be included in the model */
   LSMEANS gr /E CL PDIFF=control ("ctr") ALPHA=0.05 COV AT MEANS; /* To adapt based on the specific comparison of interest */
   ODS OUTPUT fitstatistics=stat1 Diffs=diff1 LSMeans=LSM1 COVPARMS=cov tests3=param;
run;

data LSM1;
 set LSM1;
 count+1;
run;
data diff1;
 set diff1;
 count+1;
run;

/* Delta method */
/* GM estimated frequencies */
DATA gm_lsm1;
  SET lsm1;
  FORMAT gm_est LL UL 8.4;
  gm_est = 10**estimate; /* ratio ind/Bkg*/
  LL = 10**lower;
  UL = 10**upper;
  key=1;
RUN;

PROC SQL;
  CREATE TABLE gm_lsm1 AS 
  SELECT a.*, b.m_bkg, b.m_pre, b.gm_pre, b.gm_bkg
  FROM gm_lsm1 AS a 
         LEFT JOIN desc1 AS b
      ON a.key = b.key
;
QUIT;

DATA gm_lsm1;
  SET gm_lsm1;
  freq_spe = (gm_est - 1)*gm_bkg; /* (Ind/Bkg) -1 = (Ind-Bkg)/Bkg --> specific = Ind-Bkg = ((Ind/Bkg) - 1)*Bkg */
  freq_ll = (LL - 1)*gm_bkg;
  freq_ul  = (UL - 1)*gm_bkg;
RUN;
proc sort data=gm_lsm1; by count; run;


/* Difference LSMEAN */

DATA gm_diff1;
  SET diff1;
  FORMAT gm_est LL UL 8.4;
  gm_est = 10**estimate; /* ratio of (Ind/Bkg)*/
  ll = 10**lower;
  ul  = 10**upper;
RUN;

/* GM ratio */

PROC SQL;
  CREATE TABLE gmr1 AS 
  SELECT a.*, b.freq_spe as freq_spe_ref
  FROM gm_lsm1 (where=(gr = "vac")) AS a 
         LEFT JOIN gm_lsm1 (where=(gr = "ctr")) AS b
      ON a.effect = b.effect
;
QUIT;

data gmr;
set gmr1;
gmr=freq_spe/freq_spe_ref;
run;


/* Compute variance */

data gm_lsm2;
set gm_lsm1;
 if gr in ("vac") then zg2= 1+ (1/(gm_est-1)) ; /* group vacc */
 if gr in ("ctr") then zg1= -1-(1/(gm_est-1)) ; /* group ref */
run;


data gm_lsm2;
set gm_lsm2;
 if gr="ctr" then do; var=cov1; covar1=cov2; end;
 else if gr="vac" then do; var=cov2; covar1=cov1; end;
 
run;

PROC SQL;
  CREATE TABLE wk_delta1 AS 
  SELECT a.*, b.zg2, b.var AS var2, b.covar1 as covar , c.zg1, c.var AS var1, d.df as dfdif
  FROM gmr1 AS a 
         LEFT JOIN gm_lsm2 (where=(gr ne "ctr")) AS b
                 ON a.gr = b.gr
         LEFT JOIN gm_lsm2 (WHERE=(gr="ctr")) AS c                                                       
                 ON a.effect = c.effect  
         LEFT JOIN gm_diff1 (where=(_gr="ctr")) AS d                                                       
                 ON a.gr = d.gr  
;
QUIT;

DATA wk_delta1;
  SET wk_delta1;
  gmr=freq_spe/freq_spe_ref;
  X= zg2*(zg2*var2 + zg1*covar) + zg1*(zg2*covar + zg1*var1);
  std=x**0.5;
  COMP="Vac-Ctr";
  RUN;

data wk_delta;
set wk_delta1;
run;

  DATA wk_delta; 
    SET wk_delta;
    logeff=log10(gmr);
    tinv=tinv(0.975,dfdif);
    LL=logeff-tinv*std;
    UL=logeff+tinv*std;
    ratioll=10**ll;
    ratioup=10**ul;
    std_err=std;
    Tval=ABS(logeff/std_err);
    d_probt=2*(1-probt(tval,dfdif));
  RUN;
  



