/************************************************************/
/*                                                          */
/*                                                          */
/*                                                          */
/*      MotLAct7MDB.c                                       */
/*                                                          */
/*                                                          */
/*                                                          */
/*                                                          */
/*      2022/05/19                                          */
/************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#define M_PI 3.1415926535

#define N_assm    8      // the number of cell assemblies 
#define N_T       19     // the number of cell units  
#define N_P       19     // 
#define DT        1e-4   // discrete time 
#define delay     500    // delay between networks
#define cmPY    500e-12  // capacitance for PY 
#define cmSB    243e-12  // for SB
#define cmLB    115e-12  // for LB
#define cmG     30e-12   // for G
#define cmM     224e-12  // for M 
#define gmPY    25.0e-9  // conductance for PY 
#define gmSB    9.7e-9   // for SB
#define gmLB    8.2e-9   // for LB
#define gmG     90.0e-9  // for G
#define gmM     16.0e-9  // for M  
#define gGap    20.0e-9  // 20 nS */
#define gAMPA   0.5e-9   // maximal conductance for AMPA-receptors 
#define gGABA   0.7e-9   // for GABA-receptors 
#define gGABAb  1e-9     // for GABAb receptor 
#define K1      9e4      //  
#define K2      1.2      //  
#define K3      180      //  
#define K4      34       //  
#define Kd      100      //  
#define nBS     4        // binding sites 

#define UPYact    -0.01   // action potential of PYC 
#define UPYres    -0.065  // reset membrane pot. of PYC  
#define USBact    -0.01   // for SB
#define USBres    -0.07   // for SB  
#define ULBact    -0.01   // for LB
#define ULBres    -0.07   // for LB 
#define UGres     -0.07   // for G  
#define UMact     -0.01   // for M
#define UMres     -0.057  // for M 

#define u_AMPA     0.0    // reversal potential for AMPA 
#define u_GABA     -0.08  // for GABA 
#define u_GABAb    -0.095 // for GABAb  
#define steep_PY   220.0  // steepness of sigmoid for PY 
#define thres_PY   -0.03 //  threshold of sigmoid for PY
#define steep_PY2  280.0  // for PY2 
#define thres_PY2  -0.033 // for PY2 
#define steep_SB   260.0  // for SB 
#define thres_SB   -0.033 // for SB    
#define steep_SB2  260.0  // for SB2
#define thres_SB2  -0.033 // for SB2    
#define steep_LB   300.0  // for LB 
#define thres_LB   -0.031 // for LB 
#define steep_LB2  300.0  // for LB2 
#define thres_LB2  -0.031 // for LB2 
#define steep_M    300.0  // for M
#define thres_M    -0.014 // for M

#define alph_AMPA    1.1e6    // open AMPA-channels 
#define beta_AMPA    190.0    // 
#define alph_GABA    5.0e6    // open GABA-channels 
#define beta_GABA    180.0    // 
#define Glut_c       0.001    // glutamate conc.    
#define GABA_cS      0.001    // GABA conc.    
#define GABA_cL1     0.001    // GABA conc. in V1
#define GABA_cL2     0.001    // GABA conc. in V2 

// GABA transporter
double m_G=4e9;               // trans. coef.
#define uG_trn       -0.07    // rev. pot. 

// feature
#define theta_inp0   0       
#define theta_inp1   1       
#define theta_inp2   2       
#define theta_inp3   3       
#define theta_inp4   4        
#define theta_inp5   5       
#define theta_inp6   6       
#define theta_inp7   7      
#define theta_inp8   8       

#define int_inp0_0     0e-12       
#define int_inp0_1     int_inp0_0 
#define int_inp0_2     int_inp0_0  
#define int_inp0_3     700e-12    
#define int_inp0_4     int_inp0_0  
#define int_inp0_5     int_inp0_0  
#define int_inp0_6     int_inp0_0 
#define int_inp0_7     int_inp0_0 
#define inp_prob       1.0      // input probability [0,1]

// amount of extrasyn. recepts.
double delta_T=8e2;             // V1
double delta_T2=8e2;            // V2


#define tau_area     1.0        // input broadness 
#define COLUMN       3          // recorded curnent in cell assemly  
#define NEURON       1          // recorded curnent in cell unit

int t;                          // time development  
#define OUT       5000          // output time  
#define PAT       "ON"          // display activity pattern 
#define INTV      99            // display interval 
#define PERIOD    OUT+30000     // time period 

#define onset_0   10000         // stim. onset  
#define period_0  20000         // stim. period  

#define wLPPV1		0.8		// weight in V1 
#define wLPPV2		0.8	    // weight in V2

#define wMM  		10.0	// weight in M
#define w_M_P       2.8     // P-to-M connection

// GABA conc.
void dfsGABA_ext(void);       
double I_GABA[N_assm+2][N_T+2];     
double GABA_ext[N_assm+2][N_T+2];  
double GABA_extw[N_assm+2][N_T+2]; 
double gamma=10.0;             
double GABA_c0=8E-07;         
double GABAamb_max=4E-06;	  
double GABAamb_min=0E-07;	 
double GABA_V2=10E-07;        

// fraction of open channels for extrasyn. GABA recepts.
double rEXT[N_assm+2][N_T+2]; 
double difrEXT(int,int);        
double rEXT2[N_assm+2][N_T+2];
double difrEXT2(int,int);    

//V1 cells
double I_MGB1[N_assm+2];      // input to tonic cell 
double I_MGBN1[N_assm+2];      
double I_SB[N_assm+2];         
double uPY1[N_assm+2][N_T+2]; // membrane pot. 
double uSB1[N_assm+2][N_T+2]; 
double uLB1[N_assm+2][N_T+2]; 
double uG[N_assm+2][N_T+2];   
double vPY1[N_assm+2][N_T+2][PERIOD+2];  // action potential 
double vSB1[N_assm+2][N_T+2][PERIOD+2]; 
double vLB1[N_assm+2][N_T+2][PERIOD+2];  
double rPY1[N_assm+2][N_T+2];            // fraction of open channels for AMPA-receptors 
double rPY1d[N_assm+2][N_T+2][PERIOD+2]; // delay  
double rSB1[N_assm+2][N_T+2];            
double sF[N_assm+2][N_T+2];    
double rLB1[N_assm+2][N_T+2];  
double GlutPY_c1[N_assm+2][N_T+2];   // glutamate conc. 
double GlutSB_c1[N_assm+2][N_T+2];   
double GlutLB_c1[N_assm+2][N_T+2];    
double duPY_leak1[N_assm+2][N_T+2];  // leak current
double duPY_rec_1[N_assm+2][N_T+2];  // recurrent excitatory current
double duPY_fed_1[N_assm+2][N_T+2];  // feedback inhibitory current
double duPY_lat_1[N_assm+2][N_T+2];  // lateral inhibitory current
double duPY_topdown[N_assm+2][N_T+2];// V2(P) to V1(Ib) current
double duPY_ext_1[N_assm+2][N_T+2];  // inhibitory current (am. GABA)
double duPY_MGB1[N_assm+2];          // excitatory current

double I_leak1;   // leak current
double I_rec_1;   // recurrent excitatory current
double I_fed_1;   // feedback inhibitory current
double I_lat_1;   // lateral inhibitory current
double I_topdown; // V2(P) to V1(Ib) current
double I_ext_1;   // inhibitory current (am. GABA)
double I_MG1;     // excitatory current
double duSB_leak1[N_assm+2][N_T+2]; 
double duSB_1[N_assm+2][N_T+2];
double duSB_topdown[N_assm+2][N_T+2];  
double duSB_ext_1;   
               

double I_leakSB1; 
double I_fed_SB1; 
double I_ext_SB1; 
double duSB_inp;
double duLB_leak1[N_assm+2][N_T+2]; 
double duLB1[N_assm+2][N_T+2];
double duLB1_topdown[N_assm+2][N_T+2]; 
double duLB_ext1;                  

double I_leakLB1; 
double I_fed_LB1; 
double I_ext_LB1; 
double duG_leak[N_assm+2][N_T+2];       
double duG_P[N_assm+2][N_T+2];           
double duG_Ia[N_assm+2][N_T+2];         
double duG_topdown[N_assm+2][N_T+2];   
double duG_G[N_assm+2][N_T+2];      
double drPY1[N_assm+2][N_T+2];   
double drSB1[N_assm+2][N_T+2];    
double dsF[N_assm+2][N_T+2];; 
double drLB1[N_assm+2][N_T+2];    

double AMPA_c1[N_assm+2][N_T+2];  
double GABASB_c1[N_assm+2][N_T+2][PERIOD+2]; 
double GABALB_c1[N_assm+2][N_T+2][PERIOD+2];  

// weights
double w_rec_1[N_assm+2][N_T+2][N_assm+2][N_T+2];  
double w_rec_10[N_assm+2][N_T+2][N_assm+2][N_T+2]; 
double w_fed_1[N_assm+2][N_T+2];                  
double w_lat_1[N_assm+2][N_T+2][N_T+2]; 
double w_v1Ia_v2[N_assm+2][N_T+2][N_assm+2][N_T+2];
double w_v1Ib_v2[N_assm+2][N_T+2][N_assm+2][N_T+2];
double w_v1P_v2[N_assm+2][N_T+2][N_assm+2][N_T+2]; 
double w_v1G_v2[N_assm+2][N_T+2][N_assm+2][N_T+2]; 
double wSB_PY1[N_assm+2][N_T+2];                   
double wLB_PY1[N_assm+2][N_assm+2][N_T+2];         
double wG_P[N_assm+2];                             
double wG_Ia[N_assm+2];                           

// difference
double difuPY1(int,int);    
double difuSB1(int,int);    
double difuLB1(int,int);   
double difrPY1(int,int);   
double difrSB1(int,int);   
double difsF(int,int);     
double difrLB1(int,int);   

//V2 cells
double I_V2[N_assm+2];        
double uPY2[N_assm+2][N_T+2]; 
double uSB2[N_assm+2][N_T+2]; 
double uLB2[N_assm+2][N_T+2]; 
double uM[N_assm+2][N_T+2]; 
double vPY2[N_assm+2][N_T+2][PERIOD+2]; 
double vSB2[N_assm+2][N_T+2][PERIOD+2]; 
double vLB2[N_assm+2][N_T+2][PERIOD+2]; 
double vM[N_assm+2][N_T+2][PERIOD+2]; 
double rPY2[N_assm+2][N_T+2]; 
double rPY2d[N_assm+2][N_T+2][PERIOD+2]; 
double rSB2[N_assm+2][N_T+2]; 
double rLB2[N_assm+2][N_T+2]; 
double rM[N_assm+2][N_T+2]; 
double rMd[N_assm+2][N_T+2][PERIOD+2]; 
double GlutPY_c2[N_assm+2][N_T+2];  
double GlutSB_c2[N_assm+2][N_T+2];   
double GlutLB_c2[N_assm+2][N_T+2]; 
double GlutM[N_assm+2][N_T+2];   
double duPY_leak2[N_assm+2][N_T+2];  
double duPY_rec_2[N_assm+2][N_T+2];  
double duPY_fed_2[N_assm+2][N_T+2];  
double duPY_lat_2[N_assm+2][N_T+2];  
double duPY_ext_2[N_assm+2][N_T+2];  
double duPY_MGB2[N_assm+2];          

double duM_leak[N_assm+2][N_T+2];  
double duM_rec[N_assm+2][N_T+2];  
double duM_fed[N_assm+2][N_T+2]; 
double duM_lat[N_assm+2][N_T+2];  
double duM_ext[N_assm+2][N_T+2];  
 
double I_leak2; 
double I_rec_2; 
double I_fed_2; 
double I_lat_2; 
double I_ext_2; 
double duPY_inp;
double duSB_leak2[N_assm+2][N_T+2]; 
double duSB_2[N_assm+2][N_T+2];
double duSB_ext_2;    
               
double I_leakSB2; 
double I_fed_SB2; 
double I_ext_SB2; 
double duLB_leak2[N_assm+2][N_T+2];     
double duLB2[N_assm+2][N_T+2];
double duLB2_bottom_up[N_assm+2][N_T+2];
double duLB_ext2;                        

double I_leakLB2; 
double I_fed_LB2; 
double I_ext_LB2; 
double drPY2[N_assm+2][N_T+2]; 
double drSB2[N_assm+2][N_T+2];
double drLB2[N_assm+2][N_T+2]; 
double drM[N_assm+2][N_T+2]; 

double AMPA_c2[N_assm+2][N_T+2];  
double GABASB_c2[N_assm+2][N_T+2][PERIOD+2]; 
double GABALB_c2[N_assm+2][N_T+2][PERIOD+2];  
double w_rec_2[N_assm+2][N_T+2][N_assm+2][N_T+2];            
double w_fed_2[N_assm+2][N_T+2];                  
double w_lat_2[N_assm+2][N_T+2][N_T+2];
double w_v2_v1[N_assm+2][N_T+2][N_assm+2][N_T+2] ;
double wSB_PY2[N_assm+2][N_T+2];                 
double wLB_PY2[N_assm+2][N_assm+2][N_T+2];       
double difuPY2(int,int);  
double difuSB2(int,int);    
double difuLB2(int,int); 
double difuM(int,int);  
double difrPY2(int,int);   
double difrSB2(int,int);  
double difrLB2(int,int);   
double difrM(int,int);   

double sigmoidPY(double);  // sigmoid for PY in V1
double sigmoidPY2(double); // sigmoid for PY in V2
double sigmoidSB(double); 
double sigmoidSB2(double); 
double sigmoidLB(double); 
double sigmoidLB2(double); 
double sigmoidM(double); 
double OFSWCA=0.01;        
double OFSBCA=0.12;        
int    xi[N_assm+2][N_assm+2][N_T+2];
double rand01(long int *); // radom number [0, 1.0]
void display(void);        

void init(void);           // initialization  
long int SEEDMP=7000;      // a seed for random number generation

//V1 cells
FILE *uPY0_01,*uPY0_11,*uPY0_21,*uPY0_31,*uPY0_41,*uPY0_51,*uPY0_61,*uPY0_71,*uPY0_81,*uPY0_91;
FILE *uPY1_01,*uPY1_11,*uPY1_21,*uPY1_31,*uPY1_41,*uPY1_51,*uPY1_61,*uPY1_71,*uPY1_81,*uPY1_91;
FILE *uPY2_01,*uPY2_11,*uPY2_21,*uPY2_31,*uPY2_41,*uPY2_51,*uPY2_61,*uPY2_71,*uPY2_81,*uPY2_91;
FILE *uPY3_01,*uPY3_11,*uPY3_21,*uPY3_31,*uPY3_41,*uPY3_51,*uPY3_61,*uPY3_71,*uPY3_81,*uPY3_91;
FILE *uPY4_01,*uPY4_11,*uPY4_21,*uPY4_31,*uPY4_41,*uPY4_51,*uPY4_61,*uPY4_71,*uPY4_81,*uPY4_91;
FILE *uPY5_01,*uPY5_11,*uPY5_21,*uPY5_31,*uPY5_41,*uPY5_51,*uPY5_61,*uPY5_71,*uPY5_81,*uPY5_91;
FILE *uPY6_01,*uPY6_11,*uPY6_21,*uPY6_31,*uPY6_41,*uPY6_51,*uPY6_61,*uPY6_71,*uPY6_81,*uPY6_91;
FILE *uPY7_01,*uPY7_11,*uPY7_21,*uPY7_31,*uPY7_41,*uPY7_51,*uPY7_61,*uPY7_71,*uPY7_81,*uPY7_91;
FILE *leak1,*rec1,*fed1,*lat1,*ext1,*topdown,*MG1; 
FILE *extSB1,*extLB1;
FILE *vPY0_01,*vPY0_11,*vPY0_21,*vPY0_31,*vPY0_41,*vPY0_51,*vPY0_61,*vPY0_71,*vPY0_81,*vPY0_91;
FILE *vPY1_01,*vPY1_11,*vPY1_21,*vPY1_31,*vPY1_41,*vPY1_51,*vPY1_61,*vPY1_71,*vPY1_81,*vPY1_91;
FILE *vPY2_01,*vPY2_11,*vPY2_21,*vPY2_31,*vPY2_41,*vPY2_51,*vPY2_61,*vPY2_71,*vPY2_81,*vPY2_91;
FILE *vPY3_01,*vPY3_11,*vPY3_21,*vPY3_31,*vPY3_41,*vPY3_51,*vPY3_61,*vPY3_71,*vPY3_81,*vPY3_91;
FILE *vPY4_01,*vPY4_11,*vPY4_21,*vPY4_31,*vPY4_41,*vPY4_51,*vPY4_61,*vPY4_71,*vPY4_81,*vPY4_91;
FILE *vPY5_01,*vPY5_11,*vPY5_21,*vPY5_31,*vPY5_41,*vPY5_51,*vPY5_61,*vPY5_71,*vPY5_81,*vPY5_91;
FILE *vPY6_01,*vPY6_11,*vPY6_21,*vPY6_31,*vPY6_41,*vPY6_51,*vPY6_61,*vPY6_71,*vPY6_81,*vPY6_91;
FILE *vPY7_01,*vPY7_11,*vPY7_21,*vPY7_31,*vPY7_41,*vPY7_51,*vPY7_61,*vPY7_71,*vPY7_81,*vPY7_91;

FILE *uSB0_01,*uSB1_01,*uSB2_01,*uSB3_01,*uSB4_01,*uSB5_01,*uSB6_01,*uSB7_01;

FILE *vSB0_01,*vSB0_11,*vSB0_21,*vSB0_31,*vSB0_41,*vSB0_51,*vSB0_61,*vSB0_71,*vSB0_81,*vSB0_91;
FILE *vSB1_01,*vSB1_11,*vSB1_21,*vSB1_31,*vSB1_41,*vSB1_51,*vSB1_61,*vSB1_71,*vSB1_81,*vSB1_91;
FILE *vSB2_01,*vSB2_11,*vSB2_21,*vSB2_31,*vSB2_41,*vSB2_51,*vSB2_61,*vSB2_71,*vSB2_81,*vSB2_91;
FILE *vSB3_01,*vSB3_11,*vSB3_21,*vSB3_31,*vSB3_41,*vSB3_51,*vSB3_61,*vSB3_71,*vSB3_81,*vSB3_91;
FILE *vSB4_01,*vSB4_11,*vSB4_21,*vSB4_31,*vSB4_41,*vSB4_51,*vSB4_61,*vSB4_71,*vSB4_81,*vSB4_91;
FILE *vSB5_01,*vSB5_11,*vSB5_21,*vSB5_31,*vSB5_41,*vSB5_51,*vSB5_61,*vSB5_71,*vSB5_81,*vSB5_91;
FILE *vSB6_01,*vSB6_11,*vSB6_21,*vSB6_31,*vSB6_41,*vSB6_51,*vSB6_61,*vSB6_71,*vSB6_81,*vSB6_91;
FILE *vSB7_01,*vSB7_11,*vSB7_21,*vSB7_31,*vSB7_41,*vSB7_51,*vSB7_61,*vSB7_71,*vSB7_81,*vSB7_91;

FILE *uLB0_01,*uLB1_01,*uLB2_01,*uLB3_01,*uLB4_01,*uLB5_01,*uLB6_01,*uLB7_01;

FILE *vLB0_01,*vLB0_11,*vLB0_21,*vLB0_31,*vLB0_41,*vLB0_51,*vLB0_61,*vLB0_71,*vLB0_81,*vLB0_91;
FILE *vLB1_01,*vLB1_11,*vLB1_21,*vLB1_31,*vLB1_41,*vLB1_51,*vLB1_61,*vLB1_71,*vLB1_81,*vLB1_91;
FILE *vLB2_01,*vLB2_11,*vLB2_21,*vLB2_31,*vLB2_41,*vLB2_51,*vLB2_61,*vLB2_71,*vLB2_81,*vLB2_91;
FILE *vLB3_01,*vLB3_11,*vLB3_21,*vLB3_31,*vLB3_41,*vLB3_51,*vLB3_61,*vLB3_71,*vLB3_81,*vLB3_91;
FILE *vLB4_01,*vLB4_11,*vLB4_21,*vLB4_31,*vLB4_41,*vLB4_51,*vLB4_61,*vLB4_71,*vLB4_81,*vLB4_91;
FILE *vLB5_01,*vLB5_11,*vLB5_21,*vLB5_31,*vLB5_41,*vLB5_51,*vLB5_61,*vLB5_71,*vLB5_81,*vLB5_91;
FILE *vLB6_01,*vLB6_11,*vLB6_21,*vLB6_31,*vLB6_41,*vLB6_51,*vLB6_61,*vLB6_71,*vLB6_81,*vLB6_91;
FILE *vLB7_01,*vLB7_11,*vLB7_21,*vLB7_31,*vLB7_41,*vLB7_51,*vLB7_61,*vLB7_71,*vLB7_81,*vLB7_91;

//V2 cells
FILE *uPY0_02,*uPY0_12,*uPY0_22,*uPY0_32,*uPY0_42,*uPY0_52,*uPY0_62,*uPY0_72,*uPY0_82,*uPY0_92;
FILE *uPY1_02,*uPY1_12,*uPY1_22,*uPY1_32,*uPY1_42,*uPY1_52,*uPY1_62,*uPY1_72,*uPY1_82,*uPY1_92;
FILE *uPY2_02,*uPY2_12,*uPY2_22,*uPY2_32,*uPY2_42,*uPY2_52,*uPY2_62,*uPY2_72,*uPY2_82,*uPY2_92;
FILE *uPY3_02,*uPY3_12,*uPY3_22,*uPY3_32,*uPY3_42,*uPY3_52,*uPY3_62,*uPY3_72,*uPY3_82,*uPY3_92;
FILE *uPY4_02,*uPY4_12,*uPY4_22,*uPY4_32,*uPY4_42,*uPY4_52,*uPY4_62,*uPY4_72,*uPY4_82,*uPY4_92;
FILE *uPY5_02,*uPY5_12,*uPY5_22,*uPY5_32,*uPY5_42,*uPY5_52,*uPY5_62,*uPY5_72,*uPY5_82,*uPY5_92;
FILE *uPY6_02,*uPY6_12,*uPY6_22,*uPY6_32,*uPY6_42,*uPY6_52,*uPY6_62,*uPY6_72,*uPY6_82,*uPY6_92;
FILE *uPY7_02,*uPY7_12,*uPY7_22,*uPY7_32,*uPY7_42,*uPY7_52,*uPY7_62,*uPY7_72,*uPY7_82,*uPY7_92;
FILE *leak2,*rec2,*fed2,*lat2,*ext2;  
FILE *extSB2,*extLB2;
FILE *vPY0_02,*vPY0_12,*vPY0_22,*vPY0_32,*vPY0_42,*vPY0_52,*vPY0_62,*vPY0_72,*vPY0_82,*vPY0_92;
FILE *vPY1_02,*vPY1_12,*vPY1_22,*vPY1_32,*vPY1_42,*vPY1_52,*vPY1_62,*vPY1_72,*vPY1_82,*vPY1_92;
FILE *vPY2_02,*vPY2_12,*vPY2_22,*vPY2_32,*vPY2_42,*vPY2_52,*vPY2_62,*vPY2_72,*vPY2_82,*vPY2_92;
FILE *vPY3_02,*vPY3_12,*vPY3_22,*vPY3_32,*vPY3_42,*vPY3_52,*vPY3_62,*vPY3_72,*vPY3_82,*vPY3_92;
FILE *vPY4_02,*vPY4_12,*vPY4_22,*vPY4_32,*vPY4_42,*vPY4_52,*vPY4_62,*vPY4_72,*vPY4_82,*vPY4_92;
FILE *vPY5_02,*vPY5_12,*vPY5_22,*vPY5_32,*vPY5_42,*vPY5_52,*vPY5_62,*vPY5_72,*vPY5_82,*vPY5_92;
FILE *vPY6_02,*vPY6_12,*vPY6_22,*vPY6_32,*vPY6_42,*vPY6_52,*vPY6_62,*vPY6_72,*vPY6_82,*vPY6_92;
FILE *vPY7_02,*vPY7_12,*vPY7_22,*vPY7_32,*vPY7_42,*vPY7_52,*vPY7_62,*vPY7_72,*vPY7_82,*vPY7_92;
FILE *uLB0_02,*uLB1_02,*uLB2_02,*uLB3_02,*uLB4_02,*uLB5_02,*uLB6_02,*uLB7_02;

FILE *vLB0_02,*vLB0_12,*vLB0_22,*vLB0_32,*vLB0_42,*vLB0_52,*vLB0_62,*vLB0_72,*vLB0_82,*vLB0_92;
FILE *vLB1_02,*vLB1_12,*vLB1_22,*vLB1_32,*vLB1_42,*vLB1_52,*vLB1_62,*vLB1_72,*vLB1_82,*vLB1_92;
FILE *vLB2_02,*vLB2_12,*vLB2_22,*vLB2_32,*vLB2_42,*vLB2_52,*vLB2_62,*vLB2_72,*vLB2_82,*vLB2_92;
FILE *vLB3_02,*vLB3_12,*vLB3_22,*vLB3_32,*vLB3_42,*vLB3_52,*vLB3_62,*vLB3_72,*vLB3_82,*vLB3_92;
FILE *vLB4_02,*vLB4_12,*vLB4_22,*vLB4_32,*vLB4_42,*vLB4_52,*vLB4_62,*vLB4_72,*vLB4_82,*vLB4_92;
FILE *vLB5_02,*vLB5_12,*vLB5_22,*vLB5_32,*vLB5_42,*vLB5_52,*vLB5_62,*vLB5_72,*vLB5_82,*vLB5_92;
FILE *vLB6_02,*vLB6_12,*vLB6_22,*vLB6_32,*vLB6_42,*vLB6_52,*vLB6_62,*vLB6_72,*vLB6_82,*vLB6_92;
FILE *vLB7_02,*vLB7_12,*vLB7_22,*vLB7_32,*vLB7_42,*vLB7_52,*vLB7_62,*vLB7_72,*vLB7_82,*vLB7_92;

FILE *uM0_02,*uM0_12,*uM0_22,*uM0_32,*uM0_42,*uM0_52,*uM0_62,*uM0_72,*uM0_82,*uM0_92;
FILE *uM1_02,*uM1_12,*uM1_22,*uM1_32,*uM1_42,*uM1_52,*uM1_62,*uM1_72,*uM1_82,*uM1_92;
FILE *uM2_02,*uM2_12,*uM2_22,*uM2_32,*uM2_42,*uM2_52,*uM2_62,*uM2_72,*uM2_82,*uM2_92;
FILE *uM3_02,*uM3_12,*uM3_22,*uM3_32,*uM3_42,*uM3_52,*uM3_62,*uM3_72,*uM3_82,*uM3_92;
FILE *uM4_02,*uM4_12,*uM4_22,*uM4_32,*uM4_42,*uM4_52,*uM4_62,*uM4_72,*uM4_82,*uM4_92;
FILE *uM5_02,*uM5_12,*uM5_22,*uM5_32,*uM5_42,*uM5_52,*uM5_62,*uM5_72,*uM5_82,*uM5_92;
FILE *uM6_02,*uM6_12,*uM6_22,*uM6_32,*uM6_42,*uM6_52,*uM6_62,*uM6_72,*uM6_82,*uM6_92;
FILE *uM7_02,*uM7_12,*uM7_22,*uM7_32,*uM7_42,*uM7_52,*uM7_62,*uM7_72,*uM7_82,*uM7_92;

FILE *vM0_02,*vM0_12,*vM0_22,*vM0_32,*vM0_42,*vM0_52,*vM0_62,*vM0_72,*vM0_82,*vM0_92;
FILE *vM1_02,*vM1_12,*vM1_22,*vM1_32,*vM1_42,*vM1_52,*vM1_62,*vM1_72,*vM1_82,*vM1_92;
FILE *vM2_02,*vM2_12,*vM2_22,*vM2_32,*vM2_42,*vM2_52,*vM2_62,*vM2_72,*vM2_82,*vM2_92;
FILE *vM3_02,*vM3_12,*vM3_22,*vM3_32,*vM3_42,*vM3_52,*vM3_62,*vM3_72,*vM3_82,*vM3_92;
FILE *vM4_02,*vM4_12,*vM4_22,*vM4_32,*vM4_42,*vM4_52,*vM4_62,*vM4_72,*vM4_82,*vM4_92;
FILE *vM5_02,*vM5_12,*vM5_22,*vM5_32,*vM5_42,*vM5_52,*vM5_62,*vM5_72,*vM5_82,*vM5_92;
FILE *vM6_02,*vM6_12,*vM6_22,*vM6_32,*vM6_42,*vM6_52,*vM6_62,*vM6_72,*vM6_82,*vM6_92;
FILE *vM7_02,*vM7_12,*vM7_22,*vM7_32,*vM7_42,*vM7_52,*vM7_62,*vM7_72,*vM7_82,*vM7_92;

void main(void){

	int theta,i,ii;
	double sigPYT,sigSBT,sigLBT,sigMT;
	int ActPY1_ON[N_assm+2][N_T+2],ActSB1_ON[N_assm+2][N_T+2],ActLB1_ON[N_assm+2][N_T+2],duration;
	int ActPY2_ON[N_assm+2][N_T+2],ActSB2_ON[N_assm+2][N_T+2],ActLB2_ON[N_assm+2][N_T+2],ActM_ON[N_assm+2][N_T+2];
	int ActPY1_OF[N_assm+2][N_T+2],ActSB1_OF[N_assm+2][N_T+2],ActLB1_OF[N_assm+2][N_T+2];
	int ActPY2_OF[N_assm+2][N_T+2],ActSB2_OF[N_assm+2][N_T+2],ActLB2_OF[N_assm+2][N_T+2],ActM_OF[N_assm+2][N_T+2];

	init(); 
    
	// V1 cells
	uPY0_01=fopen("uPY0_01.dat","w");
	uPY1_01=fopen("uPY1_01.dat","w");
	uPY2_01=fopen("uPY2_01.dat","w");
	uPY3_01=fopen("uPY3_01.dat","w");
	uPY4_01=fopen("uPY4_01.dat","w");
	uPY5_01=fopen("uPY5_01.dat","w");
	uPY6_01=fopen("uPY6_01.dat","w");
	uPY7_01=fopen("uPY7_01.dat","w");


	leak1  =fopen("leak1c.dat","w"); // a given PYC‚Ö‚Ìinput currents 
	rec1   =fopen("rec1c.dat","w");
	fed1   =fopen("fed1c.dat","w");
	lat1  =fopen("lat1c.dat","w");
	ext1   =fopen("ext1c.dat","w");
	topdown   =fopen("topdownc.dat","w");
	MG1    =fopen("MG1c.dat","w");

	extSB1  =fopen("extSBc1.dat","w");
	extLB1  =fopen("extLBc1.dat","w");

	vPY0_01=fopen("vPY0_01.dat","w");
	vPY0_11=fopen("vPY0_11.dat","w");
	vPY0_21=fopen("vPY0_21.dat","w");
	vPY0_31=fopen("vPY0_31.dat","w");
	vPY0_41=fopen("vPY0_41.dat","w");
	vPY0_51=fopen("vPY0_51.dat","w");
	vPY0_61=fopen("vPY0_61.dat","w");
	vPY0_71=fopen("vPY0_71.dat","w");
	vPY0_81=fopen("vPY0_81.dat","w");
	vPY0_91=fopen("vPY0_91.dat","w");
	vPY1_01=fopen("vPY1_01.dat","w");
	vPY1_11=fopen("vPY1_11.dat","w");
	vPY1_21=fopen("vPY1_21.dat","w");
	vPY1_31=fopen("vPY1_31.dat","w");
	vPY1_41=fopen("vPY1_41.dat","w");
	vPY1_51=fopen("vPY1_51.dat","w");
	vPY1_61=fopen("vPY1_61.dat","w");
	vPY1_71=fopen("vPY1_71.dat","w");
	vPY1_81=fopen("vPY1_81.dat","w");
	vPY1_91=fopen("vPY1_91.dat","w");
	vPY2_01=fopen("vPY2_01.dat","w");
	vPY2_11=fopen("vPY2_11.dat","w");
	vPY2_21=fopen("vPY2_21.dat","w");
	vPY2_31=fopen("vPY2_31.dat","w");
	vPY2_41=fopen("vPY2_41.dat","w");
	vPY2_51=fopen("vPY2_51.dat","w");
	vPY2_61=fopen("vPY2_61.dat","w");
	vPY2_71=fopen("vPY2_71.dat","w");
	vPY2_81=fopen("vPY2_81.dat","w");
	vPY2_91=fopen("vPY2_91.dat","w");
	vPY3_01=fopen("vPY3_01.dat","w");
	vPY3_11=fopen("vPY3_11.dat","w");
	vPY3_21=fopen("vPY3_21.dat","w");
	vPY3_31=fopen("vPY3_31.dat","w");
	vPY3_41=fopen("vPY3_41.dat","w");
	vPY3_51=fopen("vPY3_51.dat","w");
	vPY3_61=fopen("vPY3_61.dat","w");
	vPY3_71=fopen("vPY3_71.dat","w");
	vPY3_81=fopen("vPY3_81.dat","w");
	vPY3_91=fopen("vPY3_91.dat","w");
	vPY4_01=fopen("vPY4_01.dat","w");
	vPY4_11=fopen("vPY4_11.dat","w");
	vPY4_21=fopen("vPY4_21.dat","w");
	vPY4_31=fopen("vPY4_31.dat","w");
	vPY4_41=fopen("vPY4_41.dat","w");
	vPY4_51=fopen("vPY4_51.dat","w");
	vPY4_61=fopen("vPY4_61.dat","w");
	vPY4_71=fopen("vPY4_71.dat","w");
	vPY4_81=fopen("vPY4_81.dat","w");
	vPY4_91=fopen("vPY4_91.dat","w");
	vPY5_01=fopen("vPY5_01.dat","w");
	vPY5_11=fopen("vPY5_11.dat","w");
	vPY5_21=fopen("vPY5_21.dat","w");
	vPY5_31=fopen("vPY5_31.dat","w");
	vPY5_41=fopen("vPY5_41.dat","w");
	vPY5_51=fopen("vPY5_51.dat","w");
	vPY5_61=fopen("vPY5_61.dat","w");
	vPY5_71=fopen("vPY5_71.dat","w");
	vPY5_81=fopen("vPY5_81.dat","w");
	vPY5_91=fopen("vPY5_91.dat","w");
	vPY6_01=fopen("vPY6_01.dat","w");
	vPY6_11=fopen("vPY6_11.dat","w");
	vPY6_21=fopen("vPY6_21.dat","w");
	vPY6_31=fopen("vPY6_31.dat","w");
	vPY6_41=fopen("vPY6_41.dat","w");
	vPY6_51=fopen("vPY6_51.dat","w");
	vPY6_61=fopen("vPY6_61.dat","w");
	vPY6_71=fopen("vPY6_71.dat","w");
	vPY6_81=fopen("vPY6_81.dat","w");
	vPY6_91=fopen("vPY6_91.dat","w");
	vPY7_01=fopen("vPY7_01.dat","w");
	vPY7_11=fopen("vPY7_11.dat","w");
	vPY7_21=fopen("vPY7_21.dat","w");
	vPY7_31=fopen("vPY7_31.dat","w");
	vPY7_41=fopen("vPY7_41.dat","w");
	vPY7_51=fopen("vPY7_51.dat","w");
	vPY7_61=fopen("vPY7_61.dat","w");
	vPY7_71=fopen("vPY7_71.dat","w");
	vPY7_81=fopen("vPY7_81.dat","w");
	vPY7_91=fopen("vPY7_91.dat","w");

	uSB0_01=fopen("uSB0_01.dat","w");
	uSB1_01=fopen("uSB1_01.dat","w");
	uSB2_01=fopen("uSB2_01.dat","w");
	uSB3_01=fopen("uSB3_01.dat","w");
	uSB4_01=fopen("uSB4_01.dat","w");
	uSB5_01=fopen("uSB5_01.dat","w");
	uSB6_01=fopen("uSB6_01.dat","w");
	uSB7_01=fopen("uSB7_01.dat","w");

	uLB0_01=fopen("uLB0_01.dat","w");
	uLB1_01=fopen("uLB1_01.dat","w");
	uLB2_01=fopen("uLB2_01.dat","w");
	uLB3_01=fopen("uLB3_01.dat","w");
	uLB4_01=fopen("uLB4_01.dat","w");
	uLB5_01=fopen("uLB5_01.dat","w");
	uLB6_01=fopen("uLB6_01.dat","w");
	uLB7_01=fopen("uLB7_01.dat","w");

	// V2 cells
	uPY0_02=fopen("uPY0_02.dat","w");
	uPY1_02=fopen("uPY1_02.dat","w");
	uPY2_02=fopen("uPY2_02.dat","w");
	uPY3_02=fopen("uPY3_02.dat","w");
	uPY4_02=fopen("uPY4_02.dat","w");
	uPY5_02=fopen("uPY5_02.dat","w");
	uPY6_02=fopen("uPY6_02.dat","w");
	uPY7_02=fopen("uPY7_02.dat","w");

	vPY0_02=fopen("vPY0_02.dat","w");
	vPY0_12=fopen("vPY0_12.dat","w");
	vPY0_22=fopen("vPY0_22.dat","w");
	vPY0_32=fopen("vPY0_32.dat","w");
	vPY0_42=fopen("vPY0_42.dat","w");
	vPY0_52=fopen("vPY0_52.dat","w");
	vPY0_62=fopen("vPY0_62.dat","w");
	vPY0_72=fopen("vPY0_72.dat","w");
	vPY0_82=fopen("vPY0_82.dat","w");
	vPY0_92=fopen("vPY0_92.dat","w");
	vPY1_02=fopen("vPY1_02.dat","w");
	vPY1_12=fopen("vPY1_12.dat","w");
	vPY1_22=fopen("vPY1_22.dat","w");
	vPY1_32=fopen("vPY1_32.dat","w");
	vPY1_42=fopen("vPY1_42.dat","w");
	vPY1_52=fopen("vPY1_52.dat","w");
	vPY1_62=fopen("vPY1_62.dat","w");
	vPY1_72=fopen("vPY1_72.dat","w");
	vPY1_82=fopen("vPY1_82.dat","w");
	vPY1_92=fopen("vPY1_92.dat","w");
	vPY2_02=fopen("vPY2_02.dat","w");
	vPY2_12=fopen("vPY2_12.dat","w");
	vPY2_22=fopen("vPY2_22.dat","w");
	vPY2_32=fopen("vPY2_32.dat","w");
	vPY2_42=fopen("vPY2_42.dat","w");
	vPY2_52=fopen("vPY2_52.dat","w");
	vPY2_62=fopen("vPY2_62.dat","w");
	vPY2_72=fopen("vPY2_72.dat","w");
	vPY2_82=fopen("vPY2_82.dat","w");
	vPY2_92=fopen("vPY2_92.dat","w");
	vPY3_02=fopen("vPY3_02.dat","w");
	vPY3_12=fopen("vPY3_12.dat","w");
	vPY3_22=fopen("vPY3_22.dat","w");
	vPY3_32=fopen("vPY3_32.dat","w");
	vPY3_42=fopen("vPY3_42.dat","w");
	vPY3_52=fopen("vPY3_52.dat","w");
	vPY3_62=fopen("vPY3_62.dat","w");
	vPY3_72=fopen("vPY3_72.dat","w");
	vPY3_82=fopen("vPY3_82.dat","w");
	vPY3_92=fopen("vPY3_92.dat","w");
	
	vPY4_02=fopen("vPY4_02.dat","w");
	vPY4_12=fopen("vPY4_12.dat","w");
	vPY4_22=fopen("vPY4_22.dat","w");
	vPY4_32=fopen("vPY4_32.dat","w");
	vPY4_42=fopen("vPY4_42.dat","w");
	vPY4_52=fopen("vPY4_52.dat","w");
	vPY4_62=fopen("vPY4_62.dat","w");
	vPY4_72=fopen("vPY4_72.dat","w");
	vPY4_82=fopen("vPY4_82.dat","w");
	vPY4_92=fopen("vPY4_92.dat","w");
	
	vPY5_02=fopen("vPY5_02.dat","w");
	vPY5_12=fopen("vPY5_12.dat","w");
	vPY5_22=fopen("vPY5_22.dat","w");
	vPY5_32=fopen("vPY5_32.dat","w");
	vPY5_42=fopen("vPY5_42.dat","w");
	vPY5_52=fopen("vPY5_52.dat","w");
	vPY5_62=fopen("vPY5_62.dat","w");
	vPY5_72=fopen("vPY5_72.dat","w");
	vPY5_82=fopen("vPY5_82.dat","w");
	vPY5_92=fopen("vPY5_92.dat","w");
	
	vPY6_02=fopen("vPY6_02.dat","w");
	vPY6_12=fopen("vPY6_12.dat","w");
	vPY6_22=fopen("vPY6_22.dat","w");
	vPY6_32=fopen("vPY6_32.dat","w");
	vPY6_42=fopen("vPY6_42.dat","w");
	vPY6_52=fopen("vPY6_52.dat","w");
	vPY6_62=fopen("vPY6_62.dat","w");
	vPY6_72=fopen("vPY6_72.dat","w");
	vPY6_82=fopen("vPY6_82.dat","w");
	vPY6_92=fopen("vPY6_92.dat","w");
	vPY7_02=fopen("vPY7_02.dat","w");
	vPY7_12=fopen("vPY7_12.dat","w");
	vPY7_22=fopen("vPY7_22.dat","w");
	vPY7_32=fopen("vPY7_32.dat","w");
	vPY7_42=fopen("vPY7_42.dat","w");
	vPY7_52=fopen("vPY7_52.dat","w");
	vPY7_62=fopen("vPY7_62.dat","w");
	vPY7_72=fopen("vPY7_72.dat","w");
	vPY7_82=fopen("vPY7_82.dat","w");
	vPY7_92=fopen("vPY7_92.dat","w");

	uLB0_02=fopen("uLB0_02.dat","w");
	uLB1_02=fopen("uLB1_02.dat","w");
	uLB2_02=fopen("uLB2_02.dat","w");
	uLB3_02=fopen("uLB3_02.dat","w");
	uLB4_02=fopen("uLB4_02.dat","w");
	uLB5_02=fopen("uLB5_02.dat","w");
	uLB6_02=fopen("uLB6_02.dat","w");
	uLB7_02=fopen("uLB7_02.dat","w");

	uM0_02=fopen("uM0_02.dat","w");
	uM1_02=fopen("uM1_02.dat","w");
	uM2_02=fopen("uM2_02.dat","w");
	uM3_02=fopen("uM3_02.dat","w");
	uM4_02=fopen("uM4_02.dat","w");
	uM5_02=fopen("uM5_02.dat","w");
	uM6_02=fopen("uM6_02.dat","w");
	uM7_02=fopen("uM7_02.dat","w");

	vM3_02=fopen("vM3_02.dat","w");
	vM3_12=fopen("vM3_12.dat","w");
	vM3_22=fopen("vM3_22.dat","w");
	vM3_32=fopen("vM3_32.dat","w");
	vM3_42=fopen("vM3_42.dat","w");
	vM3_52=fopen("vM3_52.dat","w");
	vM3_62=fopen("vM3_62.dat","w");
	vM3_72=fopen("vM3_72.dat","w");
	vM3_82=fopen("vM3_82.dat","w");
	vM3_92=fopen("vM3_92.dat","w");

	for (theta=0; theta<=N_assm; ++theta){  
		for (ii=0; ii<=N_T; ++ii){
			ActPY1_ON[theta][ii] = 0; 
			ActSB1_ON[theta][ii] = 0; 
			ActLB1_ON[theta][ii] = 0; 
			ActPY1_OF[theta][ii] = 0; 
			ActSB1_OF[theta][ii] = 0; 
			ActLB1_OF[theta][ii] = 0; 
			uPY1[theta][ii]=UPYres;
			uSB1[theta][ii]=USBres;
			uLB1[theta][ii]=ULBres;
		}
	}
	for (theta=0; theta<=N_assm; ++theta){  
		for (ii=0; ii<=N_T; ++ii){
			ActPY2_ON[theta][ii] = 0; 
			ActM_ON[theta][ii] = 0; 
			ActSB2_ON[theta][ii] = 0; 
			ActLB2_ON[theta][ii] = 0; 
			ActPY2_OF[theta][ii] = 0; 
			ActM_OF[theta][ii] = 0; 
			ActSB2_OF[theta][ii] = 0; 
			ActLB2_OF[theta][ii] = 0; 
			uPY2[theta][ii]=UPYres;
			uM[theta][ii]=UMres;
			uSB2[theta][ii]=USBres;
			uLB2[theta][ii]=ULBres;
		}
	}
	duration  = (int)(0.001/DT);
	for (t=0; t<PERIOD; ++t){    
		display();
		for (theta=0; theta<=N_assm-1; ++theta){
			for (i=0; i<=N_T; ++i){	
				if (theta==COLUMN && i==NEURON){  
					I_leak1  = duPY_leak1[theta][i]*(cmPY/DT);
					I_rec_1  = duPY_rec_1[theta][i]*(cmPY/DT);
					I_fed_1  = duPY_fed_1[theta][i]*(cmPY/DT);
					I_lat_1  = duPY_lat_1[theta][i]*(cmPY/DT);
					I_topdown  = duPY_topdown[theta][i]*(cmPY/DT);
					I_ext_1  = duPY_ext_1[theta][i]*(cmPY/DT);
					I_MG1    = duPY_MGB1[theta]*(cmPY/DT)*(int)(rand()/32768.0+inp_prob);
				}
				uPY1[theta][i] += difuPY1(theta,i); 
				sigPYT = sigmoidPY(uPY1[theta][i]);
				if (ActPY1_ON[theta][i]==0 && ActPY1_OF[theta][i]==0) vPY1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigPYT);
				if (vPY1[theta][i][t]==1.0){
					uPY1[theta][i]=UPYact;
					GlutPY_c1[theta][i]=Glut_c;
					ActPY1_ON[theta][i]=1;
					ActPY1_OF[theta][i]=1;
				}else{
					GlutPY_c1[theta][i]=0.0;
				}
				if (ActPY1_ON[theta][i] != 0){
					++ActPY1_ON[theta][i];
					uPY1[theta][i]=UPYact;
					GlutPY_c1[theta][i]=Glut_c;
				}
				if (ActPY1_OF[theta][i] != 0){
					++ActPY1_OF[theta][i];
				}
				if (ActPY1_ON[theta][i] > duration){
					uPY1[theta][i]=UPYres;
					ActPY1_ON[theta][i] = 0;
				}
				if (ActPY1_OF[theta][i] > 10*duration){
					ActPY1_OF[theta][i] = 0;
				}
				rPY1[theta][i] += difrPY1(theta,i); 
				rPY1d[theta][i][t] = rPY1[theta][i]; 
				if (t%INTV == 0){
				}

				uSB1[theta][i] += difuSB1(theta,i); 
				sigSBT = sigmoidSB(uSB1[theta][i]);
				if (ActSB1_ON[theta][i]==0 && ActSB1_OF[theta][i]==0) vSB1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigSBT);
				if (vSB1[theta][i][t]==1.0){
					uSB1[theta][i]=USBact;
					GABASB_c1[theta][i][t]=GABA_cS;
					ActSB1_ON[theta][i]=1;
					ActSB1_OF[theta][i]=1;
				}else{
					GABASB_c1[theta][i][t]=0.0;
				}
				if (ActSB1_ON[theta][i] != 0){
					++ActSB1_ON[theta][i];
					uSB1[theta][i]=USBact;
					GABASB_c1[theta][i][t]=GABA_cS;
				}
				if (ActSB1_OF[theta][i] != 0){
					++ActSB1_OF[theta][i];
				}
				if (ActSB1_ON[theta][i] > duration){
					uSB1[theta][i]=USBres;
					ActSB1_ON[theta][i] = 0;
				}
				if (ActSB1_OF[theta][i] > 10*duration){
					ActSB1_OF[theta][i] = 0;
				}	
				rSB1[theta][i] += difrSB1(theta,i); 
				sF[theta][i]   += difsF(theta,i);     
				if (t%INTV == 0){
				}
				uLB1[theta][i] += difuLB1(theta,i); 
				sigLBT = sigmoidLB(uLB1[theta][i]);
				if (ActLB1_ON[theta][i]==0 && ActLB1_OF[theta][i]==0) vLB1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigLBT);
				if (vLB1[theta][i][t]==1.0){
					uLB1[theta][i]=ULBact;
					GABALB_c1[theta][i][t]=GABA_cL1;
					ActLB1_ON[theta][i]=1;
					ActLB1_OF[theta][i]=1;
				}else{
					GABALB_c1[theta][i][t]=0.0;
				}
				if (ActLB1_ON[theta][i] != 0){
					++ActLB1_ON[theta][i];
					uLB1[theta][i]=ULBact;
					GABALB_c1[theta][i][t]=GABA_cL1;
				}
				if (ActLB1_OF[theta][i] != 0){
					++ActLB1_OF[theta][i];
				}				
				if (ActLB1_ON[theta][i] > duration){
					uLB1[theta][i]=ULBres;
					ActLB1_ON[theta][i] = 0;
				}
				if (ActLB1_OF[theta][i] > 10*duration){
					ActLB1_OF[theta][i] = 0;
				}
				rLB1[theta][i] += difrLB1(theta,i); 
				if (t%INTV == 0){
				}

				if (theta==COLUMN && i==NEURON){  
					I_leak2  = duPY_leak2[theta][i]*(cmPY/DT);
					I_rec_2  = duPY_rec_2[theta][i]*(cmPY/DT);
					I_fed_2  = duPY_fed_2[theta][i]*(cmPY/DT);
					I_lat_2 = duPY_lat_2[theta][i]*(cmPY/DT);
					I_ext_2  = duPY_ext_2[theta][i]*(cmPY/DT);
				}
				uPY2[theta][i] += difuPY2(theta,i); 
				sigPYT = sigmoidPY2(uPY2[theta][i]);
				if (ActPY2_ON[theta][i]==0 && ActPY2_OF[theta][i]==0) vPY2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigPYT);
				if (vPY2[theta][i][t]==1.0){
					uPY2[theta][i]=UPYact;
					GlutPY_c2[theta][i]=Glut_c;
					ActPY2_ON[theta][i]=1;
					ActPY2_OF[theta][i]=1;
				}else{
					GlutPY_c2[theta][i]=0.0;
				}
				if (ActPY2_ON[theta][i] != 0){
					++ActPY2_ON[theta][i];
					uPY2[theta][i]=UPYact;
					GlutPY_c2[theta][i]=Glut_c;
				}
				if (ActPY2_OF[theta][i] != 0){
					++ActPY2_OF[theta][i];
				}
				if (ActPY2_ON[theta][i] > duration){
					uPY2[theta][i]=UPYres;
					ActPY2_ON[theta][i] = 0;
				}
				if (ActPY2_OF[theta][i] > 10*duration){
					ActPY2_OF[theta][i] = 0;
				}
				rPY2[theta][i] += difrPY2(theta,i); 
				rPY2d[theta][i][t] = rPY2[theta][i]; 
				if (t%INTV == 0){
				}
				uM[theta][i] += difuM(theta,i); 
				sigMT = sigmoidM(uM[theta][i]);
				if (ActM_ON[theta][i]==0 && ActM_OF[theta][i]==0) vM[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigMT);
				if (vM[theta][i][t]==1.0){
					uM[theta][i]=UMact;
					GlutM[theta][i]=Glut_c;
					ActM_ON[theta][i]=1;
					ActM_OF[theta][i]=1;
				}else{
					GlutM[theta][i]=0.0;
				}
				if (ActM_ON[theta][i] != 0){
					++ActM_ON[theta][i];
					uM[theta][i]=UMact;
					GlutM[theta][i]=Glut_c;
				}
				if (ActM_OF[theta][i] != 0){
					++ActM_OF[theta][i];
				}
				if (ActM_ON[theta][i] > duration){
					uM[theta][i]=UMres;
					ActM_ON[theta][i] = 0;
				}
				if (ActM_OF[theta][i] > 10*duration){
					ActM_OF[theta][i] = 0;
				}
				rM[theta][i] += difrM(theta,i); 
				rMd[theta][i][t] = rM[theta][i]; 
				if (t%INTV == 0){
				}
				uSB2[theta][i] += difuSB2(theta,i); 
				sigSBT = sigmoidSB2(uSB2[theta][i]);
				if (ActSB2_ON[theta][i]==0 && ActSB2_OF[theta][i]==0) vSB2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigSBT);
				if (vSB2[theta][i][t]==1.0){
					uSB2[theta][i]=USBact;
					GABASB_c2[theta][i][t]=GABA_cS;
					ActSB2_ON[theta][i]=1;
					ActSB2_OF[theta][i]=1;
				}else{
					GABASB_c2[theta][i][t]=0.0;
				}
				if (ActSB2_ON[theta][i] != 0){
					++ActSB2_ON[theta][i];
					uSB2[theta][i]=USBact;
					GABASB_c2[theta][i][t]=GABA_cS;
				}
				if (ActSB2_OF[theta][i] != 0){
					++ActSB2_OF[theta][i];
				}
				if (ActSB2_ON[theta][i] > duration){
					uSB2[theta][i]=USBres;
					ActSB2_ON[theta][i] = 0;
				}
				if (ActSB2_OF[theta][i] > 10*duration){
					ActSB2_OF[theta][i] = 0;
				}	
				rSB2[theta][i] += difrSB2(theta,i); 
				if (t%INTV == 0){
				}
				uLB2[theta][i] += difuLB2(theta,i); 
				sigLBT = sigmoidLB2(uLB2[theta][i]);
				if (ActLB2_ON[theta][i]==0 && ActLB2_OF[theta][i]==0) vLB2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigLBT);
				if (vLB2[theta][i][t]==1.0){
					uLB2[theta][i]=ULBact;
					GABALB_c2[theta][i][t]=GABA_cL2;
					ActLB2_ON[theta][i]=1;
					ActLB2_OF[theta][i]=1;
				}else{
					GABALB_c2[theta][i][t]=0.0;
				}
				if (ActLB2_ON[theta][i] != 0){
					++ActLB2_ON[theta][i];
					uLB2[theta][i]=ULBact;
					GABALB_c2[theta][i][t]=GABA_cL2;
				}
				if (ActLB2_OF[theta][i] != 0){
					++ActLB2_OF[theta][i];
				}				
				if (ActLB2_ON[theta][i] > duration){
					uLB2[theta][i]=ULBres;
					ActLB2_ON[theta][i] = 0;
				}
				if (ActLB2_OF[theta][i] > 10*duration){
					ActLB2_OF[theta][i] = 0;
				}
				rLB2[theta][i] += difrLB2(theta,i); 
				if (t%INTV == 0){
				}

				if (t>=OUT){  
					if (theta==0 && i==0){									
						fprintf(uPY0_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM0_02,"%f\n",uM[theta][i]); 
						fprintf(uSB0_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB0_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB0_02,"%f\n",uLB2[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_01,"%f\n",-10.0);
						}else{
							fprintf(vPY0_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_02,"%f\n",-10.0);
						}else{
							fprintf(vPY0_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					} 
					if (theta==0 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_11,"%f\n",-10.0);
						}else{
							fprintf(vPY0_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_12,"%f\n",-10.0);
						}else{
							fprintf(vPY0_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_21,"%f\n",-10.0);
						}else{
							fprintf(vPY0_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_22,"%f\n",-10.0);
						}else{
							fprintf(vPY0_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_31,"%f\n",-10.0);
						}else{
							fprintf(vPY0_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_32,"%f\n",-10.0);
						}else{
							fprintf(vPY0_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_41,"%f\n",-10.0);
						}else{
							fprintf(vPY0_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_42,"%f\n",-10.0);
						}else{
							fprintf(vPY0_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==0 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_51,"%f\n",-10.0);
						}else{
							fprintf(vPY0_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_52,"%f\n",-10.0);
						}else{
							fprintf(vPY0_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_61,"%f\n",-10.0);
						}else{
							fprintf(vPY0_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_62,"%f\n",-10.0);
						}else{
							fprintf(vPY0_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_71,"%f\n",-10.0);
						}else{
							fprintf(vPY0_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_72,"%f\n",-10.0);
						}else{
							fprintf(vPY0_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_81,"%f\n",-10.0);
						}else{
							fprintf(vPY0_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_82,"%f\n",-10.0);
						}else{
							fprintf(vPY0_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==0 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY0_91,"%f\n",-10.0);
						}else{
							fprintf(vPY0_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY0_92,"%f\n",-10.0);
						}else{
							fprintf(vPY0_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==0){									
						fprintf(uPY1_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM1_02,"%f\n",uM[theta][i]); 
						fprintf(uSB1_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB1_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB1_02,"%f\n",uLB2[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_01,"%f\n",-10.0);
						}else{
							fprintf(vPY1_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_02,"%f\n",-10.0);
						}else{
							fprintf(vPY1_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
				}
					if (theta==1 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_11,"%f\n",-10.0);
						}else{
							fprintf(vPY1_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_12,"%f\n",-10.0);
						}else{
							fprintf(vPY1_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_21,"%f\n",-10.0);
						}else{
							fprintf(vPY1_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_22,"%f\n",-10.0);
						}else{
							fprintf(vPY1_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_31,"%f\n",-10.0);
						}else{
							fprintf(vPY1_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_32,"%f\n",-10.0);
						}else{
							fprintf(vPY1_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
				}
					if (theta==1 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_41,"%f\n",-10.0);
						}else{
							fprintf(vPY1_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_42,"%f\n",-10.0);
						}else{
							fprintf(vPY1_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_51,"%f\n",-10.0);
						}else{
							fprintf(vPY1_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_52,"%f\n",-10.0);
						}else{
							fprintf(vPY1_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_61,"%f\n",-10.0);
						}else{
							fprintf(vPY1_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_62,"%f\n",-10.0);
						}else{
							fprintf(vPY1_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_71,"%f\n",-10.0);
						}else{
							fprintf(vPY1_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_72,"%f\n",-10.0);
						}else{
							fprintf(vPY1_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_81,"%f\n",-10.0);
						}else{
							fprintf(vPY1_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_82,"%f\n",-10.0);
						}else{
							fprintf(vPY1_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==1 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY1_91,"%f\n",-10.0);
						}else{
							fprintf(vPY1_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY1_92,"%f\n",-10.0);
						}else{
							fprintf(vPY1_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==0){									
						fprintf(uPY2_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM2_02,"%f\n",uM[theta][i]); 
						fprintf(uSB2_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB2_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB2_02,"%f\n",uLB2[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_01,"%f\n",-10.0);
						}else{
							fprintf(vPY2_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_02,"%f\n",-10.0);
						}else{
							fprintf(vPY2_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_11,"%f\n",-10.0);
						}else{
							fprintf(vPY2_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_12,"%f\n",-10.0);
						}else{
							fprintf(vPY2_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_21,"%f\n",-10.0);
						}else{
							fprintf(vPY2_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_22,"%f\n",-10.0);
						}else{
							fprintf(vPY2_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_31,"%f\n",-10.0);
						}else{
							fprintf(vPY2_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_32,"%f\n",-10.0);
						}else{
							fprintf(vPY2_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_41,"%f\n",-10.0);
						}else{
							fprintf(vPY2_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_42,"%f\n",-10.0);
						}else{
							fprintf(vPY2_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_51,"%f\n",-10.0);
						}else{
							fprintf(vPY2_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_52,"%f\n",-10.0);
						}else{
							fprintf(vPY2_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_61,"%f\n",-10.0);
						}else{
							fprintf(vPY2_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_62,"%f\n",-10.0);
						}else{
							fprintf(vPY2_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_71,"%f\n",-10.0);
						}else{
							fprintf(vPY2_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_72,"%f\n",-10.0);
						}else{
							fprintf(vPY2_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_81,"%f\n",-10.0);
						}else{
							fprintf(vPY2_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_82,"%f\n",-10.0);
						}else{
							fprintf(vPY2_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==2 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY2_91,"%f\n",-10.0);
						}else{
							fprintf(vPY2_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY2_92,"%f\n",-10.0);
						}else{
							fprintf(vPY2_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==3 && i==0){									
						fprintf(uPY3_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM3_02,"%f\n",uM[theta][i]); 
						fprintf(uSB3_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB3_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB3_02,"%f\n",uLB2[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_01,"%f\n",-10.0);
						}else{
							fprintf(vPY3_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_02,"%f\n",-10.0);
						}else{
							fprintf(vPY3_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_02,"%f\n",-10.0);
						}else{
							fprintf(vM3_02,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_11,"%f\n",-10.0);
						}else{
							fprintf(vPY3_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_12,"%f\n",-10.0);
						}else{
							fprintf(vPY3_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_12,"%f\n",-10.0);
						}else{
							fprintf(vM3_12,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_21,"%f\n",-10.0);
						}else{
							fprintf(vPY3_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_22,"%f\n",-10.0);
						}else{
							fprintf(vPY3_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_22,"%f\n",-10.0);
						}else{
							fprintf(vM3_22,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_31,"%f\n",-10.0);
						}else{
							fprintf(vPY3_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_32,"%f\n",-10.0);
						}else{
							fprintf(vPY3_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_32,"%f\n",-10.0);
						}else{
							fprintf(vM3_32,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_41,"%f\n",-10.0);
						}else{
							fprintf(vPY3_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_42,"%f\n",-10.0);
						}else{
							fprintf(vPY3_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_42,"%f\n",-10.0);
						}else{
							fprintf(vM3_42,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_51,"%f\n",-10.0);
						}else{
							fprintf(vPY3_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_52,"%f\n",-10.0);
						}else{
							fprintf(vPY3_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_52,"%f\n",-10.0);
						}else{
							fprintf(vM3_52,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_61,"%f\n",-10.0);
						}else{
							fprintf(vPY3_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_62,"%f\n",-10.0);
						}else{
							fprintf(vPY3_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_62,"%f\n",-10.0);
						}else{
							fprintf(vM3_62,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_71,"%f\n",-10.0);
						}else{
							fprintf(vPY3_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_72,"%f\n",-10.0);
						}else{
							fprintf(vPY3_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_72,"%f\n",-10.0);
						}else{
							fprintf(vM3_72,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_81,"%f\n",-10.0);
						}else{
							fprintf(vPY3_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_82,"%f\n",-10.0);
						}else{
							fprintf(vPY3_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_82,"%f\n",-10.0);
						}else{
							fprintf(vM3_82,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==3 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY3_91,"%f\n",-10.0);
						}else{
							fprintf(vPY3_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY3_92,"%f\n",-10.0);
						}else{
							fprintf(vPY3_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vM[theta][i][t]==0.0){
							fprintf(vM3_92,"%f\n",-10.0);
						}else{
							fprintf(vM3_92,"%f\n",vM[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==4 && i==0){									
						fprintf(uPY4_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM4_02,"%f\n",uM[theta][i]); 
						fprintf(uSB4_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB4_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB4_02,"%f\n",uLB2[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_01,"%f\n",-10.0);
						}else{
							fprintf(vPY4_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_02,"%f\n",-10.0);
						}else{
							fprintf(vPY4_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_11,"%f\n",-10.0);
						}else{
							fprintf(vPY4_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_12,"%f\n",-10.0);
						}else{
							fprintf(vPY4_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_21,"%f\n",-10.0);
						}else{
							fprintf(vPY4_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_22,"%f\n",-10.0);
						}else{
							fprintf(vPY4_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_31,"%f\n",-10.0);
						}else{
							fprintf(vPY4_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_32,"%f\n",-10.0);
						}else{
							fprintf(vPY4_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_41,"%f\n",-10.0);
						}else{
							fprintf(vPY4_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_42,"%f\n",-10.0);
						}else{
							fprintf(vPY4_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_51,"%f\n",-10.0);
						}else{
							fprintf(vPY4_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_52,"%f\n",-10.0);
						}else{
							fprintf(vPY4_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_61,"%f\n",-10.0);
						}else{
							fprintf(vPY4_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_62,"%f\n",-10.0);
						}else{
							fprintf(vPY4_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_71,"%f\n",-10.0);
						}else{
							fprintf(vPY4_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_72,"%f\n",-10.0);
						}else{
							fprintf(vPY4_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_81,"%f\n",-10.0);
						}else{
							fprintf(vPY4_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_82,"%f\n",-10.0);
						}else{
							fprintf(vPY4_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==4 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY4_91,"%f\n",-10.0);
						}else{
							fprintf(vPY4_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY4_92,"%f\n",-10.0);
						}else{
							fprintf(vPY4_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}					
					}
					if (theta==5 && i==0){									
						fprintf(uPY5_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM5_02,"%f\n",uM[theta][i]); 
						fprintf(uSB5_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB5_01,"%f\n",uLB1[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_01,"%f\n",-10.0);
						}else{
							fprintf(vPY5_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_02,"%f\n",-10.0);
						}else{
							fprintf(vPY5_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}						
					}
					if (theta==5 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_11,"%f\n",-10.0);
						}else{
							fprintf(vPY5_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_12,"%f\n",-10.0);
						}else{
							fprintf(vPY5_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_21,"%f\n",-10.0);
						}else{
							fprintf(vPY5_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_22,"%f\n",-10.0);
						}else{
							fprintf(vPY5_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_31,"%f\n",-10.0);
						}else{
							fprintf(vPY5_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_32,"%f\n",-10.0);
						}else{
							fprintf(vPY5_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_41,"%f\n",-10.0);
						}else{
							fprintf(vPY5_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_42,"%f\n",-10.0);
						}else{
							fprintf(vPY5_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_51,"%f\n",-10.0);
						}else{
							fprintf(vPY5_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_52,"%f\n",-10.0);
						}else{
							fprintf(vPY5_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_61,"%f\n",-10.0);
						}else{
							fprintf(vPY5_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_62,"%f\n",-10.0);
						}else{
							fprintf(vPY5_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_71,"%f\n",-10.0);
						}else{
							fprintf(vPY5_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_72,"%f\n",-10.0);
						}else{
							fprintf(vPY5_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_81,"%f\n",-10.0);
						}else{
							fprintf(vPY5_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_82,"%f\n",-10.0);
						}else{
							fprintf(vPY5_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==5 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY5_91,"%f\n",-10.0);
						}else{
							fprintf(vPY5_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY5_92,"%f\n",-10.0);
						}else{
							fprintf(vPY5_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==0){									
						fprintf(uPY6_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM6_02,"%f\n",uM[theta][i]); 
						fprintf(uSB6_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB6_01,"%f\n",uLB1[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_01,"%f\n",-10.0);
						}else{
							fprintf(vPY6_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_02,"%f\n",-10.0);
						}else{
							fprintf(vPY6_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_11,"%f\n",-10.0);
						}else{
							fprintf(vPY6_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_12,"%f\n",-10.0);
						}else{
							fprintf(vPY6_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_21,"%f\n",-10.0);
						}else{
							fprintf(vPY6_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_22,"%f\n",-10.0);
						}else{
							fprintf(vPY6_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_31,"%f\n",-10.0);
						}else{
							fprintf(vPY6_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_32,"%f\n",-10.0);
						}else{
							fprintf(vPY6_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_41,"%f\n",-10.0);
						}else{
							fprintf(vPY6_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_42,"%f\n",-10.0);
						}else{
							fprintf(vPY6_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_51,"%f\n",-10.0);
						}else{
							fprintf(vPY6_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_52,"%f\n",-10.0);
						}else{
							fprintf(vPY6_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_61,"%f\n",-10.0);
						}else{
							fprintf(vPY6_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_62,"%f\n",-10.0);
						}else{
							fprintf(vPY6_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_71,"%f\n",-10.0);
						}else{
							fprintf(vPY6_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_72,"%f\n",-10.0);
						}else{
							fprintf(vPY6_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_81,"%f\n",-10.0);
						}else{
							fprintf(vPY6_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_82,"%f\n",-10.0);
						}else{
							fprintf(vPY6_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==6 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY6_91,"%f\n",-10.0);
						}else{
							fprintf(vPY6_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY6_92,"%f\n",-10.0);
						}else{
							fprintf(vPY6_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==0){									
						fprintf(uPY7_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_02,"%f\n",uPY2[theta][i]); 
						fprintf(uM7_02,"%f\n",uM[theta][i]); 
						fprintf(uSB7_01,"%f\n",uSB1[theta][i]); 
						fprintf(uLB7_01,"%f\n",uLB1[theta][i]); 
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_01,"%f\n",-10.0);
						}else{
							fprintf(vPY7_01,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_02,"%f\n",-10.0);
						}else{
							fprintf(vPY7_02,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==1){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_11,"%f\n",-10.0);
						}else{
							fprintf(vPY7_11,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_12,"%f\n",-10.0);
						}else{
							fprintf(vPY7_12,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==2){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_21,"%f\n",-10.0);
						}else{
							fprintf(vPY7_21,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_22,"%f\n",-10.0);
						}else{
							fprintf(vPY7_22,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==3){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_31,"%f\n",-10.0);
						}else{
							fprintf(vPY7_31,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_32,"%f\n",-10.0);
						}else{
							fprintf(vPY7_32,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==4){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_41,"%f\n",-10.0);
						}else{
							fprintf(vPY7_41,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_42,"%f\n",-10.0);
						}else{
							fprintf(vPY7_42,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==5){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_51,"%f\n",-10.0);
						}else{
							fprintf(vPY7_51,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_52,"%f\n",-10.0);
						}else{
							fprintf(vPY7_52,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==6){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_61,"%f\n",-10.0);
						}else{
							fprintf(vPY7_61,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_62,"%f\n",-10.0);
						}else{
							fprintf(vPY7_62,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==7){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_71,"%f\n",-10.0);
						}else{
							fprintf(vPY7_71,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_72,"%f\n",-10.0);
						}else{
							fprintf(vPY7_72,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==8){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_81,"%f\n",-10.0);
						}else{
							fprintf(vPY7_81,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_82,"%f\n",-10.0);
						}else{
							fprintf(vPY7_82,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==7 && i==9){									
						if (vPY1[theta][i][t]==0.0){
							fprintf(vPY7_91,"%f\n",-10.0);
						}else{
							fprintf(vPY7_91,"%f\n",vPY1[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
						if (vPY2[theta][i][t]==0.0){
							fprintf(vPY7_92,"%f\n",-10.0);
						}else{
							fprintf(vPY7_92,"%f\n",vPY2[theta][i][t]+OFSBCA*theta+i*OFSWCA); 
						}
					}
					if (theta==COLUMN && i==NEURON){  
						fprintf(leak1,"%e\n",I_leak1); 
						fprintf(rec1,"%e\n",I_rec_1); 
						fprintf(fed1,"%e\n",I_fed_1); 
						fprintf(lat1,"%e\n",I_lat_1); 
						fprintf(ext1,"%e\n",I_ext_1); 
						fprintf(topdown,"%e\n",I_topdown); 
						fprintf(MG1,"%e\n",I_MG1);
					}
					if (theta==COLUMN && i==NEURON){  
						fprintf(extSB1,"%e\n",I_ext_SB1); 
					}
					if (theta==COLUMN && i==NEURON){  
						fprintf(extLB1,"%e\n",I_ext_LB1); 
					}
				} 
			}
		} 
		for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
			rEXT[theta][i] += difrEXT(theta,i); 
			}
		}
		for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
			rEXT2[theta][i] += difrEXT2(theta,i); 
			}
		}
	} 
    
	fclose(uPY0_01);
	fclose(uPY1_01);
	fclose(uPY2_01);
	fclose(uPY3_01);
	fclose(uPY4_01);
	fclose(uPY5_01);
	fclose(uPY6_01);
	fclose(uPY7_01);

	fclose(vPY0_01);
	fclose(vPY0_11);
	fclose(vPY0_21);
	fclose(vPY0_31);
	fclose(vPY0_41);
	fclose(vPY0_51);
	fclose(vPY0_61);
	fclose(vPY0_71);
	fclose(vPY0_81);
	fclose(vPY0_91);
	fclose(vPY1_01);
	fclose(vPY1_11);
	fclose(vPY1_21);
	fclose(vPY1_31);
	fclose(vPY1_41);
	fclose(vPY1_51);
	fclose(vPY1_61);
	fclose(vPY1_71);
	fclose(vPY1_81);
	fclose(vPY1_91);
	fclose(vPY2_01);
	fclose(vPY2_11);
	fclose(vPY2_21);
	fclose(vPY2_31);
	fclose(vPY2_41);
	fclose(vPY2_51);
	fclose(vPY2_61);
	fclose(vPY2_71);
	fclose(vPY2_81);
	fclose(vPY2_91);
	fclose(vPY3_01);
	fclose(vPY3_11);
	fclose(vPY3_21);
	fclose(vPY3_31);
	fclose(vPY3_41);
	fclose(vPY3_51);
	fclose(vPY3_61);
	fclose(vPY3_71);
	fclose(vPY3_81);
	fclose(vPY3_91);
	fclose(vPY4_01);
	fclose(vPY4_11);
	fclose(vPY4_21);
	fclose(vPY4_31);
	fclose(vPY4_41);
	fclose(vPY4_51);
	fclose(vPY4_61);
	fclose(vPY4_71);
	fclose(vPY4_81);
	fclose(vPY4_91);
	fclose(vPY5_01);
	fclose(vPY5_11);
	fclose(vPY5_21);
	fclose(vPY5_31);
	fclose(vPY5_41);
	fclose(vPY5_51);
	fclose(vPY5_61);
	fclose(vPY5_71);
	fclose(vPY5_81);
	fclose(vPY5_91);
	fclose(vPY6_01);
	fclose(vPY6_11);
	fclose(vPY6_21);
	fclose(vPY6_31);
	fclose(vPY6_41);
	fclose(vPY6_51);
	fclose(vPY6_61);
	fclose(vPY6_71);
	fclose(vPY6_81);
	fclose(vPY6_91);
	fclose(vPY7_01);
	fclose(vPY7_11);
	fclose(vPY7_21);
	fclose(vPY7_31);
	fclose(vPY7_41);
	fclose(vPY7_51);
	fclose(vPY7_61);
	fclose(vPY7_71);
	fclose(vPY7_81);
	fclose(vPY7_91);

	fclose(leak1);
	fclose(rec1);
	fclose(fed1);
	fclose(lat1);
	fclose(ext1);
	fclose(topdown);
	fclose(MG1);

	fclose(extSB1);
	fclose(extLB1);

	fclose(uSB0_01);
	fclose(uSB1_01);
	fclose(uSB2_01);
	fclose(uSB3_01);
	fclose(uSB4_01);
	fclose(uSB5_01);
	fclose(uSB6_01);
	fclose(uSB7_01);

	fclose(uLB0_01);
	fclose(uLB1_01);
	fclose(uLB2_01);
	fclose(uLB3_01);
	fclose(uLB4_01);
	fclose(uLB5_01);
	fclose(uLB6_01);
	fclose(uLB7_01);
	
	fclose(uPY0_02);
	fclose(uPY1_02);
	fclose(uPY2_02);
	fclose(uPY3_02);
	fclose(uPY4_02);
	fclose(uPY5_02);
	fclose(uPY6_02);
	fclose(uPY7_02);
	
	fclose(uLB0_02);
	fclose(uLB1_02);
	fclose(uLB2_02);
	fclose(uLB3_02);
	fclose(uLB4_02);
	fclose(uLB5_02);
	fclose(uLB6_02);
	fclose(uLB7_02);

	fclose(uM0_02);
	fclose(uM1_02);
	fclose(uM2_02);
	fclose(uM3_02);
	fclose(uM4_02);
	fclose(uM5_02);
	fclose(uM6_02);
	fclose(uM7_02);

	fclose(vM3_02);
	fclose(vM3_12);
	fclose(vM3_22);
	fclose(vM3_32);
	fclose(vM3_42);
	fclose(vM3_52);
	fclose(vM3_62);
	fclose(vM3_72);
	fclose(vM3_82);
	fclose(vM3_92);

	printf("\a");
    printf("\a");
    printf("\a");
} 

double difuPY1(int thetaa, int ii){
	int jj,thetdash;
	double duPY1=0.0,rec_exT1=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double top_down = 0.0;
	double alpha_T=1.0; 
	double tau_T0=0.1;   

	duPY_leak1[thetaa][ii] = -(DT*gmPY/cmPY)*(uPY1[thetaa][ii]-UPYres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT1 += w_rec_1[thetaa][ii][thetaa][jj]*rPY1[thetaa][jj]; 
	}
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT1 += w_rec_1[thetaa][ii][N_assm][jj]*rPY1[N_assm][jj]; 
	}
   		for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
			for (jj=0; jj<=N_T; ++jj){
				if (thetaa!=thetdash) rec_exT1 += w_rec_1[thetaa][ii][thetdash][jj]*rPY1[thetdash][jj];
			}
		}
 	duPY_rec_1[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY1[thetaa][ii]-u_AMPA)*rec_exT1; 
	duPY_fed_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*w_fed_1[thetaa][ii]*rSB1[thetaa][ii]; 
	for (jj=0; jj<N_T; ++jj){ // 
		lat_ihTT += w_lat_1[thetaa][ii][jj]*rLB1[thetaa][jj];
	}
	duPY_lat_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*lat_ihTT; 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			top_down += w_v1P_v2[thetaa][ii][thetdash][jj]*rPY2d[thetdash][jj][t-delay];
		}
	}
	duPY_topdown[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY1[thetaa][ii]-u_AMPA)*top_down; 
	duPY_ext_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*delta_T*rEXT[thetaa][ii]; 
	duPY1 = duPY_leak1[thetaa][ii]                  
	      + duPY_rec_1[thetaa][ii] 
		  + duPY_fed_1[thetaa][ii] 
		  + duPY_lat_1[thetaa][ii]  
	      + duPY_topdown[thetaa][ii]
		  + duPY_ext_1[thetaa][ii]; 
	return(duPY1);
}
double difuPY2(int thetaa, int ii){
	int jj,thetdash;
	double duPY2=0.0,rec_exT2=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double alpha_T=1.0; // 1.0 
	double tau_T0=4; 

	duPY_leak2[thetaa][ii] = -(DT*gmPY/cmPY)*(uPY2[thetaa][ii]-UPYres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT2 += w_rec_2[thetaa][ii][thetaa][jj]*rPY2[thetaa][jj]; 
	}
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			rec_exT2 += w_v2_v1[thetaa][ii][thetdash][jj]*rPY1d[thetdash][jj][t-delay];
		}
	}
	duPY_rec_2[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY2[thetaa][ii]-u_AMPA)*rec_exT2; 

	for (jj=0; jj<N_T; ++jj){ 
		lat_ihTT += w_lat_2[thetaa][ii][jj]*rLB2[thetaa][jj];
	}
	duPY_lat_2[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY2[thetaa][ii]-u_GABA)*lat_ihTT; 
	duPY_ext_2[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY2[thetaa][ii]-u_GABA)*delta_T2*rEXT2[thetaa][ii]; 
	if(t>=onset_0 && t<onset_0+period_0 && thetaa<N_assm){    
		I_MGB1[thetaa]  = int_inp0_0*exp(-(pow((thetaa-theta_inp0)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_1*exp(-(pow((thetaa-theta_inp1)/tau_T0,2.0)));
		I_MGB1[thetaa] += int_inp0_2*exp(-(pow((thetaa-theta_inp2)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_3*exp(-(pow((thetaa-theta_inp3)/tau_T0,2.0)));
		I_MGB1[thetaa] += int_inp0_4*exp(-(pow((thetaa-theta_inp4)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_5*exp(-(pow((thetaa-theta_inp5)/tau_T0,2.0))); 
		I_MGB1[thetaa] += int_inp0_6*exp(-(pow((thetaa-theta_inp6)/tau_T0,2.0)));
		I_MGB1[thetaa] += int_inp0_7*exp(-(pow((thetaa-theta_inp7)/tau_T0,2.0)));
	}else{
		I_MGB1[thetaa] = 0.0;
	}
	duPY_MGB1[thetaa] = (DT/cmPY)*(I_MGB1[thetaa]+I_MGBN1[thetaa]); 
	if (ii>N_T) duPY_MGB1[thetaa] = 0.0; 
	duPY2 = duPY_leak2[thetaa][ii]                  
		  + duPY_rec_2[thetaa][ii] 
		  + duPY_fed_2[thetaa][ii] 
		  + duPY_lat_2[thetaa][ii]  
		  + duPY_ext_2[thetaa][ii]
		  + duPY_MGB1[thetaa]*(int)(rand()/32768.0+inp_prob); 
	return(duPY2);
}

double difuM(int thetaa, int ii){
	int jj;
	double duM=0.0,rec_exT2=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double alpha_T=1.0; // 1.0 
	double tau_T0=0.1,tau_T1=0.1,tau_T2=0.1,tau_T3=0.1;    

	duM_leak[thetaa][ii] = -(DT*gmM/cmM)*(uM[thetaa][ii]-UMres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT2 += wMM*rM[thetaa][jj]; 
	}
	for (jj=0; jj<=N_T; ++jj){
		rec_exT2 += w_M_P*rPY1d[thetaa][jj][t-delay];
	}
	duM_rec[thetaa][ii] = -(DT/cmM)*gAMPA*(uM[thetaa][ii]-u_AMPA)*rec_exT2; 
	duM = duM_leak[thetaa][ii]                  
		 + duM_rec[thetaa][ii] 
		 + duM_fed[thetaa][ii] 
		 + duM_lat[thetaa][ii]  
		 + duM_ext[thetaa][ii]; 
	return(duM);
}

double difuSB1(int thetaa, int ii){
	double duSB1; 
	int thetdash,jj;
	double top_down = 0.0;

	duSB_leak1[thetaa][ii] = -(DT*gmSB/cmSB)*(uSB1[thetaa][ii]-USBres); 
	if (thetaa<N_assm){
		duSB_1[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)*(wSB_PY1[thetaa][ii]*rPY1[thetaa][ii]); 
	}
	if (thetaa==N_assm){
		duSB_1[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)*(wSB_PY1[thetaa][ii]*rPY1[thetaa][ii]); 
	}
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			top_down += w_v1Ia_v2[thetaa][ii][thetdash][jj]*rPY2[thetdash][jj];
		}
	}
	duSB_topdown[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)*top_down; 
	duSB1 = duSB_leak1[thetaa][ii]                  
		  + duSB_1[thetaa][ii]
		  + duSB_topdown[thetaa][ii];
	return(duSB1);
}

double difuSB2(int thetaa, int ii){
	double duSB2; 
	duSB_leak2[thetaa][ii] = -(DT*gmSB/cmSB)*(uSB2[thetaa][ii]-USBres); 
	if (thetaa<N_assm){
		duSB_2[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB2[thetaa][ii]-u_AMPA)
						   *(wSB_PY1[thetaa][ii]*rPY2[thetaa][ii]); 
	}
	if (thetaa==N_assm){
		duSB_2[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB2[thetaa][ii]-u_AMPA)*(wSB_PY1[thetaa][ii]*rPY2[thetaa][ii]); 
	}
	duSB2 = duSB_leak2[thetaa][ii]                  
		  + duSB_2[thetaa][ii];
	return(duSB2);
}

double difuLB1(int thetaa, int ii){
	int jj,thetdash;
	double duLBB,lat_ih=0.0;
	double top_down = 0.0;
	
	duLB_leak1[thetaa][ii] = -(DT*gmLB/cmLB)*(uLB1[thetaa][ii]-ULBres); 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){ 
		for (jj=0; jj<=N_T; ++jj){
			lat_ih += (wLB_PY1[thetaa][thetdash][ii]*rPY1[thetdash][jj] + wLB_PY1[thetaa][N_assm][ii]*rPY1[N_assm][jj]); 
		}
	}
	duLB1[thetaa][ii] = -(DT/cmLB)*gAMPA*(uLB1[thetaa][ii]-u_AMPA)*lat_ih; 
	duLBB = duLB_leak1[thetaa][ii]                  
		  + duLB1[thetaa][ii]   
		  + duLB1_topdown[thetaa][ii]
		  + duLB_ext1;
	return(duLBB);
}
double difuLB2(int thetaa, int ii){
	int jj,thetdash;
	double duLBB,lat_ih=0.0,bottom_up=0.0;
	
	duLB_leak2[thetaa][ii] = -(DT*gmLB/cmLB)*(uLB2[thetaa][ii]-ULBres); 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){ 
		for (jj=0; jj<=N_T; ++jj){
			lat_ih += (wLB_PY2[thetaa][thetdash][ii]*rPY2[thetdash][jj] + wLB_PY2[thetaa][N_assm][ii]*rPY2[N_assm][jj]); 
		}
	}
	duLB2[thetaa][ii] = -(DT/cmLB)*gAMPA*(uLB2[thetaa][ii]-u_AMPA)*lat_ih; 
	duLBB = duLB_leak2[thetaa][ii]                  
		  + duLB2[thetaa][ii]   
		  + duLB2_bottom_up[thetaa][ii]   
		  + duLB_ext2;
	return(duLBB);
}

double difrPY1(int thetaa, int ii){
	double drPYT;
	drPYT = DT*(alph_AMPA*GlutPY_c1[thetaa][ii]*(1.0-rPY1[thetaa][ii])-beta_AMPA*rPY1[thetaa][ii]); 
	return(drPYT);
}

double difrPY2(int thetaa, int ii){
	double drPYT;
	drPYT = DT*(alph_AMPA*GlutPY_c2[thetaa][ii]*(1.0-rPY2[thetaa][ii])-beta_AMPA*rPY2[thetaa][ii]); 
	return(drPYT);
}

double difrM(int thetaa, int ii){
	double drMT;
	drMT = DT*(alph_AMPA*GlutM[thetaa][ii]*(1.0-rM[thetaa][ii])-beta_AMPA*rM[thetaa][ii]); 
	return(drMT);
}

double difrSB1(int thetaa, int ii){
	double drSBT;
	drSBT = DT*(K1*GABASB_c1[thetaa][ii][t]*(1.0-rSB1[thetaa][ii])-K2*rSB1[thetaa][ii]); 
	return(drSBT);
}

double difsF(int thetaa, int ii){
	double dsF;
	dsF = DT*(K3*rSB1[thetaa][ii]-K4*sF[thetaa][ii]); 
	return(dsF);
}

double difrSB2(int thetaa, int ii){
	double drSBT;
	drSBT = DT*(alph_GABA*GABASB_c2[thetaa][ii][t]*(1.0-rSB2[thetaa][ii])-beta_GABA*rSB2[thetaa][ii]); 
	return(drSBT);

}double difrLB1(int thetaa, int ii){
	double drLB;
	drLB = DT*(alph_GABA*GABALB_c1[thetaa][ii][t]*(1.0-rLB1[thetaa][ii])-beta_GABA*rLB1[thetaa][ii]); 
	return(drLB);
}

double difrLB2(int thetaa, int ii){
	double drLB;
	drLB = DT*(alph_GABA*GABALB_c2[thetaa][ii][t]*(1.0-rLB2[thetaa][ii])-beta_GABA*rLB2[thetaa][ii]); 
	return(drLB);
}

double difrEXT(int thetaa, int ii){
	double drEXT;
	drEXT = DT*(alph_GABA*GABA_ext[thetaa][ii]*(1.0-rEXT[thetaa][ii])-beta_GABA*rEXT[thetaa][ii]); 
	return(drEXT);
}

double difrEXT2(int thetaa, int ii){
	double drEXT2;
	drEXT2 = DT*(alph_GABA*GABA_V2*(1.0-rEXT2[thetaa][ii])-beta_GABA*rEXT2[thetaa][ii]); 
	return(drEXT2);
}

void dfsGABA_ext(void){
	int theta,i;

	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				I_GABA[theta][i] = 0.0;
			}
	}
	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				I_GABA[theta][i] += -gamma*(GABA_ext[theta][i]-GABA_c0)*DT + m_G*(uG[theta][i]-uG_trn)*(GABAamb_max-GABA_ext[theta][i])*(GABA_ext[theta][i]-GABAamb_min)*DT; 
				GABA_extw[theta][i] = GABA_ext[theta][i] + I_GABA[theta][i];
			}
	}
	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				GABA_ext[theta][i] = GABA_extw[theta][i];
			}
	}
}

double sigmoidPY(double uPYY){
	return(1/(1+exp(-steep_PY*(uPYY-thres_PY))));
}

double sigmoidPY2(double uPYY){
	return(1/(1+exp(-steep_PY2*(uPYY-thres_PY2))));
}

double sigmoidM(double uMM){
	return(1/(1+exp(-steep_M*(uMM-thres_M))));
}

double sigmoidSB(double uSBB){
	return(1/(1+exp(-steep_SB*(uSBB-thres_SB))));
}
double sigmoidSB2(double uSBB){
	return(1/(1+exp(-steep_SB2*(uSBB-thres_SB2))));
}

double sigmoidLB(double uLBB){
	return(1/(1+exp(-steep_LB*(uLBB-thres_LB))));
}
double sigmoidLB2(double uLBB){
	return(1/(1+exp(-steep_LB2*(uLBB-thres_LB2))));
}

double rand01(long int *ix){
   double x;
   long int of;
   *ix=(*ix)*48828125;
   if(*ix<0){
	   of=2147483647;
       *ix=(*ix+of)+1;
   }
   x=(double)(*ix)*0.4656613e-9;
   return(x);
 }

void init(void){
	int theta,i,thetaa,ii,thetdash,jj;
	double tau_lat=1000.0;   // no use 
	double w_lat_excit=0.0;  // no use 
    double alph_w_L_TV1=1.6; // V1(P)-to-LB weight    
    double alph_w_L_TV2=1.2; // V2(P)-to-LB weight
	double w_v2_1=4.0;       // V2(P)-to-V2(Ib) weight
	double w_v1Ia_2=0.0;     // V2(P)-to-V1(Ia) weight
	double w_v1Ib_2=0.0;     // V2(P)-to-V1(Ib) weight
	double w_v1P_2=10.0;     // V2(P)-to-V1(P) weight
	double tauV2_1=1000;     // diffusive feedback

	srand(17);
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_rec_1[N_assm][ii][thetdash][jj] = 0.0; // 
						w_rec_2[N_assm][ii][thetdash][jj] = 0.0; // 
					}
				}
			}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v2_v1[thetaa][ii][thetdash][jj] = 0.0; // 
					}
				}
			}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v1Ia_v2[thetaa][ii][thetdash][jj] = 0.0; // 
					}
				}
			}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v1G_v2[thetaa][ii][thetdash][jj] = 0.0; // 
					}
				}
			}
	}
	for (thetaa=4; thetaa<=N_assm-1; ++thetaa){
		for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=3; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						if (abs(thetaa-thetdash) == 4) 
						{
							w_rec_1[thetaa][ii][thetdash][jj] = w_lat_excit;
						}
					}
				}
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				if (ii!=jj)	w_rec_1[thetaa][ii][thetaa][jj]  = wLPPV1; 
				if (ii!=jj)	w_rec_10[thetaa][ii][thetaa][jj] = wLPPV1;   
				if (ii!=jj)	w_rec_2[thetaa][ii][thetaa][jj]  = wLPPV2; 
			}
		}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			w_fed_1[thetaa][ii] = 0.0; 
			w_fed_2[thetaa][ii] = 0.0;  
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				w_lat_1[thetaa][ii][jj] = 6.0;    
				w_lat_2[thetaa][ii][jj] = 1.0;  		
			}
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			wSB_PY1[thetaa][ii] = 0.0;  
			wSB_PY2[thetaa][ii] = 0.0;   
		}
	}
	for (ii=0; ii<=N_T; ++ii){
		wLB_PY1[0][5][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[0][6][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[0][7][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[0][1][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[0][2][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[0][3][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[0][4][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[1][6][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[1][7][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[1][0][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[1][2][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[1][3][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[1][4][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0)));  
		wLB_PY1[1][5][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));   

		wLB_PY1[2][7][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[2][0][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[2][1][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[2][3][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[2][4][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[2][5][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0)));  
		wLB_PY1[2][6][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));  

		wLB_PY1[3][0][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[3][1][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[3][2][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[3][4][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0)));  
		wLB_PY1[3][5][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[3][6][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[3][7][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));  

		wLB_PY1[4][1][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[4][2][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[4][3][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[4][5][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[4][6][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[4][7][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[4][0][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[5][2][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[5][3][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[5][4][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[5][6][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[5][7][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[5][0][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[5][1][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[6][3][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[6][4][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[6][5][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[6][7][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[6][0][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[6][1][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[6][2][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[7][4][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[7][5][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[7][6][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[7][0][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[7][1][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[7][2][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0)));  
		wLB_PY1[7][3][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));  

		wLB_PY2[0][5][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[0][6][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[0][7][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[0][1][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[0][2][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[0][3][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[0][4][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[1][6][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[1][7][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[1][0][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[1][2][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[1][3][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[1][4][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[1][5][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[2][7][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[2][0][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[2][1][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[2][3][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[2][4][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[2][5][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[2][6][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[3][0][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[3][1][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[3][2][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[3][4][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[3][5][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[3][6][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[3][7][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[4][1][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[4][2][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[4][3][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[4][5][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[4][6][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[4][7][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[4][0][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[5][2][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[5][3][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[5][4][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[5][6][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[5][7][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[5][0][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[5][1][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[6][3][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[6][4][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[6][5][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[6][7][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[6][0][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[6][1][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[6][2][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[7][4][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[7][5][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[7][6][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[7][0][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[7][1][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[7][2][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[7][3][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 
	}
	for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				w_v1Ia_v2[0][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-0)/tauV2_1,2.0))); 
				w_v1Ia_v2[1][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-1)/tauV2_1,2.0))); 
				w_v1Ia_v2[2][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-2)/tauV2_1,2.0))); 
				w_v1Ia_v2[3][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-3)/tauV2_1,2.0))); 
				w_v1Ia_v2[4][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-0)/tauV2_1,2.0))); 
				w_v1Ia_v2[5][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-1)/tauV2_1,2.0))); 
				w_v1Ia_v2[6][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-2)/tauV2_1,2.0))); 
				w_v1Ia_v2[7][ii][0][jj] = w_v1Ia_2*exp(-(pow((0-3)/tauV2_1,2.0))); 
			}
	}
	for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				w_v1Ib_v2[0][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-0)/tauV2_1,2.0))); 
				w_v1Ib_v2[1][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-1)/tauV2_1,2.0))); 
				w_v1Ib_v2[2][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-2)/tauV2_1,2.0))); 
				w_v1Ib_v2[3][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-3)/tauV2_1,2.0))); 
				w_v1Ib_v2[4][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-0)/tauV2_1,2.0))); 
				w_v1Ib_v2[5][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-1)/tauV2_1,2.0))); 
				w_v1Ib_v2[6][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-2)/tauV2_1,2.0))); 
				w_v1Ib_v2[7][ii][0][jj] = w_v1Ib_2*exp(-(pow((0-3)/tauV2_1,2.0))); 
			}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj) w_v1P_v2[thetaa][ii][thetaa][jj] = w_v1P_2;
		}
	}

	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj) w_v2_v1[thetaa][ii][thetaa][jj] = w_v2_1;
		}
	}
		for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				GABA_ext[theta][i]  = GABA_c0; 
				GABA_extw[theta][i] = GABA_c0; 
			}
		}
}	

void display(void){
	printf("t=%d \n",t);
}