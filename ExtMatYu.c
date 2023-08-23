/** External material model to linearly hardening plastic solid with Mises Yielding and damaged-induced shrinkage **/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif


EXPORT int eval(double e[6],         // Input: Green-Lagrange strain tensor components in Voigt order (xx,yy,zz,yz,zx,xy)
                double s[6],         // Output: Second Piola-Kirchhoff stress components in Voigt order (xx,yy,zz,yz,zx,xy)
                double Jac[6][6],    // Output: Jacobian of stress with respect to strain, 6-by-6 matrix in row-major order
                int *nPar,           // Input: Number of material model parameters, scalar
                double *par,         // Input: Parameters: par[0] = K, par[1] = G, par[2] = sigyYs0, par[3] = Eiso, Par[4] = comp.d, Par[5] = Fna, Par[6] = Fnb
                
                int *nStates1,       // Input: Number of states, scalar        
                double *states1,     // Internal state: Plastic strain tensor components in Voigt order (xx,yy,zz,yz,zx,xy)
				int *nStates2,       // Input: Number of states, scalar        
                double *states2) {   // Internal state: Accumulated effective plastic strain

  int i,j;                       // Iteration indicies
  double Kb,Mu,sigYs0,Eiso,d;    // Material parameters: Bulk modulus, Shear modulus, Yield stress, Hardening modulus, & phase field
  double K,G;                    // Bulk and shear modulus with phase field induced degradation
  double Fna, Fnb;               // Control parameter for yield surface degradation function
  double Fgy;                    // Yield surface degradation multiplier & degraded hardening modulus
  double kp,dd,cc;               // perturbation coefficient and damage multipliers for volumetric/deviatoric parts
  double ep[6];                  // Plastic strain tensor
  double trE;                    // First invariants of trial elastic strain tensor
  double eTrial[6];              // Trial strain deviator
  double sTrial[6];              // Trial stress deviator
  double N[6];                   // Plastic flow direction at trial state
  double sigTrialEff,sigYs,Fn;   // Effective trial stress, current yield stress, & trial yield function
  double dep,epEff;              // Plastic strain increment & effective plastic strain
  double q1,q2;                  // Two coefficients for the elastic-plastic Jacobian
  double theta1,theta2;          // Elastoplasticity related terms

  // Check inputs
  if (nPar[0] != 7)              // space for input parameters
    return 1;                    // error code 1 = "Wrong number of parameters"
  if (nStates1[0] !=6)           // space for the plastic strain tensor (six components)
    return 2;                    // error code 2 = "Wrong number of states"
  if (nStates2[0] !=1)           // space for the effective plastic strain
    return 2;                    // error code 2 = "Wrong number of states" 

  // Read input parameters from the parameter vector      
  Kb      = par[0];
  Mu      = par[1];
  sigYs0  = par[2];
  Eiso    = par[3];
  d       = par[4];
  Fna     = par[5];
  Fnb     = par[6];


  // Check values of input parameters
  if (Kb < 0.0) return -1;
  if (Mu < 0.0) return -1;
  if (sigYs0 < 0.0) return -1;
  if (Eiso < 0.0) return -1;
  if (d < 0.0) return -1;

  // Read the plastic strain tensor, effective plastic strain from the state variables
  for (i=0;i<6;i++) ep[i] = states1[i];
  epEff = states2[0];

  // Initialize the Jacobian
  for (i = 0; i < 6; i++){
    for (j = 0; j < 6; j++) {
      Jac[i][j] = 0.0;
      }
  }
  
  /*
   Stress update algorithm & Jacobian computation
  */

  // 1. Compute (deviatoric) trial strain tensor
  for (i=0;i<6;i++) eTrial[i] = e[i]-ep[i];
  trE = eTrial[0] + eTrial[1] + eTrial[2];
  for (i=0;i<3;i++) eTrial[i] = eTrial[i]-trE/3.0;

  // 2. compute damage multiplier and update material properties
  kp = pow(10.0, -9.0);
  cc = (1.0-kp)*pow((1.0-d), 2.0)+kp;
  if (trE <= 0.0) dd = 1.0;
  if (trE > 0.0) dd = cc;
  G = Mu*cc;
  K = Kb*dd;

  // 3. Compute (deviatoric) trial stress tensor and its effective value
  for (i = 0; i < 6; i++) sTrial[i] = 2.0*G*eTrial[i];
  sigTrialEff = 0.0;
  for (i = 0; i < 3; i++) sigTrialEff += sTrial[i]*sTrial[i] + 2.0*sTrial[i+3]*sTrial[i+3];
  sigTrialEff = sqrt(sigTrialEff);

  // 4. Re-construct the yield stress from the previous increment
  Fgy = (1.0-Fnb)*pow((1.0-d), Fna)+Fnb;
  sigYs = Fgy*sqrt(2.0/3.0)*(sigYs0 + Eiso*epEff);

  // 5. Construct trial yield function
  Fn = sigTrialEff - sigYs;

  // 6. Check yield condition
 
  if (Fn <= 0.0) {

     // Elastic step

     // E1. Compute the stress tensor (trial evaluation is the solution)
     for (i=0;i<6;i++) s[i] = sTrial[i];
     for (i=0;i<3;i++) s[i] += K*trE;

     // E2. Compute the elastic Jacobian
	 Jac[0][0] = Jac[1][1] = Jac[2][2] = 4.0/3.0*G + K;                        
     Jac[3][3] = Jac[4][4] = Jac[5][5] = 2.0*G;                                                      
     Jac[0][1] = Jac[0][2] = Jac[1][0] = Jac[1][2] = Jac[2][0] = Jac[2][1] = K - 2.0/3.0*G;

     return 0; // Return value 0 if success, any other value triggers an exception
  }
  else {

     // Elastic-Plastic step
     
     // EP1. Compute the effective plastic strain increment
     dep = Fn/(2.0*G + 2.0/3.0*Eiso);

     // EP2. Compute plastic flow direction
     for (i=0;i<6;i++) N[i] = sTrial[i]/sigTrialEff;

	 // EP3. Compute the stress deviator, and update plastic strains
     for (i=0;i<6;i++) s[i] = sTrial[i] - 2.0*G*dep*N[i];

     for (i=0;i<6;i++) states1[i] = ep[i] + dep*N[i]; 
	 epEff += sqrt(2.0/3.0)*dep;
	 states2[0] = epEff; 

     // EP4. Compute Jacobian related parameters
     theta1 = 3.0*G/(3.0*G + Eiso);
     theta2 = 2.0*G/sigTrialEff;
     q1 = 1.0 - theta2*dep;
     q2 = theta2*dep - theta1;
	
     // EP5. Compute Jacobian (Part I with Id in Voigt order)
	 Jac[0][0] = Jac[1][1] = Jac[2][2] = 2.0/3.0*q1*2.0*G; 
     Jac[3][3] = Jac[4][4] = Jac[5][5] = 1.0*q1*2.0*G; 
     Jac[0][1] = Jac[0][2] = Jac[1][0] = Jac[1][2] = Jac[2][0] = Jac[2][1] = -1.0/3.0*q1*2.0*G; 

     for (i=0;i<3;i++) {
       for (j=0;j<3;j++) {
         Jac[i][j] += K;
         }
     }

	 // EP6. Compute Jacobian (Part II with plastic flow direction in Voigt order)

  	 for (i=0;i<6;i++) {
       for (j=0;j<6;j++) {
         Jac[i][j] += 2.0*G*q2*N[i]*N[j];
		 }
     }

     // EP7. Add hydrostatic part to the second Piola-Kirchhoff stress tensor
	 for (i=0;i<3;i++) s[i] += K*trE;

     return 0;  // Return value 0 if success, any other value triggers an exception
  }

}