//! \file elastic_correspondence.cxx

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "elastic_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>

namespace CORRESPONDENCE {



template<typename ScalarT>
void updateElasticCauchyStress
(
const ScalarT* unrotatedRateOfDeformation, 
const ScalarT* unrotatedCauchyStressN, 
ScalarT* unrotatedCauchyStressNP1, 
const int numPoints, 
const double bulkMod,
const double shearMod,
const double dt
)
{
  // Hooke's law
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  for(int iID=0 ; iID<numPoints ; ++iID, 
        rateOfDef+=9, sigmaN+=9, sigmaNP1+=9){

      //strainInc = dt * rateOfDef
      for (int i = 0; i < 9; i++) {
          strainInc[i] = *(rateOfDef+i)*dt;
          deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //update stress
      for (int i = 0; i < 9; i++) {
          *(sigmaNP1+i) = *(sigmaN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(sigmaNP1) += bulkMod*dilatationInc;
      *(sigmaNP1+4) += bulkMod*dilatationInc;
      *(sigmaNP1+8) += bulkMod*dilatationInc;

  }
}

// Explicit template instantiation for double
template void updateElasticCauchyStress<double>
(
const double* unrotatedRateOfDeformation, 
const double* unrotatedCauchyStressN, 
double* unrotatedCauchyStressNP1, 
int numPoints, 
double bulkMod,
double shearMod,
double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */




template<typename ScalarT>
void updateElasticCauchyStressSmallDef
(
const ScalarT* DeformationGradient, 
const ScalarT* unrotatedCauchyStressN, 
ScalarT* unrotatedCauchyStressNP1, 
const int numPoints, 
const double bulkMod,
const double shearMod,
const double* C,
const int type,
const double dt
)
{
  // Hooke's law
  const ScalarT* defGrad = DeformationGradient;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  ScalarT strain[9];


  // 0 -> xx,  1 -> xy, 2 -> xz
  // 3 -> yx,  4 -> yy, 5 -> yz
  // 6 -> zx,  7 -> zy, 8 -> zz

//  if (type == 0){
//    C[0]  = 2*shearMod + bulkMod;C[1] = bulkMod;             C[2] = bulkMod;              C[3] = 0.0;      C[4]  = 0.0;      C[5]  = 0.0;
//    C[6]  = bulkMod;             C[7] = 2*shearMod + bulkMod;C[8] = bulkMod;              C[9] = 0.0;      C[10] = 0.0;      C[11] = 0.0;
//    C[12] = bulkMod;             C[13] = bulkMod;            C[14] = 2*shearMod + bulkMod;C[15] = 0.0;     C[16] = 0.0;      C[17] = 0.0;
//    C[18] = 0.0;                 C[19] = 0.0;                C[20] = 0.0;                 C[21] = shearMod;C[22] = 0.0;      C[23] = 0.0;
//    C[24] = 0.0;                 C[25] = 0.0;                C[26] = 0.0;                 C[27] = 0.0;     C[28] = shearMod; C[29]  = 0.0;
//    C[30] = 0.0;                 C[31] = 0.0;                C[32] = 0.0;                 C[33] = 0.0;     C[34] = 0.0;      C[35] = shearMod;
//  }  
//  if (type == 1){//plane strain  
//        double E = 9*bulkMod*shearMod/(3*bulkMod+shearMod);
//        double nu = (3*bulkMod-2*shearMod)/(6*bulkMod-2*shearMod);
//        double temp = E / (1-nu*nu);
//        C[0]  = temp;C[1] = nu*temp;   C[1]  = nu*temp;           C[2] = 0;              C[3] = 0.0;      C[4]  = 0.0;      C[5]  = 0.0;
//        C[6]  = nu*temp;             C[7] = temp;C[8] = 0;              C[9] = 0.0;      C[10] = 0.0;      C[11] = 0.0;
//        C[12] = 0;             C[13] = 0;            C[14] = 0;C[15] = 0.0;     C[16] = 0.0;      C[17] = 0.0;
//        C[18] = 0.0;                 C[19] = 0.0;                C[20] = 0.0;                 C[21] = 0;C[22] = 0.0;      C[23] = 0.0;
//        C[24] = 0.0;                 C[25] = 0.0;                C[26] = 0.0;                 C[27] = 0.0;     C[28] = 0; C[29]  = 0.0;
//        C[30] = 0.0;                 C[31] = 0.0;                C[32] = 0.0;                 C[33] = 0.0;     C[34] = 0.0;      C[35] = (1-nu)/2*temp;
//  }
//  
//  if (type == 2){//plane stress
//        double E = 9*bulkMod*shearMod/(3*bulkMod+shearMod);
//        double nu = (3*bulkMod-2*shearMod)/(6*bulkMod-2*shearMod);
//        double temp = E / (1-nu*nu);
//        C[0]  = temp;C[1] = nu*temp;  C[1]  = nu*temp;            C[2] = 0;              C[3] = 0.0;      C[4]  = 0.0;      C[5]  = 0.0;
//        C[6]  = nu*temp;             C[7] = temp;C[8] = 0;              C[9] = 0.0;      C[10] = 0.0;      C[11] = 0.0;
//        C[12] = 0;             C[13] = 0;            C[14] = 0;C[15] = 0.0;     C[16] = 0.0;      C[17] = 0.0;
//        C[18] = 0.0;                 C[19] = 0.0;                C[20] = 0.0;                 C[21] = 0;C[22] = 0.0;      C[23] = 0.0;
//        C[24] = 0.0;                 C[25] = 0.0;                C[26] = 0.0;                 C[27] = 0.0;     C[28] = 0; C[29]  = 0.0;
//        C[30] = 0.0;                 C[31] = 0.0;                C[32] = 0.0;                 C[33] = 0.0;     C[34] = 0.0;      C[35] = (1-nu)/2*temp;
//  }

  for(int iID=0 ; iID<numPoints ; ++iID, 
        defGrad+=9, sigmaN+=9, sigmaNP1+=9){
            
            // df muss genutzt werden EQ 34 Silling Stability --> d.h. ich muss das Delta bestimmen
            // Schadensmodell muss angepasst werden
            // Referenzkonfiguration muss angepasst werden
      
      if (type==0){
        strain[0] = 0.5 * ( *(defGrad)*   *(defGrad)   + *(defGrad+1)* *(defGrad+3) + *(defGrad+2) * *(defGrad+6)  - 1.0 );
        strain[1] = 0.5 * ( *(defGrad)*   *(defGrad+1) + *(defGrad+1)* *(defGrad+4) + *(defGrad+2) * *(defGrad+7)  );
        strain[2] = 0.5 * ( *(defGrad)*   *(defGrad+2) + *(defGrad+1)* *(defGrad+5) + *(defGrad+2) * *(defGrad+8)  );
        strain[3] = 0.5 * ( *(defGrad+3)* *(defGrad)   + *(defGrad+4)* *(defGrad+3) + *(defGrad+5) * *(defGrad+6)  );
        strain[4] = 0.5 * ( *(defGrad+3)* *(defGrad+1) + *(defGrad+4)* *(defGrad+4) + *(defGrad+5) * *(defGrad+7)  - 1.0 );
        strain[5] = 0.5 * ( *(defGrad+3)* *(defGrad+2) + *(defGrad+4)* *(defGrad+5) + *(defGrad+5) * *(defGrad+8)  );
        strain[6] = 0.5 * ( *(defGrad+6)* *(defGrad)   + *(defGrad+7)* *(defGrad+3) + *(defGrad+8) * *(defGrad+6)  );
        strain[7] = 0.5 * ( *(defGrad+6)* *(defGrad+1) + *(defGrad+7)* *(defGrad+4) + *(defGrad+8) * *(defGrad+7)  );
        strain[8] = 0.5 * ( *(defGrad+6)* *(defGrad+2) + *(defGrad+7)* *(defGrad+5) + *(defGrad+8) * *(defGrad+8)  - 1.0 );
        
        *(sigmaNP1)   = C[0] *strain[0] + C[1]* strain[4] + C[2]* strain[8] + C[3]*(strain[5]+strain[7]) + C[4]*(strain[2]+strain[6]) + C[5]*(strain[1]+strain[3]);
        *(sigmaNP1+1) = C[30]*strain[0] + C[31]*strain[4] + C[32]*strain[8] + C[33]*(strain[5]+strain[7])+ C[34]*(strain[2]+strain[6])+ C[35]*(strain[1]+strain[3]);
        *(sigmaNP1+2) = C[24]*strain[0] + C[25]*strain[4] + C[26]*strain[8] + C[27]*(strain[5]+strain[7])+ C[28]*(strain[2]+strain[6]) + C[29]*(strain[1]+strain[3]);
        *(sigmaNP1+3) = *(sigmaNP1+1);
        *(sigmaNP1+4) = C[6] *strain[0] + C[7]* strain[4] + C[8]* strain[8] + C[9]*(strain[5]+strain[7]) + C[10]*(strain[2]+strain[6]) + C[11]*(strain[1]+strain[3]);
        *(sigmaNP1+5) = C[18]*strain[0] + C[19]*strain[4] + C[20]*strain[8] + C[21]*(strain[5]+strain[7])+ C[22]*(strain[2]+strain[6]) + C[23]*(strain[1]+strain[3]);
        *(sigmaNP1+6) = *(sigmaNP1+2);
        *(sigmaNP1+7) = *(sigmaNP1+5);
        *(sigmaNP1+8) = C[12]*strain[0] + C[13]*strain[4] + C[14]*strain[8] + C[15]*(strain[5]+strain[7])+ C[16]*(strain[2]+strain[6]) + C[17]*(strain[1]+strain[3]);
        
      }
      if (type!=0){
        strain[0] = 0.5 * ( *(defGrad)*   *(defGrad)   + *(defGrad+1)* *(defGrad+3) - 1.0 );
        strain[1] = 0.5 * ( *(defGrad)*   *(defGrad+1) + *(defGrad+1)* *(defGrad+4) );
        strain[3] = 0.5 * ( *(defGrad+3)* *(defGrad)   + *(defGrad+4)* *(defGrad+3) );
        strain[4] = 0.5 * ( *(defGrad+3)* *(defGrad+1) + *(defGrad+4)* *(defGrad+4) +   - 1.0 );
        
        *(sigmaNP1)   = C[0] *strain[0] + C[1]* strain[4] + C[5] *(strain[1]+strain[3]);
        *(sigmaNP1+1) = C[30]*strain[0] + C[31]*strain[4] + C[35]*(strain[1]+strain[3]);

        *(sigmaNP1+2) = 0.0;
        
        *(sigmaNP1+3) = *(sigmaNP1+1);
        *(sigmaNP1+4) = C[6] *strain[0] + C[7]* strain[4] + C[11]*(strain[1]+strain[3]);
        *(sigmaNP1+5) = 0.0;
        *(sigmaNP1+6) = 0.0;
        *(sigmaNP1+7) = 0.0;
        *(sigmaNP1+8) = 0.0;

      }

  }
}

// Explicit template instantiation for double
template void updateElasticCauchyStressSmallDef<double>
(
const double* DeformationGradient, 
const double* unrotatedCauchyStressN, 
double* unrotatedCauchyStressNP1, 
int numPoints, 
double bulkMod,
double shearMod,
const double* C,
const int type,
double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */

}
