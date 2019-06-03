/*! \file Peridigm_ElasticLinearCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticLinearCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticLinearCorrespondenceMaterial::ElasticLinearCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_unrotatedRateOfDeformationFieldId(-1),
    m_unrotatedCauchyStressFieldId(-1),
    m_deformationGradientFieldId(-1)
{
  bool m_planeStrain = false, m_planeStress = false;
  m_type = 0;
  m_density = params.get<double>("Density");
  if (params.isParameter("Plane Strain")){
    m_planeStrain = params.get<bool>("Plane Strain");
    }
  if (params.isParameter("Plane Stress")){
    m_planeStress = params.get<bool>("Plane Stress");
    }
  if (m_planeStrain==true)m_type=1;
  if (m_planeStress==true)m_type=2;
  
  double C11, C44, C55, C66, C12, C13, C14, C15, C16, C22, C33, C23, C24, C25, C26, C34, C35, C36, C45, C46, C56;   
  if (params.isParameter("Material Symmetry")){
    if (params.get<string>("Material Symmetry")=="Isotropic"){
       C11 = params.get<double>("C11");
       C44 = params.get<double>("C44");
       C55 = params.get<double>("C44");
       C66 = params.get<double>("C44");
       C12 = C11 - 2*C55;
       C13 = C12;
       C14 = 0.0;
       C15 = 0.0;
       C16 = 0.0;
       C22 = params.get<double>("C11");
       C33 = params.get<double>("C11");
       C23 = C12;
       C24 = 0.0;
       C25 = 0.0;
       C26 = 0.0;
       C34 = 0.0;
       C35 = 0.0;
       C36 = 0.0;
       C45 = 0.0;
       C46 = 0.0;
       C56 = 0.0; 
    }
   if (params.get<string>("Material Symmetry")=="Anisotropic"){
       C11 = params.get<double>("C11");
       C12 = params.get<double>("C12");
       C13 = params.get<double>("C13");
       C14 = params.get<double>("C14");
       C15 = params.get<double>("C15");
       C16 = params.get<double>("C16");
       C22 = params.get<double>("C22");
       C23 = params.get<double>("C23");
       C24 = params.get<double>("C24");
       C25 = params.get<double>("C25");
       C26 = params.get<double>("C26");
       C33 = params.get<double>("C33");
       C34 = params.get<double>("C34");
       C35 = params.get<double>("C35");
       C36 = params.get<double>("C36");
       C44 = params.get<double>("C44");
       C45 = params.get<double>("C45");
       C46 = params.get<double>("C46");
       C55 = params.get<double>("C55");
       C56 = params.get<double>("C56");
       C66 = params.get<double>("C66");
   }
  }
  else{
       m_bulkModulus = calculateBulkModulus(params);
       m_shearModulus = calculateShearModulus(params);
  
       C11 = 2*m_shearModulus + m_bulkModulus;
       C44 = m_shearModulus;
       C55 = m_shearModulus;
       C66 = m_shearModulus;
       C12 = C11 - 2*C55;
       C13 = C12;
       C14 = 0.0;
       C15 = 0.0;
       C16 = 0.0;
       C22 = 2*m_shearModulus + m_bulkModulus;
       C33 = 2*m_shearModulus + m_bulkModulus;
       C23 = C12;
       C24 = 0.0;
       C25 = 0.0;
       C26 = 0.0;
       C34 = 0.0;
       C35 = 0.0;
       C36 = 0.0;
       C45 = 0.0;
       C46 = 0.0;
       C56 = 0.0;
  }
// Equation (8) Dipasquale, D., Sarego, G., Zaccariotto, M., Galvanetto, U., A discussion on failure criteria
  // for ordinary state-based Peridynamics, Engineering Fracture Mechanics (2017), doi: https://doi.org/10.1016/
  // j.engfracmech.2017.10.011
  if (m_planeStrain==true)m_plane=true;
  if (m_planeStress==true)m_plane=true;
  if (m_plane==false){
   C[0]  = C11;C[1]  = C12;C[2]  = C13; C[3]  = C14; C[4]  = C15; C[5]  = C16;
   C[6]  = C12;C[7]  = C22;C[8]  = C23; C[9]  = C24; C[10] = C25; C[11] = C26;
   C[12] = C13;C[13] = C23;C[14] = C33; C[15] = C34; C[16] = C35; C[17] = C36;
   C[18] = C14;C[19] = C24;C[20] = C34; C[21] = C44; C[22] = C45; C[23] = C46;
   C[24] = C15;C[25] = C25;C[26] = C35; C[27] = C45; C[28] = C55; C[29] = C56;
   C[30] = C16;C[31] = C26;C[32] = C36; C[33] = C46; C[34] = C56; C[35] = C66;
  }
  // tbd in future
  if (m_planeStress==true){
      //only transversal isotropic in the moment
   C[0]  = C11-C13*C13/C22;C[1]  = C12-C13*C23/C22;C[2]  = 0.0; C[3]  = 0.0; C[4]  = 0.0; C[5]  = 0.0;
   C[6]  = C12-C13*C23/C22;C[7]  = C22-C13*C23/C22;C[8]  = 0.0; C[9]  = 0.0; C[10] = 0.0; C[11] = 0.0;
   C[12] = 0.0;C[13] = 0.0;C[14] = 0.0; C[15] = 0.0; C[16] = 0.0; C[17] = 0.0;
   C[18] = 0.0;C[19] = 0.0;C[20] = 0.0; C[21] = 0.0; C[22] = 0.0; C[23] = 0.0;
   C[24] = 0.0;C[25] = 0.0;C[26] = 0.0; C[27] = 0.0; C[28] = 0.0; C[29] = 0.0;
   C[30] = 0.0;C[31] = 0.0;C[32] = 0.0; C[33] = 0.0; C[34] = 0.0; C[35] = C66;
   

  }
  // not correct for plane stress!!!
  if (m_planeStrain==true | m_planeStress==true){
   C[0]  = C11;C[1]  = C12;C[2]  = 0.0; C[3]  = 0.0; C[4]  = 0.0; C[5]  = C16;
   C[6]  = C12;C[7]  = C22;C[8]  = 0.0; C[9]  = 0.0; C[10] = 0.0; C[11] = C26;
   C[12] = 0.0;C[13] = 0.0;C[14] = 0.0; C[15] = 0.0; C[16] = 0.0; C[17] = 0.0;
   C[18] = 0.0;C[19] = 0.0;C[20] = 0.0; C[21] = 0.0; C[22] = 0.0; C[23] = 0.0;
   C[24] = 0.0;C[25] = 0.0;C[26] = 0.0; C[27] = 0.0; C[28] = 0.0; C[29] = 0.0;
   C[30] = C16;C[31] = C26;C[32] = 0.0; C[33] = 0.0; C[34] = 0.0; C[35] = C66;
  }

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");

  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
}

PeridigmNS::ElasticLinearCorrespondenceMaterial::~ElasticLinearCorrespondenceMaterial()
{
}

void
PeridigmNS::ElasticLinearCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               PeridigmNS::DataManager& dataManager) const
{
  double *unrotatedCauchyStressN;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);

  double *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);

 // double *unrotatedRateOfDeformation;
 // dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);
  
  double *defGrad;
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&defGrad);
  
  CORRESPONDENCE::updateElasticCauchyStressSmallDef(defGrad, 
                                            unrotatedCauchyStressN, 
                                            unrotatedCauchyStressNP1,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_shearModulus,
                                            C,
                                            m_type,
                                            dt);
}
