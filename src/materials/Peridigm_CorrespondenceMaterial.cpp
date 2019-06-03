/*! \file Peridigm_CorrespondenceMaterial.cpp */

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

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::CorrespondenceMaterial::CorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0), m_hourglassCoefficient(0.0),
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_horizonFieldId(-1), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_velocitiesFieldId(-1), 
    m_hourglassForceDensityFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
    m_deformationGradientFieldId(-1),
    m_shapeTensorInverseFieldId(-1),
    m_leftStretchTensorFieldId(-1),
    m_rotationTensorFieldId(-1), 
    m_unrotatedCauchyStressFieldId(-1),
    m_cauchyStressFieldId(-1), 
    m_unrotatedRateOfDeformationFieldId(-1),
    m_partialStressFieldId(-1)
{
     
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
  m_stabilizationType = 2;
  D = 0.0;
  if (params.isParameter("Velocity Stabilization")){
    D = params.get<double>("Velocity Stabilization");
  }
  if (params.isParameter("Stabilizaton Type")){
      
    if (params.get<string>("Stabilizaton Type")=="Bond Based") m_stabilizationType = 1;
    if (params.get<string>("Stabilizaton Type")=="State Based") m_stabilizationType = 2;
    }
  
  if (params.isParameter("Elastic Correspondence")){
     // std::cout<<"Use Material: Elastic Correspondence"<<std::endl;
	  nonLin = true;
  }
  m_plane = false;
  //************************************
  // wie komme ich an den Namen??
  //************************************
  //if (params.isParameter("Linear Elastic Correspondence")){
    nonLin = false;
    //std::cout<<"Use Material: Linear Elastic Correspondence"<<std::endl;
    bool m_planeStrain = false, m_planeStress = false;
    if (params.isParameter("Plane Strain"))
        m_planeStrain = params.get<bool>("Plane Strain");
        
    if (params.isParameter("Plane Stress"))
        m_planeStress = params.get<bool>("Plane Stress");
        
    if (m_planeStrain==true){
        m_plane=true;
        std::cout<<"Method 2D Plane Strain"<<std::endl;
    }
    if (m_planeStress==true){
        m_plane=true;
        std::cout<<"WRN: Method 2D Plane Stress --> not fully implemented yet"<<std::endl;
    }
    
  //}

  
  
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_horizonFieldId                    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_hourglassForceDensityFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Hourglass_Force_Density");
  m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
  m_leftStretchTensorFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_shapeTensorInverseFieldId         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Shape_Tensor_Inverse");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_partialStressFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  m_detachedNodesFieldId              = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached Nodes");

  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_hourglassForceDensityFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_partialStressFieldId);
  m_fieldIds.push_back(m_detachedNodesFieldId);
}

PeridigmNS::CorrespondenceMaterial::~CorrespondenceMaterial()
{
}

void
PeridigmNS::CorrespondenceMaterial::initialize(const double dt,
                                               const int numOwnedPoints,
                                               const int* ownedIDs,
                                               const int* neighborhoodList,
                                               PeridigmNS::DataManager& dataManager)
{
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *leftStretchTensorN;
  double *leftStretchTensorNP1;
  double *rotationTensorN;
  double *rotationTensorNP1;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);

  //Initialize the left stretch and rotation tenor to the identity matrix
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorNP1, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorNP1, numOwnedPoints);

  //Initialize the inverse of the shape tensor and the deformation gradient
  double *volume;
  double *horizon;
  double *modelCoordinates;
  double *coordinates;
  double *shapeTensorInverse;
  double *deformationGradient;
  double *bondDamage;
  double *detachedNodes;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  
  int shapeTensorReturnCode = 
    CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                               horizon,
                                                                               modelCoordinates,
                                                                               coordinates,
                                                                               shapeTensorInverse,
                                                                               deformationGradient,
                                                                               bondDamage,
                                                                               neighborhoodList,
                                                                               numOwnedPoints,
                                                                               detachedNodes,
                                                                               m_plane);

  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::initialize() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);

}

void
PeridigmNS::CorrespondenceMaterial::computeForce(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces and partial stress
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *horizon, *volume, *modelCoordinates, *coordinates, *velocities, *shapeTensorInverse, *deformationGradient, *bondDamage;
  double *detachedNodes;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  //*(dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_N));
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
 
  // Compute the inverse of the shape tensor and the approximate deformation gradient
  // The approximate deformation gradient will be used by the derived class (specific correspondence material model)
  // to compute the Cauchy stress.
  // The inverse of the shape tensor is stored for later use after the Cauchy stress calculation
   int shapeTensorReturnCode = 
    CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                               horizon,
                                                                               modelCoordinates,
                                                                               coordinates,
                                                                               shapeTensorInverse,
                                                                               deformationGradient,
                                                                               bondDamage,
                                                                               neighborhoodList,
                                                                               numOwnedPoints,
                                                                               detachedNodes,
                                                                               m_plane);
  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);

  double *leftStretchTensorN, *leftStretchTensorNP1, *rotationTensorN, *rotationTensorNP1, *unrotatedRateOfDeformation;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);
  
  // Compute left stretch tensor, rotation tensor, and unrotated rate-of-deformation.
  // Performs a polar decomposition via Flanagan & Taylor (1987) algorithm.
  
  int rotationTensorReturnCode = 0;
  if (nonLin==true){  
     rotationTensorReturnCode = CORRESPONDENCE::computeUnrotatedRateOfDeformationAndRotationTensor(volume,
                                                                                                      horizon,
                                                                                                      modelCoordinates, 
                                                                                                      velocities, 
                                                                                                      deformationGradient,
                                                                                                      shapeTensorInverse,
                                                                                                      leftStretchTensorN,
                                                                                                      rotationTensorN,
                                                                                                      leftStretchTensorNP1,
                                                                                                      rotationTensorNP1,
                                                                                                      unrotatedRateOfDeformation,
                                                                                                      neighborhoodList, 
                                                                                                      numOwnedPoints, 
                                                                                                      dt,
                                                                                                      bondDamage,
                                                                                                      detachedNodes, 
                                                                                                      m_plane);

      string rotationTensorErrorMessage =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
      rotationTensorErrorMessage +=
        "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
      string rotationTensorErrorMessage2 =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to invert deformation gradient tensor.\n";
        
      string rotationTensorErrorMessage3 =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to invert temp.\n";

      TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode == 1, rotationTensorErrorMessage);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode == 2, rotationTensorErrorMessage2);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode == 3, rotationTensorErrorMessage3);
 
      CORRESPONDENCE::setValuesForDetachedNodes(deformationGradient,
                                                leftStretchTensorNP1,
                                                rotationTensorNP1,
                                                unrotatedRateOfDeformation,
                                                shapeTensorInverse,
                                                detachedNodes,
                                                numOwnedPoints);
     }  

  // Evaluate the Cauchy stress using the routine implemented in the derived class (specific correspondence material model)
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  // computeCauchyStress() typically uses the following fields which are accessed via the DataManager:
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1
  

  computeCauchyStress(dt, numOwnedPoints, dataManager);

  // rotate back to the Eulerian frame
  double *unrotatedCauchyStressNP1, *cauchyStressNP1;
  
  if (nonLin == true) {
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);
  
  CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                     unrotatedCauchyStressNP1,
                                     cauchyStressNP1,
                                     numOwnedPoints);
  }
  else
  {
    dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
    dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);
    cauchyStressNP1 = unrotatedCauchyStressNP1;
  }
  // Cauchy stress is now updated and in the rotated state.  Proceed with
  // conversion to Piola-Kirchoff and force-vector states.

  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double *partialStress;
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);

  double *delta = horizon;
  double* stress = cauchyStressNP1;
  double* shapeTensorInv = shapeTensorInverse;
  double* defGrad = deformationGradient;

  double *modelCoordinatesPtr, *neighborModelCoordinatesPtr, *forceDensityPtr, *neighborForceDensityPtr, *partialStressPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double TX, TY, TZ, omega, vol, neighborVol;
  int numNeighbors, neighborIndex;

  string matrixInversionErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to invert deformation gradient.\n";
  matrixInversionErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  matrixInversionErrorMessage +=
    "****         Force calculation";
     // Compute hourglass forces for stabilization of low-energy and/or zero-energy modes
  dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  
  // Missbrauch des Feldes zum Testen
  double *hourglassForceDensity;
  dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&hourglassForceDensity);

    
  vector<double>  piolaStressVector(9), tempVector(9), defGradInvVector(9);
  double* defGradInv = &defGradInvVector[0];
  double* piolaStress = &piolaStressVector[0];
  double* temp = &tempVector[0];
  double jacobianDeterminant;

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  const int *neighborListPtr = neighborhoodList;
  int bondIndex = 0;
  int matrixInversionReturnCode(0);
  for(int iID=0 ; iID<numOwnedPoints ; ++iID, 
          ++delta, defGrad+=9, stress+=9, shapeTensorInv+=9){

    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient and store the determinant
   
    //std::cout<<defGrad[0]<<" "<<defGrad[1]<<" "<<defGrad[2]<<" "<<defGrad[3]<<" "<<defGrad[4]<<" "<<defGrad[5]<<" "<<defGrad[6]<<" "<<defGrad[7]<<" "<<defGrad[8]<<std::endl;
    if (detachedNodes[iID]==0){
		
		if (m_plane==true){
			matrixInversionReturnCode =
			CORRESPONDENCE::Invert2by2Matrix(defGrad, jacobianDeterminant, defGradInv);
			}
        else{
			matrixInversionReturnCode =
			CORRESPONDENCE::Invert3by3Matrix(defGrad, jacobianDeterminant, defGradInv);
			}
		

        TEUCHOS_TEST_FOR_EXCEPT_MSG(matrixInversionReturnCode != 0, matrixInversionErrorMessage);
        
        //P = J * \sigma * F^(-T)
        CORRESPONDENCE::MatrixMultiply(false, true, jacobianDeterminant, stress, defGradInv, piolaStress);

        // Inner product of Piola stress and the inverse of the shape tensor
        CORRESPONDENCE::MatrixMultiply(false, false, 1.0, piolaStress, shapeTensorInv, temp);
    }
    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;
    numNeighbors = *neighborListPtr; neighborListPtr++;
    forceDensityPtr = forceDensity + 3*iID;
    int countNeighbors = 0;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborModelCoordinatesPtr = modelCoordinates + 3*neighborIndex;
      if (detachedNodes[neighborIndex]==0){
        undeformedBondX = *(neighborModelCoordinatesPtr)   - *(modelCoordinatesPtr);
        undeformedBondY = *(neighborModelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
        undeformedBondZ = *(neighborModelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);
        undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                    undeformedBondY*undeformedBondY +
                                    undeformedBondZ*undeformedBondZ);
    
        omega = m_OMEGA(undeformedBondLength, *delta);
        
        TX =  (1-bondDamage[bondIndex]) * omega * ( *(temp)   * undeformedBondX + *(temp+1) * undeformedBondY + *(temp+2) * undeformedBondZ );
        TY =  (1-bondDamage[bondIndex]) * omega * ( *(temp+3) * undeformedBondX + *(temp+4) * undeformedBondY + *(temp+5) * undeformedBondZ );
        TZ =  (1-bondDamage[bondIndex]) * omega * ( *(temp+6) * undeformedBondX + *(temp+7) * undeformedBondY + *(temp+8) * undeformedBondZ );
       
        neighborVol = volume[neighborIndex];
        vol = volume[iID];
        
        neighborForceDensityPtr = forceDensity + 3*neighborIndex;
    
        *(forceDensityPtr)   += TX * neighborVol;
        *(forceDensityPtr+1) += TY * neighborVol;
        *(forceDensityPtr+2) += TZ * neighborVol;
        *(neighborForceDensityPtr)   -= TX * vol;
        *(neighborForceDensityPtr+1) -= TY * vol;
        *(neighborForceDensityPtr+2) -= TZ * vol;
    
        partialStressPtr = partialStress + 9*iID;
        *(partialStressPtr)   += TX*undeformedBondX*neighborVol;
        *(partialStressPtr+1) += TX*undeformedBondY*neighborVol;
        *(partialStressPtr+2) += TX*undeformedBondZ*neighborVol;
        *(partialStressPtr+3) += TY*undeformedBondX*neighborVol;
        *(partialStressPtr+4) += TY*undeformedBondY*neighborVol;
        *(partialStressPtr+5) += TY*undeformedBondZ*neighborVol;
        *(partialStressPtr+6) += TZ*undeformedBondX*neighborVol;
        *(partialStressPtr+7) += TZ*undeformedBondY*neighborVol;
        *(partialStressPtr+8) += TZ*undeformedBondZ*neighborVol;
        
        
      //  hourglassForceDensity[iID]   += (1-bondDamage[bondIndex]) *velocities[neighborIndex]  ;
      //  hourglassForceDensity[iID+1] += (1-bondDamage[bondIndex]) *velocities[neighborIndex+1];
      //  hourglassForceDensity[iID+2] += (1-bondDamage[bondIndex]) *velocities[neighborIndex+2];
      //  countNeighbors += bondDamage[bondIndex];
      }
      bondIndex += 1;
    //  hourglassForceDensity[iID]  /= (numNeighbors-countNeighbors+1);
    //  hourglassForceDensity[iID+1]/= (numNeighbors-countNeighbors+1);
    //  hourglassForceDensity[iID+2]/= (numNeighbors-countNeighbors+1);
    }
  }
  
  
  //for(int iID=0 ; iID<numOwnedPoints ; ++iID){
  //  velocities[iID] = D*hourglassForceDensity[iID] + (1-D)*velocities[iID]; 
  //  velocities[iID+1] = D*hourglassForceDensity[iID+1] + (1-D)*velocities[iID+1];
  //  velocities[iID+2] = D*hourglassForceDensity[iID+2] + (1-D)*velocities[iID+2];
  //}  
 
  //double *hourglassForceDensity;
  //dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&hourglassForceDensity);
  //dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
   // \todo HOURGLASS FORCES ARE NOT OUTPUT TO EXODUS CORRECTLY BECAUSE THEY ARE NOT ASSEMBLED ACROSS PROCESSORS.
  //       They are summed into the force vector below, and the force vector is assembled across processors,
  //       so the calculation runs correctly, but the hourglass output is off.

  if (m_stabilizationType == 1){
  CORRESPONDENCE::computeHourglassForce(volume,
                                        horizon,
                                        modelCoordinates,
                                        coordinates,
                                        deformationGradient,
                                        forceDensity,
                                        neighborhoodList,
                                        numOwnedPoints,
                                        m_bulkModulus,
                                        m_hourglassCoefficient,
                                        bondDamage
                                        );
  }
  
  if (m_stabilizationType == 2){
  
  CORRESPONDENCE::computeCorrespondenceStabilityForce(volume,
                                        horizon,
                                        modelCoordinates,
                                        coordinates,
                                        deformationGradient,
                                        forceDensity,
                                        neighborhoodList,
                                        numOwnedPoints,
                                        m_bulkModulus,
                                        m_hourglassCoefficient,
                                        bondDamage
                                        
                                        );
  }
  // Sum the hourglass force densities into the force densities
  Teuchos::RCP<Epetra_Vector> forceDensityVector = dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1);
  Teuchos::RCP<Epetra_Vector> hourglassForceDensityVector = dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1);
  //forceDensityVector->Update(1.0, *hourglassForceDensityVector, 1.0);
}
