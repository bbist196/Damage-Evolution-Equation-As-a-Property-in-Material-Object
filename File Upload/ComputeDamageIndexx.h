/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "TimeIntegrator.h"
#include "FEProblem.h"
#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "StrainEnergyDensity.h"


/**
 * Scalar damage model that defines the damage parameter using a material property
 */
 template <bool is_ad>
class ComputeDamageIndexxTempl :  public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  ComputeDamageIndexxTempl(const InputParameters & parameters);
  
 

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  const Real & dt() const { return _dt; }
  
protected:


  
  // Parameer defining damage evolution equaiton
  const MaterialProperty<Real> & _Yin;
  const MaterialProperty<Real> & _P1;
  const MaterialProperty<Real> & _P2;
  const MaterialProperty<Real> & _Mu;
  Real & _dt;

  /// The damage index as material property
  MaterialProperty<Real> & _damage_indexx;
  const MaterialProperty<Real> & _damage_indexx_old;
  /// Strain energy density as material property


  /// for computing strain energy density
  /// Base name of the material system
  const std::string _base_name;

  /// stress name, for different stresses on the same material system
  const std::string _stress_name;



  ///{@ Current and old values of stress
 // const GenericMaterialProperty<RankTwoTensor, is_ad> & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  ///@}

  /// Current value of mechanical strain which includes elastic and
  /// inelastic components of the strain
  const GenericMaterialProperty<RankTwoTensor, is_ad> & _mechanical_strain;

  
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  
 
 
};
typedef ComputeDamageIndexxTempl<false> ComputeDamageIndexx;
typedef ComputeDamageIndexxTempl<true> ADComputeDamageIndexx;
