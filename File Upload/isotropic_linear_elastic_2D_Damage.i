[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 1
  xmin = 0
  ymin = 0
  xmax = 1.0
  ymax = 1.0
  
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[AuxVariables]
  
   [w]
    order = CONSTANT
    family = MONOMIAL
  []
  
 
[]

[Variables]
  [./disp_x]
     order = first
     family = lagrange
  [../]
  [./disp_y]
     order = first
     family = lagrange
  [../]
[]


[Physics/SolidMechanics/QuasiStatic]
  [./all]
    incremental = true
 
    generate_output = 'strain_xx'
    planar_formulation = PLANE_STRAIN
  [../]
[]
[AuxKernels]
 
  [w]
    type = MaterialRealAux
    variable = w
    property = damage_indexx
    execute_on = timestep_end
  []
 
   
  
[]

[BCs]
  [./bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = 't * 0.00001'
  [../]
[]
[Materials]
  
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.0
    youngs_modulus = 2.5e9
  [../]
 

  [Yin]
    type = GenericConstantMaterial
    prop_names = 'Yin'
    prop_values = 300.0
  []
  [P2]
   type = GenericConstantMaterial
   prop_names = 'P2'
   prop_values = 1.0
  []
  [P1]
   type = GenericConstantMaterial
   prop_names = 'P1'
   prop_values = 100.0
  []
  [damage_indexx]
    type = ComputeDamageIndexx
    
    Yin = Yin
    P1 = P1
    P2 = P2
    Mu = 20.0
  []
   [damage]
    type = ScalarMaterialDamage
    damage_index = damage_indexx
  []
   [stress]
  
    type = ComputeDamageStress
    damage_model = damage
  []
  
   

  
[]

[Executioner]
  type = Transient
  end_time = 100
  dt = 1

[]
[Postprocessors]
  [./disp_x1]
    type = PointValue
    point = '0.3 0.1 0'
    variable = disp_x
  [../]
  [./disp_x2]
    type = PointValue
    point = '0.3 0.0 0'
    variable = disp_x
  [../]
  
 
  [w]
    type = ElementAverageValue
    variable = w
  []
  
  
[]
  

[Outputs]

  exodus = true
  csv = true
  console = true

[]
