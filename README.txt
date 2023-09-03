SurfactantInterFoamPlus

# Author

- Olivier Huet
- Email: olivier.hueto@gmail.com

# Environment

- OpenFOAM-v1912

# Project Description

This project was developed as part of a PhD thesis titled "Development of improved mathematical models for droplet impaction on plant leaf surfaces." It enables the simulation of surfactant-laden drops impacting a solid surface within an axisymmetric domain.

# Key Components

- surfactantInterFoam
- transportModels


# surfactantInterFoam

SurfactantInterFoam is a modified version of the InterFoam solver for simulating two incompressible, isothermal immiscible fluids using a VOF method. The major modifications include:
- Coupling the VOF method with an LS method using the A-VOFCLS to improve interface accuracy.
- Implementing an adaptive coefficient of compression to enhance interface accuracy.
- Adding governing equations for surfactants.
- Splitting the sink-source terms of the surfactant equations into explicit and implicit terms for improved stability.
- Incorporating the interaction between surfactants and the substrate.

## Modified Files from interFoam

- surfactantInterFoam.C
- alphaEqn.H
- createFields.H
- pEqn.H
- UEqn.H

## Additional Files

- cellInfDist.H
- cellWidth.H
- createFieldsCA.H
- createFieldsDict.H
- createFieldsIP.H
- createFieldsLS.H
- createFieldsMesh.H
- createFieldsOptionsDrop.H
- createFieldsOptionsLS.H
- createFieldsOptionsSaa.H
- createFieldsReadLS.H
- createFieldsReadSaa.H
- createFieldsSaa.H
- Dirac.H
- gradPsi
- gradPsiDummy.H (not in use)
- gradPsiWall.H
- Gtot.H
- interfacePropertiesFields.H
- LSCorr.H
- meshLS.H
- psi0.H
- psiInf.H
- psiSeeding
- saaConcEqn.H
- saaConcInfArea.H
- saaSourceSink.H
- saaSourceSinkLink.H
- setInitialAlpha.H
- setInitialLS.H
- setInitialSaa.H
- setInitialU.H
- signFunc.H
- sorptionCourantNo.H

# transportModels

The transportModels module has been modified to incorporate Marangoni stresses due to the presence of surfactants.

## Modified Files from transportModels

- interfaceProperties.C
- interfaceProperties.H

## Additional Files

- saaConcDependentSurfaceTension.C (influence of surfactant on surface tension)
- saaConcDependentSurfaceTension.H (influence of surfactant on surface tension)

## Additional Files from source

- saaDynamicKistlerAlphaContactAngleFvPatchScalarField (Kistler contact angle model)
- saaDynamicKistlerAlphaContactAngleFvPatchScalarField.H (Kistler contact angle model)

# Others

The new solver needs material develloped by other OpenFOAM users: dynamicMesh, dynamicFvMesh and MakeAxialMesh that are available online.

#References and Links

- A-VOFCLS:  Haghshenas, M., Wilson, J.A., Kumar, R., 2017. Algebraic coupled level set-volume of fluid method for surface
- Adaptive coefficient of compression: Mehmani, Y., 2018. Wrinkle-free interface compression for two-fluid flows. arXiv preprint arXiv:1811.09744 .
- Governing equations for the surfactants: Antritter, T., Hachmann, P., Gambaryan-Roisman, T., Buck, B., Stephan, P., 2019. Spreading of micrometer-sized
droplets under the influence of insoluble and soluble surfactants: A numerical study. Colloids and Interfaces 3, 56.
- Dynamic mesh: Rettenmaier D, Deising D, Ouedraogo Y, Gjonaj E, De Gersem H, Bothe D, Tropea C, Marschall H. 2019 Load balanced 2D and 3D
adaptive mesh refinement in OpenFOAM. SoftwareX 10, 100317
- Kistler model: https://github.com/Swagga5aur/interThermalPhaseFoam/tree/master/Libraries/DynamicKistlerContactAngle
- MakeAxialMesh: https://openfoamwiki.net/index.php/Contrib/MakeAxialMesh