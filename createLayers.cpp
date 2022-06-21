#include "BPPE.h"
#include "Structure.h"
#include "createLayers.h"
#include "physicalConstants.h"

void generateLayers(MaterialDB& theMaterialDB, Structure& theStructure)
{
	generatePlasmaTestMaterialsAndStructure(theMaterialDB, theStructure);
	//generateSilicaMaterialsAndStructure(theMaterialDB, theStructure);
	//generateDefectMaterialsAndStructure(theMaterialDB, theStructure);
	//generateLayerTestMaterialsAndStructure(theMaterialDB,  theStructure);
}

void generateLayerTestMaterialsAndStructure(MaterialDB &theMaterialDB,  Structure &theStructure)
{

	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	//Material argon("Argon", n0_Argon, n2_Argon, chi2_Argon, chi3_Argon);
	//theMaterialDB.addMaterial(argon);
	Material mat1("dieMat1", n0_Material1, n2_Material1, chi2_Material1, chi3_Material1);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);
	

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
    theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), sampleLayerThickness, zStepMaterial1);
	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}

void generateSilicaMaterialsAndStructure(MaterialDB &theMaterialDB,  Structure &theStructure)
{
	// Indices of refraction (quoted in 10.1103/PhysRevLett.110.115003)
	const double n0_Silica = 1.4395;
	const double n2_Silica = 2.3e-20; // [m^2 / W]
	const double chi3_Silica = (4 / 3) * epsilon_0 * cLight * pow(n0_Silica, 2) * n2_Silica;
	const double chi2_Silica = 0.0;

	// Parameters for ionization in Durand paper
	//const double mpi_sigmaK = 1.0e-224; // [m^(28) W^(-14) s^(-1)] MPI cross section
	//const double mpi_k = 14.0; // Photons for MPI
	//const double sigmaBremsstrahlung = 4.7e-25; // [m^2] Bremsstrahlung cross section
	//const double silicaUi = 1.44e-18; // [J] = 9 [eV] Ionization potential of silica
	//const double recombTime = 150e-15; // [s] Recombination time

	// -- (FAKE) MPI parameters
	const double mpi_sigmaK = 1.0e-23;
	const double mpi_k = 3.0;

	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);

	Material silica("Silica", n0_Silica, n2_Silica, chi2_Silica, chi3_Silica);
	//silica.setAsPlasmaMaterial(3, mpi_sigmaK, mpi_k, sigmaBremsstrahlung, silicaUi, recombTime); 
	silica.setAsPlasmaMaterial(2, mpi_sigmaK, mpi_k); // MPI ONLY
	theMaterialDB.addMaterial(silica);

	Material silicaLin("SilicaLinear", n0_Silica, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(silicaLin);
	

	//theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
    //theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
	theStructure.addLayer(theMaterialDB.getMaterialByName("SilicaLinear"), LHSsourceLayerThickness, zStepMaterial1);
	theStructure.addLayer(theMaterialDB.getMaterialByName("Silica"), sampleLayerThickness, zStepMaterial1);
	theStructure.addLayer(theMaterialDB.getMaterialByName("SilicaLinear"), RHSbufferLayerThickness, zStepMaterial1);
	//theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}

void generateApp1MaterialsAndStructure(MaterialDB &theMaterialDB,  Structure &theStructure)
{

	int numLayersInSample = 20;

	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	Material mat1("dieMat1", n0_Material1, n2_Material1, chi2_Material1, chi3_Material1);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);
	Material mat3("linMieMat3", n0_Material2, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(mat3);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
	for (int i = 0; i < (numLayersInSample/2); i++)
	{
		//theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), sampleLayerThickness, zStepMaterial1);
		//theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
	}
	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}

void generateDefectMaterialsAndStructure(MaterialDB& theMaterialDB, Structure& theStructure)
{
	int numLayersInSample = 20;

	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	Material mat1("dieMat1", n0_Material1, n2_Material1, chi2_Material1, chi3_Material1);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);
	Material argon("Argon", n0_Argon, n2_Argon, chi2_Argon, chi3_Argon);
	theMaterialDB.addMaterial(argon);

	const double lamb0 = 630e-9; // Wavelength (MAKE SURE CORRESPONDS TO ACTUAL WAVELENGTH)
	const double thickness1 = lamb0 / (4.0 * n0_Material1);
	const double thickness2 = lamb0 / (4.0 * n0_Material2);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
	for (int i = 0; i < (numLayersInSample / 4); i++)
	{
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), thickness1, zStepMaterial1);
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), thickness2, zStepMaterial1);
	}

	//theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), sampleLayerThickness + sampleLayerThickness/3.0, zStepMaterial1);

	for (int i = 0; i < (numLayersInSample / 4); i++)
	{
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), thickness1, zStepMaterial1);
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), thickness2, zStepMaterial1);
		//theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), thickness2, zStepMaterial1);
	}
	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}

void generatePlasmaTestMaterialsAndStructure(MaterialDB& theMaterialDB, Structure& theStructure)
{
	const double mpi_sigmaK = 3.4e-128;
	const double mpi_k = 8.0;

	double U_Ar = 15.759; // Ionization potential of Argon [eV]
	//cout << sampleLayerThickness << endl;
	Material vacuum("Vacuum", 1.0, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);

	//Material argon("Argon", n0_Argon, n2_Argon, chi2_Argon, chi3_Argon);
	Material argon("Argon", n0_Argon, 0.0, 0.0, 0.0);
	//argon.setAsPlasmaMaterial(2, mpi_sigmaK, mpi_k);
	//argon.setAsPlasmaMaterial(1, U_Ar);
	theMaterialDB.addMaterial(argon);

	Material plasmaMat("PlasmaMat", 1.0, 0.0, 0.0, 0.0);
	//plasmaMat.setAsPlasmaMaterial(2, mpi_sigmaK, mpi_k);
	//plasmaMat.setAsPlasmaMaterial(1, U_Ar);
	theMaterialDB.addMaterial(plasmaMat);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
	
	theStructure.addLayer(theMaterialDB.getMaterialByName("Argon"), sampleLayerThickness, zStepMaterial1);
	//theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), sampleLayerThickness, zStepMaterial1);
	//theStructure.addLayer(theMaterialDB.getMaterialByName("PlasmaMat"), sampleLayerThickness, zStepMaterial1);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}
