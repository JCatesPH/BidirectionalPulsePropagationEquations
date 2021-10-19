#include "BPPE.h"
#include "Structure.h"
#include "createLayers.h"
#include "physicalConstants.h"

//double sampleLayerThickness = 1.4e-6;
int numLayersInSample = 20;


void generateLayers(MaterialDB& theMaterialDB, Structure& theStructure)
{
	generatePlasmaTestMaterialsAndStructure(theMaterialDB, theStructure);
	//generateDefectMaterialsAndStructure(theMaterialDB, theStructure);
}

void generateLayerTestMaterialsAndStructure(MaterialDB &theMaterialDB,  Structure &theStructure)
{
	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	Material mat1("dieMat1", n0_Material1, n2_Material1, chi2_Material1, chi3_Material1);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);

	double aLayerThickness = 10*microns;

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
    theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), aLayerThickness, zStepMaterial1);
	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}

void generateApp1MaterialsAndStructure(MaterialDB &theMaterialDB,  Structure &theStructure)
{
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
	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	Material mat1("dieMat1", n0_Material1, n2_Material1, chi2_Material1, chi3_Material1);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);

	const double lamb0 = 630.0e-9; // Wavelength (MAKE SURE CORRESPONDS TO ACTUAL WAVELENGTH)
	const double thickness1 = n0_Material1 * lamb0 / 4.0;
	const double thickness2 = n0_Material2 * lamb0 / 4.0;

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
	double U_Ar = 15.759; // Ionization potential of Argon [eV]
	cout << sampleLayerThickness << endl;
	Material vacuum("Vacuum", 1.0, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);

	Material argon("Argon", n0_Argon, n2_Argon, chi2_Argon, chi3_Argon);
	//argon.setAsPlasmaMaterial(2, mpi_sigmaK, mpi_k);
	argon.setAsPlasmaMaterial(1, U_Ar);
	theMaterialDB.addMaterial(argon);

	Material plasmaMat("PlasmaMat", 1.0, 0.0, 0.0, 0.0);
	//plasmaMat.setAsPlasmaMaterial(2, mpi_sigmaK, mpi_k);
	//plasmaMat.setAsPlasmaMaterial(1, U_Ar);
	theMaterialDB.addMaterial(plasmaMat);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
	
	//theStructure.addLayer(theMaterialDB.getMaterialByName("Argon"), sampleLayerThickness, zStepMaterial1);
	//theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), sampleLayerThickness, zStepMaterial1);
	theStructure.addLayer(theMaterialDB.getMaterialByName("PlasmaMat"), sampleLayerThickness, zStepMaterial1);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}
