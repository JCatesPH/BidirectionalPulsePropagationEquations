#include "BPPE.h"
#include "Structure.h"

double sampleLayerThickness = 5e-6;
int numLayersInSample = 5;

void generateLayerTestMaterialsAndStructure(MaterialDB &theMaterialDB,  Structure &theStructure)
{
	Material vacuum("Vacuum", n0_Vacuum, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	Material mat1("dieMat1", n0_Material1, n2_Material1, chi2_Material1, chi3_Material1);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
	for (int i = 0; i < numLayersInSample; i++)
	{
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
	}
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

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);
	for (int i = 0; i < (numLayersInSample / 4); i++)
	{
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), sampleLayerThickness, zStepMaterial1);
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
	}

	theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), sampleLayerThickness + sampleLayerThickness/3.0, zStepMaterial1);

	for (int i = 0; i < (numLayersInSample / 4); i++)
	{
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat2"), sampleLayerThickness, zStepMaterial1);
		theStructure.addLayer(theMaterialDB.getMaterialByName("dieMat1"), sampleLayerThickness, zStepMaterial1);
	}
	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}

void generatePlasmaTestMaterialsAndStructure(MaterialDB& theMaterialDB, Structure& theStructure)
{
	Material vacuum("Vacuum", 1.0, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(vacuum);
	Material mat1("dieMat1", 1.75, 0.0, 0.0, 0.0);
	theMaterialDB.addMaterial(mat1);
	Material mat2("dieMat2", n0_Material2, n2_Material2, chi2_Material2, chi3_Material2);
	theMaterialDB.addMaterial(mat2);

	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), LHSsourceLayerThickness, zStepMaterial1);

	//theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), sampleLayerThickness * numLayersInSample/4, zStepMaterial1);
	for (int i = 0; i < numLayersInSample; i++)
	{
		theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), sampleLayerThickness, zStepMaterial1);
	}
	theStructure.addLayer(theMaterialDB.getMaterialByName("Vacuum"), RHSbufferLayerThickness, zStepMaterial1);
}