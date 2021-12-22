#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "gsl/gsl_math.h"
#include "Materials.h"

//#define VERBOSE_BOUNDARY 

// global functions defined in main.cpp file
void boundary(double z, complex<double>* k_0, complex<double>* k_1, double* y);

//Forward declarations
class Boundary;

//typedef struct {
//
//    double chi_2;
//    double chi_3;
//    double* omega, * kx, * rho;
//    complex<double>* k, * ee_p, * ee_m, * nl_k, * nl_p, * j_e;
//    fftw_plan nk_f, ep_b, em_b, np_f;
//
//}param_type;

class Layer
{
private:
    Material* m_material=NULL;
    size_t m_layerIDnum = 0;
    double m_thickness = 0.0; 
    double m_startZ = 0.0; 
    double m_endZ = 0.0; 
    double m_stepSize = 0.0; //propagation step in material
    int m_numStepsInLayer = 0;
    Boundary* lowSideBoundary = NULL;
    Boundary* hiSideBoundary = NULL;
    
public:
    Layer() {}

    Layer(Material* aMaterial, double aThickness, double aStepSize) : m_material(aMaterial), m_thickness(aThickness), m_stepSize(aStepSize) {
        m_numStepsInLayer = int(m_thickness / m_stepSize);
    }

    string getMaterialName() { return m_material->getName(); }
    Material& getMaterial() { return *m_material; }

    void setThickness(double aThickness) { m_thickness = aThickness; }
    double getThickness() { return m_thickness; }
    double getNumStepsInLayer() { return m_numStepsInLayer; }

    void setlayerIDnum(size_t aLayerIDnum) { m_layerIDnum = aLayerIDnum; }
    size_t getlayerIDnum() { return m_layerIDnum; }

    void setLowSideBoundary(Boundary* aBoundary) { lowSideBoundary = aBoundary; }
    Boundary* getLowSideBoundary() {
        #ifdef VERBOSE_BOUNDARY
        if (lowSideBoundary == NULL) {
            std::cout << " WARNING: the requested lowSideBoundary for layer#" << m_layerIDnum << " does not exist" << endl;
            //exit(-1);
        }
        #endif
        return lowSideBoundary;
    }
    void setHiSideBoundary(Boundary* aBoundary) { hiSideBoundary = aBoundary; }
    Boundary* getHiSideBoundary() {
        #ifdef VERBOSE_BOUNDARY
        if (hiSideBoundary == NULL) {
            std::cout << " WARNING: the requested hiSideBoundary for layer#" << m_layerIDnum << " does not exist" << endl;
            //exit(-1);
        }
        #endif
        return hiSideBoundary;
    }
   
    void setStartZpos(double aZpos) { m_startZ = aZpos; }
    double getStartZpos() { return m_startZ; }
    double getEndZpos() { return m_startZ+m_thickness; }

    void setStepSize(double aZstep) { m_stepSize = aZstep; }
    double getStepSize() { return m_stepSize; }
    
    // void setMaterial(Material aMaterial) { m_material = aMaterial; }

};


class Boundary
{
public:
    size_t m_BoundaryIDnum = 0;
    Layer* lowSideLayer;
    Layer* hiSideLayer;
    double m_zPos;

public:
    Boundary(Layer* aLowSideLayer, Layer* aHiSideLayer) : lowSideLayer(aLowSideLayer), hiSideLayer(aHiSideLayer) {
        if ((lowSideLayer->getStartZpos() + lowSideLayer->getThickness()) != hiSideLayer->getStartZpos()) {
            std::cout << " ERROR: invalid boundary number ::" << endl;
            exit(-1);
        }
        m_zPos = hiSideLayer->getStartZpos();
    };

    void setBoundaryIDnum(size_t aBoundaryIDnum) { m_BoundaryIDnum = aBoundaryIDnum; }
    size_t getBoundaryIDnum() { return m_BoundaryIDnum; }
};

class Structure
{
private:
    string m_structureName;

    list<Boundary> m_boundaries;
    int m_numLayers = 0;
    double m_totalThickness = 0;
public:
    list<Layer> m_layers;

public:

    // Structure() : m_structureName("Undefined Structure") {}
    Structure() {}

    Structure(string theStructureName) : m_structureName(theStructureName) {}

    size_t numLayers() { return m_layers.size(); }
    double getThickness() { return m_totalThickness; }

    double getZstepSizeAtZpos(double aZpos)
    {
        double foundZstepSize = std::numeric_limits<double>::quiet_NaN(); //  -666.666;
        for (Layer& aLayer : m_layers) {
            if (aZpos >= aLayer.getStartZpos() && aZpos <= aLayer.getEndZpos())
                foundZstepSize =   aLayer.getStepSize();
        }
        if (isnan(foundZstepSize)) {
            std::cout << " ERROR: the requested ZstepSize at zPosition " << aZpos << " does not exist (is aZpos outside domain?)" << endl;
            exit(-1);
        }
        return foundZstepSize;
    }


    void addLayer(Material* layerMaterial, double layerThickness, double aStepSize) {
        Layer aLayer(layerMaterial, layerThickness, aStepSize);
        addLayer(aLayer);
    }

    void addLayer(Layer aLayer) {
        aLayer.setStartZpos(m_totalThickness);
        aLayer.setlayerIDnum(m_layers.size());
        m_layers.push_back(aLayer);
        m_totalThickness += aLayer.getThickness();
        if (m_layers.size() > 1) {
            createNewBoundary();
        }
    }

    void createNewBoundary() {
        list<Layer>::reverse_iterator rit = m_layers.rbegin();
        Layer* RHSlayer, * LHSLayer;
        RHSlayer = &*rit;
        rit++;
        LHSLayer = &*rit;
        Boundary aBoundary(LHSLayer, RHSlayer);
        aBoundary.setBoundaryIDnum(m_boundaries.size());
        m_boundaries.push_back(aBoundary);
        
        // set the boundry pointers in the layers on either 
        // side of the new boundary
        list<Boundary>::reverse_iterator brit = m_boundaries.rbegin();
        Boundary* theNewBoundary;
        theNewBoundary = &*brit;
        LHSLayer->setHiSideBoundary(theNewBoundary);
        RHSlayer->setLowSideBoundary(theNewBoundary);
    }

    void doForwardPassThroughAllBoundaries(double* y) {
        if (m_boundaries.size() > 0) {
            for (std::list<Boundary>::iterator it = m_boundaries.begin(); it != m_boundaries.end(); ++it) {
                boundary(it->m_zPos, it->lowSideLayer->getMaterial().getK(), it->hiSideLayer->getMaterial().getK(), y);
            }
        }
    }

    void doBackwardPassThroughAllBoundaries(double* y) {
        if (m_boundaries.size() > 0) {
            for (std::list<Boundary>::reverse_iterator it = m_boundaries.rbegin(); it != m_boundaries.rend(); ++it) {
                boundary(it->m_zPos, it->hiSideLayer->getMaterial().getK(), it->lowSideLayer->getMaterial().getK(), y);
            }
        }
    }

    void doBoundaryUpdate(double aBoundaryzPos, double* y) {
        if (m_boundaries.size() > 0) {
            for (std::list<Boundary>::reverse_iterator it = m_boundaries.rbegin(); it != m_boundaries.rend(); ++it) {
                if(it->m_zPos == aBoundaryzPos)
                    boundary(it->m_zPos, it->hiSideLayer->getMaterial().getK(), it->lowSideLayer->getMaterial().getK(), y);
            }
        }
    }

    void writeStructureLayoutToASCIIFile(string aFilePath) {
        if (true) { printf("Writing Structure Layout file %s\n", aFilePath.c_str()); }
        ofstream outFile(aFilePath);
        outFile << "Structure Name:\t" << m_structureName << endl;
        outFile << "Total Thickness:\t" << m_totalThickness * m2um << " [um]" << endl;
        outFile << "Number of layers:\t" << m_layers.size() << endl;
        outFile << endl;
        outFile << setw(8) << "Layer#" << setw(12) << "zStart[um]" << setw(15) << "Thickness[um]" << setw(15) << "stepSize[um]" << setw(15) << "zEnd[um]" << "    " << "Material" << endl;
        outFile << "-------------------------------------------------------------------------------------------------" << endl;
        for (Layer& aLayer : m_layers) {
            outFile << setw(5) << aLayer.getlayerIDnum();
            outFile << setw(15) << aLayer.getStartZpos() * m2um;
            outFile << setw(15) << aLayer.getThickness() * m2um;
            outFile << setw(15) << aLayer.getStepSize() * m2um;
            outFile << setw(15) << aLayer.getEndZpos() * m2um;
            outFile << "    " << aLayer.getMaterial().serializedToString();
            if (aLayer.getlayerIDnum() != 0)
                outFile << " LHSBoundaryID= " << aLayer.getLowSideBoundary()->getBoundaryIDnum();
            if (aLayer.getlayerIDnum() != (m_layers.size()-1))
                outFile << " RHSBoundaryID= " << aLayer.getHiSideBoundary()->getBoundaryIDnum();
            outFile << endl;
        }
        outFile.close();
    }

    void writeStructureToDATFile(string aFilePath) {
        if (true) { printf("Writing Structure to DATA file %s\n", aFilePath.c_str()); }
        ofstream outFile(aFilePath);
        outFile << "# z_pos\tMaterialID\tn0 []\tn2 []\tchi_2 []\tchi_3 []" << endl;
        outFile << std::scientific;
        for (Layer& aLayer : m_layers) {
           double theZpos = aLayer.getStartZpos();
           for (int i = 0; i < aLayer.getNumStepsInLayer(); i++) {
                outFile << theZpos << "\t";
                outFile << aLayer.getMaterial().getidNum() << "\t";
                outFile << aLayer.getMaterial().getn0() << "\t";
                outFile << aLayer.getMaterial().getn2() << "\t";
                outFile << aLayer.getMaterial().getChi2() << "\t";
                outFile << aLayer.getMaterial().getChi3() << "\t";
                outFile << endl;
                theZpos += aLayer.getStepSize();              
            }
        }
        outFile.close();
    }

    void writeBoundaryLayoutToASCIIFile(string aFilePath) {
        if (true) { printf("Writing Structure Layout file %s\n", aFilePath.c_str()); }
        ofstream outFile(aFilePath);
        outFile << "Structure Name:\t" << m_structureName << endl;
        outFile << "Number of Boundaries:\t" << (m_boundaries.size()) << endl;
        outFile << endl;
        outFile << setw(12) << "boundary#" << setw(12) << "LHS Layer" << setw(12) << "zPos" << setw(12) << "RHS Layer" << endl;
        outFile << "-------------------------------------------------------------------" << endl;
        for (Boundary& aBoundary : m_boundaries) {
            outFile << setw(12) << aBoundary.getBoundaryIDnum();
            outFile << setw(12) << aBoundary.lowSideLayer->getMaterialName();
            outFile << setw(12) << aBoundary.m_zPos;
            outFile << setw(12) << aBoundary.hiSideLayer->getMaterialName();
            outFile << endl;
        }
        outFile.close();
    }
};

class Simulation
{
public:
    string m_simulationName;
    MaterialDB* m_materialDatabase = NULL;
    Structure* m_structure = NULL;
    const double m_domain_t = 1500e-15; //time domain
    const double m_domain_x = 125.0e-6; //x domain
    const int m_numt = int(pow(2, 14)); //number of time points
    const int m_numActiveOmega = m_numt / 2 + 1; //num_x * num_t / 2 + 2;
    complex<double>* eFieldPlus = NULL;
    complex<double>* eFieldMinus = NULL;


public:
    Simulation() { }
    Simulation(string theSimulationName) : m_simulationName(theSimulationName) { }

    void setMaterialDatabase(MaterialDB* aMaterialDatabase) { m_materialDatabase = aMaterialDatabase; }
    MaterialDB* getMaterialDatabase() { return m_materialDatabase; }
    void setStructure(Structure* aStructure) { m_structure = aStructure; }
    Structure* getStructure() { return m_structure; }
};
