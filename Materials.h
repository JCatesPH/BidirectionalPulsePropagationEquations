#pragma once
#include <string>
#include <complex>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm> // Necessary for find_if(...) function
#include "physicalConstants.h"
//#include "Structure.h"

using namespace std;

class Material
{
private:
	string m_name;
	size_t m_idNum;
	string m_matInfoAsString;
	double m_n0 = 0.0;
	double m_n2;
	double m_chi_2;
	double m_chi_3;
	int m_doPlasmaCalc = 0;
	double m_mpi_sigmaK = 0.0, m_mpi_k =0.0, m_ionE = 0.0;

	//complex<double>* m_k=nullptr;
	int m_numActiveOmega = 0;
public:
	Material() : m_name("Undefined Material"), m_n0(0.0), m_n2(0.0), m_chi_2(666.666), m_chi_3(333.333) {}
	Material(string aName, double n0, double n2, double chi_2, double chi_3) : m_name(aName), m_n0(n0), m_n2(n2), m_chi_2(chi_2), m_chi_3(chi_3) {
		std::stringstream ss;
		ss << m_name;
		ss << "  [ n0=" << m_n0;
		ss << " n2=" << m_n2;
		ss << " Chi2=" << m_chi_2;
		ss << " Chi3=" << m_chi_3;
		ss << " ]";
		m_matInfoAsString.append(ss.str());
	}

	complex<double>* m_k=nullptr;

	void setName(string aName) { m_name = aName; }
	string getName() { return m_name; }
	void setidNum(size_t aidNum) { m_idNum = aidNum; }
	size_t getidNum() { return m_idNum; }
	string serializedToString() { return m_matInfoAsString; }

	double getn0() { return m_n0; }
	double getn2() { return m_n2; }
	double getChi2() { return m_chi_2; }
	double getChi3() { return m_chi_3; }
	double getdoPlasmaCalc() { return m_doPlasmaCalc; }
	double getmpi_sigmaK() { return m_mpi_sigmaK; }
	double getmpi_k() { return m_mpi_k; }
	double getIonizationEnergy() {return m_ionE; }

	void setAsPlasmaMaterial(int doPlasmaCalc, double mpi_sigmaK, double mpi_k) {
		m_doPlasmaCalc = doPlasmaCalc;
		m_mpi_sigmaK = mpi_sigmaK;
		m_mpi_k = mpi_k;
		std::stringstream ss;
		ss << m_name;
		ss << "  [ doPlasmaCalc=" << m_doPlasmaCalc;
		ss << " mpi_sigmaK=" << m_mpi_sigmaK;
		ss << " mpi_k=" << m_mpi_k;
		ss << " ]";
		m_matInfoAsString.append(ss.str());
	}
	void setAsPlasmaMaterial(int doPlasmaCalc, double ionEnergy) {
		m_doPlasmaCalc = doPlasmaCalc;
		m_ionE = ionEnergy;
		std::stringstream ss;
		ss << m_name;
		ss << "  [ doPlasmaCalc=" << m_doPlasmaCalc;
		ss << " Ionization Energy =" << m_ionE;
		ss << " ]";
		m_matInfoAsString.append(ss.str());
	}
	//void getMaterialInformationAsString(string& aString) {
	//	char materialAsCharacters[200];
	//	sprintf_s(materialAsCharacters, "%s\t[n0=%g\tn2=%g\tchi2=%g\tk[1]=%g", m_name.c_str(), m_n0, m_n2, m_chi_2, real(m_k[1]));
	//	// materialAsCharacters << m_name << "     [n0=" << m_n0 << " n2=" << m_n2 << " chi2=" << m_chi_2 << "k[1]=" << real(m_k[1]);
	//	aString.append(materialAsCharacters);
	//	return;
	//}

	complex<double> getIndex(double omg, double kx) {
		return m_n0;
	}

	complex<double>* getK() {return m_k; }
		
	complex<double>* mallocK(int theNumActiveOmega) {
		m_numActiveOmega = theNumActiveOmega;
		if (m_k == nullptr) {
			m_k = (complex<double>*)malloc(sizeof(complex<double>) * m_numActiveOmega);
			if (!m_k) {
				cout << "Memory Allocation of K-array Failed for material " << m_name << endl;
				exit(1);
			}
		}
		else {
			cout << "Attempt to re-allocate memory for K-array for material " << m_name << endl;
			// exit(1);
		}
		return m_k;
	}

	void fillK(double* omg) {
		for (int i = 0; i <m_numActiveOmega; i++)
		{
			if (i == 0) 
				m_k[i] = 0.0;
			else 
				m_k[i] = omg[i] * getIndex(omg[i], 0.0) / cLight;
		}
		//m_matInfoAsString.append("    k[1] = "); m_matInfoAsString.append(to_string(real(m_k[1])));
	}

};

class MaterialDB
{
private:
	string m_databaseName;
	list<Material> m_materials;

public:
	MaterialDB() : m_databaseName("Undefined Material") {}
	MaterialDB(string theDatabaseName) { m_databaseName = theDatabaseName; }

	void addMaterial(Material aMaterial) {
		aMaterial.setidNum(m_materials.size());
		m_materials.push_back(aMaterial);
	}

	Material* getMaterialByName(string aMatName) {
		auto it = find_if(m_materials.begin(), m_materials.end(), [&aMatName](Material& obj) {return obj.getName() == aMatName; });
		if (it == m_materials.end()) {
			std::cout << " ERROR: Item with name ::" << aMatName << " not Found" << std::endl;
			exit(-1);
		}
		return &*it;
	}

	void writeMaterialDBToASCIIFile(string aFilePath) {
		if (true) { cout << "Writing MaterialDB to file " << aFilePath << endl; }
		ofstream outFile(aFilePath);
		outFile << "Material Database Name::\t" << m_databaseName << endl;
		outFile << "Number of Materials:\t" << m_materials.size() << endl;
		outFile << endl;
		outFile << setw(12) << "MaterialName" << setw(25) << "Material Info" << endl;
		outFile << "-------------------------------------------------------------------" << endl;
		for (Material& aMaterial : m_materials) {	
			outFile << setw(25) << aMaterial.serializedToString() << endl;
		}
		outFile.close();
}

	void initAllMaterialKs(double* omg, int numOmega) {
		for (Material& aMaterial : m_materials) {
			aMaterial.mallocK(numOmega);
			aMaterial.fillK(omg);
			if (false) 	cout << "   Filled k in material:" << aMaterial.serializedToString() << endl;
		}
	}
		
};

