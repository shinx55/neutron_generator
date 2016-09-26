/*
	simulate.c : This is a simulation program for the reactor of neutron generator.
    Copyright (C) 2016 by SHINEUKE ONO, JPN (given sur, country) <shinx55@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
	You can see the idea of reactor of neutron generator in
	https://docs.google.com/document/d/1ONZLusSrrTWbh4wBKiO8K6TwLG2wExtH-jtinQRzBt4/edit?usp=sharing
	https://docs.google.com/document/d/1ZmPn4N57MOAG2C02d_nFATG-t7OUhTOJm-7tQdR2cZA/edit?usp=sharing
	
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include "hashtable.h"
#include "timeformat.h"
//---------------------------------------------------------------------
//[Physics constants]
#define NAvogadro 6.022140857E23 //Avogadro constant 
double e_energyJouleOfeV = 1.6021766208E-19; //[J/eV]
double e_energyJouleOfCalorie = 4.184; //[J/cal]
double e_energyJouleOfKWH = 3.6E6; //[J/KWH]
double e_elementaryCharge = 1.6021766208E-19; //[C]The elementary charge is the electric charge of an electron.
double e_vacuumPermittivity = 8.854187817e-12;//[F/m] The physical constant ε0, commonly called the vacuum permittivity
double e_unifiedAtomicMassUnitKg = 1.660538921E-27; //[kg], unified Atomic Mass Unit 1[u]/[dalton, Da], the mass of 1/12 of mass of an atom of 12C(carbon 12) 
double e_massNeutronKg = 1.674927471E-27; //[kg], the mass of a neutron
double e_massNeutronMeV = 939.5654133;//[MeV/C^2], the mass of a neutron
double e_massProtonKg = 1.672621898E-27;//[kg], the mass of a proton
double e_massProtonMeV = 938.2720813;//[MeV/C^2], the mass of a proton
double e_massElectronKg = 9.10938356E-31; //[kg], the mass of an electron
double e_massElectronMeV = 0.5109989461;//[MeV/C^2], the mass of an electron
//[calcuated Physics constants]
double e_coefElectroPotentialMeV;//(e_elementaryCharge * e_elementaryCharge) / (4.0 * M_PI * e_vacuumPermittivity * e_energyJouleOfeV * 1.0E6)
double e_coefMassUToMeV;
double e_coefMeVtoMassU;
double e_massNeutronMassU;//[u], the mass of a neutron
double e_massProtonMassU;//[u], the mass of a proton
double e_massElectronMassU;//[u], the mass of an electron
double e_betaEnergyMassU;//
double e_betaEnergyMeV;// the kinetic energy that is released in beta decay, (It is random split into electrons and the anti-electron neutrino)
//---------------------------------------------------------------------
//[User conditions]
int e_rangeDecimalDigitPrecision;//2 - 6;//The precision expressed by the number of dicimal digit, 6 means that the tolelance rate is 10^-6 = 0.000001, 2 means 10^-2 = 0.01. The function "registOutput" will treat that a gamma ray is almost equal to another gamma ray with in the tolelance specified by "e_rangeDecimalDigitPrecision". Why 2 and not 6? If you use 6, there are too many grades of strength of gamma, the program need very large memories and take long time of the computation of simulation.
double e_detectLimitRateForIsotope; 

//The negative electrode will emmit electrons, 
//The positive electrode will emmit protons, becuase positive electrode absorbs hydrogen.
//The following variables show the numbers of moles of the elements in the electrodes.
#define TEXT_LEN_OF_MOLS 256
char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS];
char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS];
//There is a negative electrode over the positive electrode in the neutron generator.
//upper[negative Electrode]==>(electron)==>|collide|<==(proton)<==[positive Electrode]lower

double e_appliedVoltageScale;
double e_emittedElectronMol;//The number of moles of electrons emitted per one(1) second from the negative electrode by high-voltage pulse like the Cockcroft-Walton voltage multiplier circuit.

double e_collideElectronRateOnElectrode;//Probability to collide between (the electron emitted from the negative electrode) and (the nucleus in the positive electrode), the electron has large kinetic energy than the energy of beta decay, the electron emitted from the negative electrode, that reaches on the positive electrode
double e_collideElectronRateForMidiMeV;//This is paired with 'e_collideElectronMidiMeV'
double e_collideElectronRateForMiniMeV;//This is paired with 'e_collideElectronMiniMeV'
double e_emittedProtonMol;//The number of moles of protons emitted per one(1) second from the by positive electrode by high-voltage pulse like the Cockcroft-Walton voltage multiplier circuit.
double e_neutronGenInSpaceProtonRate;//The probability applied for the proton that is emitted from the positive electrode and will collide with the electron emitted from the negative electrode in the intermediate space of electrodes and will transform into a neutron with enough kinetic energy greater than the energy of beta decay.
double e_neutronGenInSpaceFallRate;//Almost all of the neutrons generated in the intermediate space fall on the positive electrode. The neutrons have smaller energy(= e_appliedDefectMeV).
double e_hydrogenGenInSpaceProtonRate;//The probability applied for the proton that is emitted from the positive electrode and will collide with the electron emitted from the negative electrode in the intermediate space of electrodes and will transform into a hydrogen with small kinetic energy less than the energy of beta decay.
double e_collideProtonRateOnElectrode;//Probability to collide between the proton and ((the nucleus or the electron) in the negative electrode), the proton has large kinetic energy than the energy of beta decay, the proton emitted from the positive electrode reach the negative electrode
double e_collideElectronMidiMeV;//This is almost half of e_betaEnergyMeV = 0.78 [MeV]
double e_collideElectronMiniMeV; //10[keV] is 1.16e+8 [K]. It is almost as same as D-T fusion reactor. It is about 100 times of the heat of the corona of the Sun, 1.0e+6 [K]. It is about 1,000 times of the energy to ionize a hydrogen, 13.6 [eV].
int e_usebulletCrossSection;// 1 (dafault) or 0
int e_useComptonEffect;// 1 (dafault) or 0

extern void initUserConditionsByDefault()
{
	e_rangeDecimalDigitPrecision = 2;
	e_detectLimitRateForIsotope = 1.0E-15; 
	strcpy(e_negativeElectrodeAtomicMols, "Ni=1.0");
	strcpy(e_positiveElectrodeAtomicMols, "H=0.6, Ni=1.0");
	e_appliedVoltageScale = 1.01;
	e_emittedElectronMol = 1.0E-10;
	e_collideElectronRateOnElectrode = 0.2;
	e_collideElectronRateForMidiMeV = 0.1;
	e_collideElectronRateForMiniMeV = 0.01;
	e_emittedProtonMol = 0.4E-10;
	e_neutronGenInSpaceProtonRate = 0.01;
	e_neutronGenInSpaceFallRate = 0.95;
	e_hydrogenGenInSpaceProtonRate = 0.01;
	e_collideProtonRateOnElectrode = 0.2;
	e_collideElectronMidiMeV = 0.36;
	e_collideElectronMiniMeV = 0.01;
	e_usebulletCrossSection = 1;
	e_useComptonEffect = 1;
}
//---------------------------------------------------------------------
//[calcuated User conditions]
int e_stepsOfLostingEnergy;//When un-perfect collision, the particles will lost their energy step by step. We simulate it only two steps.
double e_appliedVoltageMassU;
double e_appliedVoltageMeV;
double e_appliedDefectMeV;

double e_stepByStepLostingEnergyMeV;

double e_neutronGenInSpaceFlyRate;
double e_genedNeutronInSpaceMol;
double e_genedHydrogenInSpaceMol;
double e_arrivedProtonMol;
double e_arrivedProtonImperfectCollideMol;
double e_arrivedProtonRate;

double e_arrivedElectronMol;
double e_arrivedElectronRate;
double e_neutronGenInSpaceElectronRate;
double e_hydrogenGenInSpaceElectronRate;

double e_emittedPulseEnergyMeV;
double e_emittedPulseEnergyJ;
double e_emittedPulseEnergyCal;
double e_emittedPulseEnergyKW;
double e_electricCurrent;

//double e_lostHeatOfGenedNeutronnInSpaceMeV;
//double e_lostHeatOfGenedHydrogenInSpaceMeV;
//double e_lostHeatInSpaceMeV;
//double e_lostHeatByGenedNeutronnIna_rwSpaceRate;
//double e_lostHeatByGenedHydrogenInSpaceRate;
//double e_lostHeatInSpaceRate;
//---------------------------------------------------------------------
#define VERSION_CHECK(VERSION_NUMBER) \
		if(version != VERSION_NUMBER){\
			fprintf(stderr, "FATAL ERROR:%s:version mismatch\n", __FUNCTION__);\
			exit(1);\
		}\

#define VERSION_CHECK2(VERSION_NUMBER, OLD_VERSION_NUMBER) \
		if(version != VERSION_NUMBER && version != OLD_VERSION_NUMBER){\
			fprintf(stderr, "FATAL ERROR:%s:version mismatch\n", __FUNCTION__);\
			exit(1);\
		}\

#define FREAD_VALUE_NULL(VALUE, FP) \
	if(fread(&VALUE, sizeof(VALUE), 1, FP) != 1){return NULL;}

#define FREAD_MEMBER_FREE(PTR, MEMBER, FP) \
	if(fread(&PTR->MEMBER, sizeof(PTR->MEMBER), 1, FP) != 1){free(PTR); return NULL;}

#define FREAD_VERSION_CHECK(VERSION_NUMBER, FP)\
		int version;\
		FREAD_VALUE_NULL(version, FP)\
		VERSION_CHECK(VERSION_NUMBER)\

#define FWRITE_VALUE(VALUE, FP) \
	if(fwrite(&VALUE, sizeof(VALUE), 1, FP) != 1){return 0;}

#define FWRITE_VERSION(VERSION_NUMBER, FP)\
		int version = VERSION_NUMBER;\
		FWRITE_VALUE(version, FP)\

#define SERIALIZE_VALUE(RW, VALUE, FP) \
	if(RW(&VALUE, sizeof(VALUE), 1, FP) != 1){return 0;}

#define SERIALIZE_ARRAY(RW, ARRAY_PTR, TYPE, COUNT, FP) \
	if(COUNT > 0){\
		if(RW(ARRAY_PTR, sizeof(TYPE), COUNT, FP) != COUNT){return 0;}\
	}

#define SERIALIZE_VERSION_CHECK(RW, VERSION_NUMBER, FP) \
		int version = VERSION_NUMBER;\
		SERIALIZE_VALUE(RW, version, FP)\
		if((void *)RW == (void *)fread){\
			VERSION_CHECK(VERSION_NUMBER) \
		}

#define SERIALIZE_VERSION_CHECK2(RW, VERSION_NUMBER, OLD_VERSION_NUMBER, FP) \
		int version = VERSION_NUMBER;\
		SERIALIZE_VALUE(RW, version, FP)\
		if((void *)RW == (void *)fread){\
			VERSION_CHECK2(VERSION_NUMBER, OLD_VERSION_NUMBER) \
		}

//---------------------------------------------------------------------
#define OLD_VERSION_OF_serializeUserConditions 1
#define VERSION_OF_serializeUserConditions 2

#define SERIALIZE_USER_CONDITIONS(RW) \
extern int RW ## UserConditions(FILE * a_fp) \
{\
	SERIALIZE_VERSION_CHECK2(RW, VERSION_OF_serializeUserConditions, OLD_VERSION_OF_serializeUserConditions, a_fp) \
	SERIALIZE_VALUE(RW, e_rangeDecimalDigitPrecision, a_fp)\
	SERIALIZE_VALUE(RW, e_detectLimitRateForIsotope, a_fp)\
	if((void *)RW == (void *)fread && version == OLD_VERSION_OF_serializeUserConditions){\
		double negativeElectrodeHydrogenMol;\
		double negativeElectrodeNickelMol;\
		double positiveElectrodeHydrogenMol;\
		double positiveElectrodeNickelMol;\
		SERIALIZE_VALUE(RW, negativeElectrodeHydrogenMol, a_fp)\
		SERIALIZE_VALUE(RW, negativeElectrodeNickelMol, a_fp)\
		SERIALIZE_VALUE(RW, positiveElectrodeHydrogenMol, a_fp)\
		SERIALIZE_VALUE(RW, positiveElectrodeNickelMol, a_fp)\
		sprintf(e_negativeElectrodeAtomicMols, "H=%lg, Ni=%lg", negativeElectrodeHydrogenMol, negativeElectrodeNickelMol);\
		sprintf(e_positiveElectrodeAtomicMols, "H=%lg, Ni=%lg", positiveElectrodeHydrogenMol, positiveElectrodeNickelMol);\
	}else{\
		SERIALIZE_ARRAY(RW, e_negativeElectrodeAtomicMols, char, TEXT_LEN_OF_MOLS, a_fp)\
		SERIALIZE_ARRAY(RW, e_positiveElectrodeAtomicMols, char, TEXT_LEN_OF_MOLS, a_fp)\
	}\
	SERIALIZE_VALUE(RW, e_appliedVoltageScale, a_fp)\
	SERIALIZE_VALUE(RW, e_emittedElectronMol, a_fp)\
	SERIALIZE_VALUE(RW, e_collideElectronRateOnElectrode, a_fp)\
	SERIALIZE_VALUE(RW, e_collideElectronRateForMidiMeV, a_fp)\
	SERIALIZE_VALUE(RW, e_collideElectronRateForMiniMeV, a_fp)\
	SERIALIZE_VALUE(RW, e_emittedProtonMol, a_fp)\
	SERIALIZE_VALUE(RW, e_neutronGenInSpaceProtonRate, a_fp)\
	SERIALIZE_VALUE(RW, e_neutronGenInSpaceFallRate, a_fp)\
	SERIALIZE_VALUE(RW, e_hydrogenGenInSpaceProtonRate, a_fp)\
	SERIALIZE_VALUE(RW, e_collideProtonRateOnElectrode, a_fp)\
	SERIALIZE_VALUE(RW, e_collideElectronMidiMeV, a_fp)\
	SERIALIZE_VALUE(RW, e_collideElectronMiniMeV, a_fp)\
	SERIALIZE_VALUE(RW, e_usebulletCrossSection, a_fp)\
	SERIALIZE_VALUE(RW, e_useComptonEffect, a_fp)\
	return 1;\
}
#define CALL_SERIALIZE_USER_CONDITIONS(RW, FP) \
	RW ## UserConditions(FP)
SERIALIZE_USER_CONDITIONS(fread);//freadUserConditions
SERIALIZE_USER_CONDITIONS(fwrite);//fwriteUserConditions

extern void initiallyCalculatePhysicsContants()
{
	//[calculated Physics constants]
	e_coefElectroPotentialMeV = (e_elementaryCharge * e_elementaryCharge) / (4.0 * M_PI * e_vacuumPermittivity * e_energyJouleOfeV * 1.0E6);
	e_coefMassUToMeV = e_massNeutronMeV * ( e_unifiedAtomicMassUnitKg / e_massNeutronKg);
	e_coefMeVtoMassU = 1 / e_coefMassUToMeV;
	e_massNeutronMassU = e_massNeutronMeV * e_coefMeVtoMassU;
	e_massProtonMassU = e_massProtonMeV * e_coefMeVtoMassU;
	e_massElectronMassU = e_massElectronMeV * e_coefMeVtoMassU;
	e_betaEnergyMassU = e_massNeutronMassU - (e_massProtonMassU + e_massElectronMassU);
	e_betaEnergyMeV = e_massNeutronMeV - (e_massProtonMeV + e_massElectronMeV);	
}

extern void initiallyCalculateUserContants()
{
	//[calcuated User conditions]
	if(e_appliedVoltageScale >= 30.0){
		fprintf(stderr,"FATAL ERROR:%s:too large e_appliedVoltageScale:%lg\n", __FUNCTION__, e_appliedVoltageScale);
		exit(1);
	}
	if(e_appliedVoltageScale < 0.05){
		fprintf(stderr,"FATAL ERROR:%s:too small e_appliedVoltageScale:%lg\n", __FUNCTION__, e_appliedVoltageScale);
		exit(1);
	}
	e_stepsOfLostingEnergy = (int)(e_appliedVoltageScale * 2.0);
	e_appliedVoltageMassU = e_betaEnergyMassU * e_appliedVoltageScale;
	e_appliedVoltageMeV = e_betaEnergyMeV * e_appliedVoltageScale;
	if(e_appliedVoltageMeV < e_collideElectronMidiMeV){
		fprintf(stderr,"FATAL ERROR:%s:e_appliedVoltageMeV:%lg < e_collideElectronMidiMeV:%lg\n", __FUNCTION__, e_appliedVoltageMeV, e_collideElectronMidiMeV);
		exit(1);
	}
	if(e_collideElectronMiniMeV > e_collideElectronMidiMeV){
		fprintf(stderr,"FATAL ERROR:%s:e_collideElectronMiniMeV:%lg < e_collideElectronMidiMeV:%lg\n", __FUNCTION__, e_collideElectronMiniMeV, e_collideElectronMidiMeV);
		exit(1);
	}
	e_appliedDefectMeV = e_appliedVoltageMeV - e_betaEnergyMeV;
	e_stepByStepLostingEnergyMeV = e_appliedVoltageMeV / e_stepsOfLostingEnergy;
	e_electricCurrent = e_elementaryCharge * (e_emittedElectronMol + e_emittedProtonMol) * NAvogadro;
	
	e_neutronGenInSpaceFlyRate = 1.0 - e_neutronGenInSpaceFallRate;
	e_genedNeutronInSpaceMol = e_emittedProtonMol * e_neutronGenInSpaceProtonRate;
	e_genedHydrogenInSpaceMol = e_emittedProtonMol * e_hydrogenGenInSpaceProtonRate;
	e_arrivedProtonMol = e_emittedProtonMol - (e_genedNeutronInSpaceMol + e_genedHydrogenInSpaceMol);
	e_arrivedProtonImperfectCollideMol = e_arrivedProtonMol * (1.0 - e_collideProtonRateOnElectrode);
	e_arrivedProtonRate = 1.0 - e_neutronGenInSpaceProtonRate - e_hydrogenGenInSpaceProtonRate;
	if(e_arrivedProtonMol >= 0.0
	&& 0.0 <= e_arrivedProtonRate && e_arrivedProtonRate <= 1.0
	&& 0.0 <= e_neutronGenInSpaceProtonRate && e_neutronGenInSpaceProtonRate <= 1.0
	&& 0.0 <= e_hydrogenGenInSpaceProtonRate && e_hydrogenGenInSpaceProtonRate <= 1.0){
		
		e_arrivedElectronMol = e_emittedElectronMol - (e_genedNeutronInSpaceMol + e_genedHydrogenInSpaceMol);
		if(e_emittedElectronMol > 0.0){
			e_arrivedElectronRate = e_arrivedElectronMol / e_emittedElectronMol;
			e_neutronGenInSpaceElectronRate = e_genedNeutronInSpaceMol / e_emittedElectronMol;
			e_hydrogenGenInSpaceElectronRate = e_genedHydrogenInSpaceMol / e_emittedElectronMol;
		}else{
			e_arrivedElectronRate = 0.0;
			e_neutronGenInSpaceElectronRate = 0.0;
			e_hydrogenGenInSpaceElectronRate = 0.0;
		}
		
		if(e_arrivedElectronMol >= 0.0
		&& 0.0 <= e_arrivedElectronRate && e_arrivedElectronRate <= 1.0 
		&& 0.0 <= e_neutronGenInSpaceElectronRate && e_neutronGenInSpaceElectronRate <= 1.0
		&& 0.0 <= e_hydrogenGenInSpaceElectronRate && e_hydrogenGenInSpaceElectronRate <= 1.0){

			e_emittedPulseEnergyMeV = e_appliedVoltageMeV
				* (e_genedNeutronInSpaceMol + e_genedHydrogenInSpaceMol + e_arrivedProtonMol + e_arrivedElectronMol)
				* NAvogadro;
			e_emittedPulseEnergyJ = e_emittedPulseEnergyMeV * 1.0E6 * e_energyJouleOfeV;
			e_emittedPulseEnergyCal = e_emittedPulseEnergyJ / e_energyJouleOfCalorie;
			e_emittedPulseEnergyKW = 3600 * e_emittedPulseEnergyJ / e_energyJouleOfKWH;

			/*
			e_lostHeatOfGenedNeutronnInSpaceMeV = e_appliedDefectMeV * e_genedNeutronInSpaceMol * NAvogadro;
			e_lostHeatOfGenedHydrogenInSpaceMeV = e_appliedVoltageMeV * e_genedHydrogenInSpaceMol * NAvogadro;
			e_lostHeatInSpaceMeV = e_lostHeatOfGenedNeutronnInSpaceMeV + e_lostHeatOfGenedHydrogenInSpaceMeV;
			
			e_lostHeatByGenedNeutronnInSpaceRate = e_lostHeatOfGenedNeutronnInSpaceMeV / e_emittedPulseEnergyMeV;
			e_lostHeatByGenedHydrogenInSpaceRate = e_lostHeatOfGenedHydrogenInSpaceMeV / e_emittedPulseEnergyMeV;
			e_lostHeatInSpaceRate = e_lostHeatInSpaceMeV / e_emittedPulseEnergyMeV;
			*/
		}else{
			fprintf(stderr,"FATAL ERROR:%s:too small e_arrivedElectronMol:%lg)\n", __FUNCTION__,  e_arrivedElectronMol);
			exit(1);
		}
	}else{
		fprintf(stderr,"FATAL ERROR:%s:too small e_arrivedProtonMol:%lg\n", __FUNCTION__, e_arrivedProtonMol);
		exit(1);
	}
}

extern void printContants(FILE * a_fp)
{
	int i = 1;
	fprintf(a_fp, "[Physics constants]\n");
	fprintf(a_fp, "P%d e_energyJouleOfeV      %lg [J/eV]\n", i++, e_energyJouleOfeV);
	fprintf(a_fp, "P%d e_energyJouleOfCalorie %lg [J/cal]\n", i++, e_energyJouleOfCalorie);
	fprintf(a_fp, "P%d e_energyJouleOfKWH %lg [J/KWH]\n", i++, e_energyJouleOfKWH);
	
	fprintf(a_fp, "P%d e_elementaryCharge %lg [C]\n", i++, e_elementaryCharge);
	fprintf(a_fp, "P%d e_vacuumPermittivity %lg [C]\n", i++, e_vacuumPermittivity);
	
	fprintf(a_fp, "P%d e_unifiedAtomicMassUnitKg %lg [kg]\n", i++, e_unifiedAtomicMassUnitKg);
	fprintf(a_fp, "P%d e_massNeutronKg  %lg [kg]\n", i++, e_massNeutronKg);
	fprintf(a_fp, "P%d e_massNeutronMeV %lg [MeV/C^2]\n", i++, e_massNeutronMeV);
	fprintf(a_fp, "P%d e_massProtonKg  %lg [kg]\n", i++, e_massProtonKg);
	fprintf(a_fp, "P%d e_massProtonMeV %lg [MeV/C^2]\n", i++, e_massProtonMeV);
	fprintf(a_fp, "P%d e_massElectronKg  %lg [kg]\n", i++, e_massElectronKg);
	fprintf(a_fp, "P%d e_massElectronMeV %lg [MeV/C^2]\n", i++, e_massElectronMeV);
	
	fprintf(a_fp, "[calculated Physics constants]\n");
	i = 1;
	fprintf(a_fp, "Q%d coefMeVtoKg(by Neutron)  %lg [kg]/[MeV/C^2]\n", i++, e_massNeutronKg / e_massNeutronMeV);
	fprintf(a_fp, "Q%d coefMeVtoKg(by Proton)   %lg [kg]/[MeV/C^2]\n", i++, e_massProtonKg / e_massProtonMeV);
	fprintf(a_fp, "Q%d coefMeVtoKg(by Electron) %lg [kg]/[MeV/C^2]\n", i++, e_massElectronKg / e_massElectronMeV);
	
	fprintf(a_fp, "Q%d e_coefElectroPotentialMeV %lg\n", i++, e_coefElectroPotentialMeV);

	fprintf(a_fp, "Q%d e_coefMassUToMeV %lg\n", i++, e_coefMassUToMeV);
	fprintf(a_fp, "Q%d e_coefMeVtoMassU %lg\n", i++, e_coefMeVtoMassU);
	fprintf(a_fp, "Q%d e_massNeutronMassU  %lg [u]\n", i++, e_massNeutronMassU);
	fprintf(a_fp, "Q%d e_massProtonMassU   %lg [u]\n", i++, e_massProtonMassU);
	fprintf(a_fp, "Q%d e_massElectronMassU %lg [u]\n", i++, e_massElectronMassU);
	fprintf(a_fp, "Q%d e_betaEnergyMassU %lg [u]\n", i++, e_betaEnergyMassU);
	fprintf(a_fp, "Q%d e_betaEnergyMeV   %lg [MeV/C^2]\n", i++, e_betaEnergyMeV);

	fprintf(a_fp, "[User conditions]\n");
	i = 1;
	fprintf(a_fp, "U%d e_rangeDecimalDigitPrecision %d\n", i++, e_rangeDecimalDigitPrecision);
	fprintf(a_fp, "U%d e_detectLimitRateForIsotope %lg\n", i++, e_detectLimitRateForIsotope);

	
	fprintf(a_fp, "U%d e_negativeElectrodeAtomicMols %s [mol]\n", i++, e_negativeElectrodeAtomicMols);
	fprintf(a_fp, "U%d e_positiveElectrodeAtomicMols %s [mol]\n", i++, e_positiveElectrodeAtomicMols);
	
	fprintf(a_fp, "U%d e_appliedVoltageScale %lg \n", i++, e_appliedVoltageScale);
	fprintf(a_fp, "U%d e_emittedElectronMol %lg [mol]\n", i++, e_emittedElectronMol);
	
	fprintf(a_fp, "U%d e_collideElectronRateOnElectrode %lg\n", i++, e_collideElectronRateOnElectrode);
	fprintf(a_fp, "U%d e_collideElectronRateForMidiMeV %lg\n", i++, e_collideElectronRateForMidiMeV);
	fprintf(a_fp, "U%d e_collideElectronRateForMiniMeV %lg\n", i++, e_collideElectronRateForMiniMeV);
	fprintf(a_fp, "U%d e_emittedProtonMol %lg [mol]\n", i++, e_emittedProtonMol);

	fprintf(a_fp, "U%d e_neutronGenInSpaceProtonRate %lg\n", i++, e_neutronGenInSpaceProtonRate);
	fprintf(a_fp, "U%d e_neutronGenInSpaceFallRate   %lg\n", i++, e_neutronGenInSpaceFallRate);
	
	fprintf(a_fp, "U%d e_hydrogenGenInSpaceProtonRate %lg\n", i++, e_hydrogenGenInSpaceProtonRate);
	fprintf(a_fp, "U%d e_collideProtonRateOnElectrode %lg\n", i++, e_collideProtonRateOnElectrode);
	fprintf(a_fp, "U%d e_collideElectronMidiMeV %lg\n", i++, e_collideElectronMidiMeV);
	fprintf(a_fp, "U%d e_collideElectronMiniMeV %lg\n", i++, e_collideElectronMiniMeV);
	fprintf(a_fp, "U%d e_usebulletCrossSection %d\n", i++, e_usebulletCrossSection);
	fprintf(a_fp, "U%d e_useComptonEffect %d\n", i++, e_useComptonEffect);
	
	fprintf(a_fp, "U%d e_stepsOfLostingEnergy %d\n", i++, e_stepsOfLostingEnergy);
	
	fprintf(a_fp, "[Calcuated User conditions]\n");
	i = 1;
	fprintf(a_fp, "V%d e_appliedVoltageMassU %lg [u]\n", i++, e_appliedVoltageMassU);
	fprintf(a_fp, "V%d e_appliedVoltageMeV   %lg [MeV/C^2]\n", i++, e_appliedVoltageMeV);
	fprintf(a_fp, "V%d e_appliedDefectMeV    %lg [MeV/C^2]\n", i++, e_appliedDefectMeV);
	fprintf(a_fp, "V%d e_stepByStepLostingEnergyMeV %lg [MeV/C^2]\n", i++, e_stepByStepLostingEnergyMeV);

	fprintf(a_fp, "V%d e_neutronGenInSpaceFlyRate %lg\n", i++, e_neutronGenInSpaceFlyRate);
	fprintf(a_fp, "V%d e_genedNeutronInSpaceMol %lg [mol]\n", i++, e_genedNeutronInSpaceMol);
	
	fprintf(a_fp, "V%d e_genedHydrogenInSpaceMol %lg [mol]\n", i++, e_genedHydrogenInSpaceMol);
	fprintf(a_fp, "V%d e_arrivedProtonMol %lg [mol]\n", i++, e_arrivedProtonMol);
	fprintf(a_fp, "V%d e_arrivedProtonImperfectCollideMol %lg [mol]\n", i++, e_arrivedProtonImperfectCollideMol);
	fprintf(a_fp, "V%d e_arrivedProtonRate %lg\n", i++, e_arrivedProtonRate);

	fprintf(a_fp, "V%d e_arrivedElectronMol %lg [mol]\n", i++, e_arrivedElectronMol);
	
	fprintf(a_fp, "V%d e_arrivedElectronRate %lg\n", i++, e_arrivedElectronRate);
	fprintf(a_fp, "V%d e_neutronGenInSpaceElectronRate %lg\n", i++, e_neutronGenInSpaceElectronRate);
	fprintf(a_fp, "V%d e_hydrogenGenInSpaceElectronRate %lg\n", i++, e_hydrogenGenInSpaceElectronRate);

	fprintf(a_fp, "V%d e_emittedPulseEnergyMeV %lg [MeV/C^2/sec]\n", i++, e_emittedPulseEnergyMeV);
	fprintf(a_fp, "V%d e_emittedPulseEnergyJ   %lg [J/sec]\n", i++, e_emittedPulseEnergyJ);
	fprintf(a_fp, "V%d e_emittedPulseEnergyCal %lg [cal/sec]\n", i++, e_emittedPulseEnergyCal);
	fprintf(a_fp, "V%d e_emittedPulseEnergyKW %lg [KWH/h]\n", i++, e_emittedPulseEnergyKW);
	
	fprintf(a_fp, "V%d e_electricCurrent %lg [C]\n", i++, e_electricCurrent);
	
	/*
	fprintf(a_fp, "V%d e_lostHeatOfGenedNeutronnInSpaceMeV %lg [MeV/C^2]\n", i++, e_lostHeatOfGenedNeutronnInSpaceMeV);
	fprintf(a_fp, "V%d e_lostHeatOfGenedHydrogenInSpaceMeV %lg [MeV/C^2]\n", i++, e_lostHeatOfGenedHydrogenInSpaceMeV);
	fprintf(a_fp, "V%d e_lostHeatInSpaceMeV %lg [MeV/C^2]\n", i++, e_lostHeatInSpaceMeV);
	fprintf(a_fp, "V%d e_lostHeatByGenedNeutronnInSpaceRate %lg\n", i++, e_lostHeatByGenedNeutronnInSpaceRate);
	fprintf(a_fp, "V%d e_lostHeatByGenedHydrogenInSpaceRate %lg\n", i++, e_lostHeatByGenedHydrogenInSpaceRate);
	fprintf(a_fp, "V%d e_lostHeatInSpaceRate %lg\n", i++, e_lostHeatInSpaceRate);
	*/
	
	fputs("\n", a_fp);
}

FILE * e_logFp = NULL;

#define SET_VALUE(FORMAT, TYPE, NAME) \
if(strcmp(type, #TYPE) == 0 && strcmp(name, #NAME) == 0){\
	sscanf(value, FORMAT, &NAME); \
	fprintf(stderr, "[INFO]%s %s = " FORMAT ";\n", type, name, NAME); \
	fprintf(e_logFp, "[INFO]%s %s = " FORMAT ";\n", type, name, NAME); \
}

extern void readUserConditions(char * a_fname)
{
	if(strcmp(a_fname, "default") == 0){
		initUserConditionsByDefault();
		fprintf(stderr, "[INFO]All parameter is set to default values.\n");	
		fprintf(e_logFp, "[INFO]All parameter is set to default values.\n");	
	}else{
		FILE * fp;
		fp = fopen(a_fname, "r");
		if(fp){
			char buff[4096];
			char buff2[4096], type[100], name[100], eq[100], value[100], sc[100];
			char * from;
			char * to;
			fprintf(stderr, "{ Open parameter file: %s\n", a_fname);
			fprintf(e_logFp, "{ Open parameter file: %s\n", a_fname);
			while(fgets(buff, 4096, fp)){
				fprintf(stderr, "%s", buff);
				fprintf(e_logFp, "%s", buff);
				{//trim comment
					char * slsl;
					slsl = strstr(buff, "//");
					if(slsl){
						slsl[0] = 0;
					}
				}
				//fprintf(stderr, "DEBUG:%s:buff %s\n", __FUNCTION__, buff);
				{//trim spaces
					int i, j, found;
					for(i = j = found = 0; buff[i]; ++i){
						//Change a tab code to a space code.
						if(buff[i] == '\t'){
							buff[i] = ' ';
						}
						//trim leading blank.
						if(buff[i] != ' '){
							found = 1;
						}
						if(found){
							if(j > 0 && (buff2[j - 1] == '=' || buff2[j - 1] == ';')){
								if(buff[i] != ' '){
									//insert a space code after '=' or ';'
									buff2[j] = ' ';
									j++;
								}
							}
							
							if(buff[i] == '=' || buff[i] == ';' || buff[i] == ',' || buff[i] == '"'){
								if(j > 0 && buff2[j - 1] != ' '){
									//insert a space code before '=', ';', ',' or '"'
									buff2[j] = ' ';
									j++;
								}
								//copy '=', ';', ',' or '"'
								buff2[j] = buff[i];
								j++;
							}else if(buff[i] == ' '){
								if(j == 0 || buff2[j - 1] != ' '){
									buff2[j] = buff[i];
									j++;
								}else{
									;//skip multi space
								}
							}else{
								//copy a character code
								buff2[j] = buff[i];
								j++;
							}
						}
						if(j == sizeof(buff2) - 1){
							break;//check the buffer over flow.
						}
					}
					buff2[j] = 0;
				}
				//fprintf(stderr, "DEBUG:%s:buff2 %s\n", __FUNCTION__, buff2);
#define NEGATIVE_ELECTRODE_HEAD "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \""
#define POSITIVE_ELECTRODE_HEAD "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \""
				if(memcmp(buff2, NEGATIVE_ELECTRODE_HEAD, strlen(NEGATIVE_ELECTRODE_HEAD)) == 0){
					from = buff2 + strlen(NEGATIVE_ELECTRODE_HEAD);
					to = strstr(from, "\"");
					if(to && to - from < TEXT_LEN_OF_MOLS){
						*to = 0;
						strcpy(e_negativeElectrodeAtomicMols, from);
					}else{
						fprintf(stderr, "ERROR:%s:Can't find last quotation. %s \n", __FUNCTION__, buff2);
						exit(1);
					}
				}else if(memcmp(buff2, POSITIVE_ELECTRODE_HEAD, strlen(POSITIVE_ELECTRODE_HEAD)) == 0){
					from = buff2 + strlen(POSITIVE_ELECTRODE_HEAD);
					to = strstr(from, "\"");
					if(to && to - from < TEXT_LEN_OF_MOLS){
						*to = 0;
						strcpy(e_positiveElectrodeAtomicMols, from);
					}else{
						fprintf(stderr, "ERROR:%s:Can't find last quotation. %s \n", __FUNCTION__, buff2);
						exit(1);
					}
				}else{
					sscanf(buff2, "%s %s %s %s %s", type, name, eq, value, sc);
					//fprintf(stderr, "DEBUG:%s:%s %s %s %s %s\n", __FUNCTION__, type, name, eq, value, sc);
					if(strcmp(eq, "=") == 0){
						if(strcmp(sc, ";") == 0){
							SET_VALUE("%d", int, e_rangeDecimalDigitPrecision)
							SET_VALUE("%lg", double, e_detectLimitRateForIsotope)
							SET_VALUE("%lg", double, e_appliedVoltageScale)
							SET_VALUE("%lg", double, e_emittedElectronMol)
							SET_VALUE("%lg", double, e_collideElectronRateOnElectrode)
							SET_VALUE("%lg", double, e_collideElectronRateForMidiMeV)
							SET_VALUE("%lg", double, e_collideElectronRateForMiniMeV)
							SET_VALUE("%lg", double, e_emittedProtonMol)
							SET_VALUE("%lg", double, e_neutronGenInSpaceProtonRate)
							SET_VALUE("%lg", double, e_neutronGenInSpaceFallRate)
							SET_VALUE("%lg", double, e_hydrogenGenInSpaceProtonRate)
							SET_VALUE("%lg", double, e_collideProtonRateOnElectrode)
							SET_VALUE("%lg", double, e_collideElectronMidiMeV)
							SET_VALUE("%lg", double, e_collideElectronMiniMeV)
							SET_VALUE("%d", int, e_usebulletCrossSection)
							SET_VALUE("%d", int, e_useComptonEffect)
						}else{
							//fprintf(stderr, "DEBUG:%s:sc %s\n", __FUNCTION__, sc);
						}
					}else{
						//fprintf(stderr, "DEBUG:%s:eq %s\n", __FUNCTION__, eq);
					}
				}
			}
			fclose(fp);
			fprintf(stderr, "\n} Close parameter file\n");
			fprintf(e_logFp, "\n} Close parameter file\n");
		}else{
			fprintf(stderr, "ERROR:%s:can't fopen(%s, r)\n", __FUNCTION__, a_fname);
			exit(1);
		}
	}
}
//---------------------------------------------------------------------

#define DECAY_MODE_ALPHA 1
#define DECAY_MODE_BETA_PLUS 2 
#define DECAY_MODE_BETA_PLUS_PROTON 3
#define DECAY_MODE_BETA_PLUS_ALPHA 4
#define DECAY_MODE_ELECTRON_CAPTURE 5
#define DECAY_MODE_ELECTRON_CAPTURE_PROTON 6
#define DECAY_MODE_ELECTRON_CAPTURE_ALPHA 7
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE 8 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_PROTON 9 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA 10 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_BETA_MINUS 11 
#define DECAY_MODE_BETA_MINUS_AND_NEWTRON 12
#define DECAY_MODE_PROTON_EMISSION 13
#define DECAY_MODE_2PROTON_EMISSION 14
#define DECAY_MODE_BETA_PLUS_2PROTON 15
#define DECAY_MODE_BETA_PLUS_3PROTON 16
#define DECAY_MODE_DOUBLE_BETA_MINUS 17 
#define DECAY_MODE_BETA_MINUS_AND_2NEWTRON 18
#define DECAY_MODE_BETA_MINUS_AND_3NEWTRON 19
#define DECAY_MODE_BETA_MINUS_AND_4NEWTRON 20
#define DECAY_MODE_BETA_MINUS_AND_ALPHA 21 
#define DECAY_MODE_NEWTRON_EMISSION 22 
#define DECAY_MODE_SELF_FISSION_80KR 23 //Don't use it! It is not programmed! // 180Tl -> 100Ru +  80Kr
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON 24 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON 25 // only use in the situation, can't use in isotopes table
#define MAX_DECAY_MODE DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON
extern const char * getCollideName(int a_collideType);
extern const char * getDecayModeText(int a_decayMode)
{
	const char * modeText = "unkown";
	switch(a_decayMode){
		case DECAY_MODE_ALPHA: modeText = "ALPHA"; break;
		case DECAY_MODE_BETA_PLUS: modeText = "BETA_PLUS"; break;
		case DECAY_MODE_BETA_PLUS_PROTON: modeText = "BETA_PLUS_PROTON"; break;
		case DECAY_MODE_BETA_PLUS_ALPHA: modeText = "BETA_PLUS_ALPHA"; break;
		case DECAY_MODE_ELECTRON_CAPTURE: modeText = "ELECTRON_CAPTURE"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_PROTON: modeText = "ELECTRON_CAPTURE_PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_ALPHA: modeText = "ELECTRON_CAPTURE_ALPHA"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE: modeText = "ELECTRON_CAPTURE_DEGRADE"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_PROTON: modeText = "ELECTRON_CAPTURE_DEGRADE_PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA: modeText = "ELECTRON_CAPTURE_DEGRADE_ALPHA"; break;
		case DECAY_MODE_BETA_MINUS: modeText = "BETA_MINUS"; break;
		case DECAY_MODE_BETA_MINUS_AND_NEWTRON: modeText = "BETA_MINUS_AND_NEWTRON"; break;
		case DECAY_MODE_PROTON_EMISSION: modeText = "PROTON_EMISSION"; break;
		case DECAY_MODE_2PROTON_EMISSION: modeText = "2PROTON_EMISSION"; break;
		case DECAY_MODE_BETA_PLUS_2PROTON: modeText = "BETA_PLUS_2PROTON"; break;
		case DECAY_MODE_BETA_PLUS_3PROTON: modeText = "BETA_PLUS_3PROTON"; break;
		case DECAY_MODE_DOUBLE_BETA_MINUS: modeText = "DOUBLE_BETA_MINUS"; break;
		case DECAY_MODE_BETA_MINUS_AND_2NEWTRON: modeText = "BETA_MINUS_AND_2NEWTRON"; break;
		case DECAY_MODE_BETA_MINUS_AND_3NEWTRON: modeText = "BETA_MINUS_AND_3NEWTRON"; break;
		case DECAY_MODE_BETA_MINUS_AND_4NEWTRON: modeText = "BETA_MINUS_AND_4NEWTRON"; break;
		case DECAY_MODE_BETA_MINUS_AND_ALPHA: modeText = "BETA_MINUS_AND_ALPHA"; break;
		case DECAY_MODE_NEWTRON_EMISSION: modeText = "NEWTRON_EMISSION"; break;
		case DECAY_MODE_SELF_FISSION_80KR: modeText = "SELF_FISSION_80KR"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON: modeText = "ELECTRON_CAPTURE_DEGRADE_2PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON: modeText = "ELECTRON_CAPTURE_DEGRADE_3PROTON"; break;
		default:
			if(a_decayMode > MAX_DECAY_MODE){
				modeText = getCollideName(a_decayMode);
			}
			break;
	}
	return modeText;
}

//Unit of the half-life.
#define HLU_STABLE 1
#define HLU_YEAR 2
#define HLU_DAY 3
#define HLU_HOUR 4
#define HLU_MINUTE 5
#define HLU_SECOND 6

extern const char * getHalfLifeUnitString(int a_halfLifeUnit)
{
	char * p = "???";
	switch(a_halfLifeUnit){
		case HLU_STABLE: p = "STABLE"; break;
		case HLU_YEAR: p = "year"; break;
		case HLU_DAY: p = "day"; break;
		case HLU_HOUR: p = "hour"; break;
		case HLU_MINUTE: p = "mimute"; break;
		case HLU_SECOND: p = "sec"; break;
		default: break;
	}
	return p;
}
#define SIZE_OF_DECAY_MODE 6

struct atomProperty;
struct isotopeProperty {
	struct atomProperty * atomPropertyPtr;
	int atomicNumber;
	int massNumber;
	double massU;//mass of isotope with electrons by the unified atomic mass unit
	double massMeV;
	double halfLife;//Numerical value of the half-life. //if(halfLifeUnit == HLU_STABLE){ "you can ignore this."}
	int halfLifeUnit;//Unit of the half-life. (HLU_*)
	int decayModeSize;
	int decayMode[SIZE_OF_DECAY_MODE];
	double decayModeRate[SIZE_OF_DECAY_MODE];
	double relativeIsotopicAbundance;//Natural abundance ratio by mole
	double nucleusRadius;//
	double relativeCollideCrossSection;//Relative collide-cross-sectional area, pow(massNumber, 2/3)
    char * symbol;
    char * name;
    char * japaneseName;
    char * nucleusName;
	//Nuclear isomer is not treated!
 };
extern double calcHalfLifeSec(const struct isotopeProperty * a_isotopePropertyPtr)
{
	double halfLifeSec;
	switch(a_isotopePropertyPtr->halfLifeUnit){
		case HLU_STABLE: halfLifeSec = (100.0 * 10000.0 * 10000.0 * 365.0 * 24.0 * 3600.0); break;
		case HLU_YEAR: halfLifeSec = a_isotopePropertyPtr->halfLife * (365.0 * 24.0 * 3600.0); break;
		case HLU_DAY: halfLifeSec = a_isotopePropertyPtr->halfLife * (24.0 * 3600.0); break;
		case HLU_HOUR: halfLifeSec = a_isotopePropertyPtr->halfLife * 3600.0; break;
		case HLU_MINUTE: halfLifeSec = a_isotopePropertyPtr->halfLife * 60.0; break;
		case HLU_SECOND: halfLifeSec = a_isotopePropertyPtr->halfLife; break;
	}
	return halfLifeSec;
}
struct atomProperty { 
	int atomicNumber;
    char * symbol;
    char * name;
    char * japaneseName;
	int minMassNumber;// The minimum mass number of isotope
	int maxMassNumber;// The maximum mass number of isotope
	int min1secMassNumber;// The mimimum mass number of isotope that has a long half life greater than 1 second.
	int max1secMassNumber;// The maximum mass number of isotope that has a long half life greater than 1 second.
	double heatCapacity;// [J/(mol*k)]
    struct isotopeProperty * isotopePtr;
} * e_atomPropertyPtr;

#define MIN_ATOM_NUMBER 0
#define MAX_ATOM_NUMBER 140
#define UNDEF_ATOMIC_NUMBER (-1)
#define UNDEF_MASS_NUMBER (-1)

//#define ATOMICNUMBER_PHOTON 0
//#define MASSNUMBER_PHOTON -1
#define ATOMICNUMBER_ELECTRON 0
#define MASSNUMBER_ELECTRON 0
#define ATOMICNUMBER_NEUTRON 0
#define MASSNUMBER_NEUTRON 1
#define MASSNUMBER_DINEUTRON 2
#define MASSNUMBER_TRINEUTRON 3
#define ATOMICNUMBER_HYDROGEN 1
#define MASSNUMBER_HYDROGEN 1
#define ATOMICNUMBER_DEUTERIUM ATOMICNUMBER_HYDROGEN
#define MASSNUMBER_DEUTERIUM 2
#define ATOMICNUMBER_TRITIUM ATOMICNUMBER_HYDROGEN
#define MASSNUMBER_TRITIUM 3
#define ATOMICNUMBER_HELIUM 2
#define ATOMICNUMBER_HELIUM2 ATOMICNUMBER_HELIUM
#define ATOMICNUMBER_HELIUM3 ATOMICNUMBER_HELIUM
#define MASSNUMBER_HELIUM2 2
#define MASSNUMBER_HELIUM3 3
#define MASSNUMBER_HELIUM 4
#define ATOMICNUMBER_Li 3
#define MASSNUMBER_Li4 4
#define ATOMICNUMBER_Fe 26
#define ATOMICNUMBER_Co 27
#define ATOMICNUMBER_Ni 28
#define ATOMICNUMBER_Cu 29
#define ATOMICNUMBER_Zn 30
#define ATOMICNUMBER_Os 76
#define ATOMICNUMBER_Ir 77
#define ATOMICNUMBER_Pt 78
#define ATOMICNUMBER_Au 79
#define ATOMICNUMBER_Hg 80
#define ATOMICNUMBER_Tl 81
#define ATOMICNUMBER_Pb 82

int e_logErrorUnregistedAtom = 0;

int e_unregistedAtomCnt = 0;
int e_unregistedAtomNumbers[MAX_ATOM_NUMBER];

extern int compareInt(const void * a, const void* b)
{
	return *((int *)a) - *((int *)b);
}
extern void printUnregistedAtomNumbers(FILE * fp)
{
	qsort(e_unregistedAtomNumbers, e_unregistedAtomCnt, sizeof(int), compareInt);
	if(e_unregistedAtomCnt > 0){
		int i;
		fprintf(fp, "\nUnregisted Atoms\n");
		for(i = 0; i < e_unregistedAtomCnt; ++i){
			fprintf(fp, "%d:AtomNumber=%d\n", i + 1, e_unregistedAtomNumbers[i]);
		}
	}else{
		fprintf(fp, "\nThere is none of unregisted atom.\n");
	}
}

#define VERSION_OF_serializeunregistedAtom 1
#define SERIALIZE_UNREGISTED_ATOM(RW)\
extern int RW ## UnregistedAtom(FILE * a_fp)\
{\
	SERIALIZE_VERSION_CHECK(RW, VERSION_OF_serializeunregistedAtom, a_fp) \
	SERIALIZE_VALUE(RW, e_unregistedAtomCnt, a_fp)\
	SERIALIZE_ARRAY(RW, e_unregistedAtomNumbers, int, e_unregistedAtomCnt, a_fp)\
	if((void*)RW == (void*)fread){\
		if(e_unregistedAtomCnt > 0){\
			fprintf(stderr, "WARNING:%s:There are unregisted Atoms.\n", __FUNCTION__);\
			printUnregistedAtomNumbers(stderr);\
		}\
	}\
	return 1;\
}
#define CALL_SERIALIZE_UNREGISTED_ATOM(RW, FP) \
	RW ## UnregistedAtom(FP)
SERIALIZE_UNREGISTED_ATOM(fread)//freadUnregistedAtom
SERIALIZE_UNREGISTED_ATOM(fwrite)//fwriteUnregistedAtom

extern struct atomProperty * getAtomPropertyPtr(int a_atomicNumber)
{
	//fprintf(stderr, "debug:%s[%d]a_atomicNumber=%d\n", __FUNCTION__, __LINE__, a_atomicNumber);
	if(MIN_ATOM_NUMBER <= a_atomicNumber && a_atomicNumber <= MAX_ATOM_NUMBER){
		if(e_atomPropertyPtr[a_atomicNumber].atomicNumber != UNDEF_ATOMIC_NUMBER){
			return (e_atomPropertyPtr + a_atomicNumber);
		}else{
			if(e_logErrorUnregistedAtom){
				int i;
				for(i = 0; i < e_unregistedAtomCnt; ++i){
					if(a_atomicNumber == e_unregistedAtomNumbers[i]){
						break;//found in e_unregistedAtomNumbers[]
					}
				}
				if(i == e_unregistedAtomCnt){
					if(e_unregistedAtomCnt < MAX_ATOM_NUMBER){
						e_unregistedAtomNumbers[i] = a_atomicNumber;
						++e_unregistedAtomCnt;
					}else{
						fprintf(stderr, "ERROR:%s:too many unregisted atom MAX_ATOM_NUMBER:%d a_atomicNumber:%d\n", __FUNCTION__, MAX_ATOM_NUMBER, a_atomicNumber);
					}
				}
				//fprintf(stderr, "ERROR:%s:unregisted atom:%d\n", __FUNCTION__, a_atomicNumber);
			}
		}
	}else{
		fprintf(stderr, "ERROR:%s:wrong a_atomicNumber:%d\n", __FUNCTION__, a_atomicNumber);
	}
	return NULL;
}
extern struct atomProperty * findAtomPropertyPtr(char * a_elementSymbolPtr)
{
	int atomicNumber;
	struct atomProperty * ret = NULL;
	for(atomicNumber = MIN_ATOM_NUMBER; atomicNumber <= MAX_ATOM_NUMBER; ++atomicNumber){
		if(e_atomPropertyPtr[atomicNumber].atomicNumber != UNDEF_ATOMIC_NUMBER){
			if(strcmp(e_atomPropertyPtr[atomicNumber].symbol, a_elementSymbolPtr) == 0){
				ret = e_atomPropertyPtr + atomicNumber;
				break;
			}
		}
	}
	return ret;
}

#define MAX_UNREGISTED_ISOTOPE_CNT 500
int e_unregistedIsotopeCnt = 0;
int e_unregistedIsotopeNumbers[MAX_UNREGISTED_ISOTOPE_CNT];

extern void printUnregistedIsotopeNumbers(FILE * fp)
{
	qsort(e_unregistedIsotopeNumbers, e_unregistedIsotopeCnt, sizeof(int), compareInt);
	if(e_unregistedIsotopeCnt > 0){
		int i;
		fprintf(fp, "\nUnregisted Isotopes\n");
		for(i = 0; i < e_unregistedIsotopeCnt; ++i){
			int atomicNumber = (e_unregistedIsotopeNumbers[i] >> 12) & 0x0fff;
			int massNumber = e_unregistedIsotopeNumbers[i] & 0x0fff;
			struct atomProperty * atomPropertyPtr = getAtomPropertyPtr(atomicNumber);
			if(atomPropertyPtr){
				fprintf(fp, "%d:%d%s AtomNumber=%d MassNumber:%d\n", i + 1, 
					massNumber, atomPropertyPtr->symbol, atomicNumber, massNumber);
			}else{
				fprintf(fp, "%d: AtomNumber=%d MassNumber:%d\n", i + 1, 
					atomicNumber, massNumber);
			}
		}
	}else{
		fprintf(fp, "\nThere is none of unregisted Isotope.\n");
	}
}
#define VERSION_OF_serializeunregistedIsotope 1
#define SERIALIZE_UNREGISTED_ISOTOPE(RW)\
extern int RW ## UnregistedIsotope(FILE * a_fp)\
{\
	SERIALIZE_VERSION_CHECK(RW, VERSION_OF_serializeunregistedIsotope, a_fp) \
	SERIALIZE_VALUE(RW, e_unregistedIsotopeCnt, a_fp)\
	SERIALIZE_ARRAY(RW, e_unregistedIsotopeNumbers, int, e_unregistedIsotopeCnt, a_fp)\
	if((void *)RW == (void *)fread){\
		if(e_unregistedIsotopeCnt > 0){\
			fprintf(stderr, "WARNING:%s:There are unregisted Isotopes.\n", __FUNCTION__);\
			printUnregistedIsotopeNumbers(stderr);\
		}\
	}\
	return 1;\
}
#define CALL_SERIALIZE_UNREGISTED_ISOTOPE(RW, FP) \
	RW ## UnregistedIsotope(FP)
SERIALIZE_UNREGISTED_ISOTOPE(fread);//freadUnregistedIsotope
SERIALIZE_UNREGISTED_ISOTOPE(fwrite);//fwriteUnregistedIsotope

extern struct isotopeProperty * getIsotopePropertyPtr(int a_atomicNumber, int a_massNumber)
{
	struct atomProperty * atomPropertyPtr = getAtomPropertyPtr(a_atomicNumber);
	struct isotopeProperty * ret = NULL;
	if(atomPropertyPtr){
		if(atomPropertyPtr->minMassNumber <= a_massNumber && a_massNumber <= atomPropertyPtr->maxMassNumber){
			int index = a_massNumber - atomPropertyPtr->minMassNumber;
			ret = (atomPropertyPtr->isotopePtr + index);
		}else{
			if(atomPropertyPtr->min1secMassNumber <= a_massNumber && a_massNumber <= atomPropertyPtr->max1secMassNumber){
				if(e_logErrorUnregistedAtom){
					//fprintf(stderr, "DEBUG:%s:The isotope %d~%s has a longer half life than 1 sec.\n", 
					//	__FUNCTION__, a_massNumber, atomPropertyPtr->symbol);

					int atomMass = (a_atomicNumber << 12) + a_massNumber;
					int i;
					for(i = 0; i < e_unregistedIsotopeCnt; ++i){
						if(atomMass == e_unregistedIsotopeNumbers[i]){
							break;//found in e_unregistedIsotopeNumbers[]
						}
					}
					if(i == e_unregistedIsotopeCnt){
						if(e_unregistedIsotopeCnt < MAX_UNREGISTED_ISOTOPE_CNT){
							e_unregistedIsotopeNumbers[i] = atomMass;
							++e_unregistedIsotopeCnt;
						}else{
							fprintf(stderr, "ERROR:%s:too many unregisted isotopes MAX_UNREGISTED_ISOTOPE_CNT:%d a_atomicNumber:%d a_massNumber:%d\n", __FUNCTION__, MAX_UNREGISTED_ISOTOPE_CNT, a_atomicNumber, a_massNumber);
						}
					}
					//fprintf(stderr, "ERROR:%s:unregisted isotope %d~%s\n", __FUNCTION__, a_massNumber, atomPropertyPtr->symbol);
				}
			}else{
				//fprintf(stderr, "DEBUG:%s:The isotope %d~%s has a shorter half life than 1 sec.\n", 
				//	__FUNCTION__, a_massNumber, atomPropertyPtr->symbol);
			}
		}
	}
	return ret;
}
extern void initIsotopeProperty(int a_atomicNumber, int a_massNumber, double a_massU, double a_halfLife, int a_halfLifeUnit, double a_relativeIsotopicAbundance)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		int n;
		double upLimitMassU = (e_massProtonMassU + e_massElectronMassU) * a_atomicNumber + e_massNeutronMassU * (a_massNumber - a_atomicNumber);
		double lowLimitMassU = a_massNumber * (55.9349375/56);//55.9349375 is massU of 56Fe
		if(a_massU > upLimitMassU && a_massNumber > 0){
			fprintf(stderr, "ERROR:%s:wrong a_massU:%lg > %lg, a_atomicNumber:%d a_massNumber:%d\n", __FUNCTION__, a_massU, upLimitMassU, a_atomicNumber, a_massNumber);
			exit(1);
		}
		if(a_massU < lowLimitMassU){
			fprintf(stderr, "ERROR:%s:wrong a_massU:%lg < %lg, a_atomicNumber:%d a_massNumber:%d\n", __FUNCTION__, a_massU, lowLimitMassU, a_atomicNumber, a_massNumber);
			exit(1);
		}
		if(a_halfLife < 0.0){
			fprintf(stderr, "ERROR:%s:wrong a_halfLife:%lg, a_atomicNumber:%d a_massNumber:%d\n", __FUNCTION__, a_halfLife, a_atomicNumber, a_massNumber);
			exit(1);
		}
		isotopePropertyPtr->atomicNumber = a_atomicNumber;
		isotopePropertyPtr->massNumber = a_massNumber;
		isotopePropertyPtr->massU = a_massU;
		isotopePropertyPtr->massMeV = a_massU * e_coefMassUToMeV;
		isotopePropertyPtr->halfLife = a_halfLife;
		isotopePropertyPtr->halfLifeUnit = a_halfLifeUnit;
		isotopePropertyPtr->relativeIsotopicAbundance = a_relativeIsotopicAbundance;
		if(a_massNumber > 0){
			isotopePropertyPtr->nucleusRadius = pow(a_massNumber, 1.0 / 3.0) * 1E-15;
			isotopePropertyPtr->relativeCollideCrossSection = pow(a_massNumber, 2.0 / 3.0);
			//pow(1, 2.0 / 3.0); = 1
			//pow(2, 2.0 / 3.0); = 1.5874010519681994
			//pow(3, 2.0 / 3.0); = 2.080083823051904
			//pow(4, 2.0 / 3.0); = 2.5198420997897464
		}else if(a_massNumber == MASSNUMBER_ELECTRON && a_atomicNumber == ATOMICNUMBER_ELECTRON){
			isotopePropertyPtr->nucleusRadius = 2.5e-18;
			isotopePropertyPtr->relativeCollideCrossSection = pow(2.5e-18 / 0.8775E-15, 2.0);// = 0.00000811681723362635
			// (2.5e-18 = electron radius by weak boson) / (0.8775E-15 = rms charge radius of proton)
		}else{
			isotopePropertyPtr->nucleusRadius = 0.0;
			isotopePropertyPtr->relativeCollideCrossSection = 0;
		}
		n = strlen(getAtomPropertyPtr(a_atomicNumber)->symbol) + 4;
		isotopePropertyPtr->symbol = clearAlloc(n, "symbol");
		snprintf(isotopePropertyPtr->symbol, n, "%d%s", a_massNumber, getAtomPropertyPtr(a_atomicNumber)->symbol);
		isotopePropertyPtr->name = "";
		isotopePropertyPtr->japaneseName = "";
		isotopePropertyPtr->nucleusName = "";
	}
}
extern void initIsotopePropertySymbols(int a_atomicNumber, int a_massNumber, char * a_symbol, char * a_name, char * a_japaneseName, char * a_nucleusName)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		isotopePropertyPtr->symbol = allocStrcpy(a_symbol);
		isotopePropertyPtr->name = allocStrcpy(a_name);
		isotopePropertyPtr->japaneseName = allocStrcpy(a_japaneseName);
		isotopePropertyPtr->nucleusName = allocStrcpy(a_nucleusName);
	}
}
extern void initDecayMode(int a_atomicNumber, int a_massNumber, int a_decayMode)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		isotopePropertyPtr->decayModeSize = 1;
		isotopePropertyPtr->decayMode[0] = a_decayMode;
		isotopePropertyPtr->decayModeRate[0] = 100.0;
	}
}
extern void initDecayMode2(int a_atomicNumber, int a_massNumber, int a_decayMode, double a_decayModePercent, int a_decayMode2)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		isotopePropertyPtr->decayModeSize = 2;
		isotopePropertyPtr->decayMode[0] = a_decayMode;
		isotopePropertyPtr->decayModeRate[0] = a_decayModePercent / 100.0;
		isotopePropertyPtr->decayMode[1] = a_decayMode2;
		isotopePropertyPtr->decayModeRate[1] = 1.0 - isotopePropertyPtr->decayModeRate[0];
	}
}
extern void initDecayMode3(int a_atomicNumber, int a_massNumber, int a_decayMode, double a_decayModePercent, int a_decayMode2, double a_decayModePercent2, int a_decayMode3)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		isotopePropertyPtr->decayModeSize = 3;
		isotopePropertyPtr->decayMode[0] = a_decayMode;
		isotopePropertyPtr->decayModeRate[0] = a_decayModePercent / 100.0;
		isotopePropertyPtr->decayMode[1] = a_decayMode2;
		isotopePropertyPtr->decayModeRate[1] = a_decayModePercent2 / 100.0;
		if(isotopePropertyPtr->decayModeRate[0] + isotopePropertyPtr->decayModeRate[1] < 1.0){
			isotopePropertyPtr->decayMode[2] = a_decayMode3;
			isotopePropertyPtr->decayModeRate[2] = 1.0 - (isotopePropertyPtr->decayModeRate[0] + isotopePropertyPtr->decayModeRate[1]);
		}else{
			fprintf(stderr, "FATAL ERROR:%s:a_atomicNumber:%d a_massNumber:%d a_decayModePercent:%lg + a_decayModePercent2:%lg >= 100.0 \n", __FUNCTION__, a_atomicNumber, a_massNumber, a_decayModePercent, a_decayModePercent2);
			exit(1);
		}
	}
}
extern void initDecayMode4(int a_atomicNumber, int a_massNumber, int a_decayMode, double a_decayModePercent, int a_decayMode2, double a_decayModePercent2, int a_decayMode3, double a_decayModePercent3, int a_decayMode4)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		isotopePropertyPtr->decayModeSize = 4;
		isotopePropertyPtr->decayMode[0] = a_decayMode;
		isotopePropertyPtr->decayModeRate[0] = a_decayModePercent / 100.0;
		isotopePropertyPtr->decayMode[1] = a_decayMode2;
		isotopePropertyPtr->decayModeRate[1] = a_decayModePercent2 / 100.0;
		isotopePropertyPtr->decayMode[2] = a_decayMode3;
		isotopePropertyPtr->decayModeRate[2] = a_decayModePercent3 / 100.0;
		if(isotopePropertyPtr->decayModeRate[0] + isotopePropertyPtr->decayModeRate[1] + isotopePropertyPtr->decayModeRate[2] < 1.0){
			isotopePropertyPtr->decayMode[3] = a_decayMode4;
			isotopePropertyPtr->decayModeRate[3] = 1.0 - (isotopePropertyPtr->decayModeRate[0] + isotopePropertyPtr->decayModeRate[1] + isotopePropertyPtr->decayModeRate[2]);
		}else{
			fprintf(stderr, "FATAL ERROR:%s:a_decayModePercent:%lg + a_decayModePercent2:%lg + a_decayModePercent3:%lg >= 100.0 \n", __FUNCTION__, a_decayModePercent, a_decayModePercent2, a_decayModePercent3);
			exit(1);
		}
	}
}
extern void initDecayMode5(int a_atomicNumber, int a_massNumber, int a_decayMode, double a_decayModePercent, int a_decayMode2, double a_decayModePercent2, int a_decayMode3, double a_decayModePercent3, int a_decayMode4, double a_decayModePercent4, int a_decayMode5)
{
	struct isotopeProperty * isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(isotopePropertyPtr){
		isotopePropertyPtr->decayModeSize = 5;
		isotopePropertyPtr->decayMode[0] = a_decayMode;
		isotopePropertyPtr->decayModeRate[0] = a_decayModePercent / 100.0;
		isotopePropertyPtr->decayMode[1] = a_decayMode2;
		isotopePropertyPtr->decayModeRate[1] = a_decayModePercent2 / 100.0;
		isotopePropertyPtr->decayMode[2] = a_decayMode3;
		isotopePropertyPtr->decayModeRate[2] = a_decayModePercent3 / 100.0;
		isotopePropertyPtr->decayMode[3] = a_decayMode4;
		isotopePropertyPtr->decayModeRate[3] = a_decayModePercent4 / 100.0;
		if(isotopePropertyPtr->decayModeRate[0] + isotopePropertyPtr->decayModeRate[1] + isotopePropertyPtr->decayModeRate[2] + isotopePropertyPtr->decayModeRate[3] < 1.0){
			isotopePropertyPtr->decayMode[4] = a_decayMode5;
			isotopePropertyPtr->decayModeRate[4] = 1.0 - (isotopePropertyPtr->decayModeRate[0] + isotopePropertyPtr->decayModeRate[1] + isotopePropertyPtr->decayModeRate[2] + isotopePropertyPtr->decayModeRate[3]);
		}else{
			fprintf(stderr, "FATAL ERROR:%s:a_decayModePercent:%lg + a_decayModePercent2:%lg + a_decayModePercent3:%lg + a_decayModePercent4:%lg >= 100.0 \n", __FUNCTION__, a_decayModePercent, a_decayModePercent2, a_decayModePercent3, a_decayModePercent4);
			exit(1);
		}
	}
}
extern void initAtomProperty(int a_atomicNumber, char * a_symbol, char * a_name, char * a_japaneseName, int a_minMassNumber, int a_maxMassNumber, int a_min1secMassNumber, int a_max1secMassNumber, double a_heatCapacity)
{
	struct atomProperty * atomPropertyPtr = e_atomPropertyPtr + a_atomicNumber; 
	if(atomPropertyPtr){
		int i;
		int countOfIsotope = a_maxMassNumber - a_minMassNumber + 1;
		atomPropertyPtr->atomicNumber = a_atomicNumber;
		atomPropertyPtr->symbol = allocStrcpy(a_symbol);
		atomPropertyPtr->name = allocStrcpy(a_name);
		atomPropertyPtr->japaneseName = allocStrcpy(a_japaneseName);
		atomPropertyPtr->minMassNumber = a_minMassNumber;
		atomPropertyPtr->maxMassNumber = a_maxMassNumber;
		atomPropertyPtr->min1secMassNumber = a_min1secMassNumber;
		atomPropertyPtr->max1secMassNumber = a_max1secMassNumber;
		atomPropertyPtr->heatCapacity = a_heatCapacity;
		atomPropertyPtr->isotopePtr = clearAlloc(sizeof(struct isotopeProperty) * countOfIsotope, "isotopes");
		for(i = 0; i < countOfIsotope; ++ i){
			atomPropertyPtr->isotopePtr[i].atomPropertyPtr = atomPropertyPtr;
			atomPropertyPtr->isotopePtr[i].atomicNumber = a_atomicNumber;
		}
	}
}
struct isotopeProperty * e_hydrogenPtr;
struct isotopeProperty * e_helium4Ptr;
extern void initAtomProperties()
{
	e_atomPropertyPtr = (struct atomProperty *)clearAlloc(sizeof(struct atomProperty) * (MAX_ATOM_NUMBER + 1), "atomPropertyPtr");
	if(e_atomPropertyPtr){
		int i;
		//fprintf(stderr, "debug:%s[%d]\n", __FUNCTION__, __LINE__);
		for(i = MIN_ATOM_NUMBER; i <= MAX_ATOM_NUMBER; ++i){
			//fprintf(stderr, "debug:%s[%d]i=%d\n", __FUNCTION__, __LINE__, i);
			e_atomPropertyPtr[i].atomicNumber = UNDEF_ATOMIC_NUMBER;
		}
		initAtomProperty(0, "n", "elementary particles", "素粒子", 0, 3, 0, 1, 0.0);
		//initIsotopeProperty(0, -1, 0.0, 0.0, HLU_STABLE, 0.0);
		initIsotopeProperty(0, 0, e_massElectronMassU, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(0, 1, e_massNeutronMassU, 886.7, HLU_SECOND, 0.0);
		initIsotopeProperty(0, 2, e_massNeutronMassU * 2.0, 0.000001, HLU_SECOND, 0.0);
		initIsotopeProperty(0, 3, e_massNeutronMassU * 3.0, 0.000001, HLU_SECOND, 0.0);
		//initIsotopePropertySymbols(0, -1, "γ", "photon", "ガンマ光子", "γ");
		initIsotopePropertySymbols(0, 0, "e", "electron", "電子", "e");
		initIsotopePropertySymbols(0, 1, "n", "neutron", "中性子", "n");
		initIsotopePropertySymbols(0, 2, "2n", "dineutron", "仮想二重中性子", "2n");
		initIsotopePropertySymbols(0, 3, "3n", "trineutron", "仮想三重中性子", "3n");
		initDecayMode(0, 1, DECAY_MODE_BETA_MINUS);
		
		initAtomProperty(1, "H", "hydrogen", "水素", 1, 3, 1, 3, 28.836);
		initIsotopeProperty(1, 1, 1.00782503207, 0.0, HLU_STABLE, 0.999885);
		initIsotopeProperty(1, 2, 2.0141017778, 0.0, HLU_STABLE, 0.000115); 
		initIsotopeProperty(1, 3, 3.0160492777, 12.32, HLU_YEAR, 0.0);
		initIsotopePropertySymbols(1, 1, "H", "hydrogen", "水素", "p");
		initIsotopePropertySymbols(1, 2, "D", "deuterium", "重水素", "d");
		initIsotopePropertySymbols(1, 3, "T", "tritium", "三重水素", "t");
		initDecayMode(1, 3, DECAY_MODE_BETA_MINUS);
		
		initAtomProperty(2, "He", "helium", "ヘリウム", 3, 8, 3, 4, 20.786);
		initIsotopeProperty(2, 3, 3.0160293191, 0.0, HLU_STABLE, 0.00000134);
		initIsotopeProperty(2, 4, 4.00260325415, 0.0, HLU_STABLE, 0.99999866);
		initIsotopeProperty(2, 5, 5.01222, 700.0e-24, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 6, 6.0188891, 806.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 7, 7.028021, 2.9e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 8, 8.033922, 119.0e-3, HLU_SECOND, 0.0);
		initDecayMode(2, 5, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode2(2, 6, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode(2, 7, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode2(2, 8, DECAY_MODE_BETA_MINUS, 84.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);//ignore "β−, fission (0.9%, 5He + 3H)"

		initAtomProperty(3, "Li", "lithium", "リチウム", 4, 11, 6, 7, 24.860);
		initIsotopeProperty(3, 4, 4.02719, 91.0E-24, HLU_SECOND, 0.0);
		initIsotopeProperty(3, 5, 5.01254, 370.0E-24, HLU_SECOND, 0.0);
		initIsotopeProperty(3, 6, 6.015122795, 0.0, HLU_STABLE, 0.0759);
		initIsotopeProperty(3, 7, 7.01600455, 0.0, HLU_STABLE, 0.9241);
		initIsotopeProperty(3, 8, 8.02248736, 840.3E-3, HLU_SECOND, 0.0);
		initIsotopeProperty(3, 9, 9.0267895, 178.3E-3, HLU_SECOND, 0.0);
		initIsotopeProperty(3, 10, 10.035481, 2.0E-21, HLU_SECOND, 0.0);
		initIsotopeProperty(3, 11, 11.043798, 8.75E-3, HLU_SECOND, 0.0);
		initDecayMode(3, 4, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(3, 5, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(3, 8, DECAY_MODE_BETA_MINUS);
		initDecayMode2(3, 9, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 50.8, DECAY_MODE_BETA_MINUS);
		initDecayMode(3, 10, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode5(3, 11, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 84.9, DECAY_MODE_BETA_MINUS, 8.07, DECAY_MODE_BETA_MINUS_AND_2NEWTRON, 4.1, DECAY_MODE_BETA_MINUS_AND_3NEWTRON, 1.9, DECAY_MODE_BETA_MINUS_AND_ALPHA); //ignore "β−, fission (.014%) 8Li, 3H", "β−, fission (.013%), 9Li 2H"

		initAtomProperty(4, "Be", "beryllium", "ベリリウム", 6, 14, 7, 11, 16.443);
		initIsotopeProperty(4, 6, 6.019726, 5.0e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(4, 7, 7.01692983, 53.22, HLU_DAY, 0.0);
		initIsotopeProperty(4, 8, 8.00530510, 6.7e-17, HLU_SECOND, 0.0);
		initIsotopeProperty(4, 9, 9.0121822, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(4, 10, 10.0135338, 1.39e6, HLU_YEAR, 0.0);
		initIsotopeProperty(4, 11, 11.021658, 13.81, HLU_SECOND, 0.0);
		initIsotopeProperty(4, 12, 12.026921, 21.49e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(4, 13, 13.03569, 2.7e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(4, 14, 14.04289, 4.84e-3, HLU_SECOND, 0.0);
		initDecayMode(4, 6, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode(4, 7, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(4, 8, DECAY_MODE_ALPHA);
		initDecayMode(4, 10, DECAY_MODE_BETA_MINUS);
		initDecayMode2(4, 11, DECAY_MODE_BETA_MINUS, 97.1, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode2(4, 12, DECAY_MODE_BETA_MINUS, 99.48, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(4, 13, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode3(4, 14, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 81.0, DECAY_MODE_BETA_MINUS, 14.0, DECAY_MODE_BETA_MINUS_AND_2NEWTRON);

		initAtomProperty(5, "B", "boron", "ホウ素", 7, 19, 10, 11, 11.087);
		initIsotopeProperty(5, 7, 7.02992, 350.0e-24, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 8, 8.0246072, 770.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 9, 9.0133288, 800.0e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 10, 10.0129370, 0.0, HLU_STABLE, 19.9);
		initIsotopeProperty(5, 11, 11.0093054, 0.0, HLU_STABLE, 80.1);
		initIsotopeProperty(5, 12, 12.0143521, 20.20e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 13, 13.0177802, 17.33e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 14, 14.025404, 12.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 15, 15.031103, 9.87e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 16, 16.03981, 190.0e-12, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 17, 17.04699, 5.08e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 18, 18.05617, 26.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(5, 19, 19.06373, 2.92e-3, HLU_SECOND, 0.0);
		initDecayMode(5, 7, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(5, 8, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(5, 9, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(5, 12, DECAY_MODE_BETA_MINUS, 98.4, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode2(5, 13, DECAY_MODE_BETA_MINUS, 99.72, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(5, 14, DECAY_MODE_BETA_MINUS, 93.96, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode3(5, 15, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 93.6, DECAY_MODE_BETA_MINUS, 6.0, 	DECAY_MODE_BETA_MINUS_AND_2NEWTRON);
		initDecayMode(5, 16, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode5(5, 17, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 63.0, DECAY_MODE_BETA_MINUS, 22.1, DECAY_MODE_BETA_MINUS_AND_2NEWTRON, 11.0, DECAY_MODE_BETA_MINUS_AND_3NEWTRON, 3.5, DECAY_MODE_BETA_MINUS_AND_4NEWTRON);
		initDecayMode(5, 18, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode(5, 19, DECAY_MODE_BETA_MINUS);

		initAtomProperty(6, "C", "carbon", "炭素", 8, 22, 10, 16, 8.517);//8.517:graphite, 6.155:diamond
		initIsotopeProperty(6, 8, 8.037675, 2.0e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 9, 9.0310367, 126.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 10, 10.0168532, 19.290, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 11, 11.0114336, 20.334, HLU_MINUTE, 0.0);
		initIsotopeProperty(6, 12, 12.0, 0.0, HLU_STABLE, 0.9893);
		initIsotopeProperty(6, 13, 13.0033548378, 0.0, HLU_STABLE, 0.0107);
		initIsotopeProperty(6, 14, 14.003241989, 5730.0, HLU_YEAR, 0.0);
		initIsotopeProperty(6, 15, 15.0105993, 2.449, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 16, 16.014701, 0.747, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 17, 17.022586, 193.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 18, 18.02676, 92.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 19, 19.03481, 46.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 20, 20.04032, 16.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 21, 21.04934, 30.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(6, 22, 22.05720, 6.2e-3, HLU_SECOND, 0.0);
		initDecayMode(6, 8, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode3(6, 9, DECAY_MODE_BETA_PLUS, 60.0, DECAY_MODE_BETA_PLUS_PROTON, 23.0, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(6, 10, DECAY_MODE_BETA_PLUS);
		initDecayMode2(6, 11, DECAY_MODE_BETA_PLUS, 99.79, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(6, 14, DECAY_MODE_BETA_MINUS);
		initDecayMode(6, 15, DECAY_MODE_BETA_MINUS);
		initDecayMode2(6, 16, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 97.9, DECAY_MODE_BETA_MINUS);
		initDecayMode2(6, 17, DECAY_MODE_BETA_MINUS, 71.59, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(6, 18, DECAY_MODE_BETA_MINUS, 68.5, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode3(6, 19, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 47.0, DECAY_MODE_BETA_MINUS, 46.0, DECAY_MODE_BETA_MINUS_AND_2NEWTRON);
		initDecayMode2(6, 20, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 72.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(6, 21, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode(6, 22, DECAY_MODE_BETA_MINUS);

		initAtomProperty(7, "N", "nitrogen", "窒素", 10, 24, 12, 17, 29.124);
		initIsotopeProperty(7, 10, 10.04165, 200e-24, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 11, 11.02609, 590e-24, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 12, 12.0186132, 11.000, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 13, 13.00573861, 9.965, HLU_MINUTE, 0.0);
		initIsotopeProperty(7, 14, 14.0030740048, 0.0, HLU_STABLE, 0.99636);
		initIsotopeProperty(7, 15, 15.0001088982, 0.0, HLU_STABLE, 0.00364);
		initIsotopeProperty(7, 16, 16.0061017, 7.13, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 17, 17.008450, 4.173, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 18, 18.014079, 622.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 19, 19.017029, 271.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 20, 20.02337, 130.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 21, 21.02711, 87.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 22, 22.03439, 13.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 23, 23.04122, 14.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(7, 24, 24.05104, 52.0e-9, HLU_SECOND, 0.0);
		initDecayMode(7, 10, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(7, 11, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(7, 12, DECAY_MODE_BETA_PLUS, 96.5, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(7, 13, DECAY_MODE_BETA_PLUS);
		initDecayMode2(7, 16, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode3(7, 17, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 95.0, DECAY_MODE_BETA_MINUS, 4.99, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode3(7, 18, DECAY_MODE_BETA_MINUS, 76.9, DECAY_MODE_BETA_MINUS_AND_ALPHA, 12.2, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(7, 19, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 54.6, DECAY_MODE_BETA_MINUS);
		initDecayMode2(7, 20, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 56.99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(7, 21, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 80.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(7, 22, DECAY_MODE_BETA_MINUS, 65.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(7, 23, DECAY_MODE_BETA_MINUS);
		initDecayMode(7, 24, DECAY_MODE_NEWTRON_EMISSION);

		initAtomProperty(8, "O", "oxygen", "酸素", 10, 24, 14, 22, 29.378);
		initIsotopeProperty(8, 12, 12.034405, 580.0e-24, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 13, 13.024812, 8.58e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 14, 14.00859625, 70.598, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 15, 15.0030656, 122.24, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 16, 15.99491461956, 0.0, HLU_STABLE, 0.99757);
		initIsotopeProperty(8, 17, 16.99913170, 0.0, HLU_STABLE, 3.8e-4);
		initIsotopeProperty(8, 18, 17.9991610, 0.0, HLU_STABLE, 2.05e-3);
		initIsotopeProperty(8, 19, 19.003580, 26.464, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 20, 20.0040767, 13.51, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 21, 21.008656, 3.42, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 22, 22.00997, 2.25, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 23, 23.01569, 82e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(8, 24, 24.02047, 65e-3, HLU_SECOND, 0.0);
		initDecayMode2(8, 12, DECAY_MODE_2PROTON_EMISSION, 60.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(8, 13, DECAY_MODE_BETA_PLUS, 89.1, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(8, 14, DECAY_MODE_BETA_PLUS);
		initDecayMode(8, 15, DECAY_MODE_BETA_PLUS);
		initDecayMode(8, 19, DECAY_MODE_BETA_MINUS);
		initDecayMode(8, 20, DECAY_MODE_BETA_MINUS);
		initDecayMode(8, 21, DECAY_MODE_BETA_MINUS);
		initDecayMode2(8, 22, DECAY_MODE_BETA_MINUS, 78.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(8, 23, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 57.99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(8, 24, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 57.99, DECAY_MODE_BETA_MINUS);

		initAtomProperty(9, "F", "fluorine", "フッ素", 15, 29, 17, 23, 31.304);
		initIsotopeProperty(9, 15, 15.01801, 410.0e-22, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 16, 16.011466, 11.0e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 17, 17.00209524, 64.49, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 18, 18.0009380, 109.771, HLU_MINUTE, 0.0);
		initIsotopeProperty(9, 19, 18.99840322, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(9, 20, 19.99998132, 11.163, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 21, 20.9999490, 4.158, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 22, 22.002999, 4.23, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 23, 23.00357, 2.23, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 24, 24.00812, 400.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 25, 25.01210, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 26, 26.01962, 9.6e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 27, 27.02676, 4.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 28, 28.03567, 40.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(9, 29, 29.04326, 2.6e-3, HLU_SECOND, 0.0);
		initDecayMode(9, 15, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(9, 16, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(9, 17, DECAY_MODE_BETA_PLUS);
		initDecayMode(9, 18, DECAY_MODE_BETA_PLUS);
		initDecayMode(9, 20, DECAY_MODE_BETA_MINUS);
		initDecayMode(9, 21, DECAY_MODE_BETA_MINUS);
		initDecayMode2(9, 22, DECAY_MODE_BETA_MINUS, 89.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(9, 23, DECAY_MODE_BETA_MINUS, 86.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(9, 24, DECAY_MODE_BETA_MINUS, 94.1, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(9, 25, DECAY_MODE_BETA_MINUS, 76.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(9, 26, DECAY_MODE_BETA_MINUS, 68.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(9, 27, DECAY_MODE_BETA_MINUS);
		initDecayMode(9, 28, DECAY_MODE_NEWTRON_EMISSION);
		initDecayMode(9, 29, DECAY_MODE_BETA_MINUS);

		initAtomProperty(10, "Ne", "neon", "ネオン", 16, 32, 18, 24, 20.786);
		initIsotopeProperty(10, 16, 16.025761, 9.0e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 17, 17.017672, 109.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 18, 18.0057082, 1.672, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 19, 19.0018802, 17.296, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 20, 19.9924401754, 0.0, HLU_STABLE, 0.9048);
		initIsotopeProperty(10, 21, 20.99384668, 0.0, HLU_STABLE,   0.0027);
		initIsotopeProperty(10, 22, 21.991385114, 0.0, HLU_STABLE,  0.0925);
		initIsotopeProperty(10, 23, 22.99446690, 37.24, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 24, 23.9936108, 3.38, HLU_MINUTE, 0.0);
		initIsotopeProperty(10, 25, 24.997737, 602.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 26, 26.000461, 197.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 27, 27.00759, 32.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 28, 28.01207, 18.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 29, 29.01939, 15.6e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 30, 30.02480, 5.8e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 31, 31.03311, 3.4e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(10, 32, 32.04002, 3.5e-3, HLU_SECOND, 0.0);
		initDecayMode(10, 16, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode3(10, 17, DECAY_MODE_BETA_PLUS_PROTON, 96.0, DECAY_MODE_BETA_PLUS_ALPHA, 2.7, DECAY_MODE_BETA_PLUS);
		initDecayMode2(10, 18, DECAY_MODE_ELECTRON_CAPTURE, 50.0, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode(10, 19, DECAY_MODE_BETA_PLUS);
		initDecayMode(10, 23, DECAY_MODE_BETA_MINUS);
		initDecayMode(10, 24, DECAY_MODE_BETA_MINUS);
		initDecayMode(10, 25, DECAY_MODE_BETA_MINUS);
		initDecayMode2(10, 26, DECAY_MODE_BETA_MINUS, 99.87, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(10, 27, DECAY_MODE_BETA_MINUS, 98.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(10, 28, DECAY_MODE_BETA_MINUS, 78.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(10, 29, DECAY_MODE_BETA_MINUS);
		initDecayMode(10, 30, DECAY_MODE_BETA_MINUS);
		initDecayMode2(10, 31, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(10, 32, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 50.0, DECAY_MODE_BETA_MINUS);

		initAtomProperty(11, "Na", "sodium", "ナトリウム", 18, 35, 21, 26, 28.230);
		initIsotopeProperty(11, 18, 18.02597, 1.3e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 19, 19.013877, 40.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 20, 20.007351, 447.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 21, 20.9976552, 22.49, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 22, 21.9944364, 2.6027, HLU_YEAR, 0.0);
		initIsotopeProperty(11, 23, 22.9897692809, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(11, 24, 23.99096278, 14.9590, HLU_HOUR, 0.0);
		initIsotopeProperty(11, 25, 24.9899540, 59.1, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 26, 25.992633, 1.077, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 27, 26.994077, 301.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 28, 27.998938, 30.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 29, 29.002861, 44.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 30, 30.008976, 48.4e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 31, 31.01359, 17.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 32, 32.02047, 12.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 33, 33.02672, 8.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 34, 34.03517, 5.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(11, 35, 35.04249, 1.5e-3, HLU_SECOND, 0.0);
		initDecayMode2(11, 18, DECAY_MODE_PROTON_EMISSION, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode(11, 19, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(11, 20, DECAY_MODE_BETA_PLUS, 75.0, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(11, 21, DECAY_MODE_BETA_PLUS);
		initDecayMode(11, 22, DECAY_MODE_BETA_PLUS);
		initDecayMode(11, 24, DECAY_MODE_BETA_MINUS);
		initDecayMode(11, 25, DECAY_MODE_BETA_MINUS);
		initDecayMode(11, 26, DECAY_MODE_BETA_MINUS);
		initDecayMode2(11, 27, DECAY_MODE_BETA_MINUS, 99.87, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(11, 28, DECAY_MODE_BETA_MINUS, 99.421, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(11, 29, DECAY_MODE_BETA_MINUS, 74.09, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode3(11, 30, DECAY_MODE_BETA_MINUS, 68.83, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 30.0, 	DECAY_MODE_BETA_MINUS_AND_2NEWTRON); //	DECAY_MODE_BETA_MINUS_AND_ALPHA
		initDecayMode3(11, 31, DECAY_MODE_BETA_MINUS, 62.05, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 33.0, 	DECAY_MODE_BETA_MINUS_AND_2NEWTRON);// DECAY_MODE_BETA_MINUS_AND_3NEWTRON
		initDecayMode3(11, 32, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 35.0, 	DECAY_MODE_BETA_MINUS_AND_2NEWTRON);
		initDecayMode3(11, 33, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 52.0, DECAY_MODE_BETA_MINUS, 36.0, 	DECAY_MODE_BETA_MINUS_AND_2NEWTRON);
		initDecayMode3(11, 34, DECAY_MODE_BETA_MINUS_AND_2NEWTRON, 50.0, DECAY_MODE_BETA_MINUS, 35.0, 	DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(11, 35, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(12, "Mg", "magnesium", "マグネシウム", 20, 37, 22, 29, 24.869);
		initIsotopeProperty(12, 20, 20.018863, 90.8e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 21, 21.011713, 122.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 22, 21.9995738, 3.8755, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 23, 22.9941237, 11.317, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 24, 23.985041700, 0.0, HLU_STABLE, 0.7899);
		initIsotopeProperty(12, 25, 24.98583692, 0.0, HLU_STABLE, 0.1000);
		initIsotopeProperty(12, 26, 25.982592929, 0.0, HLU_STABLE, 0.1101);
		initIsotopeProperty(12, 27, 26.98434059, 9.458, HLU_MINUTE, 0.0);
		initIsotopeProperty(12, 28, 27.9838768, 20.915, HLU_HOUR, 0.0);
		initIsotopeProperty(12, 29, 28.988600, 1.30, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 30, 29.990434, 335.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 31, 30.996546, 230.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 32, 31.998975, 86.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 33, 33.005254, 90.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 34, 34.00946, 20.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 35, 35.01734, 70.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 36, 36.02300, 3.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(12, 37, 37.03140, 40.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(12, 20, DECAY_MODE_BETA_PLUS, 97.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(12, 21, DECAY_MODE_BETA_PLUS, 66.9, DECAY_MODE_BETA_PLUS_PROTON, 32.6, 	DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(12, 22, DECAY_MODE_BETA_PLUS);
		initDecayMode(12, 23, DECAY_MODE_BETA_PLUS);
		initDecayMode(12, 27, DECAY_MODE_BETA_MINUS);
		initDecayMode(12, 28, DECAY_MODE_BETA_MINUS);
		initDecayMode(12, 29, DECAY_MODE_BETA_MINUS);
		initDecayMode2(12, 30, DECAY_MODE_BETA_MINUS, 99.94, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(12, 31, DECAY_MODE_BETA_MINUS, 98.3, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(12, 32, DECAY_MODE_BETA_MINUS, 97.6, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(12, 33, DECAY_MODE_BETA_MINUS, 83.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(12, 34, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(12, 35, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 52.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(12, 36, DECAY_MODE_BETA_MINUS);
		initDecayMode2(12, 37, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(13, "Al", "aluminium", "アルミニウム", 19, 39, 24, 30, 24.200);
		initIsotopeProperty(13, 19, 19.0218, 35.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 20, 20.0194, 35.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 21, 21.02804, 35.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 22, 22.01952, 59.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 23, 23.007267, 470.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 24, 23.9999389, 2.053, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 25, 24.9904281, 7.183, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 26, 25.98689169, 7.17e5, HLU_YEAR, 0.0);
		initIsotopeProperty(13, 27, 26.98153863, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(13, 28, 27.98191031, 2.2414, HLU_MINUTE, 0.0);
		initIsotopeProperty(13, 29, 28.9804450, 6.56, HLU_MINUTE, 0.0);
		initIsotopeProperty(13, 30, 29.982960, 3.60, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 31, 30.983947, 644.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 32, 31.98812, 31.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 33, 32.99084, 41.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 34, 33.99685, 56.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 35, 34.99986, 38.6e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 36, 36.00621, 90.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 37, 37.01068, 10.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 38, 38.01723, 7.6e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(13, 39, 39.02297, 7.6e-3, HLU_SECOND, 0.0);
		initDecayMode(13, 19, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(13, 20, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(13, 21, DECAY_MODE_PROTON_EMISSION);
		initDecayMode3(13, 22, DECAY_MODE_BETA_PLUS, 96.7, DECAY_MODE_BETA_PLUS_2PROTON, 2.5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(13, 23, DECAY_MODE_BETA_PLUS, 92.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(13, 24, DECAY_MODE_BETA_PLUS, 99.95, DECAY_MODE_BETA_PLUS_ALPHA, 0.0349, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(13, 25, DECAY_MODE_BETA_PLUS);
		initDecayMode(13, 26, DECAY_MODE_BETA_PLUS);
		initDecayMode(13, 28, DECAY_MODE_BETA_MINUS);
		initDecayMode(13, 29, DECAY_MODE_BETA_MINUS);
		initDecayMode(13, 30, DECAY_MODE_BETA_MINUS);
		initDecayMode2(13, 31, DECAY_MODE_BETA_MINUS, 98.4, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(13, 32, DECAY_MODE_BETA_MINUS, 99.3, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(13, 33, DECAY_MODE_BETA_MINUS, 91.5, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(13, 34, DECAY_MODE_BETA_MINUS, 87.5, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(13, 35, DECAY_MODE_BETA_MINUS, 74.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(13, 36, DECAY_MODE_BETA_MINUS, 69.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(13, 37, DECAY_MODE_BETA_MINUS);
		initDecayMode(13, 38, DECAY_MODE_BETA_MINUS);
		initDecayMode(13, 39, DECAY_MODE_BETA_MINUS);

		initAtomProperty(14, "Si", "silicon", "珪素", 22, 42, 26, 34, 19.789 );
		initIsotopeProperty(14, 22, 22.03453, 29.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 23, 23.02552, 42.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 24, 24.011546, 140.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 25, 25.004106, 220.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 26, 25.992330, 2.234, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 27, 26.98670491, 4.16, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 28, 27.9769265325, 0.0, HLU_STABLE, 0.92223);
		initIsotopeProperty(14, 29, 28.976494700, 0.0, HLU_STABLE, 0.04685);
		initIsotopeProperty(14, 30, 29.97377017, 0.0, HLU_STABLE, 0.03092);
		initIsotopeProperty(14, 31, 30.97536323, 157.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(14, 32, 31.97414808, 153.0, HLU_YEAR, 0.0);
		initIsotopeProperty(14, 33, 32.978000, 6.18, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 34, 33.978576, 2.77, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 35, 34.98458, 780.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 36, 35.98660, 0.45, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 37, 36.99294, 90.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 38, 37.99563, 90.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 39, 39.00207, 47.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 40, 40.00587, 33.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 41, 41.01456, 20.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(14, 42, 42.01979, 13.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(14, 22, DECAY_MODE_BETA_PLUS, 68.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(14, 23, DECAY_MODE_BETA_PLUS);
		initDecayMode2(14, 24, DECAY_MODE_BETA_PLUS, 92.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(14, 25, DECAY_MODE_BETA_PLUS, 63.19, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(14, 26, DECAY_MODE_BETA_PLUS);
		initDecayMode(14, 27, DECAY_MODE_BETA_PLUS);
		initDecayMode(14, 31, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 32, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 33, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 34, DECAY_MODE_BETA_MINUS);
		initDecayMode2(14, 35, DECAY_MODE_BETA_MINUS, 94.74, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(14, 36, DECAY_MODE_BETA_MINUS, 88.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(14, 37, DECAY_MODE_BETA_MINUS, 83.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(14, 38, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 50.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 39, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 40, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 41, DECAY_MODE_BETA_MINUS);
		initDecayMode(14, 42, DECAY_MODE_BETA_MINUS);

		initAtomProperty(15, "P", "phosphorus", "燐", 24, 46, 29, 37, 23.824);
		initIsotopeProperty(15, 24, 24.03435, 30.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 25, 25.02026, 30.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 26, 26.01178, 43.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 27, 26.999230, 260.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 28, 27.992315, 270.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 29, 28.9818006, 4.142, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 30, 29.9783138, 2.498, HLU_MINUTE, 0.0);
		initIsotopeProperty(15, 31, 30.97376163, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(15, 32, 31.97390727, 14.263, HLU_DAY, 0.0);
		initIsotopeProperty(15, 33, 32.9717255, 25.34, HLU_DAY, 0.0);
		initIsotopeProperty(15, 34, 33.973636, 12.43, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 35, 34.9733141, 47.3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 36, 35.978260, 5.6, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 37, 36.97961, 2.31, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 38, 37.98416, 0.64, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 39, 38.98618, 190.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 40, 39.99130, 153.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 41, 40.99434, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 42, 42.00101, 48.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 43, 43.00619, 36.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 44, 44.01299, 18.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 45, 45.01922, 8.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(15, 46, 46.02738, 4.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(15, 24, DECAY_MODE_PROTON_EMISSION, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode(15, 25, DECAY_MODE_PROTON_EMISSION);
		initDecayMode3(15, 26, DECAY_MODE_BETA_PLUS, 98.1, DECAY_MODE_BETA_PLUS_2PROTON, 1.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(15, 27, DECAY_MODE_BETA_PLUS, 99.93, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(15, 28, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON, 0.0013, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(15, 29, DECAY_MODE_BETA_PLUS);
		initDecayMode(15, 30, DECAY_MODE_BETA_PLUS);
		initDecayMode(15, 32, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 33, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 34, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 35, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 36, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 37, DECAY_MODE_BETA_MINUS);
		initDecayMode2(15, 38, DECAY_MODE_BETA_MINUS, 88.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(15, 39, DECAY_MODE_BETA_MINUS, 74.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(15, 40, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(15, 41, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(15, 42, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(15, 43, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(15, 44, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 45, DECAY_MODE_BETA_MINUS);
		initDecayMode(15, 46, DECAY_MODE_BETA_MINUS);

		initAtomProperty(16, "S", "sulfur", "硫黄", 26, 49, 30, 42, 22.75);
		initIsotopeProperty(16, 26, 26.02788, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 27, 27.01883, 15.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 28, 28.00437, 125.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 29, 28.99661, 187.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 30, 29.984903, 1.178, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 31, 30.9795547, 2.572, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 32, 31.97207100, 0.0, HLU_STABLE, 0.9493);
		initIsotopeProperty(16, 33, 32.97145876, 0.0, HLU_STABLE, 0.0076);
		initIsotopeProperty(16, 34, 33.96786690, 0.0, HLU_STABLE, 0.0429);
		initIsotopeProperty(16, 35, 34.96903216, 87.51, HLU_DAY, 0.0);
		initIsotopeProperty(16, 36, 35.96708076, 0.0, HLU_STABLE, 2.0e-4);
		initIsotopeProperty(16, 37, 36.97112557, 5.05, HLU_MINUTE, 0.0);
		initIsotopeProperty(16, 38, 37.971163, 170.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(16, 39, 38.97513, 11.5, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 40, 39.97545, 8.8, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 41, 40.97958, 1.99, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 42, 41.98102, 1.013, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 43, 42.98715, 260.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 44, 43.99021, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 45, 44.99651, 68.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 46, 46.00075, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 47, 47.00859, 20.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 48, 48.01417, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(16, 49, 49.02362, 200.0e-9, HLU_SECOND, 0.0);
		initDecayMode(16, 26, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode3(16, 27, DECAY_MODE_BETA_PLUS, 98.0, DECAY_MODE_BETA_PLUS_2PROTON, 1.91, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(16, 28, DECAY_MODE_BETA_PLUS, 79.3, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(16, 29, DECAY_MODE_BETA_PLUS, 53.6, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(16, 30, DECAY_MODE_BETA_PLUS);
		initDecayMode(16, 31, DECAY_MODE_BETA_PLUS);
		initDecayMode(16, 35, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 37, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 38, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 39, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 40, DECAY_MODE_BETA_MINUS);
		initDecayMode2(16, 41, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(16, 42, DECAY_MODE_BETA_MINUS, 96.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(16, 43, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(16, 44, DECAY_MODE_BETA_MINUS, 82.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(16, 45, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 54.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 46, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 47, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 48, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 49, DECAY_MODE_NEWTRON_EMISSION);

		initAtomProperty(17, "Cl", "chlorine", "塩素", 29, 51, 33, 43, 33.949);
		initIsotopeProperty(17, 29, 29.01411, 20.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 30, 30.00477, 30.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 31, 30.99241, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 32, 31.985690, 298.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 33, 32.9774519, 2.511, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 34, 33.97376282, 1.5264, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 35, 34.96885268, 0.0, HLU_STABLE, 0.7576);
		initIsotopeProperty(17, 36, 35.96830698, 3.01e5, HLU_YEAR, 0.0);
		initIsotopeProperty(17, 37, 36.96590259, 0.0, HLU_STABLE, 0.2424);
		initIsotopeProperty(17, 38, 37.96801043, 37.24, HLU_MINUTE, 0.0);
		initIsotopeProperty(17, 39, 38.9680082, 55.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(17, 40, 39.97042, 1.35, HLU_MINUTE, 0.0);
		initIsotopeProperty(17, 41, 40.97068, 38.4, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 42, 41.97325, 6.8, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 43, 42.97405, 3.07, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 44, 43.97828, 0.56, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 46, 45.98421, 232.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 48, 47.99495, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 49, 49.00032, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 50, 50.00784, 20.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(17, 51, 51.01449, 2.0e-3, HLU_SECOND, 0.0);
		initDecayMode(17, 29, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(17, 30, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(17, 31, DECAY_MODE_BETA_PLUS, 99.3, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(17, 32, DECAY_MODE_BETA_PLUS, 99.92, DECAY_MODE_BETA_PLUS_ALPHA, 0.054, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(17, 33, DECAY_MODE_BETA_PLUS);
		initDecayMode(17, 34, DECAY_MODE_BETA_PLUS);
		initDecayMode2(17, 36, DECAY_MODE_BETA_MINUS, 98.1, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(17, 38, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 39, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 40, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 41, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 42, DECAY_MODE_BETA_MINUS);
		initDecayMode2(17, 43, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(17, 44, DECAY_MODE_BETA_MINUS, 92.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(17, 45, DECAY_MODE_BETA_MINUS, 76.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(17, 46, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 60.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(17, 47, DECAY_MODE_BETA_MINUS, 97.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(17, 48, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 49, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 50, DECAY_MODE_BETA_MINUS);
		initDecayMode(17, 51, DECAY_MODE_BETA_MINUS);

		initAtomProperty(18, "Ar", "argon", "アルゴン", 30, 53, 35, 47, 20.786 );
		initIsotopeProperty(18, 30, 30.02156, 20.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 31, 31.01212, 14.4e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 32, 31.9976380, 98.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 33, 32.9899257, 173.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 34, 33.9802712, 844.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 35, 34.9752576, 1.775, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 36, 35.967545106, 0.0, HLU_STABLE, 0.003336);
		initIsotopeProperty(18, 37, 36.96677632, 35.04, HLU_DAY, 0.0);
		initIsotopeProperty(18, 38, 37.9627324, 0.0, HLU_STABLE, 6.29e-4);
		initIsotopeProperty(18, 39, 38.964313, 269.0, HLU_YEAR, 0.0);
		initIsotopeProperty(18, 40, 39.9623831225, 0.0, HLU_STABLE, 0.996035);
		initIsotopeProperty(18, 41, 40.9645006, 109.61, HLU_MINUTE, 0.0);
		initIsotopeProperty(18, 42, 41.963046, 32.9, HLU_YEAR, 0.0);
		initIsotopeProperty(18, 43, 42.965636, 5.37, HLU_MINUTE, 0.0);
		initIsotopeProperty(18, 44, 43.9649240, 11.87, HLU_MINUTE, 0.0);
		initIsotopeProperty(18, 45, 44.9680400, 21.48, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 46, 45.96809, 8.4, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 47, 46.97219, 1.23, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 48, 47.97454, 0.48, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 49, 48.98052, 170.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 50, 49.98443, 85.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 51, 50.99163, 60.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 52, 51.99678, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(18, 53, 53.00494, 3.0e-3, HLU_SECOND, 0.0);
		initDecayMode(18, 30, DECAY_MODE_PROTON_EMISSION);
		initDecayMode4(18, 31, DECAY_MODE_BETA_PLUS_PROTON, 55.0, DECAY_MODE_BETA_PLUS, 40.4, DECAY_MODE_BETA_PLUS_2PROTON, 2.48, DECAY_MODE_BETA_PLUS_3PROTON);
		initDecayMode2(18, 32, DECAY_MODE_BETA_PLUS, 56.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(18, 33, DECAY_MODE_BETA_PLUS, 61.35, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(18, 34, DECAY_MODE_BETA_PLUS);
		initDecayMode(18, 35, DECAY_MODE_BETA_PLUS);
		initDecayMode(18, 37, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(18, 39, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 41, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 42, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 43, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 44, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 45, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 46, DECAY_MODE_BETA_MINUS);
		initDecayMode2(18, 47, DECAY_MODE_BETA_MINUS, 99.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(18, 48, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 49, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 50, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 51, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 52, DECAY_MODE_BETA_MINUS);
		initDecayMode2(18, 53, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(19, "K", "potassium", "カリウム", 33, 55, 37, 49, 29.6);
		initIsotopeProperty(19, 33, 33.00726, 25.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 34, 33.99841, 25.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 35, 34.988010, 178.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 36, 35.981292, 342.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 37, 36.97337589, 1.226, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 38, 37.9690812, 7.636, HLU_MINUTE, 0.0);
		initIsotopeProperty(19, 39, 38.96370668, 0.0, HLU_STABLE,   0.932581);
		initIsotopeProperty(19, 40, 39.96399848, 1.248e9, HLU_YEAR, 0.000117);
		initIsotopeProperty(19, 41, 40.96182576, 0.0, HLU_STABLE,   0.067302);
		initIsotopeProperty(19, 42, 41.96240281, 12.360, HLU_HOUR, 0.0);
		initIsotopeProperty(19, 43, 42.960716, 22.3, HLU_HOUR, 0.0);
		initIsotopeProperty(19, 44, 43.96156, 22.13, HLU_MINUTE, 0.0);
		initIsotopeProperty(19, 45, 44.960699, 17.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(19, 46, 45.961977, 105.0, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 47, 46.961678, 17.50, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 48, 47.965514, 6.8, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 49, 48.96745, 1.26, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 50, 49.97278, 472.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 51, 50.97638, 365.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 52, 51.98261, 105.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 53, 52.98712, 30.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 54, 53.99420, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(19, 55, 54.99971, 3.0e-3, HLU_SECOND, 0.0);
		initDecayMode(19, 33, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(19, 34, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(19, 35, DECAY_MODE_BETA_PLUS, 99.63, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(19, 36, DECAY_MODE_BETA_PLUS, 99.94, DECAY_MODE_BETA_PLUS_PROTON, 0.048, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(19, 37, DECAY_MODE_BETA_PLUS);
		initDecayMode(19, 38, DECAY_MODE_BETA_PLUS);
		initDecayMode3(19, 40, DECAY_MODE_BETA_MINUS, 89.28, DECAY_MODE_ELECTRON_CAPTURE, 10.719, DECAY_MODE_BETA_PLUS);
		initDecayMode(19, 42, DECAY_MODE_BETA_MINUS);
		initDecayMode(19, 43, DECAY_MODE_BETA_MINUS);
		initDecayMode(19, 44, DECAY_MODE_BETA_MINUS);
		initDecayMode(19, 45, DECAY_MODE_BETA_MINUS);
		initDecayMode(19, 46, DECAY_MODE_BETA_MINUS);
		initDecayMode(19, 47, DECAY_MODE_BETA_MINUS);
		initDecayMode2(19, 48, DECAY_MODE_BETA_MINUS, 98.86, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(19, 49, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 86.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(19, 50, DECAY_MODE_BETA_MINUS, 71.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(19, 51, DECAY_MODE_BETA_MINUS, 53.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode3(19, 52, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 64.0, DECAY_MODE_BETA_MINUS_AND_2NEWTRON, 21.0, DECAY_MODE_BETA_MINUS);
		initDecayMode3(19, 53, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 67.0, DECAY_MODE_BETA_MINUS_AND_2NEWTRON, 17.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(19, 54, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(19, 55, DECAY_MODE_BETA_MINUS, 99.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(20, "Ca", "calcium", "カルシウム", 34, 57, 40, 52, 25.929);
		initIsotopeProperty(20, 34, 34.01412, 35.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 35, 35.00494, 25.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 36, 35.99309, 102.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 37, 36.985870, 181.1e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 38, 37.976318, 440.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 39, 38.9707197, 859.6e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 40, 39.96259098, 0.0, HLU_STABLE, 0.96941);
		initIsotopeProperty(20, 41, 40.96227806, 1.02e5, HLU_YEAR, 0.0);
		initIsotopeProperty(20, 42, 41.95861801, 0.0, HLU_STABLE, 0.00647);
		initIsotopeProperty(20, 43, 42.9587666, 0.0, HLU_STABLE, 0.00135);
		initIsotopeProperty(20, 44, 43.9554818, 0.0, HLU_STABLE, 0.02086);
		initIsotopeProperty(20, 45, 44.9561866, 162.67, HLU_DAY, 0.0);
		initIsotopeProperty(20, 46, 45.9536926, 0.0, HLU_STABLE, 4.0e-5);
		initIsotopeProperty(20, 47, 46.9545460, 4.536, HLU_DAY, 0.0);
		initIsotopeProperty(20, 48, 47.952534, 43.0e18, HLU_YEAR, 0.00187);
		initIsotopeProperty(20, 49, 48.955674, 8.718, HLU_MINUTE, 0.0);
		initIsotopeProperty(20, 50, 49.957519, 13.9, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 51, 50.9615, 10.0, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 52, 51.96510, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 53, 52.97005, 90.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 54, 53.97435, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 55, 54.98055, 30.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 56, 55.98557, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(20, 57, 56.99236, 5.0e-3, HLU_SECOND, 0.0);
		initDecayMode(20, 34, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(20, 35, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(20, 36, DECAY_MODE_BETA_PLUS_PROTON, 56.8, DECAY_MODE_BETA_PLUS);
		initDecayMode2(20, 37, DECAY_MODE_BETA_PLUS_PROTON, 74.5, DECAY_MODE_BETA_PLUS);
		initDecayMode(20, 38, DECAY_MODE_BETA_PLUS);
		initDecayMode(20, 39, DECAY_MODE_BETA_PLUS);
		initDecayMode(20, 41, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(20, 45, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 47, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 48, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(20, 49, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 50, DECAY_MODE_BETA_MINUS);
		initDecayMode2(20, 51, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(20, 52, DECAY_MODE_BETA_MINUS, 98.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(20, 53, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(20, 54, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 50.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 55, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 56, DECAY_MODE_BETA_MINUS);
		initDecayMode2(20, 57, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(21, "Sc", "scandium", "スカンジウム", 38, 59, 43, 53, 25.52);
		initIsotopeProperty(21, 38, 37.99470, 300.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 39, 38.984790, 300.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 40, 39.977967, 182.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 41, 40.96925113, 596.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 42, 41.96551643, 681.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 43, 42.9611507, 3.891, HLU_HOUR, 0.0);
		initIsotopeProperty(21, 44, 43.9594028, 3.97, HLU_HOUR, 0.0);
		initIsotopeProperty(21, 45, 44.9559119, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(21, 46, 45.9551719, 83.79, HLU_DAY, 0.0);
		initIsotopeProperty(21, 47, 46.9524075, 3.3492, HLU_DAY, 0.0);
		initIsotopeProperty(21, 48, 47.952131, 43.67, HLU_HOUR, 0.0);
		initIsotopeProperty(21, 49, 48.950024, 57.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(21, 50, 49.952188, 102.5, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 51, 50.953603, 12.4, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 52, 51.95668, 8.2, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 53, 52.95961, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 54, 53.96326, 260.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 55, 54.96824, 0.115, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 56, 55.97287, 35.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 57, 56.97779, 13.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 58, 57.98371, 12.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(21, 59, 58.98922, 10.0e-3, HLU_SECOND, 0.0);
		initDecayMode(21, 38, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(21, 38, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(21, 39, DECAY_MODE_PROTON_EMISSION);
		initDecayMode3(21, 40, DECAY_MODE_BETA_PLUS, 99.543, DECAY_MODE_BETA_PLUS_PROTON, 0.44, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(21, 41, DECAY_MODE_BETA_PLUS);
		initDecayMode(21, 42, DECAY_MODE_BETA_PLUS);
		initDecayMode(21, 43, DECAY_MODE_BETA_PLUS);
		initDecayMode(21, 44, DECAY_MODE_BETA_PLUS);
		initDecayMode(21, 46, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 47, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 48, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 49, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 50, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 51, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 52, DECAY_MODE_BETA_MINUS);
		initDecayMode2(21, 53, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(21,54, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(21,55, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(21, 56, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 57, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 58, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 59, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(22, "Ti", "titanium", "チタン", 38, 61, 44, 54, 25.060);
		initIsotopeProperty(22, 38, 38.00977, 120.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 39, 39.00161, 31.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 40, 39.99050, 53.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 41, 40.98315, 80.4e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 42, 41.973031, 199.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 43, 42.968522, 509.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 44, 43.9596901, 60.0, HLU_YEAR, 0.0);
		initIsotopeProperty(22, 45, 44.9581256, 184.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(22, 46, 45.9526316, 0.0, HLU_STABLE, 0.0825);
		initIsotopeProperty(22, 47, 46.9517631, 0.0, HLU_STABLE, 0.0744);
		initIsotopeProperty(22, 48, 47.9479463, 0.0, HLU_STABLE, 0.7372);
		initIsotopeProperty(22, 49, 48.9478700, 0.0, HLU_STABLE, 0.0541);
		initIsotopeProperty(22, 50, 49.9447912, 0.0, HLU_STABLE, 0.0518);
		initIsotopeProperty(22, 51, 50.946615, 5.76, HLU_MINUTE, 0.0);
		initIsotopeProperty(22, 52, 51.946897, 1.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(22, 53, 52.94973, 32.7, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 54, 53.95105, 1.5, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 55, 54.95527, 490.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 56, 55.95820, 164.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 57, 56.96399, 60.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 58, 57.96697, 54.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 59, 58.97293, 30.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 60, 59.97676, 22.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(22, 61, 60.98320, 10.0e-3, HLU_SECOND, 0.0);
		initDecayMode(22, 38, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode3(22, 39, DECAY_MODE_BETA_PLUS_PROTON, 85.0, DECAY_MODE_BETA_PLUS, 14.9, DECAY_MODE_BETA_PLUS_2PROTON);
		initDecayMode2(22, 40, DECAY_MODE_BETA_PLUS, 56.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(22, 41, DECAY_MODE_BETA_PLUS_PROTON, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode(22, 42, DECAY_MODE_BETA_PLUS);
		initDecayMode(22, 43, DECAY_MODE_BETA_PLUS);
		initDecayMode(22, 44, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(22, 45, DECAY_MODE_BETA_PLUS);
		initDecayMode(22, 51, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 52, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 53, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 54, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 55, DECAY_MODE_BETA_MINUS);
		initDecayMode2(22, 56, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(22, 57, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(22, 58, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 59, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 60, DECAY_MODE_BETA_MINUS);
		initDecayMode2(22, 61, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(23, "V", "vanadium", "バナジウム", 42, 63, 47, 55, 24.89);
		initIsotopeProperty(23, 42, 41.99123, 55.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 43, 42.98065, 80.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 44, 43.97411, 111.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 45, 44.965776, 547.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 46, 45.9602005, 422.50e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 47, 46.9549089, 32.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(23, 48, 47.9522537, 15.9735, HLU_DAY, 0.0);
		initIsotopeProperty(23, 49, 48.9485161, 329.0, HLU_DAY, 0.0);
		initIsotopeProperty(23, 50, 49.9471585, 1.4e17, HLU_YEAR, 0.0);
		initIsotopeProperty(23, 51, 50.9439595, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(23, 52, 51.9447755, 3.743, HLU_MINUTE, 0.0);
		initIsotopeProperty(23, 53, 52.944338, 1.60, HLU_MINUTE, 0.0);
		initIsotopeProperty(23, 54, 53.946440, 49.8, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 55, 54.94723, 6.54, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 56, 55.95053, 216.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 57, 56.95256, 0.35, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 58, 57.95683, 191.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 59, 58.96021, 75.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 60, 59.96503, 122.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 61, 60.96848, 47.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 62, 61.97378, 33.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 63, 62.97755, 17.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(23, 59, 58.94859, 0.460, HLU_SECOND, 0.0);
		initDecayMode(23, 42, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(23, 43, DECAY_MODE_BETA_PLUS);
		initDecayMode2(23, 44, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode(23, 45, DECAY_MODE_BETA_PLUS);
		initDecayMode(23, 46, DECAY_MODE_BETA_PLUS);
		initDecayMode(23, 47, DECAY_MODE_BETA_PLUS);
		initDecayMode(23, 48, DECAY_MODE_BETA_PLUS);
		initDecayMode(23, 49, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(23, 50, DECAY_MODE_BETA_PLUS, 83.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(23, 52, DECAY_MODE_BETA_MINUS);
		initDecayMode(23, 53, DECAY_MODE_BETA_MINUS);
		initDecayMode(23, 54, DECAY_MODE_BETA_MINUS);
		initDecayMode(23, 55, DECAY_MODE_BETA_MINUS);
		initDecayMode2(23, 56, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(23, 57, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(23, 58, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(23, 59, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(23, 60, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(23, 61, DECAY_MODE_BETA_MINUS);
		initDecayMode(23, 62, DECAY_MODE_BETA_MINUS);
		initDecayMode(23, 63, DECAY_MODE_BETA_MINUS);

		initAtomProperty(24, "Cr", "chromium", "クロム", 46, 59, 46, 58, 23.35);
		initIsotopeProperty(24, 46, 45.968359, 0.26, HLU_SECOND, 0.0);
		initIsotopeProperty(24, 47, 46.962900, 0.5, HLU_SECOND, 0.0);
		initIsotopeProperty(24, 48, 47.954032, 21.56, HLU_HOUR, 0.0);
		initIsotopeProperty(24, 49, 48.9513357, 42.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(24, 50, 49.9460442, 0.0, HLU_STABLE, 0.04345);
		initIsotopeProperty(24, 51, 50.9447674, 27.7025, HLU_DAY, 0.0);
		initIsotopeProperty(24, 52, 51.9405075, 0.0, HLU_STABLE, 0.83789);
		initIsotopeProperty(24, 53, 52.9406494, 0.0, HLU_STABLE, 0.09501);
		initIsotopeProperty(24, 54, 53.9388804, 0.0, HLU_STABLE, 0.02365);
		initIsotopeProperty(24, 55, 54.9408397, 3.497, HLU_MINUTE, 0.0);
		initIsotopeProperty(24, 56, 55.9406531, 5.94, HLU_MINUTE, 0.0);
		initIsotopeProperty(24, 57, 56.943613, 21.1, HLU_SECOND, 0.0);
		initIsotopeProperty(24, 58, 57.94435, 7.0, HLU_SECOND, 0.0);
		initIsotopeProperty(24, 59, 58.94859, 0.460, HLU_SECOND, 0.0);
		initDecayMode(24, 46, DECAY_MODE_BETA_PLUS);
		initDecayMode(24, 47, DECAY_MODE_BETA_PLUS);
		initDecayMode(24, 48, DECAY_MODE_BETA_PLUS);
		initDecayMode(24, 49, DECAY_MODE_BETA_PLUS);
		initDecayMode(24, 51, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(24, 55, DECAY_MODE_BETA_MINUS);
		initDecayMode(24, 56, DECAY_MODE_BETA_MINUS);
		initDecayMode(24, 57, DECAY_MODE_BETA_MINUS);
		initDecayMode(24, 58, DECAY_MODE_BETA_MINUS);
		initDecayMode(24, 59, DECAY_MODE_BETA_MINUS);

		initAtomProperty(25, "Mn", "manganese", "マンガン", 50, 61, 50, 61, 26.32);
		initIsotopeProperty(25, 50, 49.9542382, 0.28329, HLU_SECOND, 0.0);
		initIsotopeProperty(25, 51, 50.9482108, 46.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(25, 52, 51.9455655, 5.591, HLU_DAY, 0.0);
		initIsotopeProperty(25, 53, 52.9412901, 3.7e6, HLU_YEAR, 0.0);
		initIsotopeProperty(25, 54, 53.9403589, 312.03, HLU_DAY, 0.0);
		initIsotopeProperty(25, 55, 54.9380451, 0.0, HLU_STABLE, 1.0000);
		initIsotopeProperty(25, 56, 55.9389049, 2.5789, HLU_HOUR, 0.0);
		initIsotopeProperty(25, 57, 56.9382854, 85.4, HLU_SECOND, 0.0);
		initIsotopeProperty(25, 58, 57.93998, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(25, 59, 58.94044, 4.59, HLU_SECOND, 0.0);
		initIsotopeProperty(25, 60, 59.94291, 51.0, HLU_SECOND, 0.0);
		initIsotopeProperty(25, 61, 60.94465, 0.67, HLU_SECOND, 0.0);
		initDecayMode(25, 50, DECAY_MODE_BETA_PLUS);
		initDecayMode(25, 51, DECAY_MODE_BETA_PLUS);
		initDecayMode(25, 52, DECAY_MODE_BETA_PLUS);
		initDecayMode(25, 53, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode3(25, 54, DECAY_MODE_ELECTRON_CAPTURE, 99.99, DECAY_MODE_BETA_MINUS, 2.9e-4, DECAY_MODE_BETA_PLUS);
		initDecayMode(25, 56, DECAY_MODE_BETA_MINUS);
		initDecayMode(25, 57, DECAY_MODE_BETA_MINUS);
		initDecayMode(25, 58, DECAY_MODE_BETA_MINUS);
		initDecayMode(25, 59, DECAY_MODE_BETA_MINUS);
		initDecayMode(25, 60, DECAY_MODE_BETA_MINUS);
		initDecayMode(25, 61, DECAY_MODE_BETA_MINUS);
		
		initAtomProperty(26, "Fe", "iron", "鉄", 52, 65, 52, 65, 25.10);
		initIsotopeProperty(26, 52, 51.948114, 8.275, HLU_HOUR, 0.0);
		initIsotopeProperty(26, 53, 52.9453079, 8.51, HLU_MINUTE, 0.0);
		initIsotopeProperty(26, 54, 53.9396105, 0.0, HLU_STABLE, 0.05845);
		initIsotopeProperty(26, 55, 54.9382934, 2.737, HLU_YEAR, 0.0);
		initIsotopeProperty(26, 56, 55.9349375, 0.0, HLU_STABLE, 0.91754);
		initIsotopeProperty(26, 57, 56.9353940, 0.0, HLU_STABLE, 0.02119);
		initIsotopeProperty(26, 58, 57.9332756, 0.0, HLU_STABLE, 0.00282);
		initIsotopeProperty(26, 59, 58.9348755, 44.495, HLU_DAY, 0.0);
		initIsotopeProperty(26, 60, 59.934072, 1.5, HLU_YEAR, 0.0);
		initIsotopeProperty(26, 61, 60.936745, 5.98, HLU_MINUTE, 0.0);
		initIsotopeProperty(26, 62, 61.936767, 68, HLU_SECOND, 0.0);
		initIsotopeProperty(26, 63, 62.94037, 6.1, HLU_SECOND, 0.0);
		initIsotopeProperty(26, 64, 63.9412, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(26, 65, 64.94538, 1.3, HLU_SECOND, 0.0);
		initDecayMode(26, 52, DECAY_MODE_BETA_PLUS);
		initDecayMode(26, 53, DECAY_MODE_BETA_PLUS);
		initDecayMode(26, 55, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(26, 59, DECAY_MODE_BETA_MINUS);
		initDecayMode(26, 60, DECAY_MODE_BETA_MINUS);
		initDecayMode(26, 61, DECAY_MODE_BETA_MINUS);
		initDecayMode(26, 62, DECAY_MODE_BETA_MINUS);
		initDecayMode(26, 63, DECAY_MODE_BETA_MINUS);
		initDecayMode(26, 64, DECAY_MODE_BETA_MINUS);
		initDecayMode(26, 65, DECAY_MODE_BETA_MINUS);

		initAtomProperty(27, "Co", "cobalt", "コバルト", 54, 65, 54, 65, 24.81);
		initIsotopeProperty(27, 54, 53.9484596, 0.193, HLU_SECOND, 0.0);
		initIsotopeProperty(27, 55, 54.9419990, 17.53, HLU_HOUR, 0.0);
		initIsotopeProperty(27, 56, 55.9398393, 77.233, HLU_DAY, 0.0);
		initIsotopeProperty(27, 57, 56.9362914, 271.74, HLU_DAY, 0.0);
		initIsotopeProperty(27, 58, 57.9357528, 70.86, HLU_DAY, 0.0);
		initIsotopeProperty(27, 59, 58.9331950, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(27, 60, 59.9338171, 5.2713, HLU_YEAR, 0.0);
		initIsotopeProperty(27, 61, 60.9324758, 1.650, HLU_HOUR, 0.0);
		initIsotopeProperty(27, 62, 61.934051, 1.50, HLU_MINUTE, 0.0);
		initIsotopeProperty(27, 63, 62.933612, 26.9, HLU_SECOND, 0.0);
		initIsotopeProperty(27, 64, 63.935810, 0.30, HLU_SECOND, 0.0);
		initIsotopeProperty(27, 65, 64.936478, 1.20, HLU_SECOND, 0.0);
		initDecayMode(27, 54, DECAY_MODE_BETA_PLUS);
		initDecayMode(27, 55, DECAY_MODE_BETA_PLUS);
		initDecayMode(27, 56, DECAY_MODE_BETA_PLUS);
		initDecayMode(27, 57, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(27, 58, DECAY_MODE_BETA_PLUS);
		initDecayMode(27, 60, DECAY_MODE_BETA_MINUS);
		initDecayMode(27, 61, DECAY_MODE_BETA_MINUS);
		initDecayMode(27, 62, DECAY_MODE_BETA_MINUS);
		initDecayMode(27, 63, DECAY_MODE_BETA_MINUS);
		initDecayMode(27, 64, DECAY_MODE_BETA_MINUS);
		initDecayMode(27, 65, DECAY_MODE_BETA_MINUS);
		
		initAtomProperty(28, "Ni", "nickel", "ニッケル", 56, 72, 56, 72, 26.07);
		initIsotopeProperty(28, 56, 55.942132, 6.075, HLU_DAY, 0.0);
		initIsotopeProperty(28, 57, 56.9397935, 35.60, HLU_HOUR, 0.0);
		initIsotopeProperty(28, 58, 57.9353429, 0.0, HLU_STABLE, 0.680769);
		initIsotopeProperty(28, 59, 58.9343467, 7.6, HLU_YEAR, 0.0);
		initIsotopeProperty(28, 60, 59.9307864, 0.0, HLU_STABLE, 0.262231);
		initIsotopeProperty(28, 61, 60.9310560, 0.0, HLU_STABLE, 0.011399);
		initIsotopeProperty(28, 62, 61.9283451, 0.0, HLU_STABLE, 0.036345);
		initIsotopeProperty(28, 63, 62.9296694, 100.1, HLU_YEAR, 0.0);
		initIsotopeProperty(28, 64, 63.9279660, 0.0, HLU_STABLE, 0.009256);
		initIsotopeProperty(28, 65, 64.9300843, 2.5172, HLU_HOUR, 0.0);
		initIsotopeProperty(28, 66, 65.9291393, 54.6, HLU_HOUR, 0.0);
		initIsotopeProperty(28, 67, 66.931569, 21.0, HLU_SECOND, 0.0);
		initIsotopeProperty(28, 68, 67.931869, 29.0, HLU_SECOND, 0.0);
		initIsotopeProperty(28, 69, 68.935610, 11.5, HLU_SECOND, 0.0);
		initIsotopeProperty(28, 70, 69.93650, 6.0, HLU_SECOND, 0.0);
		initIsotopeProperty(28, 71, 70.94074, 2.56, HLU_SECOND, 0.0);
		initIsotopeProperty(28, 72, 71.94209, 1.57, HLU_SECOND, 0.0);
		initDecayMode(28, 56, DECAY_MODE_BETA_PLUS);
		initDecayMode(28, 57, DECAY_MODE_BETA_PLUS);
		initDecayMode(28, 59, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(28, 63, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 65, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 66, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 67, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 68, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 69, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 70, DECAY_MODE_BETA_MINUS);
		initDecayMode(28, 71, DECAY_MODE_BETA_MINUS);
		initDecayMode2(28, 72, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(29, "Cu", "copper", "銅", 58, 76, 58, 76, 24.440);
		initIsotopeProperty(29, 58, 57.9445385, 3.204, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 59, 58.9394980, 81.5, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 60, 59.9373650, 23.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(29, 61, 60.9334578, 3.333, HLU_HOUR, 0.0);
		initIsotopeProperty(29, 62, 61.932584, 9.673, HLU_MINUTE, 0.0);
		initIsotopeProperty(29, 63, 62.9295975, 0.0, HLU_STABLE, 0.6915);
		initIsotopeProperty(29, 64, 63.9297642, 12.7, HLU_HOUR, 0.0);
		initIsotopeProperty(29, 65, 64.9277895, 0.0, HLU_STABLE, 0.3085);
		initIsotopeProperty(29, 66, 65.9288688, 5.120, HLU_MINUTE, 0.0);
		initIsotopeProperty(29, 67, 66.9277303, 61.83, HLU_HOUR, 0.0);
		initIsotopeProperty(29, 68, 67.9296109, 31.13, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 69, 68.9294293, 2.85, HLU_MINUTE, 0.0);
		initIsotopeProperty(29, 70, 69.9323923, 44.5, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 71, 70.9326768, 19.4, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 72, 71.9358203, 6.6, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 73, 72.936675, 4.2, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 74, 73.939875, 1.594, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 75, 74.94190, 1.224, HLU_SECOND, 0.0);
		initIsotopeProperty(29, 76, 75.945275, 0.641, HLU_SECOND, 0.0);
		initDecayMode(29, 58, DECAY_MODE_BETA_PLUS);
		initDecayMode(29, 59, DECAY_MODE_BETA_PLUS);
		initDecayMode(29, 60, DECAY_MODE_BETA_PLUS);
		initDecayMode(29, 61, DECAY_MODE_BETA_PLUS);
		initDecayMode(29, 62, DECAY_MODE_BETA_PLUS);
		initDecayMode2(29, 64, DECAY_MODE_BETA_PLUS, 61.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 66, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 67, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 68, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 69, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 70, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 71, DECAY_MODE_BETA_MINUS);
		initDecayMode(29, 72, DECAY_MODE_BETA_MINUS);
		initDecayMode2(29, 73, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(29, 74, DECAY_MODE_BETA_MINUS);
		initDecayMode2(29, 75, DECAY_MODE_BETA_MINUS, 96.5, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(29, 76, DECAY_MODE_BETA_MINUS, 97.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		
		initAtomProperty(30, "Zn", "zinc", "亜鉛", 60, 70, 60, 70, 25.470);
		initIsotopeProperty(30, 60, 59.941827, 2.38, HLU_MINUTE, 0.0);
		initIsotopeProperty(30, 61, 60.939511, 89.1, HLU_SECOND, 0.0);
		initIsotopeProperty(30, 62, 61.934330, 9.186, HLU_HOUR, 0.0);
		initIsotopeProperty(30, 63, 62.9332116, 38.47, HLU_MINUTE, 0.0);
		initIsotopeProperty(30, 64, 63.9291422, 0.0, HLU_STABLE, 0.48268);
		initIsotopeProperty(30, 65, 64.9292410, 243.66, HLU_DAY, 0.0);
		initIsotopeProperty(30, 66, 65.9260334, 0.0, HLU_STABLE, 0.27975);
		initIsotopeProperty(30, 67, 66.9271273, 0.0, HLU_STABLE, 0.04102);
		initIsotopeProperty(30, 68, 67.9248442, 0.0, HLU_STABLE, 0.19024);
		initIsotopeProperty(30, 69, 68.9265503, 56.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(30, 70, 69.9253193, 0.0, HLU_STABLE, 0.00631);
		initDecayMode(30, 60, DECAY_MODE_BETA_PLUS);
		initDecayMode(30, 61, DECAY_MODE_BETA_PLUS);
		initDecayMode(30, 62, DECAY_MODE_BETA_PLUS);
		initDecayMode(30, 63, DECAY_MODE_BETA_PLUS);
		initDecayMode(30, 65, DECAY_MODE_BETA_PLUS);
		initDecayMode(30, 69, DECAY_MODE_BETA_MINUS);

		initAtomProperty(31, "Ga", "gallium", "ガリウム", 60, 84, 60, 84, 25.86);
		initIsotopeProperty(31, 60, 59.95706, 70.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 61, 60.94945, 168.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 62, 61.944175, 116.18e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 63, 62.9392942, 32.4, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 64, 63.9368387, 2.627, HLU_MINUTE, 0.0);
		initIsotopeProperty(31, 65, 64.9327348, 15.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(31, 66, 65.931589, 9.49, HLU_HOUR, 0.0);
		initIsotopeProperty(31, 67, 66.9282017, 3.2612, HLU_DAY, 0.0);
		initIsotopeProperty(31, 68, 67.9279801, 67.71, HLU_MINUTE, 0.0);
		initIsotopeProperty(31, 69, 68.9255736, 0.0, HLU_STABLE, 0.60108);
		initIsotopeProperty(31, 70, 69.9260220, 21.14, HLU_MINUTE, 0.0);	
		initIsotopeProperty(31, 71, 70.9247013, 0.0, HLU_STABLE, 0.39892);
		initIsotopeProperty(31, 72, 71.9263663, 14.095, HLU_HOUR, 0.0);
		initIsotopeProperty(31, 73, 72.9251747, 4.86, HLU_HOUR, 0.0);
		initIsotopeProperty(31, 74, 73.926946, 8.12, HLU_MINUTE, 0.0);
		initIsotopeProperty(31, 75, 74.9265002, 126.0, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 76, 75.9288276, 32.6, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 77, 76.9291543, 13.2, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 78, 77.9316082, 5.09, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 79, 78.93289, 2.847, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 80, 79.93652, 1.697, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 81, 80.93775, 1.217, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 82, 81.94299, 0.599, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 83, 82.94698, 308.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(31, 84, 83.95265, 0.085, HLU_SECOND, 0.0);
		initDecayMode(31, 60, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 61, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 62, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 63, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 64, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 65, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 66, DECAY_MODE_BETA_PLUS);
		initDecayMode(31, 67, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(31, 68, DECAY_MODE_BETA_PLUS);
		initDecayMode2(31, 70, DECAY_MODE_BETA_MINUS, 99.59, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(31, 72, DECAY_MODE_BETA_MINUS);
		initDecayMode(31, 73, DECAY_MODE_BETA_MINUS);
		initDecayMode(31, 74, DECAY_MODE_BETA_MINUS);
		initDecayMode(31, 75, DECAY_MODE_BETA_MINUS);
		initDecayMode(31, 76, DECAY_MODE_BETA_MINUS);
		initDecayMode(31, 77, DECAY_MODE_BETA_MINUS);
		initDecayMode(31, 78, DECAY_MODE_BETA_MINUS);
		initDecayMode2(31, 79, DECAY_MODE_BETA_MINUS, 99.911, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(31, 80, DECAY_MODE_BETA_MINUS, 99.11, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(31, 81, DECAY_MODE_BETA_MINUS, 88.11, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(31, 82, DECAY_MODE_BETA_MINUS, 78.5, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(31, 83, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(31, 84, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 70.0, DECAY_MODE_BETA_MINUS);

		initAtomProperty(32, "Ge", "germanium", "ゲルマニウム", 60, 86, 60, 86, 23.222);
		initIsotopeProperty(32, 60, 59.97019, 30.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(32, 61, 60.96379, 39.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(32, 62, 61.95465, 129.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(32, 63, 62.94964, 142.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(32, 64, 63.94165, 63.7, HLU_SECOND, 0.0);	
		initIsotopeProperty(32, 65, 64.93944, 30.9, HLU_SECOND, 0.0);
		initIsotopeProperty(32, 66, 65.93384, 2.26, HLU_HOUR, 0.0);
		initIsotopeProperty(32, 67, 66.932734, 18.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(32, 68, 67.928094, 270.95, HLU_DAY, 0.0);
		initIsotopeProperty(32, 69, 68.9279645, 39.05, HLU_HOUR, 0.0);
		initIsotopeProperty(32, 70, 69.9242474, 0.0, HLU_STABLE, 0.2038);
		initIsotopeProperty(32, 71, 70.9249510, 11.43, HLU_DAY, 0.0);
		initIsotopeProperty(32, 72, 71.9220758, 0.0, HLU_STABLE, 0.2731);
		initIsotopeProperty(32, 73, 72.9234589, 0.0, HLU_STABLE, 0.0776);
		initIsotopeProperty(32, 74, 73.9211778, 0.0, HLU_STABLE, 0.3672);
		initIsotopeProperty(32, 75, 74.9228589, 82.78, HLU_MINUTE, 0.0);
		initIsotopeProperty(32, 76, 75.9214026, 1.78e21, HLU_YEAR, 0.0783);
		initIsotopeProperty(32, 77, 76.9235486, 11.30, HLU_HOUR, 0.0);
		initIsotopeProperty(32, 78, 77.922853, 88, HLU_MINUTE, 0.0);
		initIsotopeProperty(32, 79, 78.9254, 18.98, HLU_SECOND, 0.0);
		initIsotopeProperty(32, 80, 79.92537, 29.5, HLU_SECOND, 0.0);
		initIsotopeProperty(32, 81, 80.92882, 7.6, HLU_SECOND, 0.0);
		initIsotopeProperty(32, 82, 81.92955, 4.55, HLU_SECOND, 0.0);
		initIsotopeProperty(32, 83, 82.93462, 1.85, HLU_SECOND, 0.0);
		initIsotopeProperty(32, 84, 83.93747, 0.947, HLU_SECOND, 0.0);	
		initIsotopeProperty(32, 85, 84.94303, 535.0e-3, HLU_SECOND, 0.0);		
		initIsotopeProperty(32, 86, 85.94649, 150.0e-6, HLU_SECOND, 0.0);
		initDecayMode2(32, 60, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_2PROTON_EMISSION);
		initDecayMode2(32, 61, DECAY_MODE_BETA_PLUS_PROTON, 80.0, DECAY_MODE_BETA_PLUS);
		initDecayMode(32, 62, DECAY_MODE_BETA_PLUS);
		initDecayMode(32, 63, DECAY_MODE_BETA_PLUS);
		initDecayMode(32, 64, DECAY_MODE_BETA_PLUS);
		initDecayMode2(32, 65, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(32, 66, DECAY_MODE_BETA_PLUS);
		initDecayMode(32, 67, DECAY_MODE_BETA_PLUS);
		initDecayMode(32, 68, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(32, 69, DECAY_MODE_BETA_PLUS);
		initDecayMode(32, 71, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(32, 75, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 76, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(32, 77, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 78, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 79, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 80, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 81, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 82, DECAY_MODE_BETA_MINUS);
		initDecayMode(32, 83, DECAY_MODE_BETA_MINUS);
		initDecayMode2(32, 84, DECAY_MODE_BETA_MINUS, 89.2, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(32, 85, DECAY_MODE_BETA_MINUS, 86.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(32, 86, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 50.0, DECAY_MODE_BETA_MINUS);

		initAtomProperty(33, "As", "arsenic", "ヒ素", 64, 89, 64, 89, 24.64);
		initIsotopeProperty(33, 64, 63.95757, 40.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(33, 65, 64.94956, 170.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 66, 65.94471, 95.77e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 67, 66.93919, 42.5, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 68, 67.93677, 151.6, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 69, 68.93227, 15.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(33, 70, 69.93092, 52.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(33, 71, 70.927112, 65.28, HLU_HOUR, 0.0);
		initIsotopeProperty(33, 72, 71.926752, 26.0, HLU_HOUR, 0.0);
		initIsotopeProperty(33, 73, 72.923825, 80.30, HLU_DAY, 0.0);
		initIsotopeProperty(33, 74, 73.9239287, 17.77, HLU_DAY, 0.0);
		initIsotopeProperty(33, 75, 74.9215965, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(33, 76, 75.922394, 1.0942, HLU_DAY, 0.0);
		initIsotopeProperty(33, 77, 76.9206473, 38.83, HLU_HOUR, 0.0);
		initIsotopeProperty(33, 78, 77.921827, 90.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(33, 79, 78.920948, 9.01, HLU_MINUTE, 0.0);
		initIsotopeProperty(33, 80, 79.922534, 15.2, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 81, 80.922132, 33.3, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 82, 81.92450, 19.1, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 83, 82.92498, 13.4, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 84, 83.92906, 4.02, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 85, 84.93202, 2.021, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 86, 85.93650, 0.945, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 87, 86.93990, 0.56, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 88, 87.94494, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(33, 89, 88.94939, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode(33, 67, DECAY_MODE_BETA_PLUS);
		initDecayMode(33, 68, DECAY_MODE_BETA_PLUS);
		initDecayMode(33, 69, DECAY_MODE_BETA_PLUS);
		initDecayMode(33, 70, DECAY_MODE_BETA_PLUS);
		initDecayMode(33, 71, DECAY_MODE_BETA_PLUS);
		initDecayMode(33, 72, DECAY_MODE_BETA_PLUS);
		initDecayMode(33, 73, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(33, 74, DECAY_MODE_BETA_PLUS, 66.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(33, 76, DECAY_MODE_BETA_MINUS, 99.98, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(33, 77, DECAY_MODE_BETA_MINUS);
		initDecayMode(33, 78, DECAY_MODE_BETA_MINUS);
		initDecayMode(33, 79, DECAY_MODE_BETA_MINUS);
		initDecayMode(33, 80, DECAY_MODE_BETA_MINUS);
		initDecayMode(33, 81, DECAY_MODE_BETA_MINUS);
		initDecayMode(33, 82, DECAY_MODE_BETA_MINUS);
		initDecayMode(33, 83, DECAY_MODE_BETA_MINUS);
		initDecayMode2(33, 84, DECAY_MODE_BETA_MINUS, 99.721, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(33, 85, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 59.4, DECAY_MODE_BETA_MINUS);
		initDecayMode2(33, 86, DECAY_MODE_BETA_MINUS, 67.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(33, 87, DECAY_MODE_BETA_MINUS, 84.6, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(33, 88, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(33, 89, DECAY_MODE_BETA_MINUS);

		initAtomProperty(34, "Se", "selenium", "セレン", 65, 92, 65, 92, 25.363);
		initIsotopeProperty(34, 65, 64.96466, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 66, 65.95521, 33.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 67, 66.95009, 133.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 68, 67.94180, 35.5, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 69, 68.93956, 27.4, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 70, 69.93339, 41.1, HLU_MINUTE, 0.0);	
		initIsotopeProperty(34, 71, 70.93224, 4.74, HLU_MINUTE, 0.0);	
		initIsotopeProperty(34, 72, 71.927112, 8.40, HLU_DAY, 0.0);	
		initIsotopeProperty(34, 73, 72.926765, 7.15, HLU_HOUR, 0.0);	
		initIsotopeProperty(34, 74, 73.9224764, 0.0, HLU_STABLE, 0.0089);
		initIsotopeProperty(34, 75, 74.9225234, 119.779, HLU_DAY, 0.0);
		initIsotopeProperty(34, 76, 75.9192136, 0.0, HLU_STABLE, 0.0937);
		initIsotopeProperty(34, 77, 76.9199140, 0.0, HLU_STABLE, 0.0763);
		initIsotopeProperty(34, 78, 77.9173091, 0.0, HLU_STABLE, 0.2377);
		initIsotopeProperty(34, 79, 78.9184991, 3.27e5, HLU_YEAR, 0.0);
		initIsotopeProperty(34, 80, 79.9165213, 0.0, HLU_STABLE, 0.4961);
		initIsotopeProperty(34, 81, 80.9179925, 18.45, HLU_MINUTE, 0.0);
		initIsotopeProperty(34, 82, 81.9166994, 0.97e20, HLU_YEAR,	0.0873);
		initIsotopeProperty(34, 83, 82.919118, 22.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(34, 84, 83.918462, 3.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(34, 85, 84.92225, 31.7, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 86, 85.924272, 15.3, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 87, 86.92852, 5.50, HLU_SECOND, 0.0);	
		initIsotopeProperty(34, 88, 87.93142, 1.53, HLU_SECOND, 0.0);	
		initIsotopeProperty(34, 89, 88.93645, 0.41, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 90, 89.93996, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 91, 90.94596, 270.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(34, 92, 91.94992, 100.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(34, 65, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(34, 66, DECAY_MODE_BETA_PLUS);
		initDecayMode2(34, 67, DECAY_MODE_BETA_PLUS, 99.5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(34, 68, DECAY_MODE_BETA_PLUS);
		initDecayMode2(34, 69, DECAY_MODE_BETA_PLUS, 99.955, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(34, 70, DECAY_MODE_BETA_PLUS);
		initDecayMode(34, 71, DECAY_MODE_BETA_PLUS);
		initDecayMode(34, 72, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(34, 73, DECAY_MODE_BETA_PLUS);
		initDecayMode(34, 75, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(34, 79, DECAY_MODE_BETA_MINUS);
		initDecayMode(34, 81, DECAY_MODE_BETA_MINUS);
		initDecayMode(34, 82, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(34, 83, DECAY_MODE_BETA_MINUS);
		initDecayMode(34, 84, DECAY_MODE_BETA_MINUS);
		initDecayMode(34, 85, DECAY_MODE_BETA_MINUS);
		initDecayMode(34, 86, DECAY_MODE_BETA_MINUS);
		initDecayMode2(34, 87, DECAY_MODE_BETA_MINUS, 99.64, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(34, 88, DECAY_MODE_BETA_MINUS, 99.01, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(34, 89, DECAY_MODE_BETA_MINUS, 92.2, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(34, 90, DECAY_MODE_BETA_MINUS_AND_NEWTRON, 50.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(34, 91, DECAY_MODE_BETA_MINUS, 79.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(34, 92, DECAY_MODE_BETA_MINUS);

		initAtomProperty(35, "Br", "bromine", "臭素", 68, 94, 68, 94, 75.69);
		initIsotopeProperty(35, 68, 67.95852, 1.2e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 69, 68.95011, 24.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 70, 69.94479, 79.1e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 71, 70.93874, 21.4, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 72, 71.93664, 78.6, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 73, 72.93169, 3.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 74, 73.929891, 25.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 75, 74.925776, 96.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 76, 75.924541, 16.2, HLU_HOUR, 0.0);
		initIsotopeProperty(35, 77, 76.921379, 57.036, HLU_HOUR, 0.0);
		initIsotopeProperty(35, 78, 77.921146, 6.46, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 79, 78.9183371, 0.0, HLU_STABLE, 0.5069);
		initIsotopeProperty(35, 80, 79.9185293, 17.68, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 81, 80.9162906, 0.0, HLU_STABLE, 0.4931);
		initIsotopeProperty(35, 82, 81.9168041, 35.282, HLU_HOUR, 0.0);
		initIsotopeProperty(35, 83, 82.915180, 2.40, HLU_HOUR, 0.0);
		initIsotopeProperty(35, 84, 83.916479, 31.80, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 85, 84.915608, 2.90, HLU_MINUTE, 0.0);
		initIsotopeProperty(35, 86, 85.918798, 55.1, HLU_SECOND, 0.0);	
		initIsotopeProperty(35, 87, 86.920711, 55.65, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 88, 87.92407, 16.29, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 89, 88.92639, 4.40, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 90, 89.93063, 1.91, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 91, 90.93397, 541.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 92, 91.93926, 0.343, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 93, 92.94305, 102.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(35, 94, 93.94868, 70.0e-3, HLU_SECOND, 0.0);
		initDecayMode(35, 68, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(35, 69, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(35, 70, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 71, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 72, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 73, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 74, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 75, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 76, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 77, DECAY_MODE_BETA_PLUS);
		initDecayMode2(35, 78, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(35, 80, DECAY_MODE_BETA_MINUS, 91.7, DECAY_MODE_BETA_PLUS);
		initDecayMode(35, 82, DECAY_MODE_BETA_MINUS);
		initDecayMode(35, 83, DECAY_MODE_BETA_MINUS);
		initDecayMode(35, 84, DECAY_MODE_BETA_MINUS);
		initDecayMode(35, 85, DECAY_MODE_BETA_MINUS);
		initDecayMode(35, 86, DECAY_MODE_BETA_MINUS);
		initDecayMode2(35, 87, DECAY_MODE_BETA_MINUS, 97.48, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 88, DECAY_MODE_BETA_MINUS, 93.42, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 89, DECAY_MODE_BETA_MINUS, 86.2, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 90, DECAY_MODE_BETA_MINUS, 74.8, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 91, DECAY_MODE_BETA_MINUS, 80.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 92, DECAY_MODE_BETA_MINUS, 66.9, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 93, DECAY_MODE_BETA_MINUS, 89.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(35, 94, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(36, "Kr", "krypton", "クリプトン", 69, 97, 69, 97, 20.786);
		initIsotopeProperty(36, 69, 68.96518, 32.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 70, 69.95526, 52.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 71, 70.94963, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 72, 71.942092, 17.16, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 73, 72.939289, 28.6, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 74, 73.9330844, 11.50, HLU_MINUTE, 0.0);
		initIsotopeProperty(36, 75, 74.930946, 4.29, HLU_MINUTE, 0.0);
		initIsotopeProperty(36, 76, 75.925910, 14.8, HLU_HOUR, 0.0);
		initIsotopeProperty(36, 77, 76.9246700, 74.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(36, 78, 77.9203648, 0.0, HLU_STABLE, 0.00355);
		initIsotopeProperty(36, 79, 78.920082, 35.04, HLU_HOUR, 0.0);
		initIsotopeProperty(36, 80, 79.9163790, 0.0, HLU_STABLE, 0.02286);
		initIsotopeProperty(36, 81, 80.9165920, 2.29e5, HLU_YEAR, 0.0);
		initIsotopeProperty(36, 82, 81.9134836, 0.0, HLU_STABLE, 0.11593);
		initIsotopeProperty(36, 83, 82.914136, 0.0, HLU_STABLE, 0.11500);
		initIsotopeProperty(36, 84, 83.911507, 0.0, HLU_STABLE, 0.56987);
		initIsotopeProperty(36, 85, 84.9125273, 10.776, HLU_YEAR, 0.0);
		initIsotopeProperty(36, 86, 85.91061073, 0.0, HLU_STABLE, 0.17279);
		initIsotopeProperty(36, 87, 86.91335486, 76.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(36, 88, 87.914447, 2.84, HLU_HOUR, 0.0);
		initIsotopeProperty(36, 89, 88.91763, 3.15, HLU_MINUTE, 0.0);
		initIsotopeProperty(36, 90, 89.919517, 32.32, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 91, 90.92345, 8.57, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 92, 91.926156, 1.840, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 93, 92.93127, 1.286, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 94, 93.93436, 210.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 95, 94.93984, 114.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 96, 95.94307, 80.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(36, 97, 96.94856, 63.0e-3, HLU_SECOND, 0.0);
		initDecayMode(36, 69, DECAY_MODE_BETA_PLUS);
		initDecayMode(36, 70, DECAY_MODE_BETA_PLUS);
		initDecayMode2(36, 71, DECAY_MODE_BETA_PLUS, 94.8, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(36, 72, DECAY_MODE_BETA_PLUS);
		initDecayMode2(36, 73, DECAY_MODE_BETA_PLUS, 99.32, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(36, 74, DECAY_MODE_BETA_PLUS);
		initDecayMode(36, 75, DECAY_MODE_BETA_PLUS);
		initDecayMode(36, 76, DECAY_MODE_BETA_PLUS);
		initDecayMode(36, 77, DECAY_MODE_BETA_PLUS);
		initDecayMode(36, 79, DECAY_MODE_BETA_PLUS);
		initDecayMode(36, 81, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(36, 85, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 87, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 88, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 89, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 90, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 91, DECAY_MODE_BETA_MINUS);
		initDecayMode2(36, 92, DECAY_MODE_BETA_MINUS, 99.96, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(36, 93, DECAY_MODE_BETA_MINUS, 98.05, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode2(36, 94, DECAY_MODE_BETA_MINUS, 94.3, DECAY_MODE_BETA_MINUS_AND_NEWTRON);
		initDecayMode(36, 95, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 96, DECAY_MODE_BETA_MINUS);
		initDecayMode2(36, 97, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(76, "Os", "osmium", "オスミウム", 162, 196, 168, 196, 24.7);
		initIsotopeProperty(76, 162, 161.98443, 1.87e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 163, 162.98269, 5.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 164, 163.97804, 21.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 165, 164.97676, 71.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 166, 165.972691, 216.e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 167, 166.97155, 810.e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 168, 167.967804, 2.06, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 169, 168.967019, 3.40, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 170, 169.963577, 7.46, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 171, 170.963185, 8.3, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 172, 171.960023, 19.2, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 173, 172.959808, 22.4, HLU_SECOND, 0.0);
		initIsotopeProperty(76, 174, 173.957062, 44., HLU_SECOND, 0.0);
		initIsotopeProperty(76, 175, 174.956946, 1.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 176, 175.95481, 3.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 177, 176.954965, 3.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 178, 177.953251, 5.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 179, 178.953816, 6.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 180, 179.952379, 21.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 181, 180.95324, 105., HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 182, 181.952110, 22.10, HLU_HOUR, 0.0);
		initIsotopeProperty(76, 183, 182.95313, 13.0, HLU_HOUR, 0.0);
		initIsotopeProperty(76, 184, 183.9524891, 0.0, HLU_STABLE, 2.0e-4);
		initIsotopeProperty(76, 185, 184.9540423, 93.6, HLU_DAY, 0.0);
		initIsotopeProperty(76, 186, 185.9538382, 0.0, HLU_STABLE, 0.0159);// 2.0e15 year > the Universe, 0.0159
		initIsotopeProperty(76, 187, 186.9557505, 0.0, HLU_STABLE, 0.0196);
		initIsotopeProperty(76, 188, 187.9558382, 0.0,HLU_STABLE, 0.1324);
		initIsotopeProperty(76, 189, 188.9581475, 0.0, HLU_STABLE, 0.1615);
		initIsotopeProperty(76, 190, 189.9584470, 0.0, HLU_STABLE, 0.2626);
		initIsotopeProperty(76, 191, 190.9609297, 15.4, HLU_DAY, 0.0);
		initIsotopeProperty(76, 192, 191.9614807, 0.0, HLU_STABLE, 0.4078);
		initIsotopeProperty(76, 193, 192.9641516, 30.11, HLU_HOUR, 0.0);
		initIsotopeProperty(76, 194, 193.9651821, 6.0, HLU_YEAR, 0.0);
		initIsotopeProperty(76, 195, 194.96813, 6.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(76, 196, 195.96964, 34.9, HLU_MINUTE, 0.0);
		initDecayMode(76, 162, DECAY_MODE_ALPHA);
		initDecayMode(76, 163, DECAY_MODE_ALPHA);
		initDecayMode2(76, 164, DECAY_MODE_ALPHA, 98., DECAY_MODE_BETA_PLUS);
		initDecayMode2(76, 165, DECAY_MODE_ALPHA, 60., DECAY_MODE_BETA_PLUS);
		initDecayMode2(76, 166, DECAY_MODE_ALPHA, 72., DECAY_MODE_BETA_PLUS);
		initDecayMode2(76, 167, DECAY_MODE_ALPHA, 67., DECAY_MODE_BETA_PLUS);
		initDecayMode2(76, 168, DECAY_MODE_BETA_PLUS, 51., DECAY_MODE_ALPHA);
		initDecayMode2(76, 169, DECAY_MODE_BETA_PLUS, 89., DECAY_MODE_ALPHA);
		initDecayMode2(76, 170, DECAY_MODE_BETA_PLUS, 91.4, DECAY_MODE_ALPHA);
		initDecayMode2(76, 171, DECAY_MODE_BETA_PLUS, 98.3, DECAY_MODE_ALPHA);
		initDecayMode2(76, 172, DECAY_MODE_BETA_PLUS, 98.9, DECAY_MODE_ALPHA);
		initDecayMode2(76, 173, DECAY_MODE_BETA_PLUS, 99.6, DECAY_MODE_ALPHA);
		initDecayMode2(76, 174, DECAY_MODE_BETA_PLUS, 99.97, DECAY_MODE_ALPHA);
		initDecayMode(76, 175, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 176, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 177, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 178, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 179, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 180, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 181, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 182, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(76, 183, DECAY_MODE_BETA_PLUS);
		initDecayMode(76, 185, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(76, 191, DECAY_MODE_BETA_MINUS);
		initDecayMode(76, 193, DECAY_MODE_BETA_MINUS);
		initDecayMode(76, 194, DECAY_MODE_BETA_MINUS);
		initDecayMode(76, 195, DECAY_MODE_BETA_MINUS);
		initDecayMode(76, 196, DECAY_MODE_BETA_MINUS);

		initAtomProperty(77, "Ir", "iridium", "イリジウム", 165, 199, 171, 199, 25.10);
		initIsotopeProperty(77, 165, 164.98752, 1.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 166, 165.98582, 10.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 167, 166.981665, 35.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 168, 167.97988, 161.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 169, 168.976295, 780.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 170, 169.97497, 910.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 171, 170.97163, 3.6, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 172, 171.97046, 4.4, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 173, 172.967502, 9.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 174, 173.966861, 7.9, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 175, 174.964113, 9.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 176, 175.963649, 8.3, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 177, 176.961302, 30.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 178, 177.961082, 12.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 179, 178.959122, 79.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 180, 179.959229, 1.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(77, 181, 180.957625, 4.90, HLU_MINUTE, 0.0);
		initIsotopeProperty(77, 182, 181.958076, 15.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(77, 183, 182.956846, 57.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(77, 184, 183.95748, 3.09, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 185, 184.95670, 14.4, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 186, 185.957946, 16.64, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 187, 186.957363, 10.5, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 188, 187.958853, 41.5, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 189, 188.958719, 13.2, HLU_DAY, 0.0);
		initIsotopeProperty(77, 190, 189.9605460, 11.78, HLU_DAY, 0.0);
		initIsotopeProperty(77, 191, 190.9605940, 0.0, HLU_STABLE, 0.373);
		initIsotopeProperty(77, 192, 191.9626050, 73.827, HLU_DAY, 0.0);
		initIsotopeProperty(77, 193, 192.9629264, 0.0, HLU_STABLE, 0.627);
		initIsotopeProperty(77, 194, 193.9650784, 19.28, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 195, 194.9659796, 2.5, HLU_HOUR, 0.0);
		initIsotopeProperty(77, 196, 195.96840, 52.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 197, 196.969653, 5.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(77, 198, 197.97228, 8.0, HLU_SECOND, 0.0);
		initIsotopeProperty(77, 199, 198.97380, 20.0, HLU_SECOND, 0.0);
		initDecayMode(77, 165, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(77, 166, DECAY_MODE_ALPHA, 93., DECAY_MODE_PROTON_EMISSION);
		initDecayMode3(77, 167, DECAY_MODE_ALPHA, 48., DECAY_MODE_PROTON_EMISSION, 32., DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 168, DECAY_MODE_ALPHA);
		initDecayMode(77, 169, DECAY_MODE_ALPHA);
		initDecayMode2(77, 170, DECAY_MODE_BETA_PLUS, 64., DECAY_MODE_ALPHA);
		initDecayMode2(77, 171, DECAY_MODE_ALPHA, 58., DECAY_MODE_BETA_PLUS);
		initDecayMode2(77, 172, DECAY_MODE_BETA_PLUS, 98., DECAY_MODE_ALPHA);
		initDecayMode2(77, 173, DECAY_MODE_BETA_PLUS, 93., DECAY_MODE_ALPHA);
		initDecayMode2(77, 174, DECAY_MODE_BETA_PLUS, 99.5, DECAY_MODE_ALPHA);
		initDecayMode2(77, 175, DECAY_MODE_BETA_PLUS, 99.15, DECAY_MODE_ALPHA);
		initDecayMode2(77, 176, DECAY_MODE_BETA_PLUS, 97.9, DECAY_MODE_ALPHA);
		initDecayMode2(77, 177, DECAY_MODE_BETA_PLUS, 99.94, DECAY_MODE_ALPHA);
		initDecayMode(77, 178, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 179, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 180, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 181, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 182, DECAY_MODE_BETA_PLUS);
		initDecayMode2(77, 183, DECAY_MODE_BETA_PLUS, 99.95, DECAY_MODE_ALPHA);
		initDecayMode(77, 184, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 185, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 186, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 187, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 188, DECAY_MODE_BETA_PLUS);
		initDecayMode(77, 189, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(77, 190, DECAY_MODE_BETA_PLUS);
		initDecayMode2(77, 192, DECAY_MODE_BETA_MINUS, 95.24, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(77, 194, DECAY_MODE_BETA_MINUS);
		initDecayMode(77, 195, DECAY_MODE_BETA_MINUS);
		initDecayMode(77, 196, DECAY_MODE_BETA_MINUS);
		initDecayMode(77, 197, DECAY_MODE_BETA_MINUS);
		initDecayMode(77, 198, DECAY_MODE_BETA_MINUS);
		initDecayMode(77, 199, DECAY_MODE_BETA_MINUS);

		initAtomProperty(78, "Pt", "platinum", "白金", 168, 202, 175, 202, 25.86);
		initIsotopeProperty(78, 168, 167.98815, 2.00e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 169, 168.98672, 3.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 170, 169.982495, 14.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 171, 170.98124, 51.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 172, 171.977347, 98.4e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 173, 172.97644, 365.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 174, 173.972819, 0.889, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 175, 174.972421, 2.53, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 176, 175.968945, 6.33, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 177, 176.968469, 10.6, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 178, 177.965649, 21.1, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 179, 178.965363, 21.2, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 180, 179.963031, 56., HLU_SECOND, 0.0);
		initIsotopeProperty(78, 181, 180.963097, 52.0, HLU_SECOND, 0.0);
		initIsotopeProperty(78, 182, 181.961171, 2.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(78, 183, 182.961597, 6.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(78, 184, 183.959922, 17.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(78, 185, 184.96062, 70.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(78, 186, 185.959351, 2.08, HLU_HOUR, 0.0);
		initIsotopeProperty(78, 187, 186.96059, 2.35, HLU_HOUR, 0.0);
		initIsotopeProperty(78, 188, 187.959395, 10.2, HLU_DAY, 0.0);
		initIsotopeProperty(78, 189, 188.960834, 10.87, HLU_HOUR, 0.0);
		initIsotopeProperty(78, 190, 189.959932, 0.0, HLU_STABLE, 1.4e-4); // 6.5e11 YEAR is longer than the universe!
		initIsotopeProperty(78, 191, 190.961677, 2.862, HLU_DAY, 0.0);
		initIsotopeProperty(78, 192, 191.9610380, 0.0, HLU_STABLE, 0.00782);
		initIsotopeProperty(78, 193, 192.9629874, 50., HLU_YEAR, 0.0);
		initIsotopeProperty(78, 194, 193.9626803, 0.0, HLU_STABLE,	0.32967);
		initIsotopeProperty(78, 195, 194.9647911, 0.0, HLU_STABLE,	0.33832);
		initIsotopeProperty(78, 196, 195.9649515, 0.0, HLU_STABLE,	0.25242);
		initIsotopeProperty(78, 197, 196.9673402, 19.8915, HLU_HOUR, 0.0);
		initIsotopeProperty(78, 198, 197.967893, 0.0, HLU_STABLE,	0.07163);
		initIsotopeProperty(78, 199, 198.970593, 30.80, HLU_MINUTE, 0.0);
		initIsotopeProperty(78, 200, 199.971441, 12.5, HLU_HOUR, 0.0);
		initIsotopeProperty(78, 201, 200.97451, 2.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(78, 202, 201.97574, 44.0, HLU_HOUR, 0.0);
		initDecayMode(78, 168, DECAY_MODE_ALPHA);
		initDecayMode(78, 169, DECAY_MODE_ALPHA);
		initDecayMode2(78, 170, DECAY_MODE_ALPHA, 98., DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 171, DECAY_MODE_ALPHA, 99., DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 172, DECAY_MODE_ALPHA, 77., DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 173, DECAY_MODE_ALPHA, 84., DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 174, DECAY_MODE_ALPHA, 83., DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 175, DECAY_MODE_ALPHA, 64., DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 176, DECAY_MODE_BETA_PLUS, 62., DECAY_MODE_ALPHA);
		initDecayMode2(78, 177, DECAY_MODE_BETA_PLUS, 94.4, DECAY_MODE_ALPHA);
		initDecayMode2(78, 178, DECAY_MODE_BETA_PLUS, 92.3, DECAY_MODE_ALPHA);
		initDecayMode2(78, 179, DECAY_MODE_BETA_PLUS, 99.76, DECAY_MODE_ALPHA);
		initDecayMode2(78, 180, DECAY_MODE_BETA_PLUS, 99.7, DECAY_MODE_ALPHA);
		initDecayMode2(78, 181, DECAY_MODE_BETA_PLUS, 100. - 0.074, DECAY_MODE_ALPHA);
		initDecayMode2(78, 182, DECAY_MODE_BETA_PLUS, 100. - 0.038, DECAY_MODE_ALPHA);
		initDecayMode2(78, 183, DECAY_MODE_BETA_PLUS, 100. - 0.0096, DECAY_MODE_ALPHA);
		initDecayMode2(78, 184, DECAY_MODE_BETA_PLUS, 100. - 0.00169, DECAY_MODE_ALPHA);
		initDecayMode2(78, 185, DECAY_MODE_BETA_PLUS, 100. - 0.005, DECAY_MODE_ALPHA);
		initDecayMode2(78, 186, DECAY_MODE_BETA_PLUS, 100. - 1.4e-4, DECAY_MODE_ALPHA);
		initDecayMode(78, 187, DECAY_MODE_BETA_PLUS);
		initDecayMode2(78, 188, DECAY_MODE_ELECTRON_CAPTURE, 100. - 2.6e-5, DECAY_MODE_ALPHA);
		initDecayMode(78, 189, DECAY_MODE_BETA_PLUS);
		initDecayMode(78, 191, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(78, 193, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(78, 197, DECAY_MODE_BETA_MINUS);
		initDecayMode(78, 199, DECAY_MODE_BETA_MINUS);
		initDecayMode(78, 200, DECAY_MODE_BETA_MINUS);
		initDecayMode(78, 201, DECAY_MODE_BETA_MINUS);
		initDecayMode(78, 202, DECAY_MODE_BETA_MINUS);
		
		initAtomProperty(79, "Au", "gold", "金", 171, 205, 176, 205, 25.418);
		initIsotopeProperty(79, 171, 170.991879, 30.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 172, 171.99004, 4.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 173, 172.986237, 25e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 174, 173.98476, 139.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 175, 174.98127, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 176, 175.98010, 1.08, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 177, 176.976865, 1.462, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 178, 177.97603, 2.6, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 179, 178.973213, 7.1, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 180, 179.972521, 8.1, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 181, 180.970079, 13.7, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 182, 181.969618, 15.5, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 183, 182.967593, 42.8, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 184, 183.967452, 20.6, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 185, 184.965789, 4.25, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 186, 185.965953, 10.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 187, 186.964568, 8.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 188, 187.965324, 8.84, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 189, 188.963948, 28.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 190, 189.964700, 42.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 191, 190.96370, 3.18, HLU_HOUR, 0.0);
		initIsotopeProperty(79, 192, 191.964813, 4.94, HLU_HOUR, 0.0);
		initIsotopeProperty(79, 193, 192.964150, 17.65, HLU_HOUR, 0.0);
		initIsotopeProperty(79, 194, 193.965365, 38.02, HLU_HOUR, 0.0);
		initIsotopeProperty(79, 195, 194.9650346, 186.098, HLU_DAY, 0.0);
		initIsotopeProperty(79, 196, 195.966570, 6.1669, HLU_DAY, 0.0);
		initIsotopeProperty(79, 197, 196.9665687, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(79, 198, 197.9682423, 2.69517, HLU_DAY, 0.0);
		initIsotopeProperty(79, 199, 198.9687652, 3.139, HLU_DAY, 0.0);
		initIsotopeProperty(79, 200, 199.97073, 48.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 201, 200.971657, 26., HLU_MINUTE, 0.0);
		initIsotopeProperty(79, 202, 201.97381, 28.8, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 203, 202.975155, 53., HLU_SECOND, 0.0);
		initIsotopeProperty(79, 204, 203.97772, 39.8, HLU_SECOND, 0.0);
		initIsotopeProperty(79, 205, 204.97987, 31.0, HLU_SECOND, 0.0);
		initDecayMode(79, 171, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(79, 172, DECAY_MODE_ALPHA, 98., DECAY_MODE_PROTON_EMISSION);
		initDecayMode(79, 173, DECAY_MODE_ALPHA);
		initDecayMode(79, 174, DECAY_MODE_ALPHA);
		initDecayMode2(79, 175, DECAY_MODE_ALPHA, 82., DECAY_MODE_BETA_PLUS);
		initDecayMode2(79, 176, DECAY_MODE_ALPHA, 60., DECAY_MODE_BETA_PLUS);
		initDecayMode2(79, 177, DECAY_MODE_BETA_PLUS, 60., DECAY_MODE_ALPHA);
		initDecayMode2(79, 178, DECAY_MODE_BETA_PLUS, 60., DECAY_MODE_ALPHA);
		initDecayMode2(79, 179, DECAY_MODE_BETA_PLUS, 78., DECAY_MODE_ALPHA);
		initDecayMode2(79, 180, DECAY_MODE_BETA_PLUS, 98.2, DECAY_MODE_ALPHA);
		initDecayMode2(79, 181, DECAY_MODE_BETA_PLUS, 97.3, DECAY_MODE_ALPHA);
		initDecayMode2(79, 182, DECAY_MODE_BETA_PLUS, 99.87, DECAY_MODE_ALPHA);
		initDecayMode2(79, 183, DECAY_MODE_BETA_PLUS, 99.2, DECAY_MODE_ALPHA);
		initDecayMode(79, 184, DECAY_MODE_BETA_PLUS);
		initDecayMode2(79, 185, DECAY_MODE_BETA_PLUS, 99.74, DECAY_MODE_ALPHA);
		initDecayMode2(79, 186, DECAY_MODE_BETA_PLUS, 99.9992, DECAY_MODE_ALPHA);
		initDecayMode2(79, 187, DECAY_MODE_BETA_PLUS, 99.997, DECAY_MODE_ALPHA);
		initDecayMode(79, 188, DECAY_MODE_BETA_PLUS);
		initDecayMode2(79, 189, DECAY_MODE_BETA_PLUS, 99.9997, DECAY_MODE_ALPHA);
		initDecayMode2(79, 190, DECAY_MODE_BETA_PLUS, 100. - 10.0e-6, DECAY_MODE_ALPHA);
		initDecayMode(79, 191, DECAY_MODE_BETA_PLUS);
		initDecayMode(79, 192, DECAY_MODE_BETA_PLUS);
		initDecayMode2(79, 193, DECAY_MODE_BETA_PLUS, 100. - 10.0e-5, DECAY_MODE_ALPHA);
		initDecayMode(79, 194, DECAY_MODE_BETA_PLUS);
		initDecayMode(79, 195, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(79, 196, DECAY_MODE_BETA_PLUS, 93.05, DECAY_MODE_BETA_MINUS);
		initDecayMode(79, 198, DECAY_MODE_BETA_MINUS); 
		initDecayMode(79, 199, DECAY_MODE_BETA_MINUS); 
		initDecayMode(79, 200, DECAY_MODE_BETA_MINUS); 
		initDecayMode(79, 201, DECAY_MODE_BETA_MINUS);
		initDecayMode(79, 202, DECAY_MODE_BETA_MINUS); 
		initDecayMode(79, 203, DECAY_MODE_BETA_MINUS); 
		initDecayMode(79, 204, DECAY_MODE_BETA_MINUS);
		initDecayMode(79, 205, DECAY_MODE_BETA_MINUS);

		initAtomProperty(80, "Hg", "mercury", "水銀", 175, 208, 179, 208, 27.983);
		initIsotopeProperty(80, 175, 174.99142, 0.0108, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 176, 175.987355, 0.0204, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 177, 176.98628, 0.1273, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 178, 177.982483, 0.269, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 179, 178.981834, 1.09, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 180, 179.978266, 2.58, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 181, 180.977819, 3.6, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 182, 181.97469, 10.83, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 183, 182.974450, 9.4, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 184, 183.971713, 30.6, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 185, 184.971899, 49.1, HLU_SECOND, 0.0);
		initIsotopeProperty(80, 186, 185.969362, 1.38, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 187, 186.969814, 1.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 188, 187.967577, 3.25, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 189, 188.96819, 7.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 190, 189.966322, 20.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 191, 190.967157, 49.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 192, 191.965634, 4.85, HLU_HOUR, 0.0);
		initIsotopeProperty(80, 193, 192.966665, 3.80, HLU_HOUR, 0.0);
		initIsotopeProperty(80, 194, 193.965439, 444.0, HLU_YEAR, 0.0);
		initIsotopeProperty(80, 195, 194.966720, 10.53, HLU_HOUR, 0.0);
		initIsotopeProperty(80, 196, 195.965833, 0.0, HLU_STABLE, 0.0015);
		initIsotopeProperty(80, 197, 196.967213, 64.14, HLU_HOUR, 0.0);
		initIsotopeProperty(80, 198, 197.9667690, 0.0, HLU_STABLE, 0.0997);
		initIsotopeProperty(80, 199, 198.9682799, 0.0, HLU_STABLE, 0.1687);
		initIsotopeProperty(80, 200, 199.9683260, 0.0, HLU_STABLE, 0.2310);
		initIsotopeProperty(80, 201, 200.9703023, 0.0, HLU_STABLE, 0.1318);
		initIsotopeProperty(80, 202, 201.9706430, 0.0, HLU_STABLE, 0.2986);
		initIsotopeProperty(80, 203, 202.9728725, 46.595, HLU_DAY, 0.0);
		initIsotopeProperty(80, 204, 203.9734939, 0.0, HLU_STABLE, 0.0687);
		initIsotopeProperty(80, 205, 204.976073, 5.14, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 206, 205.977514, 8.15, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 207, 206.98259, 2.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(80, 208, 207.98594, 42.0, HLU_MINUTE, 0.0);
		initDecayMode(80, 175, DECAY_MODE_ALPHA);
		initDecayMode2(80, 176, DECAY_MODE_ALPHA, 98.6, DECAY_MODE_BETA_PLUS);
		initDecayMode2(80, 177, DECAY_MODE_ALPHA, 85.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(80, 178, DECAY_MODE_ALPHA, 70.0, DECAY_MODE_BETA_PLUS);
		initDecayMode3(80, 179, DECAY_MODE_ALPHA, 53.0 - 0.15, DECAY_MODE_BETA_PLUS, 47.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(80, 180, DECAY_MODE_BETA_PLUS, 52.0, DECAY_MODE_ALPHA);
		initDecayMode4(80, 181, DECAY_MODE_BETA_PLUS, 64.0 - 0.14 - 9.0e-8, DECAY_MODE_ALPHA, 36.0, DECAY_MODE_BETA_PLUS_PROTON, 0.14, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode3(80, 182, DECAY_MODE_BETA_PLUS, 84.8 - 1e-7, DECAY_MODE_ALPHA, 15.2, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(80, 183, DECAY_MODE_BETA_PLUS, 74.5 - 5.6e-6, DECAY_MODE_ALPHA, 25.5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(80, 184, DECAY_MODE_BETA_PLUS, 98.89, DECAY_MODE_ALPHA);
		initDecayMode2(80, 185, DECAY_MODE_BETA_PLUS, 94.0, DECAY_MODE_ALPHA);
		initDecayMode2(80, 186, DECAY_MODE_BETA_PLUS, 99.92, DECAY_MODE_ALPHA);
		initDecayMode2(80, 187, DECAY_MODE_BETA_PLUS, 1.0 - 1.2e-6, DECAY_MODE_ALPHA);
		initDecayMode2(80, 188, DECAY_MODE_BETA_PLUS, 1.0 - 3.7e-7, DECAY_MODE_ALPHA);
		initDecayMode2(80, 189, DECAY_MODE_BETA_PLUS, 1.0 - 3.0e-7, DECAY_MODE_ALPHA);
		initDecayMode2(80, 190, DECAY_MODE_BETA_PLUS, 1.0 - 5.0e-7, DECAY_MODE_ALPHA);
		initDecayMode(80, 191, DECAY_MODE_BETA_PLUS);
		initDecayMode2(80, 192, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 4.0e-8, DECAY_MODE_ALPHA);
		initDecayMode(80, 193, DECAY_MODE_BETA_PLUS);
		initDecayMode(80, 194, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(80, 195, DECAY_MODE_BETA_PLUS);
		initDecayMode(80, 197, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(80, 203, DECAY_MODE_BETA_MINUS);
		initDecayMode(80, 205, DECAY_MODE_BETA_MINUS);
		initDecayMode(80, 206, DECAY_MODE_BETA_MINUS);
		initDecayMode(80, 207, DECAY_MODE_BETA_MINUS);
		initDecayMode(80, 208, DECAY_MODE_BETA_MINUS);
		//initDecayMode(80, 209, DECAY_MODE_BETA_MINUS);
		//initDecayMode(80, 210, DECAY_MODE_BETA_MINUS);

		initAtomProperty(81, "Tl", "thallium", "タリウム", 177, 210, 180, 210, 26.32);
		initIsotopeProperty(81, 177, 176.996427, 18.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 178, 177.99490, 255.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 179, 178.99109, 270.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 180, 179.98991, 1.5, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 181, 180.986257, 3.2, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 182, 181.98567, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 183, 182.982193, 6.9, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 184, 183.98187, 9.7, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 185, 184.97879, 19.5, HLU_SECOND, 0.0);
		initIsotopeProperty(81, 186, 185.97833, 40., HLU_SECOND, 0.0);
		initIsotopeProperty(81, 187, 186.975906, 51., HLU_SECOND, 0.0);
		initIsotopeProperty(81, 188, 187.97601, 71., HLU_SECOND, 0.0);
		initIsotopeProperty(81, 189, 188.973588, 2.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 190, 189.97388, 2.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 191, 190.971786, 20., HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 192, 191.97223, 9.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 193, 192.97067, 21.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 194, 193.97120, 33.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 195, 194.969774, 1.16, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 196, 195.970481, 1.84, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 197, 196.969575, 2.84, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 198, 197.97048, 5.3, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 199, 198.96988, 7.42, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 200, 199.970963, 26.1, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 201, 200.970819, 72.912, HLU_HOUR, 0.0);
		initIsotopeProperty(81, 202, 201.972106, 12.23, HLU_DAY, 0.0);
		initIsotopeProperty(81, 203, 202.9723442, 0.0, HLU_STABLE, 0.2952);
		initIsotopeProperty(81, 204, 203.9738635, 3.78, HLU_YEAR, 0.0);
		initIsotopeProperty(81, 205, 204.9744275, 0.0, HLU_STABLE, 0.7048);
		initIsotopeProperty(81, 206, 205.9761103, 4.200, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 207, 206.977419, 4.77, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 208, 207.9820187, 3.053, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 209, 208.985359, 2.161, HLU_MINUTE, 0.0);
		initIsotopeProperty(81, 210, 209.990074, 1.30, HLU_MINUTE, 0.0);
		initDecayMode(81, 177, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(81, 178, DECAY_MODE_ALPHA);
		initDecayMode(81, 179, DECAY_MODE_ALPHA);
		initDecayMode2(81, 180, DECAY_MODE_ALPHA, 75., DECAY_MODE_BETA_PLUS);
		//initDecayMode3(81, 180, DECAY_MODE_ALPHA, 75. - 10.0e-4, DECAY_MODE_BETA_PLUS, 25., DECAY_MODE_SELF_FISSION_80KR);
		initDecayMode2(81, 181, DECAY_MODE_ALPHA, (75. + 96.) * .5, DECAY_MODE_BETA_PLUS);
		initDecayMode2(81, 182, DECAY_MODE_BETA_PLUS, 96., DECAY_MODE_ALPHA);
		initDecayMode2(81, 183, DECAY_MODE_BETA_PLUS, 98., DECAY_MODE_ALPHA);
		initDecayMode(81, 184, DECAY_MODE_BETA_PLUS);
		initDecayMode2(81, 185, DECAY_MODE_BETA_PLUS, 99., DECAY_MODE_ALPHA);
		initDecayMode2(81, 186, DECAY_MODE_BETA_PLUS, 1. - .006, DECAY_MODE_ALPHA); 
		initDecayMode2(81, 187, DECAY_MODE_BETA_PLUS, 1. - .00006, DECAY_MODE_ALPHA);
		initDecayMode(81, 188, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 189, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 190, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 191, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 192, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 193, DECAY_MODE_BETA_PLUS);
		initDecayMode2(81, 194, DECAY_MODE_BETA_PLUS, 1. - 10.0e-7, DECAY_MODE_ALPHA);
		initDecayMode(81, 195, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 196, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 197, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 198, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 199, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 200, DECAY_MODE_BETA_PLUS);
		initDecayMode(81, 201, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(81, 202, DECAY_MODE_BETA_PLUS);
		initDecayMode2(81, 204, DECAY_MODE_BETA_MINUS, 97.1, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(81, 206, DECAY_MODE_BETA_MINUS);
		initDecayMode(81, 207, DECAY_MODE_BETA_MINUS);
		initDecayMode(81, 208, DECAY_MODE_BETA_MINUS);
		initDecayMode(81, 209, DECAY_MODE_BETA_MINUS);
		initDecayMode2(81,210, DECAY_MODE_BETA_MINUS, 99.991, DECAY_MODE_BETA_MINUS_AND_NEWTRON);

		initAtomProperty(82, "Pb", "lead", "鉛", 181, 214, 185, 214, 26.650);
		initIsotopeProperty(82, 181, 180.99662, 45.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 182, 181.992672, 60.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 183, 182.99187, 535.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 184, 183.988142, 490.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 185, 184.987610, 6.3, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 186, 185.984239, 4.82, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 187, 186.983918, 15.2, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 188, 187.980874, 25.5, HLU_SECOND, 0.0);
		initIsotopeProperty(82, 189, 188.98081, 51., HLU_SECOND, 0.0);
		initIsotopeProperty(82, 190, 189.978082, 71., HLU_SECOND, 0.0);
		initIsotopeProperty(82, 191, 190.97827, 1.33, HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 192, 191.975785, 3.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 193, 192.97617, 5., HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 194, 193.974012, 12.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 195, 194.974542, 15., HLU_MINUTE, 0.0);	
		initIsotopeProperty(82, 196, 195.972774, 37., HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 197, 196.973431, 8.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 198, 197.972034, 2.4, HLU_HOUR, 0.0);
		initIsotopeProperty(82, 199, 198.972917, 90., HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 200, 199.971827, 21.5, HLU_HOUR, 0.0);
		initIsotopeProperty(82, 201, 200.972885, 9.33, HLU_HOUR, 0.0);
		initIsotopeProperty(82, 202, 201.972159, 52.5e3, HLU_YEAR, 0.0);
		initIsotopeProperty(82, 203, 202.973391, 51.873, HLU_HOUR, 0.0);
		initIsotopeProperty(82, 204, 203.9730436, 0.0, HLU_STABLE, 0.014);
		initIsotopeProperty(82, 205, 204.9744818, 15.3e6, HLU_YEAR, 0.0);
		initIsotopeProperty(82, 206, 205.9744653, 0.0, HLU_STABLE, 0.241);
		initIsotopeProperty(82, 207, 206.9758969, 0.0, HLU_STABLE, 0.221);
		initIsotopeProperty(82, 208, 207.9766521, 0.0, HLU_STABLE, 0.524);
		initIsotopeProperty(82, 209, 208.9810901, 3.253, HLU_HOUR, 0.0);
		initIsotopeProperty(82, 210, 209.9841885, 22.20, HLU_YEAR, 0.0);
		initIsotopeProperty(82, 211, 210.9887370, 36.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 212, 211.9918975, 10.64, HLU_HOUR, 0.0);
		initIsotopeProperty(82, 213, 212.996581, 10.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(82, 214, 213.9998054, 26.8, HLU_MINUTE, 0.0);
		initDecayMode2(82, 181, DECAY_MODE_ALPHA, 98., DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 182, DECAY_MODE_ALPHA, 98., DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 183, DECAY_MODE_ALPHA, 94., DECAY_MODE_BETA_PLUS);
		initDecayMode(82, 184, DECAY_MODE_ALPHA);
		initDecayMode(82, 185, DECAY_MODE_ALPHA);
		initDecayMode2(82, 186, DECAY_MODE_ALPHA, 56., DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 187, DECAY_MODE_BETA_PLUS, (100. - 56. + 91.5) * .05, DECAY_MODE_ALPHA);
		initDecayMode2(82, 188, DECAY_MODE_BETA_PLUS, 91.5, DECAY_MODE_ALPHA);
		initDecayMode(82, 189, DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 190, DECAY_MODE_BETA_PLUS, 99.1, DECAY_MODE_ALPHA);
		initDecayMode2(82, 191, DECAY_MODE_BETA_PLUS, 99.987, DECAY_MODE_ALPHA);
		initDecayMode2(82, 192, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_ALPHA);
		initDecayMode(82, 193, DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 194, DECAY_MODE_BETA_PLUS, 100. - 7.3e-6, DECAY_MODE_ALPHA);
		initDecayMode(82, 195, DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 196, DECAY_MODE_BETA_PLUS, 100. - 3.0e-5, DECAY_MODE_ALPHA);
		initDecayMode(82, 197, DECAY_MODE_BETA_PLUS);
		initDecayMode(82, 198, DECAY_MODE_BETA_PLUS);
		initDecayMode(82, 199, DECAY_MODE_BETA_PLUS);
		initDecayMode(82, 200, DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 201, DECAY_MODE_ELECTRON_CAPTURE, 99., DECAY_MODE_BETA_PLUS);
		initDecayMode2(82, 202, DECAY_MODE_ELECTRON_CAPTURE, 99., DECAY_MODE_ALPHA);
		initDecayMode(82, 203, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(82, 205, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(82, 209, DECAY_MODE_BETA_MINUS);
		initDecayMode2(82, 210, DECAY_MODE_BETA_MINUS, 100. - 1.9e-6, DECAY_MODE_ALPHA);
		initDecayMode(82, 211, DECAY_MODE_BETA_MINUS);
		initDecayMode(82, 212, DECAY_MODE_BETA_MINUS);
		initDecayMode(82, 213, DECAY_MODE_BETA_MINUS);
		initDecayMode(82, 214, DECAY_MODE_BETA_MINUS);

		e_hydrogenPtr = getIsotopePropertyPtr(ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN);
		if(!e_hydrogenPtr){
			fprintf(stderr, "FATAL ERROR:%s:NO hydrogen\n", __FUNCTION__);
		}
		e_helium4Ptr = getIsotopePropertyPtr(ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM);
		if(!e_helium4Ptr){
			fprintf(stderr, "FATAL ERROR:%s:NO helium4\n", __FUNCTION__);
		}
	}
}

extern void printIsotopeProperty(FILE * a_fp, struct isotopeProperty * a_isotopePropertyPtr)
{
	if(a_isotopePropertyPtr->massU > 0.0){
		//fprintf(a_fp, "AN:%d, ", a_isotopePropertyPtr->atomicNumber);
		fputs(" ", a_fp);
		fprintf(a_fp, "MN %d ", a_isotopePropertyPtr->massNumber);
		fprintf(a_fp, "massU %lg ", a_isotopePropertyPtr->massU);
		if(a_isotopePropertyPtr->halfLifeUnit == HLU_STABLE){
			fprintf(a_fp, "HL %s [] ", getHalfLifeUnitString(a_isotopePropertyPtr->halfLifeUnit));
			fprintf(a_fp, "RIA %lg ", a_isotopePropertyPtr->relativeIsotopicAbundance);
		}else{
			fprintf(a_fp, "HL %lg [%s] ", a_isotopePropertyPtr->halfLife, 
			getHalfLifeUnitString(a_isotopePropertyPtr->halfLifeUnit));
			if(a_isotopePropertyPtr->relativeIsotopicAbundance > 0.0){
				fprintf(a_fp, "RIA %lg ", a_isotopePropertyPtr->relativeIsotopicAbundance);
			}
		}
		if(strlen(a_isotopePropertyPtr->name) > 0
		|| strlen(a_isotopePropertyPtr->japaneseName) > 0
		|| strlen(a_isotopePropertyPtr->nucleusName) > 0){
			fputs("\n  ", a_fp);
			if(strlen(a_isotopePropertyPtr->symbol) > 0){
				fprintf(a_fp, "symbol %s ", a_isotopePropertyPtr->symbol);
			}
			if(strlen(a_isotopePropertyPtr->name) > 0){
				fprintf(a_fp, "name %s ", a_isotopePropertyPtr->name);
			}
			if(strlen(a_isotopePropertyPtr->japaneseName) > 0){
				fprintf(a_fp, "japaneseName \"%s\" ", a_isotopePropertyPtr->japaneseName);
			}
			if(strlen(a_isotopePropertyPtr->nucleusName) > 0){
				fprintf(a_fp, "nucleusName %s", a_isotopePropertyPtr->nucleusName);
			}
			fputs("\n", a_fp);
		}else{
			fprintf(a_fp, "symbol %s\n", a_isotopePropertyPtr->symbol);
		}
		if(a_isotopePropertyPtr->decayModeSize > 0){
			int i;
			fprintf(a_fp, " decayMode ");
			for(i = 0; i < a_isotopePropertyPtr->decayModeSize; ++i){
				
				if(i > 0){
					fprintf(a_fp, " ");
				}

				fprintf(a_fp, "( %s rate %lg )", getDecayModeText(a_isotopePropertyPtr->decayMode[i]), a_isotopePropertyPtr->decayModeRate[i]);
			}
			fputs("\n", a_fp);
		}
		fprintf(a_fp, " nucleusRadius %lg RCCS %lg \n", a_isotopePropertyPtr->nucleusRadius, a_isotopePropertyPtr->relativeCollideCrossSection);
	}
}
extern void printAtomProperty(FILE * a_fp, struct atomProperty * a_atomPropertyPtr)
{
	int massNumber;
	fprintf(a_fp, "AN %d ", a_atomPropertyPtr->atomicNumber);
	fprintf(a_fp, "symbol \"%s\" ", a_atomPropertyPtr->symbol);
	fprintf(a_fp, "name \"%s\" ", a_atomPropertyPtr->name);
	fprintf(a_fp, "japanese \"%s\" \n ", a_atomPropertyPtr->japaneseName);
	fprintf(a_fp, "minMassNumber %d ", a_atomPropertyPtr->minMassNumber);
	fprintf(a_fp, "maxMassNumber %d ", a_atomPropertyPtr->maxMassNumber);
	fprintf(a_fp, "min1secMassNumber %d ", a_atomPropertyPtr->min1secMassNumber);
	fprintf(a_fp, "max1secMassNumber %d\n", a_atomPropertyPtr->max1secMassNumber);
	fprintf(a_fp, " heatCapacity %lg [J/(mol*k)]\n", a_atomPropertyPtr->heatCapacity);
	for(massNumber = a_atomPropertyPtr->minMassNumber; massNumber <= a_atomPropertyPtr->maxMassNumber; ++massNumber){
		struct isotopeProperty * isotopePropertyPtr
		= getIsotopePropertyPtr(a_atomPropertyPtr->atomicNumber, massNumber);
		if(a_atomPropertyPtr->atomicNumber == isotopePropertyPtr->atomicNumber){
			printIsotopeProperty(a_fp, isotopePropertyPtr);
		}else{
			fprintf(a_fp, "ERROR:invalid atomicNumber:%d\n", isotopePropertyPtr->atomicNumber);
		}
	}

}
extern void printAtomProperties(FILE * a_fp)
{
	int atomicNumber;
	fputs("[Atom Properties]\n", a_fp);
	for(atomicNumber = MIN_ATOM_NUMBER; atomicNumber <= MAX_ATOM_NUMBER; ++atomicNumber){
		struct atomProperty * atomPropertyPtr;
		//fprintf(a_fp, "DEBUG:atomicNumber:%d\n", atomicNumber);
		atomPropertyPtr = getAtomPropertyPtr(atomicNumber);
		if(atomPropertyPtr){
			if(atomPropertyPtr->atomicNumber != UNDEF_ATOMIC_NUMBER){
				printAtomProperty(a_fp, atomPropertyPtr);
			}
		}else{
			//fprintf(a_fp, "DEBUG:atomPropertyPtr:NULL\n");
		}
	}
	fputs("\n", a_fp);
	e_logErrorUnregistedAtom = 1;
}
//---------------------------------------------------------------------
struct atomKey {
	const struct isotopeProperty * isotopePropertyPtr;
};
extern unsigned int calcHashSeedOfAtomKey(const void * a_keyPtr)
{
	const struct atomKey * ptr = (const struct atomKey *)a_keyPtr;
	unsigned ret = 0;
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_keyPtr);
	ret = (ptr->isotopePropertyPtr->atomicNumber << 10) + ptr->isotopePropertyPtr->massNumber;
	//fprintf(stderr, "DEBUG:%s:end}ret:%u\n", __FUNCTION__, ret);
	return ret;	
}
extern void * allocCopyAtomKey(const void * a_keyPtr)
{
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_keyPtr);
	void * ret = allocCopy(a_keyPtr, sizeof(struct atomKey), __FUNCTION__);
	//fprintf(stderr, "DEBUG:%s:end}ret:%lp\n", __FUNCTION__, ret);
	return ret;
}
#define VERSION_OF_serializeAtomKey 1
extern void * allocFreadAtomKey(FILE * a_fp)
{
	struct atomKey * ret = NULL;
	int atomicNumber, massNumber;
	FREAD_VERSION_CHECK(VERSION_OF_serializeAtomKey, a_fp)
	FREAD_VALUE_NULL(atomicNumber, a_fp)
	FREAD_VALUE_NULL(massNumber, a_fp)
	ret = clearAlloc(sizeof(struct atomKey), __FUNCTION__);
	if(ret){
		ret->isotopePropertyPtr = getIsotopePropertyPtr(atomicNumber, massNumber);
		if(!ret->isotopePropertyPtr){
			fprintf(stderr, "FATAL ERROR:%s:unknown isotope atomicNumber:%d massNumber:%d\n", __FUNCTION__, atomicNumber, massNumber);
			free(ret);
			ret = NULL;
		}
	}
	return ret;
}
extern int fwriteAtomKey(const void * a_keyPtr, FILE * a_fp)
{
	const struct atomKey * ptr = (const struct atomKey *)a_keyPtr;
	FWRITE_VERSION(VERSION_OF_serializeAtomKey, a_fp)
	FWRITE_VALUE(ptr->isotopePropertyPtr->atomicNumber, a_fp)
	FWRITE_VALUE(ptr->isotopePropertyPtr->massNumber, a_fp)
	return 1;
}
extern int compareAtomKey(const void * a_keyPtr1, const void * a_keyPtr2)
{
	int ret;
	const struct atomKey * m1 = (const struct atomKey *)a_keyPtr1;
	const struct atomKey * m2 = (const struct atomKey *)a_keyPtr2;
	//fprintf(stderr, "DEBUG:%s:begin{m1 atomicNumber:%d massNumber:%d m2 atomicNumber:%d massNumber:%d\n", __FUNCTION__, 
	//	m1->isotopePropertyPtr->atomicNumber, m1->isotopePropertyPtr->massNumber, 
	//	m2->isotopePropertyPtr->atomicNumber, m2->isotopePropertyPtr->massNumber);
	if(m1->isotopePropertyPtr->atomicNumber < m2->isotopePropertyPtr->atomicNumber){
		ret = -1;
	}else if(m1->isotopePropertyPtr->atomicNumber == m2->isotopePropertyPtr->atomicNumber){
		if(m1->isotopePropertyPtr->massNumber < m2->isotopePropertyPtr->massNumber){
			ret = -1;
		}else if(m1->isotopePropertyPtr->massNumber == m2->isotopePropertyPtr->massNumber){
			ret = 0;
		}else{
			ret = 1;			
		}
	}else{
		ret = 1;
	}
	//fprintf(stderr, "DEBUG:%s:end}ret:%d\n", __FUNCTION__, ret);
	return ret;
}

struct atomValue {
	//struct electrode * electrodePtr;
	double molIni;
	double molAdd, molAddMax, molAddMin;
	long long molAddCnt;
	double molSub, molSubMax, molSubMin;//Those values must be positive.
	long long molSubCnt;
	double pastSecond;
};
extern double getMol(const struct atomValue * a_atomValuePtr)
{
	double diff, mol;

	diff = a_atomValuePtr->molAdd - a_atomValuePtr->molSub;//for precision
	mol = a_atomValuePtr->molIni + diff;
	if(mol < 0.0){
		mol = 0.0;//collect tolelance.
	}
	return mol;
}
extern void setMolSub(struct atomValue * a_atomValuePtr, double a_molSub)
{
	if(a_molSub < 0.0){
		fprintf(stderr, "ERROR:%s:a_molSub:%lg < 0.0\n", __FUNCTION__, a_molSub);
		exit(1);
	}
	a_atomValuePtr->molSub += a_molSub;
	if(a_atomValuePtr->molSubCnt == 0){
		a_atomValuePtr->molSubMax = a_molSub;
		a_atomValuePtr->molSubMin = a_molSub;
	}else{
		if(a_atomValuePtr->molSubMax < a_molSub){
			a_atomValuePtr->molSubMax = a_molSub;
		}
		if(a_atomValuePtr->molSubMin > a_molSub){
			a_atomValuePtr->molSubMin = a_molSub;
		}
	}
	a_atomValuePtr->molSubCnt++;
}
extern void * allocCopyAtomValue(const void * a_valuePtr)
{
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_valuePtr);
	void * ret = allocCopy(a_valuePtr, sizeof(struct atomValue), __FUNCTION__);
	//fprintf(stderr, "DEBUG:%s:end}ret:%lp\n", __FUNCTION__, ret);
	return ret;
}
#define VERSION_OF_serializeAtomValue 1
extern void * allocFreadAtomValue(FILE * a_fp)
{
	struct atomValue * ret = NULL;
	FREAD_VERSION_CHECK(VERSION_OF_serializeAtomValue, a_fp)
	ret = clearAlloc(sizeof(struct atomValue), __FUNCTION__);
	if(ret){
		FREAD_MEMBER_FREE(ret, molIni, a_fp)
		FREAD_MEMBER_FREE(ret, molAdd, a_fp)
		FREAD_MEMBER_FREE(ret, molAddMax, a_fp)
		FREAD_MEMBER_FREE(ret, molAddMin, a_fp)
		FREAD_MEMBER_FREE(ret, molAddCnt, a_fp)
		FREAD_MEMBER_FREE(ret, molSub, a_fp)
		FREAD_MEMBER_FREE(ret, molSubMax, a_fp)
		FREAD_MEMBER_FREE(ret, molSubMin, a_fp)
		FREAD_MEMBER_FREE(ret, molSubCnt, a_fp)
		FREAD_MEMBER_FREE(ret, pastSecond, a_fp)
	}
	return ret;
}
extern int fwriteAtomValue(const void * a_valuePtr, FILE * a_fp)
{
	const struct atomValue * ptr = (const struct atomValue *)a_valuePtr;
	FWRITE_VERSION(VERSION_OF_serializeAtomValue, a_fp)
	FWRITE_VALUE(ptr->molIni, a_fp)
	FWRITE_VALUE(ptr->molAdd, a_fp)
	FWRITE_VALUE(ptr->molAddMax, a_fp)
	FWRITE_VALUE(ptr->molAddMin, a_fp)
	FWRITE_VALUE(ptr->molAddCnt, a_fp)
	FWRITE_VALUE(ptr->molSub, a_fp)
	FWRITE_VALUE(ptr->molSubMax, a_fp)
	FWRITE_VALUE(ptr->molSubMin, a_fp)
	FWRITE_VALUE(ptr->molSubCnt, a_fp)
	FWRITE_VALUE(ptr->pastSecond, a_fp)
	return 1;
}
extern void foundActionForAtomValue(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr)
{
	struct atomValue * m1 = (struct atomValue *)a_destPtr->valuePtr;
	const struct atomValue * m2 = (const struct atomValue *)a_srcPtr->valuePtr;
	double mol1, newMol;
	double newPastSecond;
	if(m2->molIni <= 0.0
	|| m2->molAdd != 0.0 || m2->molAddMax != 0.0 || m2->molAddMin != 0.0 || m2->molAddCnt != 0 
	|| m2->molSub != 0.0 || m2->molSubMax != 0.0 || m2->molSubMin != 0.0 || m2->molSubCnt != 0){
		fprintf(stderr, "ERROR:%s:wrong values, molIni:%lg molAdd:%lg molAddMax:%lg molAddMin:%lg molAddCnt:%lld molSub:%lg molSubMax:%lg molSubMin:%lg molSubCnt:%lld\n", __FUNCTION__, m2->molIni, m2->molAdd, m2->molAddMax, m2->molAddMin, m2->molAddCnt, m2->molSub, m2->molSubMax, m2->molSubMin, m2->molSubCnt);
		exit(1);
	}
	mol1 = getMol(m1);
	newMol = mol1 + m2->molIni;
	if(newMol > 0.0){
		newPastSecond = (m1->pastSecond * mol1 + m2->pastSecond * m2->molIni) / newMol;
	}else{
		newPastSecond = 0;
	}
	m1->pastSecond = newPastSecond;

	m1->molAdd += m2->molIni;
	if(m1->molAddCnt == 0){
		m1->molAddMax = m2->molIni;
		m1->molAddMin = m2->molIni;
	}else{
		if(m1->molAddMax < m2->molIni){
			m1->molAddMax = m2->molIni;
		}
		if(m1->molAddMin > m2->molIni){
			m1->molAddMin = m2->molIni;
		}
	}
	m1->molAddCnt++;
}
extern double printAtomList(FILE * a_fp, const struct hashTable * a_atomHashTablePtr, double * a_sumOfMassUMolIniPtr, double * a_sumOfMassUMolAddPtr, double * a_sumOfMassUMolSubPtr, double * a_heatCapacityPtr)
{
	int needSortedTable = 1;
	struct objectNodeConst ** nodePtrs;
	unsigned int size, i;
	double sumOfMassUMol, sumOfMassUMolDiff;
	*a_sumOfMassUMolIniPtr = 0.0;
	*a_sumOfMassUMolAddPtr = 0.0;
	*a_sumOfMassUMolSubPtr = 0.0;
	*a_heatCapacityPtr = 0.0;
	fprintf(a_fp, "%s Atoms\n", a_atomHashTablePtr->tableName);
	nodePtrs = getFlatTable(a_atomHashTablePtr, needSortedTable, &size);
	if(nodePtrs){
		for(i = 0; i < size; ++i){
			struct atomKey * keyPtr = (struct atomKey *)nodePtrs[i]->keyPtr;
			struct atomValue * valuePtr = (struct atomValue *)nodePtrs[i]->valuePtr;
			double mol;
			if(keyPtr->isotopePropertyPtr->symbol[0]){
				fprintf(a_fp, "%d %s", i + 1, keyPtr->isotopePropertyPtr->symbol);
			}else{
				fprintf(a_fp, "%d %d%s", i + 1, keyPtr->isotopePropertyPtr->massNumber, keyPtr->isotopePropertyPtr->atomPropertyPtr->symbol);
			}
			if(keyPtr->isotopePropertyPtr->halfLifeUnit == HLU_STABLE){
				fprintf(a_fp, " %s", getHalfLifeUnitString(keyPtr->isotopePropertyPtr->halfLifeUnit));
			}else{
				fprintf(a_fp, " halfLife %lg [%s]", keyPtr->isotopePropertyPtr->halfLife, getHalfLifeUnitString(keyPtr->isotopePropertyPtr->halfLifeUnit));
			}
			mol = getMol(valuePtr);
			*a_heatCapacityPtr += (mol * keyPtr->isotopePropertyPtr->atomPropertyPtr->heatCapacity);
			fprintf(a_fp, " %lg [mol]\n", mol);
			if(valuePtr->molAdd > 0.0 && valuePtr->molSub == 0.0){
				fprintf(a_fp, " (= %lg + %lg (=[%lg , %lg] * %lld)[mol])\n", valuePtr->molIni, valuePtr->molAdd, valuePtr->molAddMin, valuePtr->molAddMax, valuePtr->molAddCnt);
			}else if(valuePtr->molAdd == 0.0 && valuePtr->molSub > 0.0){
				fprintf(a_fp, " (= %lg - %lg (=[%lg , %lg] * %lld)[mol])\n", valuePtr->molIni, valuePtr->molSub, valuePtr->molSubMin, valuePtr->molSubMax, valuePtr->molSubCnt);
			}else if(valuePtr->molAdd > 0.0 && valuePtr->molSub > 0.0){
				fprintf(a_fp, " (= %lg + %lg (=[%lg , %lg] * %lld) - %lg (=[%lg , %lg] * %lld)[mol])\n", valuePtr->molIni, valuePtr->molAdd, valuePtr->molAddMin, valuePtr->molAddMax, valuePtr->molAddCnt, valuePtr->molSub, valuePtr->molSubMin, valuePtr->molSubMax, valuePtr->molSubCnt);
			}
			*a_sumOfMassUMolIniPtr += (keyPtr->isotopePropertyPtr->massU * valuePtr->molIni);
			*a_sumOfMassUMolAddPtr += (keyPtr->isotopePropertyPtr->massU * valuePtr->molAdd);
			*a_sumOfMassUMolSubPtr += (keyPtr->isotopePropertyPtr->massU * valuePtr->molSub);
		}
		free(nodePtrs);
		sumOfMassUMolDiff = *a_sumOfMassUMolAddPtr - *a_sumOfMassUMolSubPtr;//for precision
		sumOfMassUMol = *a_sumOfMassUMolIniPtr + sumOfMassUMolDiff;
		fprintf(a_fp, "%s Sum of Mass * Mol : %lg (= %lg + %lg - %lg) [U]\n", a_atomHashTablePtr->tableName, sumOfMassUMol, *a_sumOfMassUMolIniPtr, *a_sumOfMassUMolAddPtr, *a_sumOfMassUMolSubPtr);
		fprintf(a_fp, " = %lg [MeV]\n", sumOfMassUMol * NAvogadro * e_coefMassUToMeV);
		fprintf(a_fp, " heatCapacity = %lg [J/k]\n\n", *a_heatCapacityPtr);
	}
	return sumOfMassUMol;
}
//---------------------------------------------------------------------
//[caluculated user precizions]
int e_2raisedTothePowerPrecision = 20;//2^20=1048576 > 10^6, 20 bit of binary means that we request the precision of 6 digit in decimal. 2^7=128 > 10^2, 7 bit of binary means that we request the precision of 2 digit in decimal.
#define DOUBLE_EXPONENT_VALUE_BIT 10
//Double-precision floating-point format(IEEE standard) has the 11bit exponent(1 signbit and 10 value bit, [-2^10, 2^10]). 
//And, (2^10 - e_2raisedTothePowerPrecision) / (log(10)/log(2)) = (2^10 - 20) / (log(10)/log(2)) = 302.23.
double e_positiveNearZero = 1.0e-301;
double e_negativeNearZero = -1.0e-301;

extern void initNearZero()
{
	int miniExp;
	e_2raisedTothePowerPrecision = ceil(log2(pow(10.0, (double)e_rangeDecimalDigitPrecision)));
	miniExp = 1 - (int)floor((pow(2.0, DOUBLE_EXPONENT_VALUE_BIT) - e_2raisedTothePowerPrecision) / (log(10.0)/log(2.0)));
	e_positiveNearZero = pow(10, miniExp);
	e_negativeNearZero = -e_positiveNearZero;
	//fprintf(strerr, "DEBUG:%s:e_rangeDecimalDigitPrecision:%d\n", __FUNCTION__, e_rangeDecimalDigitPrecision);
	//fprintf(strerr, "DEBUG:%s:e_2raisedTothePowerPrecision:%d\n", __FUNCTION__, e_2raisedTothePowerPrecision);
	//fprintf(strerr, "DEBUG:%s:miniExp:%d\n", __FUNCTION__, miniExp);
	//fprintf(strerr, "DEBUG:%s:e_positiveNearZero:%lg\n", __FUNCTION__, e_positiveNearZero);
	//fprintf(strerr, "DEBUG:%s:e_negativeNearZero:%lg\n", __FUNCTION__, e_negativeNearZero);
	//fprintf(strerr, "DEBUG:%s:sizeof(int):%d\n", __FUNCTION__, sizeof(int));// 4 in GCC
	//fprintf(strerr, "DEBUG:%s:sizeof(long int):%d\n", __FUNCTION__, sizeof(long int));//8 in GCC
	//fprintf(strerr, "DEBUG:%s:sizeof(long long):%d\n", __FUNCTION__, sizeof(long long));//8 in GCC
	//fprintf(strerr, "DEBUG:%s:sizeof(double):%d\n", __FUNCTION__, sizeof(double));//8 in GCC
	//fprintf(strerr, "DEBUG:%s:sizeof(long double):%d\n", __FUNCTION__, sizeof(long double));//16 in GCC
	//exit(0);//DEBUG
}
struct range {
	double high;//high is not included.
	double low;//low is included.
};
extern int compareValueAndRange(double a_valuePtr, const struct range * a_range)
{
	int ret;
	if(a_valuePtr < a_range->low){
		ret = -1;
	}else if(a_valuePtr >= a_range->high){
		ret = 1;
	}else{// a_range->low <= a_valuePtr && a_valuePtr < a_range->high
		ret = 0;
	}
	return ret;
}
extern struct range calcPrecisionRange(double a_valuePtr)
{
	double absolute, exponent, floatingScale, smallMagnitude, bigMagnitude;
	struct range ret;
	if(e_negativeNearZero <= a_valuePtr && a_valuePtr < e_positiveNearZero ){
		ret.high = e_positiveNearZero;
		ret.low = e_negativeNearZero;
	}else{
		if(a_valuePtr > 0.0){
			absolute = a_valuePtr;
		}else{
			absolute = - a_valuePtr;
		}
		exponent = logb(absolute);//logb is a very fast function.
		floatingScale = exp2(e_2raisedTothePowerPrecision - exponent);//exp2 is a very fast function.
		if(a_valuePtr > 0.0){
			smallMagnitude = floor(absolute * floatingScale);
			bigMagnitude = smallMagnitude + 1.0;
			ret.high = bigMagnitude / floatingScale;
			ret.low = smallMagnitude / floatingScale;
			if(ret.low < e_positiveNearZero){
				ret.low = e_positiveNearZero;
			}
		}else{
			bigMagnitude = ceil(absolute * floatingScale);
			smallMagnitude = bigMagnitude - 1.0;
			ret.high = - smallMagnitude / floatingScale;
			if(ret.high > e_negativeNearZero){
				ret.high = e_negativeNearZero;
			}
			ret.low = - bigMagnitude / floatingScale;
		}
	}
	return ret;
}

/*
extern void debugCalcRange()
{
	int i, comp;
	for(i = -1000; i <= 1000; ++i){
		double r1, r2, x;
		struct range rg, rc, rd;
		r1 = (rand() % 1000000) - 500000;
		r2 /= 1000000;
		r2 = (rand() % 1000000) - 500000;
		r2 /= 1000000;
		x = r1 * r2 * exp2((double)i);
		if(i == 0){
			x = 0.0;
		}
		rg = calcPrecisionRange(x);
		fprintf(stderr, "%d:", i);
		if(rg.low <= x && x < rg.high){
			fprintf(stderr, "OK:");
		}else{
			if(rg.low >= rg.high){
				fprintf(stderr, "ERROR:rg.low >= rg.high\n   ");
				break;
			}
			if(rg.low > x){
				fprintf(stderr, "ERROR:rg.low > x\n   ");
				break;
			}
			if(x >= rg.high){
				fprintf(stderr, "ERROR:x >= rg.high\n   ");
				break;
			}
		}
		fprintf(stderr, ":%20.14lg %20.14lg %20.14lg\n", rg.low, x, rg.high);
		rc = calcPrecisionRange(rg.low);
		if(rc.low != rg.low){
			fprintf(stderr, "ERROR:rc.low:%20.14lg != rg.low:%20.14lg\n", rc.low, rg.low);
			break;
		}
		rd = calcPrecisionRange(rg.high);
		if(rd.low != rg.high){
			fprintf(stderr, "ERROR:rd.low:%20.14lg != rg.high:%20.14lg\n", rd.low, rg.high);
			break;
		}
		comp = compareValueAndRange(x, &rg);
		if(comp != 0){
			fprintf(stderr, "ERROR:(comp = compareValueAndRange(x, rg)):%d != 0\n", comp);
			break;			
		}
		comp = compareValueAndRange(rg.low, &rg);
		if(comp != 0){
			fprintf(stderr, "ERROR:(comp = compareValueAndRange(rg.low, rg)):%d != 0\n", comp);
			break;			
		}
		comp = compareValueAndRange(rg.high, &rg);
		if(comp != 1){
			fprintf(stderr, "ERROR:(comp = compareValueAndRange(rg.high, rg)):%d != 1\n", comp);
			break;			
		}
		comp = compareValueAndRange(rg.low - (rg.high - rg.low), &rg);
		if(comp != -1){
			fprintf(stderr, "ERROR:(comp = compareValueAndRange(rg.low - (rg.high - rg.low), rg)):%d != -1\n", comp);
			break;			
		}
	}
	exit(0);
}
*/
extern unsigned int calcHashSeedOfRange(const void * a_keyPtr)
{
	const unsigned int * ptr; 
	unsigned int i, ret;
	ptr = (const unsigned int *)&(((const struct range *)a_keyPtr)->low);
	ret = 0;
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_keyPtr);
	for(i = 0; i < (sizeof(double) / sizeof(unsigned int)); ++i){
		ret = ((ret << 1) | (ret >> (sizeof(unsigned int) * 8 - 1))) + ptr[i];
	}
	//fprintf(stderr, "DEBUG:%s:end}ret:%u\n", __FUNCTION__, ret);
	return ret;
}
extern void * allocCopyRange(const void * a_keyPtr)
{
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_keyPtr);
	void * ret = allocCopy(a_keyPtr, sizeof(struct range), __FUNCTION__);
	//fprintf(stderr, "DEBUG:%s:end}ret:%lp\n", __FUNCTION__, ret);
	return ret;
}
#define VERSION_OF_serializeRange 1
extern void * allocFreadRange(FILE * a_fp)
{
	struct range * ret = NULL;
	FREAD_VERSION_CHECK(VERSION_OF_serializeRange, a_fp)
	ret = clearAlloc(sizeof(struct range), __FUNCTION__);
	if(ret){
		FREAD_MEMBER_FREE(ret, low, a_fp)
		FREAD_MEMBER_FREE(ret, high, a_fp)
	}
	return ret;
}
extern int fwriteRange(const void * a_keyPtr, FILE * a_fp)
{
	const struct range * ptr = (const struct range *)a_keyPtr;
	FWRITE_VERSION(VERSION_OF_serializeRange, a_fp)
	FWRITE_VALUE(ptr->low, a_fp)
	FWRITE_VALUE(ptr->high, a_fp)
	return 1;
}
extern int comparRange(const void * a_keyPtr1, const void * a_keyPtr2)
{
	int ret = 0;
	const struct range * m1 = (const struct range *)a_keyPtr1;
	const struct range * m2 = (const struct range *)a_keyPtr2;
	//fprintf(stderr, "DEBUG:%s:begin{m2->low:%lg m2->low:%lg\n", __FUNCTION__, m1->low, m2->low);
	if(m1->low < m2->low){
		ret = -1;
	}else if(m1->low > m2->low){
		ret = 1;
	}
	//fprintf(stderr, "DEBUG:%s:end}ret:%d\n", __FUNCTION__, ret);
	return ret;
}
//---------------------------------------------------------------------

struct massDefect {
	double MeV;
	double mol;
};

extern void * allocCopyMassDefect(const void * a_valuePtr)
{
	//fprintf(stderr, "DEBUG:%s:begin{a_valuePtr:%lp\n", __FUNCTION__, a_valuePtr);
	void * ret =  allocCopy(a_valuePtr, sizeof(struct massDefect), __FUNCTION__);
	//fprintf(stderr, "DEBUG:%s:end}ret:%lp\n", __FUNCTION__, ret);
	return ret;
}
#define VERSION_OF_serializeMassDefect 1
extern void * allocFreadMassDefect(FILE * a_fp)
{
	struct massDefect * ret = NULL;
	FREAD_VERSION_CHECK(VERSION_OF_serializeMassDefect, a_fp)
	ret = clearAlloc(sizeof(struct massDefect), __FUNCTION__);
	if(ret){
		FREAD_MEMBER_FREE(ret, MeV, a_fp)
		FREAD_MEMBER_FREE(ret, mol, a_fp)
	}
	return ret;
}
extern int fwriteMassDefect(const void * a_keyPtr, FILE * a_fp)
{
	const struct massDefect * ptr = (const struct massDefect *)a_keyPtr;
	FWRITE_VERSION(VERSION_OF_serializeMassDefect, a_fp)
	FWRITE_VALUE(ptr->MeV, a_fp)
	FWRITE_VALUE(ptr->mol, a_fp)
	return 1;
}
extern void foundActionForMassDefect(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr)
{
	struct massDefect * m1 = (struct massDefect *)a_destPtr->valuePtr;
	const struct massDefect * m2 = (const struct massDefect *)a_srcPtr->valuePtr;
	//fprintf(stderr, "DEBUG:%s:begin{m1->mol:%lg m2->mol:%lg\n", __FUNCTION__, m1->mol, m2->mol);
	if(m1->MeV != 0.0){
		m1->mol += (m2->mol * m2->MeV / m1->MeV);
	}else{
		fprintf(stderr, "ERROR:%s:m1->MeV == 0.0, m1->mol:%lg m2->mol:%lg m2->MeV:%lg\n", __FUNCTION__, m1->mol, m2->mol, m2->MeV);
		exit(1);
	}
	//fprintf(stderr, "DEBUG:%s:end}m1->mol:%lg m2->mol:%lg\n", __FUNCTION__, m1->mol, m2->mol);
}
//---------------------------------------------------------------------
struct sumMeVMol{
		double minMeV;
		double maxMeV;
		double minMol;
		double maxMol;
		long long cnt;
		long double MeVMol;
};
#define VERSION_OF_serializesumMeVMol 1
#define SERIALIZE_SUMMEVMOL(RW) \
extern int RW ## SumMeVMol(struct sumMeVMol * a_sumMeVMolPtr, FILE * a_fp)\
{\
	SERIALIZE_VERSION_CHECK(RW, VERSION_OF_serializesumMeVMol, a_fp)\
	SERIALIZE_VALUE(RW, a_sumMeVMolPtr->minMeV, a_fp)\
	SERIALIZE_VALUE(RW, a_sumMeVMolPtr->maxMeV, a_fp)\
	SERIALIZE_VALUE(RW, a_sumMeVMolPtr->minMol, a_fp)\
	SERIALIZE_VALUE(RW, a_sumMeVMolPtr->maxMol, a_fp)\
	SERIALIZE_VALUE(RW, a_sumMeVMolPtr->cnt, a_fp)\
	SERIALIZE_VALUE(RW, a_sumMeVMolPtr->MeVMol, a_fp)\
	return 1;\
}
SERIALIZE_SUMMEVMOL(fread)//freadSumMeVMol
SERIALIZE_SUMMEVMOL(fwrite)//fwriteSumMeVMol
#define CALL_SERIALIZE_SUMMEVMOL(RW, SUMMEVMOLPTR, FP) \
	RW ## SumMeVMol(SUMMEVMOLPTR, FP)

extern void initSumMeVMol(struct sumMeVMol * a_sumMeVMol)
{
	a_sumMeVMol->minMeV = 0.0;
	a_sumMeVMol->maxMeV = 0.0;
	a_sumMeVMol->minMol = 0.0;
	a_sumMeVMol->maxMol = 0.0;
	a_sumMeVMol->cnt = 0;
	a_sumMeVMol->MeVMol = 0.0;
}
extern void addSumMeVMol(struct sumMeVMol * a_sumMeVMol, double a_Mev, double a_mol)
{
	if(a_sumMeVMol->cnt == 0){
		a_sumMeVMol->minMeV = a_Mev;
		a_sumMeVMol->maxMeV = a_Mev;
		a_sumMeVMol->minMol = a_mol;
		a_sumMeVMol->maxMol = a_mol;
		a_sumMeVMol->cnt = 1;
		a_sumMeVMol->MeVMol = a_Mev * a_mol;
	}else{
		if(a_sumMeVMol->minMeV > a_Mev){
			a_sumMeVMol->minMeV = a_Mev;
		}
		if(a_sumMeVMol->maxMeV < a_Mev){
			a_sumMeVMol->maxMeV = a_Mev;
		}
		if(a_sumMeVMol->minMol > a_mol){
			a_sumMeVMol->minMol = a_mol;
		}
		if(a_sumMeVMol->maxMol < a_mol){
			a_sumMeVMol->maxMol = a_mol;
		}
		a_sumMeVMol->cnt++;
		a_sumMeVMol->MeVMol += (a_Mev * a_mol);
	}
}
extern void mergeSumMeVMol(struct sumMeVMol * a_merge, struct sumMeVMol * a_src1, struct sumMeVMol * a_src2)
{
	a_merge->minMeV = (a_src1->minMeV < a_src2->minMeV) ? a_src1->minMeV : a_src2->minMeV;
	a_merge->maxMeV = (a_src1->maxMeV > a_src2->maxMeV) ? a_src1->maxMeV : a_src2->maxMeV;
	a_merge->minMol = (a_src1->minMol < a_src2->minMol) ? a_src1->minMol : a_src2->minMol;
	a_merge->maxMol = (a_src1->maxMol > a_src2->maxMol) ? a_src1->maxMol : a_src2->maxMol;
	a_merge->cnt = a_src1->cnt + a_src2->cnt;
	a_merge->MeVMol = a_src1->MeVMol + a_src2->MeVMol;
}

#define SCAT_BIG_MASS_DEFECT_NOW -3 
#define SCAT_BIG_MASS_DEFECT_ALL -2//The big mass defect type is the part of energy by mass defection from either cillide of nucleus particles or decay of nuecuses and the amont of energy is eaual or bigger than "e_collideElectronMiniMeV".
#define SCAT_LOST_BY_NEUTRINO -1 //the enegey type carryed out by neutrino in beta decays.

#define SCAT_NOT_COLLIDE 0 // The not collide type is a part of enegy of input high voltage pulse current. The part is not used to collide of nucleus perticles by high voltage pulse current.
#define SCAT_IMPERFECT 1  // The IMPERFECT type is a part of enegy of input high voltage pulse current.
#define SCAT_IMPERFECT_IN_SPACE 2  // The IMPERFECT type is a part of enegy of input high voltage pulse current.
#define SCAT_SMALL_MASS_DEFECT 3 // The small mass defect type is the output energy from neclear actions. But they are less than "e_collideElectronMiniMeV". Therefore they can't cause the next neclear action anymore.
#define SCAT_SMALL_BY_COMPTON 4 // This Comton effect type is the reflection phton with less energy than "e_collideElectronMiniMeV". Therefore they can't cause the next neclear action anymore.
#define SCAT_SMALL_BY_COMPTON_E 5 // This Comton effect type is the enegy passed to an electron from a photon. But the amount of energy is less than "e_collideElectronMiniMeV". Therefore they can't cause the next neclear action anymore.


#define SCAT_DEGRADE_IN_COMPTON_B 6 //The too many SCAT_BIG_MASS_DEFECT_NOW nodes in the hash table, "massDefectHashTable", they waste the memory and power of computer. So, we do degreade them by shopping the energy into 32 greades. The SCAT_DEGRADE_IN_COMPTON_B type contains shopped small remain part of BIG energy. The grading big energy is still in the "massDefectHashTable".
#define SCAT_DEGRADE_IN_COMPTON_S 7
#define SIZE_OF_SCATTERD 8 //The array size of scatterd photon types. Those types of photon are reason of heat.

extern const char * getScatName(int a_scatType)
{
	char * retPtr = "";
	switch(a_scatType){
		case SCAT_BIG_MASS_DEFECT_NOW: retPtr = "BIG_MASS_DEFECT_NOW"; break;
		case SCAT_BIG_MASS_DEFECT_ALL: retPtr = "BIG_MASS_DEFECT_ALL"; break;
		case SCAT_LOST_BY_NEUTRINO: retPtr = "LOST_BY_NEUTRINO"; break;
		
		case SCAT_NOT_COLLIDE: retPtr = "NOT_COLLIDE"; break;
		case SCAT_IMPERFECT: retPtr = "IMPERFECT:"; break;
		case SCAT_IMPERFECT_IN_SPACE: retPtr = "IMPERFECT_IN_SPACE"; break;
		case SCAT_SMALL_MASS_DEFECT: retPtr = "SMALL_MASS_DEFECT"; break;
		case SCAT_SMALL_BY_COMPTON: retPtr = "SMALL_BY_COMPTON"; break;
		case SCAT_SMALL_BY_COMPTON_E: retPtr = "SMALL_BY_COMPTON_E"; break;
		case SCAT_DEGRADE_IN_COMPTON_B: retPtr = "DEGRADE_IN_COMPTON_B"; break;
		case SCAT_DEGRADE_IN_COMPTON_S: retPtr = "DEGRADE_IN_COMPTON_S"; break;
	}
	return retPtr;
}
extern void printSumMeVMol(FILE * a_fp, const char * a_operatorStr, int a_scatType, const struct sumMeVMol * a_sumMeVMol, long double a_grandTotal, const char * a_rateName)
{
	long double MeVMolNAvogadro, rate = 0.0;
	MeVMolNAvogadro = a_sumMeVMol->MeVMol * NAvogadro;
	fprintf(a_fp, "%s %s %Lg [MeV] %lld[cnt]",
		a_operatorStr, getScatName(a_scatType), MeVMolNAvogadro, a_sumMeVMol->cnt);
	if(a_grandTotal > 0.0){
		rate = 100.0 * MeVMolNAvogadro / a_grandTotal;
		fprintf(a_fp, " %Lg %s", rate, a_rateName);		
	}
	fprintf(a_fp, "\n (range[ %lg , %lg ][MeV] * [ %lg , %lg ][mol])\n", 
			a_sumMeVMol->minMeV, a_sumMeVMol->maxMeV,
			a_sumMeVMol->minMol, a_sumMeVMol->maxMol);
}
extern void mergeScattered(struct sumMeVMol * a_mergeAry, struct sumMeVMol * a_src1, struct sumMeVMol * a_src2)
{
	int i;
	for(i = 0; i < SIZE_OF_SCATTERD; ++i){
		mergeSumMeVMol(&a_mergeAry[i], &a_src1[i], &a_src2[i]);
	}
}
extern long double calcGrandTotalOfScattered(struct sumMeVMol * a_scatteredAry)
{
	int i;
	long double grandTotal = 0.0;
	for(i = 0; i < SIZE_OF_SCATTERD; ++i){
		grandTotal += a_scatteredAry[i].MeVMol;
	}
	grandTotal *= NAvogadro;
	return grandTotal;
}
extern long double printScattered(FILE * a_fp, const char * a_titlePtr, struct sumMeVMol * a_scatteredAry)
{
	int i;
	long double grandTotal = calcGrandTotalOfScattered(a_scatteredAry);
	fprintf(a_fp, "%s OUTPUT HEAT =\n", a_titlePtr);
	for(i = 0; i < SIZE_OF_SCATTERD; ++i){
		printSumMeVMol(a_fp, (i == 0) ? " " : "+", i, &a_scatteredAry[i], grandTotal, "[%]");
	}
	fprintf(a_fp, "= %Lg [MeV]\n",  grandTotal);
	return grandTotal;
}

//---------------------------------------------------------------------
struct electrode {
	double detectLimitMolForIsotope;
	struct hashTable atomHashTable;
	struct hashTable massDefectHashTable;
	struct sumMeVMol massDefectAll;
	struct sumMeVMol byNeutrinoAll;
	struct sumMeVMol scattered[SIZE_OF_SCATTERD];
};

struct electrode e_negativeElectrode;
struct electrode e_positiveElectrode;

#define VERSION_OF_serializeElectrode 1

#define SERIALIZE_ELECTRODE(RW) \
extern int RW ## Electrode(struct electrode * a_electrodePtr, FILE * a_fp)\
{\
	int i;\
	SERIALIZE_VERSION_CHECK(RW, VERSION_OF_serializeElectrode, a_fp)\
	SERIALIZE_VALUE(RW, a_electrodePtr->detectLimitMolForIsotope, a_fp)\
	if((void *)RW == (void *)fwrite){\
		if(!fwritedHashTable(&a_electrodePtr->atomHashTable, a_fp, fwriteAtomKey, fwriteAtomValue)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable, a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
	}else{\
		if(!freadHashTable(&a_electrodePtr->atomHashTable, a_fp,\
			allocFreadAtomKey, allocFreadAtomValue, \
			calcHashSeedOfAtomKey,\
			allocCopyAtomKey,\
			free,\
			compareAtomKey,\
			allocCopyAtomValue,\
			free,\
			foundActionForAtomValue)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable, a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
	}\
	CALL_SERIALIZE_SUMMEVMOL(RW, &a_electrodePtr->massDefectAll, a_fp);\
	CALL_SERIALIZE_SUMMEVMOL(RW, &a_electrodePtr->byNeutrinoAll, a_fp);\
	for(i = 0; i < SIZE_OF_SCATTERD; ++i){\
		CALL_SERIALIZE_SUMMEVMOL(RW, a_electrodePtr->scattered + i, a_fp);\
	}\
	return 1;\
}
#define CALL_SERIALIZE_ELECTRODE(RW, ELECTRODEPTR, FP) \
	RW ## Electrode(ELECTRODEPTR, FP)

SERIALIZE_ELECTRODE(fread)//freadElectrode
SERIALIZE_ELECTRODE(fwrite)//fwriteElectrode

//---------------------------------------------------------------------

struct atomNodeConst {
	unsigned int hashSeed;
	const struct atomKey * key;
	struct atomValue * vPtr;
};

extern struct atomNodeConst * findAtom(struct electrode * a_electrodePtr, int a_atomicNumber, int a_massNumber)
{
	struct atomNodeConst * foundPtr;
	struct atomKey ak;
	ak.isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
	if(ak.isotopePropertyPtr){
		struct objectNodeConst search;
		search.keyPtr = &ak;
		foundPtr = (struct atomNodeConst *)findNodeInHashTable(&a_electrodePtr->atomHashTable, &search);
	}else{
		foundPtr = NULL;
	}
	return foundPtr;
}
extern struct atomNodeConst * findElectronInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON);
}
extern struct atomNodeConst * findNeutronInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON);
}

extern void increaseAtom(struct electrode * a_electrodePtr, int a_atomicNumber, int a_massNumber, double a_increaseMol)
{
	//fprintf(stderr, "DEBUG:%s:BEGIN{a_atomicNumber:%d massNumber:%d a_increaseMol:%lg \n", __FUNCTION__, a_atomicNumber, a_massNumber, a_increaseMol);
	if(a_increaseMol >= 0.0){
		//struct atomNodeConst * atomPtr;
		struct atomKey akey;
		akey.isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
		if(akey.isotopePropertyPtr){
			struct atomValue aValue;
			struct objectNodeConst nodeConst;
			memset(&aValue, 0, sizeof(aValue));
			aValue.molIni = a_increaseMol;
			nodeConst.keyPtr = &akey;
			nodeConst.valuePtr = &aValue;
			if(a_increaseMol > 1.0){//DEBUG
				fprintf(stderr, "DEBUG:%s(%d):%s a_atomicNumber:%d massNumber:%d a_increaseMol:%lg \n", __FUNCTION__, __LINE__, a_electrodePtr->atomHashTable.tableName, a_atomicNumber, a_massNumber, a_increaseMol);
				exit(1);
			}
			insertObjectInHashTable(&a_electrodePtr->atomHashTable, &nodeConst);
		}
	}else{
		struct atomNodeConst * atomPtr;
		atomPtr = findAtom(a_electrodePtr, a_atomicNumber, a_massNumber);
		if(atomPtr){
			setMolSub(atomPtr->vPtr, - a_increaseMol);
		}else{
			fprintf(stderr, "ERROR:%s:Not find a_atomicNumber:%d massNumber:%d, a_increaseMol:%lg <= 0.0 \n", __FUNCTION__, a_atomicNumber, a_massNumber, a_increaseMol);
		}
	}
	//fprintf(stderr, "DEBUG:%s:END}\n", __FUNCTION__);
}
extern double assignStableIstopesOnElectrode(struct electrode * a_electrodePtr, int a_atomicNumber, double a_mol)
{
	//double molSum = 0.0;
	double electronMol = a_atomicNumber * a_mol;
	struct atomProperty * atomPropertyPtr = getAtomPropertyPtr(a_atomicNumber);
	//fprintf(stderr, "DEBUG:%s:BEGIN{ a_atomicNumber:%d a_mol:%lg \n", __FUNCTION__, a_atomicNumber, a_mol);
	if(atomPropertyPtr){
		int massNumber;
		for(massNumber = atomPropertyPtr->minMassNumber; massNumber <= atomPropertyPtr->maxMassNumber; ++massNumber){
			struct atomKey akey;
			akey.isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, massNumber);
			if(akey.isotopePropertyPtr){
				if(akey.isotopePropertyPtr->relativeIsotopicAbundance > 0.0){
					struct atomValue aValue;
					memset(&aValue, 0, sizeof(aValue));
					aValue.molIni = a_mol * akey.isotopePropertyPtr->relativeIsotopicAbundance;
					if(aValue.molIni > 0.0){
						struct objectNodeConst nodeConst;
						nodeConst.keyPtr = &akey;
						nodeConst.valuePtr = &aValue;
						//fprintf(stderr, "DEBUG:%s:a_atomicNumber:%d massNumber:%d molIni:%lg \n", __FUNCTION__, a_atomicNumber, massNumber, aValue.molIni);
						//molSum += aValue.molIni;
						insertObjectInHashTable(&a_electrodePtr->atomHashTable, &nodeConst);
					}
				}
			}
		}
	}
	//fprintf(stderr, "DEBUG:%s:END} molSum:%lg \n", __FUNCTION__, molSum);
	return electronMol;
}
#define NAME_MAX 32

extern double registAtomOnElectrodes(struct electrode * a_electrodePtr, char * a_name, char * a_val, int a_valLen, double * a_sumOfMolPtr)
{
	double electronMol = 0.0;
	double molValue;
	char atomicName[NAME_MAX];
	int massVal;
	struct atomProperty * atomP;
	struct isotopeProperty * isoP;

	*a_sumOfMolPtr = 0.0;
	a_val[a_valLen] = 0;
	//fprintf(stderr, "DEBUG:%s:{BEGIN a_name:%s a_val:%s a_valLen:%d\n", __FUNCTION__, a_name, a_val, a_valLen);
	if(a_valLen <= 0){
		fprintf(stderr, "ERROR:%s:none of value of name:%s\n", __FUNCTION__, a_name);
		exit(1);
	}
	sscanf(a_val, "%lf", &molValue);
	if('0' <= a_name[0] && a_name[0] <= '9'){
		if('0' <= a_name[1] && a_name[1] <= '9'){
			if('0' <= a_name[2] && a_name[2] <= '9'){
				massVal = (a_name[0] - '0') * 100 + (a_name[1] - '0') * 10 + (a_name[2] - '0');
				strcpy(atomicName, a_name + 3);
			}else{
				massVal = (a_name[0] - '0') * 10 + (a_name[1] - '0');
				strcpy(atomicName, a_name + 2);
			}
		}else{
			massVal = (a_name[0] - '0');
			strcpy(atomicName, a_name + 1);
		}
		atomP = findAtomPropertyPtr(atomicName);
		if(atomP){
			isoP = getIsotopePropertyPtr(atomP->atomicNumber, massVal);
			if(isoP){
				//fprintf(stderr, "DEBUG:%s:%d%s=%lg\n", __FUNCTION__, massVal, atomicName, molValue);
				increaseAtom(a_electrodePtr, atomP->atomicNumber, massVal, molValue);
				electronMol += (atomP->atomicNumber * molValue);
				*a_sumOfMolPtr += molValue;
				//heatCapacity += (molValue * atomP->heatCapacity);
			}else{
				fprintf(stderr, "ERROR:%s:BAD massVal %d of name:%s\n", __FUNCTION__, massVal, a_name);
				exit(1);
			}
		}else{
			fprintf(stderr, "ERROR:%s:BAD atomicName %s of name:%s\n", __FUNCTION__, atomicName, a_name);
			exit(1);
		}
	}else if(strcmp(a_name, "D") == 0){
		increaseAtom(a_electrodePtr, ATOMICNUMBER_DEUTERIUM, MASSNUMBER_DEUTERIUM, molValue);
		electronMol += (ATOMICNUMBER_DEUTERIUM * molValue);
		*a_sumOfMolPtr += molValue;
	}else if(strcmp(a_name, "T") == 0){
		increaseAtom(a_electrodePtr, ATOMICNUMBER_TRITIUM, MASSNUMBER_TRITIUM, molValue);
		electronMol += (ATOMICNUMBER_TRITIUM * molValue);
		*a_sumOfMolPtr += molValue;
	}else{
		strcpy(atomicName, a_name);
		atomP = findAtomPropertyPtr(atomicName);
		if(atomP){
			//fprintf(stderr, "DEBUG:%s:%s=%lg\n", __FUNCTION__, atomicName, molValue);
			electronMol += assignStableIstopesOnElectrode(a_electrodePtr, atomP->atomicNumber, molValue);
			*a_sumOfMolPtr += molValue;
			//heatCapacity += (molValue * atomP->heatCapacity);
		}else{
			fprintf(stderr, "ERROR:%s:BAD atomicName %s of name:%s\n", __FUNCTION__, atomicName, a_name);
			exit(1);
		}
	}
	//fprintf(stderr, "DEBUG:%s:END} electronMol:%lg sumOfMol:%lg\n", __FUNCTION__, electronMol, *a_sumOfMolPtr);
	return electronMol;
}
extern void initElectrodes(struct electrode * a_electrodePtr, char * a_ElectrodeAtomicMolsStr, double a_detectLimitRateForIsotope, const char * a_electrodeName)
{
	//fprintf(stderr, "DEBUG:%s:{BEGIN a_ElectrodeAtomicMolsStr:%s a_detectLimitRateForIsotope:%lg a_electrodeName:%s\n", __FUNCTION__, a_ElectrodeAtomicMolsStr, a_detectLimitRateForIsotope, a_electrodeName);
	double electronMol = 0.0;
	double sumOfMol = 0.0, tempOfMol;
	//double heatCapacity = 0.0;
	struct atomNodeConst * electronPtr;
	int i;
	char nodeName[128];
	char name[NAME_MAX];
#define VAL_MAX 48
	char val[VAL_MAX];
	int findComma, findEqual, endText, nameLen, valLen;

	initHashTable(&a_electrodePtr->atomHashTable, 
		calcHashSeedOfAtomKey,
		allocCopyAtomKey,
		free,
		compareAtomKey,
		allocCopyAtomValue,
		free,
		foundActionForAtomValue,
		strcat(strcpy(nodeName, a_electrodeName), " Electrode"),
		30);
	increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, 0.0);
	
	findComma = 1;
	findEqual = 0;
	endText = 0;
	nameLen = valLen = 0;	
	for(i = 0; a_ElectrodeAtomicMolsStr[i] != 0; ++i){
		//fprintf(stderr, "DEBUG:%s:[%d]:%c findComma:%d findEqual:%d endText:%d nameLen:%d valLen:%d\n", __FUNCTION__, i, a_ElectrodeAtomicMolsStr[i], findComma, findEqual, endText, nameLen, valLen);
		if(a_ElectrodeAtomicMolsStr[i] == '='){
			if(findComma == 1){
				findComma = 0;
				findEqual = 1;
				endText = 0;
				name[nameLen] = 0;
				valLen = 0;
			}else{
				fprintf(stderr, "ERROR:%s:BAD EQUAL in %s\n", __FUNCTION__, a_ElectrodeAtomicMolsStr);
				exit(1);
			}
		}else if(a_ElectrodeAtomicMolsStr[i] == ','){
			if(findEqual == 1){
				findComma = 1;
				findEqual = 0;
				endText = 0;
				electronMol += registAtomOnElectrodes(a_electrodePtr, name, val, valLen, &tempOfMol);
				sumOfMol += tempOfMol;
				nameLen = 0;
			}else{
				fprintf(stderr, "ERROR:%s:BAD COMMA %s\n", __FUNCTION__, a_ElectrodeAtomicMolsStr);
				exit(1);
			}
		}else if(a_ElectrodeAtomicMolsStr[i] == ' ' || a_ElectrodeAtomicMolsStr[i] == '\t'){
			if(findComma == 1){
				if(nameLen > 0){
					endText = 1;
				}
			}
			if(findEqual == 1){
				if(valLen > 0){
					endText = 1;
				}
			}
		}else{
			if(endText == 0){
				if(findComma == 1){
					name[nameLen] = a_ElectrodeAtomicMolsStr[i];
					++nameLen;
				}
				if(findEqual == 1){
					val[valLen] = a_ElectrodeAtomicMolsStr[i];
					++valLen;
				}				
			}else{
				fprintf(stderr, "ERROR:%s:BAD FORMAT in %s\n", __FUNCTION__, a_ElectrodeAtomicMolsStr);
				exit(1);
			}
		}
	}//end of for
	if(findEqual == 1){
		electronMol += registAtomOnElectrodes(a_electrodePtr, name, val, valLen, &tempOfMol);
		sumOfMol += tempOfMol;
	}
	electronPtr = findElectronInElectrode(a_electrodePtr);
	if(electronPtr){
		electronPtr->vPtr->molIni = electronMol;
	}else{
		fprintf(stderr, "ERROR:%s:electronPtr == NULL\n", __FUNCTION__);
		exit(1);
	}
	increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, 0.0);
	a_electrodePtr->detectLimitMolForIsotope = sumOfMol * a_detectLimitRateForIsotope;
	
	initHashTable(&a_electrodePtr->massDefectHashTable, 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getScatName(SCAT_BIG_MASS_DEFECT_NOW)),
		100);
	initSumMeVMol(&a_electrodePtr->massDefectAll);
	initSumMeVMol(&a_electrodePtr->byNeutrinoAll);
	for(i = 0; i < SIZE_OF_SCATTERD; ++i){
		initSumMeVMol(a_electrodePtr->scattered + i);
	}

	//fprintf(stderr, "DEBUG:%s:}END \n", __FUNCTION__);
	//return heatCapacity;
}
//---------------------------------------------------------------------

extern void registOutput(struct electrode * a_electrodePtr, int a_scatType, double a_MeV, double a_mol)
{
	//fprintf(stderr, "DEBUG:%s:{BEGIN %lp, %s a_MeV:%lg * a_mol:%lg = %lg\n\n", __FUNCTION__, a_electrodePtr, getMassChangeName(a_scatType), a_MeV, a_mol, a_MeV * a_mol);
	if(a_mol > 0.0){
		//Check data!
		if(a_scatType == SCAT_BIG_MASS_DEFECT_NOW || a_scatType == SCAT_BIG_MASS_DEFECT_ALL){
			if(a_MeV < e_collideElectronMiniMeV){
				a_scatType = SCAT_SMALL_MASS_DEFECT;
			}else{
				a_scatType = SCAT_BIG_MASS_DEFECT_NOW;
			}
		}
		if(a_scatType == SCAT_BIG_MASS_DEFECT_NOW){
			struct range key;
			struct massDefect value;
			struct objectNodeConst nodeConst;
			key = calcPrecisionRange(a_MeV);
			value.MeV = a_MeV;
			value.mol = a_mol;
			nodeConst.keyPtr = &key;
			nodeConst.valuePtr = &value;
			insertObjectInHashTable(&a_electrodePtr->massDefectHashTable, &nodeConst);
			addSumMeVMol(&a_electrodePtr->massDefectAll, a_MeV, a_mol);
		}else if(a_scatType == SCAT_LOST_BY_NEUTRINO){
			addSumMeVMol(&a_electrodePtr->byNeutrinoAll, a_MeV, a_mol);
		}else{
			addSumMeVMol(a_electrodePtr->scattered + a_scatType, a_MeV, a_mol);
		}
	}else{
		//fprintf(stderr, "DEBUG:%s:a_electrodePtr:%lp a_scatType:%d a_MeV:%lg a_mol:%lg <= 0.0\n", __FUNCTION__, a_electrodePtr, a_scatType, a_MeV, a_mol);
	}
	//fprintf(stderr, "DEBUG:%s:}END\n", __FUNCTION__);
}

const double e_rateNeutrinoMev[]   = {0.25, 0.5, 0.75};//The rate of the some part of beta decay energy that carried out by the neutrino.
const double e_probabilityForMol[] = {0.25, 0.5, 0.25};//The sum of probabilityForMol must be 1.0. 

/*
extern void checkNeutrinoTableForDBUG()
{
	int i;
	double sum = 0.0;
	if(sizeof(e_rateNeutrinoMev) != sizeof(e_probabilityForMol)){
		fprintf(stderr, "ERROR:%s:sizeof(e_rateNeutrinoMev):%d != sizeof(e_probabilityForMol):%d\n", __FUNCTION__, sizeof(e_rateNeutrinoMev), sizeof(e_probabilityForMol));
		exit(1);			
	}
	for(i = 0; i < (sizeof(e_rateNeutrinoMev)/sizeof(double)); ++i){
		if(e_rateNeutrinoMev[i] < 0.0 || 1.0 < e_rateNeutrinoMev[i]){
			fprintf(stderr, "ERROR:%s:e_rateNeutrinoMev[%d]:%lg != 1.0\n", __FUNCTION__, i, e_rateNeutrinoMev[i]);
			exit(1);
		}
	}
	for(i = 0; i < (sizeof(e_probabilityForMol)/sizeof(double)); ++i){
		sum += e_probabilityForMol[i];
	}
	if(sum < 0.9999999999 || 1.0000000001 < sum){
		fprintf(stderr, "ERROR:%s:sum of e_probabilityForMol:%lg != 1.0\n", __FUNCTION__, sum);
		exit(1);
	}
}
*/
extern void registOutputOfBetaMinus(struct electrode * a_electrodePtr, double a_massDefectMeV, double a_mol)
{
	//The beta minus emmits an electron and an anti-neutrino. The Nuclear Reaction is : n -> p + e- + ~νe
	//This function simulate the energy(mass defect) of between a neutrino and an electron in the decay of beta minus
	//The neutrino will carry out the some part of beta decay energy.
	//The reamin part of energy is holded by the electron. Immediately the electron emmits its energy.
	int i;
	//double debugGoal = a_massDefectMeV * a_mol;//DEBUG
	//double debugSum = 0.0;//DEBUG
	//fprintf(stderr, "DEBUG:%s:{BEGIN %lp, a_massDefectMeV:%lg * a_mol:%lg = debugGoal:%lg\n\n", __FUNCTION__, a_electrodePtr, a_massDefectMeV, a_mol, debugGoal);
	for(i = 0; i < (sizeof(e_rateNeutrinoMev)/sizeof(double)); ++i){
		double energyCarryedByNeutrino = a_massDefectMeV * e_rateNeutrinoMev[i];
		double remainMassDefect = a_massDefectMeV - energyCarryedByNeutrino;
		double molOfNeutrino = a_mol * e_probabilityForMol[i];
		double molOfElectoron = molOfNeutrino;
		registOutput(a_electrodePtr, SCAT_LOST_BY_NEUTRINO, energyCarryedByNeutrino, molOfNeutrino);
		//debugSum += energyCarryedByNeutrino * molOfNeutrino;//DEBUG
		registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, remainMassDefect, molOfElectoron);
		//debugSum += remainMassDefect * molOfElectoron;//DEBUG
	}
	/*
	if(debugGoal > 0.0){
		double rate = (debugGoal - debugSum) / debugGoal;
		if(rate < 0.0){
			rate = -rate;
		}
		if(rate > 0.00001){
			fprintf(stderr, "ERROR:%s:debugGoal:%lg != debugSum:%lg\n", __FUNCTION__, debugGoal, debugSum);
		}
	}
	*/
	//fprintf(stderr, "DEBUG:%s:}END debugSum:%lg\n", __FUNCTION__, debugSum);
}
extern void registOutputOfElectronCapture(struct electrode * a_electrodePtr, double a_massDefectMeV, double a_mol)
{
	//The electron capture do not emmit a positron(e+), only emmits a neutorino.
	//The electron capture is also the collide of a proton and an electoron.
	//The Nuclear Reaction is : p + e- -> n + νe
	//In the electron capture mode, the some part of energy of mass defect is carried out by the neutrino.
	//The reamin part of energy is holded by the daughter neucus that has a neutron. Immediately the daughter neucus emmits its energy.
	int i;
	//double debugGoal = a_massDefectMeV * a_mol;//DEBUG
	//double debugSum = 0.0;//DEBUG
	//fprintf(stderr, "DEBUG:%s:{BEGIN %lp, a_massDefectMeV:%lg * a_mol:%lg = debugGoal:%lg\n\n", __FUNCTION__, a_electrodePtr, a_massDefectMeV, a_mol, debugGoal);
	for(i = 0; i < (sizeof(e_rateNeutrinoMev)/sizeof(double)); ++i){
		double energyCarryedByNeutrino = a_massDefectMeV * e_rateNeutrinoMev[i];
		double remainMassDefect = a_massDefectMeV - energyCarryedByNeutrino;
		double molOfNeutrino = a_mol * e_probabilityForMol[i];
		double molOfDaughter = molOfNeutrino;
		registOutput(a_electrodePtr, SCAT_LOST_BY_NEUTRINO, energyCarryedByNeutrino, molOfNeutrino);
		//debugSum += energyCarryedByNeutrino * molOfNeutrino;//DEBUG
		registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, remainMassDefect, molOfDaughter);	
		//debugSum += remainMassDefect * molOfDaughter;//DEBUG
	}
	/*
	if(debugGoal > 0.0){
		double rate = (debugGoal - debugSum) / debugGoal;
		if(rate < 0.0){
			rate = -rate;
		}
		if(rate > 0.00001){
			fprintf(stderr, "ERROR:%s:debugGoal:%lg != debugSum:%lg\n", __FUNCTION__, debugGoal, debugSum);
		}
	}
	*/
	//fprintf(stderr, "DEBUG:%s:}END debugSum:%lg\n", __FUNCTION__, debugSum);
}
extern void registOutputOfBetaPlusWithPositron(struct electrode * a_electrodePtr, double a_massDefectMeV, double a_mol)
{
	//The beta plus is as same as the positron(e+) emission
	//The Nuclear Reaction is : p -> n + e+ + νe
	//This function simulate the energy(mass defect) of between a neutrino and a positron in the decay of beta plus.
	//A positron is an anni-electron.
	//The neutrino will carry out the some part of beta decay energy.
	//The reamin part of energy is holded by the positron. 
	//The positron will immediately collide an electron in the neighbor, its called annihilation.
	//The annihilations of a positron and an electron with their kinetic ennergy will expose two gamma.
	//The Nuclear Reaction is : e+ + e- -> gamma + gamma
	int i;
	double debugGoal = (a_massDefectMeV + 2.0 * e_massElectronMeV) * a_mol;//DEBUG
	double debugSum = 0.0;//DEBUG
	//fprintf(stderr, "DEBUG:%s:{BEGIN %p, (a_massDefectMeV:%lg  + 2.0 * e_massElectronMeV:%lg) * a_mol:%lg = debugGoal:%lg\n\n", __FUNCTION__, a_electrodePtr, a_massDefectMeV, e_massElectronMeV, a_mol, debugGoal);
	for(i = 0; i < (sizeof(e_rateNeutrinoMev)/sizeof(double)); ++i){
		double energyCarryedByNeutrino = a_massDefectMeV * e_rateNeutrinoMev[i];
		double remainMassDefect = a_massDefectMeV - energyCarryedByNeutrino;
		double energyOfGamma = e_massElectronMeV + remainMassDefect * 0.5;
		double molOfNeutrino = a_mol * e_probabilityForMol[i];
		double molOfGamma = molOfNeutrino * 2.0;
		registOutput(a_electrodePtr, SCAT_LOST_BY_NEUTRINO, energyCarryedByNeutrino, molOfNeutrino);
		debugSum += energyCarryedByNeutrino * molOfNeutrino;//DEBUG
		registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, energyOfGamma, molOfGamma);
		debugSum += energyOfGamma * molOfGamma;//DEBUG
	}
	
	if(debugGoal > 0.0){
		double rate = (debugGoal - debugSum) / debugGoal;
		if(rate < 0.0){
			rate = -rate;
		}
		if(rate > 0.00001){
			fprintf(stderr, "ERROR:%s:debugGoal:%lg != debugSum:%lg\n", __FUNCTION__, debugGoal, debugSum);
			exit(1);
		}
	}
	
	//fprintf(stderr, "DEBUG:%s:}END debugSum:%lg\n", __FUNCTION__, debugSum);
}
//---------------------------------------------------------------------
extern int iterateSumupMassDefect(void * a_total, struct objectNodeConst * a_nodePtr)
{
	struct sumMeVMol * sumPtr = (struct sumMeVMol *)a_total;
	//struct range * keyPtr = (struct range *)a_nodePtr->key;
	struct massDefect * valuePtr = (struct massDefect *)a_nodePtr->valuePtr;
	addSumMeVMol(sumPtr, valuePtr->MeV, valuePtr->mol);
	return KEEP_NODE;
}

extern void sumupBigMassDefect(struct electrode * a_electrodePtr, struct sumMeVMol * a_massDefectNowPtr)
{
	initSumMeVMol(a_massDefectNowPtr);
	iterateInHashTable(&a_electrodePtr->massDefectHashTable, a_massDefectNowPtr, iterateSumupMassDefect);
}
extern long double printMassDefectHeat(FILE * a_fp, const char * a_titlePtr, const struct sumMeVMol * a_massDefectNowPtr, const struct sumMeVMol * a_massDefectAllPtr, const struct sumMeVMol * a_byNeutrinoAllPtr, struct sumMeVMol scatteredAry[SIZE_OF_SCATTERD])
{
	long double outputHeat;
	fprintf(a_fp, "%s MASS DEFECT\n(gamma ray >= %lg [MeV], beta energy = %lg [MeV])\n", a_titlePtr, e_collideElectronMiniMeV, e_betaEnergyMeV);
	printSumMeVMol(a_fp, " ", SCAT_BIG_MASS_DEFECT_NOW, a_massDefectNowPtr, 0.0, "");
	printSumMeVMol(a_fp, " ", SCAT_BIG_MASS_DEFECT_ALL, a_massDefectAllPtr, 0.0, "");
	outputHeat = printScattered(a_fp, a_titlePtr, scatteredAry);
	fprintf(a_fp, "%s LOST HEAT =\n", a_titlePtr);
	printSumMeVMol(a_fp, " ", SCAT_LOST_BY_NEUTRINO, a_byNeutrinoAllPtr, outputHeat, "[%(/Output heat)]");
	return outputHeat;
}
extern double printElectrode(FILE * a_fp, struct electrode * a_electrodePtr, struct sumMeVMol * a_massDefectNowPtr, double * a_sumOfMassUMolIniPtr, double * a_sumOfMassUMolAddPtr, double * a_sumOfMassUMolSubPtr, double * a_heatCapacityPtr)
{
	//long double grandTotal;
	double sumOfMassUMol;

	fprintf(a_fp, "%s detectLimitMolForIsotope %lg\n", a_electrodePtr->atomHashTable.tableName, a_electrodePtr->detectLimitMolForIsotope);
	sumOfMassUMol = printAtomList(a_fp, &a_electrodePtr->atomHashTable, a_sumOfMassUMolIniPtr, a_sumOfMassUMolAddPtr, a_sumOfMassUMolSubPtr, a_heatCapacityPtr);
	sumupBigMassDefect(a_electrodePtr, a_massDefectNowPtr);
	/*grandTotal = */printMassDefectHeat(a_fp, a_electrodePtr->atomHashTable.tableName, a_massDefectNowPtr, &a_electrodePtr->massDefectAll, &a_electrodePtr->byNeutrinoAll, a_electrodePtr->scattered);
	fputs("\n", a_fp);
	return sumOfMassUMol;
}
//---------------------------------------------------------------------

#define COLLIDE_ELECTRON  (MAX_DECAY_MODE + 1) //by high volatage
#define COLLIDE_ELECTRON_C  (MAX_DECAY_MODE + 2) //by the compton Effect
#define COLLIDE_PROTON (MAX_DECAY_MODE + 3)
#define COLLIDE_DEUTERIUM (MAX_DECAY_MODE + 4)
#define COLLIDE_TRITIUM (MAX_DECAY_MODE + 5)
#define COLLIDE_NEUTRON (MAX_DECAY_MODE + 6)
extern const char * getCollideName(int a_collideType)
{
	const char * collideNamePtr = "unkown";
	switch(a_collideType){
	case COLLIDE_ELECTRON:	 collideNamePtr = "COLLIDE_ELECTRON_BY_HIGH_VOLTAGE"; break;
	case COLLIDE_ELECTRON_C: collideNamePtr = "COLLIDE_ELECTRON_BY_COMPTON_EFFECT"; break;
	case COLLIDE_PROTON:         collideNamePtr = "COLLIDE_PROTON_BY_HIGH_VOLTAGE"; break;
	case COLLIDE_DEUTERIUM:      collideNamePtr = "COLLIDE_DEUTERIUM_BY_HIGH_VOLTAGE"; break;
	case COLLIDE_TRITIUM:        collideNamePtr = "COLLIDE_TRITIUM_BY_HIGH_VOLTAGE";  break;
	case COLLIDE_NEUTRON:        collideNamePtr = "COLLIDE_NEUTRON";  break;
	default:
		if(a_collideType <= MAX_DECAY_MODE){
			collideNamePtr = getDecayModeText(a_collideType);
		}
		break;
	}
	return collideNamePtr;
}
#define REACTION_LEN 512
struct collide {
	//{ set by collideBulletToElectrode
	struct electrode * electrodePtr;
	int collideType;
	const struct isotopeProperty * bulletPropertyPtr;
	double appliedVoltageMeV;
	double arriveBulletMol;
	double totalCollideCrossSection;//incremented by iterateCrossSection
	double totalTargetMol;//incremented by iterateCrossSection
	int printNuclearReaction;
	//}

	//{ set by useElectronOrNeutron
	int useElectron;
	int useNeutron;
	//}

	//{ set by collideBulletToElectrode
	double collideBulletMol;
	double imperfectCollideMol;
	double collidedBulletMol;// incremented by iterateCollide
	double remainBulletMol;
	//}

	//{ set by iterateCollide and iterateCrossSection
	const struct isotopeProperty * targetIsotopePropertyPtr;
	struct atomValue * targetAtomValuePtr;
	double targetAtomMol;
	//}
	
	double electroPotentialMeV;//set by calcElectroPotentialMeVByCollide
	
	//{ set by getCollidedNewIsotopePropertyPtr
	int newAtomicNumber;
	int newMassNumber;
	const struct isotopeProperty * newIsotopePropertyPtr;
	double minAppliedVoltageMeV;
	double massDefectMeV;
	double minMassDefectMeV;
	char * plusAlpha;
	char nuclearReactionStr[REACTION_LEN + 1];
	int isElectronCapture;
	//}

	double newIsotopeMol;// set by collideParticleAtom
};
extern void useElectronOrNeutron(struct collide * a_pX)
{
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_C){
		a_pX->useElectron = 0;
		a_pX->useNeutron = 0;
	}else if(a_pX->collideType == COLLIDE_PROTON){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 1;
	}else if(a_pX->collideType == COLLIDE_DEUTERIUM){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 1;
	}else if(a_pX->collideType == COLLIDE_TRITIUM){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 1;
	}else if(a_pX->collideType == COLLIDE_NEUTRON){
		a_pX->useElectron = 0;
		a_pX->useNeutron = 0;
	}
}
//---------------------------------------------------------------------
#define ELECTRODE_NEGATIVE 0
#define ELECTRODE_POSITIVE 1
extern const char * getElectrodeName(int a_electrodeType)
{
	const char * electrodeTypeStr = "";
	switch(a_electrodeType){
		case ELECTRODE_NEGATIVE: electrodeTypeStr = "NEGATIVE ELECTRODE"; break;
		case ELECTRODE_POSITIVE: electrodeTypeStr = "POSITIVE ELECTRODE"; break;
	}
	return electrodeTypeStr;	
}
extern char getElectrodeChar(int a_electrodeType)
{
	char electrodeChar;
	switch(a_electrodeType){
		case ELECTRODE_NEGATIVE: electrodeChar = 'N'; break;
		case ELECTRODE_POSITIVE: electrodeChar = 'P'; break;
	}
	return electrodeChar;	
}
#define REACTION_DETECT 1
#define REACTION_CANT_DETECT 2 //The mol of isotpoe is under the detect limit.
#define REACTION_ENDOTHERMIC 3 
#define REACTION_COULOMB_BARRIER 4
#define REACTION_ERROR 5 //The ERROR may be caused by an unregisted atom.
extern const char * getNuclearReactionName(int a_reactionType)
{
	const char * reactionTypeStr = "";
	switch(a_reactionType){
		case REACTION_DETECT: reactionTypeStr = "DETECT exothermic"; break;
		case REACTION_CANT_DETECT: reactionTypeStr = "can't DETECT exothermic by little mol"; break;
		case REACTION_ENDOTHERMIC: reactionTypeStr = "Endothermic"; break;
		case REACTION_COULOMB_BARRIER: reactionTypeStr = "Coulomb Barrier"; break;
		case REACTION_ERROR: reactionTypeStr = "ERROR";  break;
	}
	return reactionTypeStr;
}
extern char getNuclearReactionChar(int a_reactionType)
{
	char reactionTypeChar;
	switch(a_reactionType){
		case REACTION_DETECT: reactionTypeChar = '+'; break;
		case REACTION_CANT_DETECT: reactionTypeChar = '0'; break;
		case REACTION_ENDOTHERMIC: reactionTypeChar = '-'; break;
		case REACTION_COULOMB_BARRIER: reactionTypeChar = '!'; break;
		case REACTION_ERROR: reactionTypeChar = 'E';  break;
	}
	return reactionTypeChar;
}

struct nuclearReactionKey {
	int electrodeType;
	int reactionType;
	int decayOrCollideType;
	char * formPtr;
};
extern unsigned int calcCheckSum(char * a_formPtr)
{
	unsigned int i, checkSum = 0;
	for(i = 0; a_formPtr[i]; ++i){
		checkSum = ((checkSum << 1) | (checkSum >> (sizeof(unsigned int) * 8 - 1))) + (unsigned int)a_formPtr[i];
	}
	return checkSum;
}
extern unsigned int calcHashSeedOfNuclearReactionKey(const void * a_keyPtr)
{
	const struct nuclearReactionKey * ptr = (const struct nuclearReactionKey *)a_keyPtr;
	unsigned int hashSeed = (unsigned int)(ptr->electrodeType << 24) + (unsigned int)(ptr->reactionType << 20) + (unsigned int)(ptr->decayOrCollideType << 10) + calcCheckSum(ptr->formPtr);
	//"massDefectMeV" is included in text of "formPtr"
	return hashSeed;
}
extern void * allocCopyNuclearReactionKey(const void * a_keyPtr)
{
	const struct nuclearReactionKey * ptr = (const struct nuclearReactionKey *)a_keyPtr;
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_keyPtr);
	struct nuclearReactionKey * ret = allocCopy(a_keyPtr, sizeof(struct nuclearReactionKey), __FUNCTION__);
	ret->formPtr = allocStrcpy(ptr->formPtr);
	//fprintf(stderr, "DEBUG:%s:end}ret:%lp\n", __FUNCTION__, ret);
	return ret;
}
extern void freeCopyNuclearReactionKey(void * a_keyPtr)
{
	const struct nuclearReactionKey * ptr = (struct nuclearReactionKey *)a_keyPtr;
	free(ptr->formPtr);
	free(a_keyPtr);
}
#define VERSION_OF_serializeNuclearReactionKey 1
extern void * allocFreadNuclearReactionKey(FILE * a_fp)
{
	struct nuclearReactionKey * ret;
	FREAD_VERSION_CHECK(VERSION_OF_serializeNuclearReactionKey, a_fp)
	ret = clearAlloc(sizeof(struct nuclearReactionKey), __FUNCTION__);
	if(ret){
		FREAD_MEMBER_FREE(ret, electrodeType, a_fp)
		FREAD_MEMBER_FREE(ret, reactionType, a_fp)
		FREAD_MEMBER_FREE(ret, decayOrCollideType, a_fp)
		if(!(ret->formPtr = allocFreadStr(a_fp))){free(ret); return NULL;}
	}
	return ret;
}
extern int fwriteNuclearReactionKey(const void * a_keyPtr, FILE * a_fp)
{
	const struct nuclearReactionKey * ptr = (const struct nuclearReactionKey *)a_keyPtr;
	FWRITE_VERSION(VERSION_OF_serializeNuclearReactionKey, a_fp)
	FWRITE_VALUE(ptr->electrodeType, a_fp)
	FWRITE_VALUE(ptr->reactionType, a_fp)
	FWRITE_VALUE(ptr->decayOrCollideType, a_fp)
	if(!fwriteStr(ptr->formPtr, a_fp)){return 0;}
	return 1;
}
extern int compareNuclearReactionByFormPtr(const void *a, const void *b )
{
	struct nuclearReactionKey * x = (struct nuclearReactionKey *)a;
	struct nuclearReactionKey * y = (struct nuclearReactionKey *)b;
	int ret = strcmp(x->formPtr, y->formPtr);
	//if(e_degug){fprintf(stderr, "DEBUG:%s:ret:%d\nx->form:%s\ny->form:%s\n", __FUNCTION__, ret, x->formPtr, y->formPtr);}
	return ret;
}
extern int compareNuclearReactionByDecayOrCollideType(const void *a, const void *b )
{
	struct nuclearReactionKey * x = (struct nuclearReactionKey *)a;
	struct nuclearReactionKey * y = (struct nuclearReactionKey *)b;
	int ret;
	if(x->decayOrCollideType < y->decayOrCollideType){
		ret = -1;
	}else if(x->decayOrCollideType == y->decayOrCollideType){
		ret = 0;
	}else{
		ret = 1;
	}
	//if(e_degug){fprintf(stderr, "DEBUG:%s:ret:%d\nx->decayOrCollideType:%s\ny->decayOrCollideType:%s\n", __FUNCTION__, ret, getCollideName(x->decayOrCollideType), getCollideName(y->decayOrCollideType));}
	return ret;
}
extern int compareNuclearReactionByReactionType(const void *a, const void *b )
{
	struct nuclearReactionKey * x = (struct nuclearReactionKey *)a;
	struct nuclearReactionKey * y = (struct nuclearReactionKey *)b;
	int ret;	
	if(x->reactionType < y->reactionType){
		ret = -1;
	}else if(x->reactionType == y->reactionType){
		ret = 0;
	}else{
		ret = 1;
	}
	//if(e_degug){fprintf(stderr, "DEBUG:%s:ret:%d\nx->reactionType:%s\ny->reactionType:%s\n", __FUNCTION__, ret, getNuclearReactionName(x->reactionType), getNuclearReactionName(y->reactionType));}	
	return ret;
}
extern int compareNuclearReactionByElectrodeType(const void *a, const void *b )
{
	struct nuclearReactionKey * x = (struct nuclearReactionKey *)a;
	struct nuclearReactionKey * y = (struct nuclearReactionKey *)b;
	int ret;	
	if(x->electrodeType < y->electrodeType){
		ret = -1;
	}else if(x->electrodeType == y->electrodeType){
		ret = 0;
	}else{
		ret = 1;
	}
	//if(e_degug){fprintf(stderr, "DEBUG:%s:ret:%d\nx->reactionType:%s\ny->reactionType:%s\n", __FUNCTION__, ret, getNuclearReactionName(x->reactionType), getNuclearReactionName(y->reactionType));}	
	return ret;
}
extern int compareNuclearReactionKey(const void *a, const void *b )
{
	int ret;
	ret = compareNuclearReactionByElectrodeType(a, b);
	if(ret == 0){
		ret = compareNuclearReactionByReactionType(a, b);
		if(ret == 0){
			ret = compareNuclearReactionByDecayOrCollideType(a, b);
			if(ret == 0){
				ret = compareNuclearReactionByFormPtr(a, b);
			}
		}
	}
	//if(e_degug){fprintf(stderr, "DEBUG:%s:ret:%d\n", __FUNCTION__, ret);}	return ret;
	return ret;
}
	
struct nuclearReactionValue {
	double maxMassDefectMeV;
	double minMassDefectMeV;
	double maxIsotopeMol;
	double minIsotopeMol;
	double maxAppliedVoltageMeV;
	double minAppliedVoltageMeV;
	double maxElectroPotentialMeV;
	double minElectroPotentialMeV;
	long long sumupCnt;
	long double sumupOfMulMeVMol;
	double maxTimeSec;
	double minTimeSec;
};
extern void * allocCopyNuclearReactionValue(const void * a_valuePtr)
{
	//const struct nuclearReactionValue * ptr = (const struct nuclearReactionValue *)a_valuePtr;
	//fprintf(stderr, "DEBUG:%s:begin{a_keyPtr:%lp\n", __FUNCTION__, a_valuePtr);
	struct nuclearReactionValue * ret = allocCopy(a_valuePtr, sizeof(struct nuclearReactionValue), __FUNCTION__);
	//fprintf(stderr, "DEBUG:%s:end}ret:%lp\n", __FUNCTION__, ret);
	return ret;
}
#define VERSION_OF_serializeNuclearReactionValue 1
extern void * allocFreadNuclearReactionValue(FILE * a_fp)
{
	struct nuclearReactionValue * ret;
	FREAD_VERSION_CHECK(VERSION_OF_serializeNuclearReactionValue, a_fp)
	ret = clearAlloc(sizeof(struct nuclearReactionValue), __FUNCTION__);
	if(ret){
		FREAD_MEMBER_FREE(ret, maxMassDefectMeV, a_fp)
		FREAD_MEMBER_FREE(ret, minMassDefectMeV, a_fp)
		FREAD_MEMBER_FREE(ret, maxIsotopeMol, a_fp)
		FREAD_MEMBER_FREE(ret, minIsotopeMol, a_fp)
		FREAD_MEMBER_FREE(ret, maxAppliedVoltageMeV, a_fp)
		FREAD_MEMBER_FREE(ret, minAppliedVoltageMeV, a_fp)
		FREAD_MEMBER_FREE(ret, maxElectroPotentialMeV, a_fp)
		FREAD_MEMBER_FREE(ret, minElectroPotentialMeV, a_fp)
		FREAD_MEMBER_FREE(ret, sumupCnt, a_fp)
		FREAD_MEMBER_FREE(ret, sumupOfMulMeVMol, a_fp)
		FREAD_MEMBER_FREE(ret, maxTimeSec, a_fp)
		FREAD_MEMBER_FREE(ret, minTimeSec, a_fp)
	}
	return ret;
}
extern int fwriteNuclearReactionValue(const void * a_valuePtr, FILE * a_fp)
{
	const struct nuclearReactionValue * ptr = (const struct nuclearReactionValue *)a_valuePtr;
	FWRITE_VERSION(VERSION_OF_serializeNuclearReactionValue, a_fp)
	FWRITE_VALUE(ptr->maxMassDefectMeV, a_fp)
	FWRITE_VALUE(ptr->minMassDefectMeV, a_fp)
	FWRITE_VALUE(ptr->maxIsotopeMol, a_fp)
	FWRITE_VALUE(ptr->minIsotopeMol, a_fp)
	FWRITE_VALUE(ptr->maxAppliedVoltageMeV, a_fp)
	FWRITE_VALUE(ptr->minAppliedVoltageMeV, a_fp)
	FWRITE_VALUE(ptr->maxElectroPotentialMeV, a_fp)
	FWRITE_VALUE(ptr->minElectroPotentialMeV, a_fp)
	FWRITE_VALUE(ptr->sumupCnt, a_fp)
	FWRITE_VALUE(ptr->sumupOfMulMeVMol, a_fp)
	FWRITE_VALUE(ptr->maxTimeSec, a_fp)
	FWRITE_VALUE(ptr->minTimeSec, a_fp)
	return 1;
}

extern void foundActionForNuclearReactionValue(struct objectNodeConst * a_destPtr, const struct objectNodeConst * a_srcPtr)
{
	//const struct nuclearReactionKey * kPtr = (const struct nuclearReactionKey *)a_destPtr->keyPtr;
	struct nuclearReactionValue * vPtr1 = (struct nuclearReactionValue *)a_destPtr->valuePtr;
	const struct nuclearReactionValue * vPtr2 = (const struct nuclearReactionValue *)a_srcPtr->valuePtr;
	
	if(vPtr1->maxMassDefectMeV < vPtr2->maxMassDefectMeV){
		vPtr1->maxMassDefectMeV = vPtr2->maxMassDefectMeV;
	}
	if(vPtr1->minMassDefectMeV > vPtr2->minMassDefectMeV){
		vPtr1->minMassDefectMeV = vPtr2->minMassDefectMeV;
	}
	if(vPtr1->maxIsotopeMol < vPtr2->maxIsotopeMol){
		vPtr1->maxIsotopeMol = vPtr2->maxIsotopeMol;
	}
	if(vPtr1->minIsotopeMol > vPtr2->minIsotopeMol){
		vPtr1->minIsotopeMol = vPtr2->minIsotopeMol;
	}
	if(vPtr1->maxAppliedVoltageMeV < vPtr2->maxAppliedVoltageMeV){
		vPtr1->maxAppliedVoltageMeV = vPtr2->maxAppliedVoltageMeV;
	}
	if(vPtr1->minAppliedVoltageMeV > vPtr2->minAppliedVoltageMeV){
		vPtr1->minAppliedVoltageMeV = vPtr2->minAppliedVoltageMeV;
	}
	if(vPtr1->maxElectroPotentialMeV < vPtr2->maxElectroPotentialMeV){
		vPtr1->maxElectroPotentialMeV = vPtr2->maxElectroPotentialMeV;
	}
	if(vPtr1->minElectroPotentialMeV > vPtr2->minElectroPotentialMeV){
		vPtr1->minElectroPotentialMeV = vPtr2->minElectroPotentialMeV;
	}
	vPtr1->sumupCnt += vPtr2->sumupCnt;
	vPtr1->sumupOfMulMeVMol += vPtr2->sumupOfMulMeVMol;
	if(vPtr1->maxTimeSec < vPtr2->maxTimeSec){
		vPtr1->maxTimeSec = vPtr2->maxTimeSec;
	}
	if(vPtr1->minTimeSec > vPtr2->minTimeSec){
		vPtr1->minTimeSec = vPtr2->minTimeSec;
	}
}
//---------------------------------------------------------------------
#define FNAME_LEN 255
struct runningCondition {
	double startTime;
	double termTime;
	double intervalTime;
	int loopBeg;
	int loopEnd;
	int logStep;
	int logNecleusReactions;
	int dataStep;
	char inputDataFname[FNAME_LEN + 1];
	char outputDataFname[FNAME_LEN + 1];
	char logFname[FNAME_LEN + 1];
	char paramFname[FNAME_LEN + 1];
	long double inputTotalEnergy;
	double simulatedLastTime;
};
struct runningCondition e_rc;




//---------------------------------------------------------------------
struct hashTable e_nuclearReaction;

extern void registNuclearReaction(struct electrode * a_electrodePtr, int a_reactionType, int a_decayOrCollideType, char * a_formPtr, double a_massDefectMeV, double a_newIsotopeMol, double a_appliedVoltageMeV, double a_electroPotentialMeV)
{
	struct nuclearReactionKey Key;
	struct nuclearReactionValue Value;
	struct objectNodeConst nodeConst;
	
	Key.electrodeType = (a_electrodePtr == &e_negativeElectrode) ? ELECTRODE_NEGATIVE : ELECTRODE_POSITIVE;
	Key.reactionType = a_reactionType;
	Key.decayOrCollideType = a_decayOrCollideType;
	Key.formPtr = a_formPtr;
	Value.maxMassDefectMeV = a_massDefectMeV;
	Value.minMassDefectMeV = a_massDefectMeV;
	Value.maxIsotopeMol = a_newIsotopeMol;
	Value.minIsotopeMol = a_newIsotopeMol;
	Value.maxAppliedVoltageMeV = a_appliedVoltageMeV;
	Value.minAppliedVoltageMeV = a_appliedVoltageMeV;
	Value.maxElectroPotentialMeV = a_electroPotentialMeV;
	Value.minElectroPotentialMeV = a_electroPotentialMeV;
	Value.sumupCnt = 1;
	Value.sumupOfMulMeVMol = a_massDefectMeV * a_newIsotopeMol;
	Value.maxTimeSec = e_rc.simulatedLastTime;
	Value.minTimeSec = e_rc.simulatedLastTime;
	nodeConst.keyPtr = &Key;
	nodeConst.valuePtr = &Value;
	insertObjectInHashTable(&e_nuclearReaction, &nodeConst);
}


/*
struct debug_info {
	unsigned int cnt;
};
extern int iterateDebug(void * a_total, const void * a_keyPtr, void * a_valuePtr)
{
	struct debug_info * dddPtr = (struct debug_info *)a_total;
	const struct nuclearReactionKey * kPtr = (const struct nuclearReactionKey *)a_keyPtr;
	struct nuclearReactionValue * vPtr = (struct nuclearReactionValue *)a_valuePtr;
	dddPtr->cnt++;
	fprintf(stderr, "DEBUG:%s:cnt=%u %s %s %s, form:%s maxMassDefectMeV:%lg\n", __FUNCTION__, dddPtr->cnt, 
		getElectrodeName(pKey->electrodeType),
		getNuclearReactionName(kPtr->reactionType),
		getCollideName(kPtr->collideType),
		kPtr->formPtr, vPtr->maxMassDefectMeV);
	return KEEP_NODE;
}
*/
extern void printAllNuclearReaction(FILE * a_fp)
{
	int needSortedTable = 1;
	unsigned int size;
	struct objectNodeConst ** nodePtrs;
	
	//{debug
	//struct debug_info ddd;
	//ddd.cnt = 0;
	//iterateInHashTable(&e_nuclearReaction, &ddd, iterateDebug);
	//}debug
	
	fputs("[[All Nuclear Reactions]]\n", a_fp);
	nodePtrs = getFlatTable(&e_nuclearReaction, needSortedTable, &size);
	//fprintf(stderr, "DEBUG:%s:usedCnt=%d nodePtrs:%lp\n", __FUNCTION__, e_nuclearReaction.usedCnt, nodePtrs);
	if(nodePtrs){
		int i = 0, j;
		int oldElectrodeType = -1;
		int oldReactionType = -1;
		long double totalMeV = 0.0;
		double invNAvogadro = (1.0 / NAvogadro);
		const char * electrodeTypeStr;
		char electrodeChar = 'X';
		const char * reactionTypeStr;
		char reactionChar = 'X';
		const char * decayModeStr;
		for(i = 0; i < size; ++i){
			struct nuclearReactionKey * keyPtr = (struct nuclearReactionKey *)nodePtrs[i]->keyPtr;
			struct nuclearReactionValue * valuePtr = (struct nuclearReactionValue *)nodePtrs[i]->valuePtr;
			char rangeMev[80];
			char rangeMol[80];
			char maxTimeSecStr[80], minTimeSecStr[80];
			int pretty = 1;
			//fprintf(stderr, "DEBUG:%s:keyPtr:%lp valuePtr:%lp\n", __FUNCTION__, keyPtr, valuePtr);
			//fprintf(stderr, "DEBUG:%s:reactionType:%d \n", __FUNCTION__, keyPtr->reactionType);
			//fprintf(stderr, "DEBUG:%s:formPtr:%s \n", __FUNCTION__, keyPtr->formPtr);
			if(oldElectrodeType != keyPtr->electrodeType || oldReactionType != keyPtr->reactionType){
				if(oldElectrodeType != -1 && oldReactionType != -1){
					fprintf(a_fp, "%c%c Total %Lg [MeV mol] -> %Lg [MeV]\n", electrodeChar, reactionChar, totalMeV, totalMeV * NAvogadro);
					fprintf(a_fp, "[ %s ReactionType %s ]>>>\n", electrodeTypeStr, reactionTypeStr);
				}
				oldElectrodeType = keyPtr->electrodeType;
				oldReactionType = keyPtr->reactionType;
				electrodeTypeStr = getElectrodeName(keyPtr->electrodeType);
				electrodeChar = getElectrodeChar(keyPtr->electrodeType);
				reactionTypeStr = getNuclearReactionName(keyPtr->reactionType);
				reactionChar = getNuclearReactionChar(keyPtr->reactionType);
				fprintf(a_fp, "\n[ %s ReactionType %s ]<<<\n", electrodeTypeStr, reactionTypeStr);
				j = 1;
				totalMeV = 0.0;
			}
			decayModeStr = getDecayModeText(keyPtr->decayOrCollideType);
			switch(keyPtr->reactionType){
				case REACTION_DETECT: 
				case REACTION_CANT_DETECT: 
				case REACTION_ENDOTHERMIC: 
				case REACTION_COULOMB_BARRIER:
					if(valuePtr->minMassDefectMeV == valuePtr->maxMassDefectMeV){
						snprintf(rangeMev, 80, "MassDefect %lg", valuePtr->minMassDefectMeV);
					}else{
						snprintf(rangeMev, 80, "MassDefect [ %lg , %lg ]", valuePtr->minMassDefectMeV, valuePtr->maxMassDefectMeV);
					}
					if(valuePtr->minIsotopeMol == valuePtr->maxIsotopeMol){
						snprintf(rangeMol, 80, "%lg", valuePtr->minIsotopeMol);
					}else{
						snprintf(rangeMol, 80, "[ %lg , %lg ]", valuePtr->minIsotopeMol, valuePtr->maxIsotopeMol);
					}
					fprintf(a_fp, "%c%c %d %s\n %s\n %s [MeV] * %s [mol] * %lld [cnt]\n = %Lg [MeV mol] -> %Lg [MeV]\n",
						electrodeChar, reactionChar, j, decayModeStr, 
						keyPtr->formPtr, rangeMev, rangeMol, valuePtr->sumupCnt, 
						valuePtr->sumupOfMulMeVMol, valuePtr->sumupOfMulMeVMol * NAvogadro);
					if(keyPtr->reactionType == REACTION_CANT_DETECT){
						double midMeV = (valuePtr->minMassDefectMeV + valuePtr->maxMassDefectMeV) * 0.5;
						double midMol = 0.0;
						char * markPtr = "==";
						if(midMeV != 0.0){
							midMol = valuePtr->sumupOfMulMeVMol / (midMeV * valuePtr->sumupCnt);
						}
						if(midMol > invNAvogadro){
							markPtr = ">";
						}else if(midMol < invNAvogadro){
							markPtr = "<";
						}
						fprintf(a_fp, " = middy( %lg [MeV] * %lg [mol] * %lld [cnt])\n (middy %lg [mol] %s (%lg = 1/(Avogadro:%lg)))\n", midMeV, midMol, valuePtr->sumupCnt, midMol, markPtr, invNAvogadro, NAvogadro);
					}		
					if(valuePtr->minAppliedVoltageMeV == valuePtr->maxAppliedVoltageMeV){
						fprintf(a_fp, " AppliedVoltage %lg [MeV]\n", valuePtr->minAppliedVoltageMeV);
					}else{
						fprintf(a_fp, " AppliedVoltage [ %lg , %lg ] [MeV]\n", valuePtr->minAppliedVoltageMeV, valuePtr->maxAppliedVoltageMeV);
					}
					if(valuePtr->minElectroPotentialMeV == valuePtr->maxElectroPotentialMeV){
						fprintf(a_fp, " Coulomb Barrier ElectroPotential %lg [MeV]\n", valuePtr->minElectroPotentialMeV);
					}else{
						fprintf(a_fp, " Coulomb Barrier ElectroPotential [ %lg , %lg ] [MeV]\n", valuePtr->minElectroPotentialMeV, valuePtr->maxElectroPotentialMeV);
					}
					break;
				case REACTION_ERROR:
					fprintf(a_fp, "%d:%s\n", j, keyPtr->formPtr);
					break;
			}
			fprintf(a_fp, " Time[ %s , %s ]\n\n", formatSecond(minTimeSecStr, 80, valuePtr->minTimeSec, pretty), formatSecond(maxTimeSecStr, 80, valuePtr->maxTimeSec, pretty));
			totalMeV += valuePtr->sumupOfMulMeVMol;
			++j;
		}
		if(size > 0){
			fprintf(a_fp, "%c%c Total %Lg [MeV mol] -> %Lg [MeV]\n", electrodeChar, reactionChar, totalMeV, totalMeV * NAvogadro);
			fprintf(a_fp, "[ %s ReactionType %s ]>>>\n", electrodeTypeStr, reactionTypeStr);
		}
		free(nodePtrs);
	}
}


//---------------------------------------------------------------------
extern void sprintNoDataOfAtom(struct collide * a_pX, struct atomProperty * a_atomPropertyPtr)
{
	if(a_pX->printNuclearReaction){
		if(a_atomPropertyPtr){
			snprintf(a_pX->nuclearReactionStr, REACTION_LEN,
				"There is no data of atom, %s: %s + %s -> %d%s\n", 
				getCollideName(a_pX->collideType), 
				a_pX->targetIsotopePropertyPtr->symbol, a_pX->bulletPropertyPtr->nucleusName,
				a_pX->newMassNumber, a_atomPropertyPtr->symbol);
		}else{
			snprintf(a_pX->nuclearReactionStr, REACTION_LEN,
				"There is no data of atom, %s: %s + %s -> %d~%d(mass~atom)\n", 
				getCollideName(a_pX->collideType), 
				a_pX->targetIsotopePropertyPtr->symbol, a_pX->bulletPropertyPtr->nucleusName,
				a_pX->newMassNumber, a_pX->newAtomicNumber);		
		}
	}
}
extern void checkNoDataOfAtom(struct collide * a_pX)
{
	struct atomProperty * atomPropertyPtr = getAtomPropertyPtr(a_pX->newAtomicNumber);
	if(atomPropertyPtr){
		if(atomPropertyPtr->min1secMassNumber <= a_pX->newMassNumber 
		&& a_pX->newMassNumber <= atomPropertyPtr->max1secMassNumber){
			sprintNoDataOfAtom(a_pX, atomPropertyPtr);
		}else{
			;//The new isotope has very short half-time. So, we ignore them.
		}
	}else{
		sprintNoDataOfAtom(a_pX, atomPropertyPtr);
	}
}
extern void calcElectroPotentialMeVByCollide(struct collide * a_pX)
{
	int sign;
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_C){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			sign = 1;
		}else if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_NEUTRON 
		&& (a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_NEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_DINEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_TRINEUTRON
			)
		){
			sign = 0;
		}else{
			sign = -1;
		}
	}else if(a_pX->collideType == COLLIDE_PROTON 
	|| a_pX->collideType == COLLIDE_DEUTERIUM
	|| a_pX->collideType == COLLIDE_TRITIUM){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			sign = -1;
		}else if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_NEUTRON 
		&& (a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_NEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_DINEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_TRINEUTRON
			)
		){
			sign = 0;
		}else{
			sign = 1;
		}
	}else{//COLLIDE_NEUTRON
		sign = 0;
	}
	if(sign == 0){
		a_pX->electroPotentialMeV = 0.0;
	}else{
		a_pX->electroPotentialMeV = sign * e_coefElectroPotentialMeV * (a_pX->targetIsotopePropertyPtr->atomicNumber) / (a_pX->targetIsotopePropertyPtr->nucleusRadius + a_pX->bulletPropertyPtr->nucleusRadius);
	}
	//fprintf(stderr, "DEBUG:%s:%s target:%s electroPotentialMeV=%lg\n", __FUNCTION__, getCollideName(a_pX->collideType), a_pX->targetIsotopePropertyPtr->symbol, a_pX->electroPotentialMeV);
}
extern void set_minMassDefectMeV_plusAlpha(struct collide * a_pX)
{
			
	if(a_pX->minAppliedVoltageMeV >= 0.0){
		a_pX->minMassDefectMeV = 0.0;
	}else{
		a_pX->minMassDefectMeV = - a_pX->minAppliedVoltageMeV;
		a_pX->minAppliedVoltageMeV = 0.0;
	}
	if(a_pX->massDefectMeV >= 0.0){
		a_pX->plusAlpha = " + x";
	}else{
		a_pX->plusAlpha = " - x";
	}
}
extern void getCollidedNewIsotopePropertyPtr(struct collide * a_pX)
{
	//fprintf(stderr, "DEBUG:%s:BEGIN{%lp\n", __FUNCTION__, a_pX);
	calcElectroPotentialMeVByCollide(a_pX);
	a_pX->newAtomicNumber = UNDEF_ATOMIC_NUMBER;//for safety
	a_pX->newMassNumber = UNDEF_MASS_NUMBER;//for safety
	a_pX->newIsotopePropertyPtr = NULL;//for safety
	a_pX->minAppliedVoltageMeV = 0.0;//for safety
	a_pX->massDefectMeV = 0.0;//for safety
	a_pX->minMassDefectMeV = 0.0;//for safety
	a_pX->plusAlpha = "";//for safety
	a_pX->nuclearReactionStr[0] = 0;//for safety
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_C){
		a_pX->isElectronCapture = 1;
		a_pX->newAtomicNumber = a_pX->targetIsotopePropertyPtr->atomicNumber - 1;
		a_pX->newMassNumber = a_pX->targetIsotopePropertyPtr->massNumber;
		if(a_pX->newAtomicNumber != ATOMICNUMBER_ELECTRON || a_pX->newMassNumber != MASSNUMBER_ELECTRON){
			a_pX->newIsotopePropertyPtr = getIsotopePropertyPtr(a_pX->newAtomicNumber, a_pX->newMassNumber);
			if(a_pX->newIsotopePropertyPtr){
				a_pX->minAppliedVoltageMeV = a_pX->newIsotopePropertyPtr->massMeV - a_pX->targetIsotopePropertyPtr->massMeV;
				a_pX->massDefectMeV = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV;
				set_minMassDefectMeV_plusAlpha(a_pX);
				if(a_pX->printNuclearReaction){
					double targetIonMassMeV = a_pX->targetIsotopePropertyPtr->massMeV - e_massElectronMeV;
					char * targetSymbolPtr;
					if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_HYDROGEN){
						targetSymbolPtr = a_pX->targetIsotopePropertyPtr->nucleusName;
					}else{
						targetSymbolPtr = a_pX->targetIsotopePropertyPtr->symbol;
					}
					snprintf(a_pX->nuclearReactionStr, REACTION_LEN,
						"%s+(%lg) + e-(%lg) + (%lg%s)[MeV] -> %s(%lg) + νe + (%lg%s)[MeV]", 
						targetSymbolPtr, targetIonMassMeV,
						e_massElectronMeV,
						a_pX->minAppliedVoltageMeV, a_pX->plusAlpha,
						a_pX->newIsotopePropertyPtr->symbol, a_pX->newIsotopePropertyPtr->massMeV,
						a_pX->minMassDefectMeV, a_pX->plusAlpha);
					//fprintf(stderr, "DEBUG:%s:snprintf %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
				}
			}else{
				checkNoDataOfAtom(a_pX);
			}
		}else{
			fprintf(stderr, "ERROR:%s(%d):COLLIDE_ELECTRON:%d with ELECTRON\n", __FUNCTION__, __LINE__, a_pX->collideType);
		}
	}else if(a_pX->collideType == COLLIDE_PROTON 
	|| a_pX->collideType == COLLIDE_DEUTERIUM 
	|| a_pX->collideType == COLLIDE_TRITIUM){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			//fprintf(stderr, "DEBUG:%s:COLLIDE_PROTON... ATOMICNUMBER_ELECTRON\n", __FUNCTION__);
			a_pX->isElectronCapture = 1;
			a_pX->newAtomicNumber = ATOMICNUMBER_NEUTRON;
			switch(a_pX->collideType){
				case COLLIDE_PROTON:    a_pX->newMassNumber = MASSNUMBER_NEUTRON;    break;
				case COLLIDE_DEUTERIUM: a_pX->newMassNumber = MASSNUMBER_DINEUTRON;  break;
				case COLLIDE_TRITIUM:   a_pX->newMassNumber = MASSNUMBER_TRINEUTRON; break;
			}
			a_pX->newIsotopePropertyPtr = getIsotopePropertyPtr(a_pX->newAtomicNumber, a_pX->newMassNumber);
			if(a_pX->newIsotopePropertyPtr){
				a_pX->minAppliedVoltageMeV = a_pX->newIsotopePropertyPtr->massMeV - a_pX->bulletPropertyPtr->massMeV;
				a_pX->massDefectMeV = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV;
				set_minMassDefectMeV_plusAlpha(a_pX);
				if(a_pX->printNuclearReaction){
					double bulletIonMassMeV = a_pX->bulletPropertyPtr->massMeV - e_massElectronMeV;
					//fprintf(stderr, "DEBUG:%s:snprintf COLLIDE_PROTON\n", __FUNCTION__);
					snprintf(a_pX->nuclearReactionStr, REACTION_LEN,
						"e-(%lg) + %s+(%lg) + (%lg%s)[MeV] -> %s(%lg) + νe + (%lg%s)[MeV]", 
						e_massElectronMeV,
						a_pX->bulletPropertyPtr->nucleusName, bulletIonMassMeV,
						a_pX->minAppliedVoltageMeV, a_pX->plusAlpha,
						a_pX->newIsotopePropertyPtr->symbol, a_pX->newIsotopePropertyPtr->massMeV,
						a_pX->minMassDefectMeV, a_pX->plusAlpha);
					//fprintf(stderr, "DEBUG:%s:snprintf COLLIDE_PROTON %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
				}
			}else{
				checkNoDataOfAtom(a_pX);
			}
		}else{
			//fprintf(stderr, "DEBUG:%s:COLLIDE_PROTON... other\n", __FUNCTION__);
			a_pX->isElectronCapture = 0;
			a_pX->newAtomicNumber = a_pX->targetIsotopePropertyPtr->atomicNumber + a_pX->bulletPropertyPtr->atomicNumber;
			a_pX->newMassNumber = a_pX->targetIsotopePropertyPtr->massNumber + a_pX->bulletPropertyPtr->massNumber;
			a_pX->newIsotopePropertyPtr = getIsotopePropertyPtr(a_pX->newAtomicNumber, a_pX->newMassNumber);
			if(a_pX->newIsotopePropertyPtr){
				a_pX->minAppliedVoltageMeV = a_pX->newIsotopePropertyPtr->massMeV - (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV);
				a_pX->massDefectMeV = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV;
				set_minMassDefectMeV_plusAlpha(a_pX);
				if(a_pX->printNuclearReaction){
					char * targetSymbol;
					char * targetIonSign;
					char * newSymbol;
					char * newIonSign;
					double targetIonMassMeV, bulletIonMassMeV, newAtomMassMeV;
					if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_NEUTRON
					&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_NEUTRON){
						targetSymbol = a_pX->targetIsotopePropertyPtr->symbol;
						targetIonSign = "";
						targetIonMassMeV = a_pX->targetIsotopePropertyPtr->massMeV;
						bulletIonMassMeV = a_pX->bulletPropertyPtr->massMeV - e_massElectronMeV;
						newAtomMassMeV = a_pX->newIsotopePropertyPtr->massMeV - e_massElectronMeV;
						newSymbol = a_pX->newIsotopePropertyPtr->nucleusName;
						newIonSign = "+";
					}else{
						targetSymbol = a_pX->targetIsotopePropertyPtr->symbol;
						targetIonSign = "-";
						targetIonMassMeV = a_pX->targetIsotopePropertyPtr->massMeV + e_massElectronMeV;
						bulletIonMassMeV = a_pX->bulletPropertyPtr->massMeV - e_massElectronMeV;
						newAtomMassMeV = a_pX->newIsotopePropertyPtr->massMeV;
						newSymbol = a_pX->newIsotopePropertyPtr->symbol;
						newIonSign = "";
					}
					//fprintf(stderr, "DEBUG:%s:snprintf other\n", __FUNCTION__);
					snprintf(a_pX->nuclearReactionStr, REACTION_LEN,
						"%s%s(%lg) + %s+(%lg) + (%lg%s)[MeV] -> %s%s(%lg) + (%lg%s)[MeV]", 
						targetSymbol, targetIonSign, targetIonMassMeV,
						a_pX->bulletPropertyPtr->nucleusName, bulletIonMassMeV,
						a_pX->minAppliedVoltageMeV, a_pX->plusAlpha,
						newSymbol, newIonSign, newAtomMassMeV,
						a_pX->minMassDefectMeV, a_pX->plusAlpha);
					//fprintf(stderr, "DEBUG:%s:snprintf other %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
				}
			}else{
				checkNoDataOfAtom(a_pX);
			}
		}
	}else if(a_pX->collideType == COLLIDE_NEUTRON){ 
		a_pX->isElectronCapture = 0;
		a_pX->newAtomicNumber = a_pX->targetIsotopePropertyPtr->atomicNumber + a_pX->bulletPropertyPtr->atomicNumber;
		a_pX->newMassNumber = a_pX->targetIsotopePropertyPtr->massNumber + a_pX->bulletPropertyPtr->massNumber;
		a_pX->newIsotopePropertyPtr = getIsotopePropertyPtr(a_pX->newAtomicNumber, a_pX->newMassNumber);
		if(a_pX->newIsotopePropertyPtr){
			a_pX->massDefectMeV = (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV) - a_pX->newIsotopePropertyPtr->massMeV;
			if(a_pX->printNuclearReaction){
				//fprintf(stderr, "DEBUG:%s:snprintf other\n", __FUNCTION__);
				snprintf(a_pX->nuclearReactionStr, REACTION_LEN,
					"%s(%lg) + %s(%lg) -> %s(%lg) + %lg[MeV]", 
					a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->massMeV,
					a_pX->bulletPropertyPtr->nucleusName, a_pX->bulletPropertyPtr->massMeV,
					a_pX->newIsotopePropertyPtr->symbol, a_pX->newIsotopePropertyPtr->massMeV,
					a_pX->massDefectMeV);
				//fprintf(stderr, "DEBUG:%s:snprintf other %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
			}
		}else{
			checkNoDataOfAtom(a_pX);
		}
	}
	//fprintf(stderr, "DEBUG:%s:END}\n", __FUNCTION__);
}
extern int checkCollidingUsage(struct collide * a_pX)
//const struct isotopeProperty * a_isotopeProperty, int a_useElectron, int a_useNeutron)
{
	int ret = 0;
	if(a_pX->targetIsotopePropertyPtr->atomicNumber >= 1
	|| 
	(a_pX->useElectron == 1 
	&& a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
	&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON)
	||
	(a_pX->useNeutron == 1 
	&& a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_NEUTRON 
	&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_NEUTRON)
	){
		ret = 1;
	}
	return ret;
}
extern double getCrossSection(struct collide * a_pX)
{
	double CrossSection;
	if(e_usebulletCrossSection){
		CrossSection = (a_pX->targetIsotopePropertyPtr->relativeCollideCrossSection + a_pX->bulletPropertyPtr->relativeCollideCrossSection) * a_pX->targetAtomMol;
	}else{
		CrossSection = a_pX->targetIsotopePropertyPtr->relativeCollideCrossSection * a_pX->targetAtomMol;
	}
	return CrossSection;
}
extern int iterateCrossSection(void * a_total, struct objectNodeConst * a_nodePtr)
{
	struct collide * pX = (struct collide *)a_total;
	//fprintf(stderr, "DEBUG:%s:{BEGIN:%lp %lp %lp\n", __FUNCTION__, a_total, a_keyPtr, a_valuePtr);
	pX->targetIsotopePropertyPtr = ((struct atomKey *)a_nodePtr->keyPtr)->isotopePropertyPtr;
	pX->targetAtomValuePtr = (struct atomValue *)a_nodePtr->valuePtr;

	if(checkCollidingUsage(pX)){
		getCollidedNewIsotopePropertyPtr(pX);
		if(pX->newIsotopePropertyPtr){
			if(pX->massDefectMeV > 0.0){
				if(pX->appliedVoltageMeV >= pX->electroPotentialMeV){
					pX->targetAtomMol = getMol(pX->targetAtomValuePtr);
					if(pX->targetAtomMol >= pX->electrodePtr->detectLimitMolForIsotope){
						pX->totalCollideCrossSection += getCrossSection(pX);
						pX->totalTargetMol += pX->targetAtomMol;
						//fprintf(stderr, "DEBUG:%s:totalCollideCrossSection:%lg\n", __FUNCTION__, pX->totalCollideCrossSection);
					}
				}
			}
		}
	}
	//fprintf(stderr, "DEBUG:%s:{END:totalCollideCrossSection:%lg\n", __FUNCTION__, pX->totalCollideCrossSection);
	return KEEP_NODE;
}

extern double calcCollideCrossSectionRate(struct collide * a_pX)
{
	double collideCrossSectionRate;
	if(a_pX->totalCollideCrossSection > 0.0){
		collideCrossSectionRate = getCrossSection(a_pX) / a_pX->totalCollideCrossSection;
	}else{
		collideCrossSectionRate = 0.0;
	}
	return collideCrossSectionRate;
}

extern double collideParticleAtom(struct collide * a_pX)
{
	//The collision of a particle like an electron, a proton, a deterium or a tritium, it may have kinetic energy,  it will be absorbed by neuclaus.
	
	double collidedMol = 0.0, collideCrossSectionRate;
	//fprintf(stderr, "DEBUG:%s:{BEGIN%lp\n", __FUNCTION__, a_pX);
	if(checkCollidingUsage(a_pX)){
		getCollidedNewIsotopePropertyPtr(a_pX);
		if(a_pX->newIsotopePropertyPtr){
			if(a_pX->massDefectMeV > 0.0){
				if(a_pX->appliedVoltageMeV >= a_pX->electroPotentialMeV){
					a_pX->targetAtomMol = getMol(a_pX->targetAtomValuePtr);
					if(a_pX->targetAtomMol >= a_pX->electrodePtr->detectLimitMolForIsotope){
						collideCrossSectionRate = calcCollideCrossSectionRate(a_pX);
						a_pX->newIsotopeMol = a_pX->collideBulletMol * collideCrossSectionRate;
						if(a_pX->newIsotopeMol < 0.0){
							a_pX->newIsotopeMol = 0.0;//collect tolelance.
						}
						if(a_pX->newIsotopeMol > a_pX->targetAtomMol){
							if(a_pX->newIsotopeMol > a_pX->targetAtomMol * 4.0){
								fprintf(stderr, "WARN:%s(%d):%s newIsotope %s Mol:%lg > targetAtom %s mol:%lg\n",
									__FUNCTION__, __LINE__, getCollideName(a_pX->collideType), 
									a_pX->newIsotopePropertyPtr->symbol, a_pX->newIsotopeMol, 
									a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetAtomMol);
							}
							if(a_pX->newIsotopeMol > a_pX->targetAtomMol * 8.0){
								fprintf(e_logFp, "WARN:%s(%d):%s newIsotope %s Mol:%lg > targetAtom %s mol:%lg\n",
									__FUNCTION__, __LINE__, getCollideName(a_pX->collideType), 
									a_pX->newIsotopePropertyPtr->symbol, a_pX->newIsotopeMol, 
									a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetAtomMol);
							}
							a_pX->newIsotopeMol = a_pX->targetAtomMol;
						}
						//fprintf(stderr, "DEBUG:%s:%s a_pX->newIsotopeMol:%lg\n", __FUNCTION__, a_pX->newIsotopePropertyPtr->symbol, a_pX->newIsotopeMol);
						if(a_pX->newIsotopeMol >= a_pX->electrodePtr->detectLimitMolForIsotope){
							registNuclearReaction(a_pX->electrodePtr, REACTION_DETECT, a_pX->collideType, a_pX->nuclearReactionStr, a_pX->massDefectMeV, a_pX->newIsotopeMol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
							if(a_pX->isElectronCapture){
								registOutputOfElectronCapture(a_pX->electrodePtr, a_pX->massDefectMeV, a_pX->newIsotopeMol);
							}else{
								registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, a_pX->massDefectMeV, a_pX->newIsotopeMol);
							}
							collidedMol = a_pX->newIsotopeMol;
							setMolSub(a_pX->targetAtomValuePtr, a_pX->newIsotopeMol);
							if(a_pX->newAtomicNumber == ATOMICNUMBER_NEUTRON 
							&& (a_pX->newMassNumber == MASSNUMBER_DINEUTRON 
							|| a_pX->newMassNumber == MASSNUMBER_TRINEUTRON)){
								// The dineutron and trineutron will imediately separate into neutrons.
								a_pX->newIsotopeMol *= a_pX->newMassNumber;
								a_pX->newMassNumber = MASSNUMBER_NEUTRON;
								a_pX->newIsotopePropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON);
							}
							increaseAtom(a_pX->electrodePtr, a_pX->newAtomicNumber, a_pX->newMassNumber, a_pX->newIsotopeMol);
						}else if(a_pX->newIsotopeMol > 0.0){
							registNuclearReaction(a_pX->electrodePtr, REACTION_CANT_DETECT, a_pX->collideType, a_pX->nuclearReactionStr, a_pX->massDefectMeV, a_pX->newIsotopeMol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
							//fprintf(stderr, "DEBUG:%s:%s a_pX->newIsotopeMol:%lg < detectLimitMolForIsotope:%lg\n %s\n",
							//__FUNCTION__, a_pX->newIsotopePropertyPtr->symbol, 
							//a_pX->newIsotopeMol, a_pX->electrodePtr->detectLimitMolForIsotope, a_pX->nuclearReactionStr);
						}else{ //a_pX->newIsotopeMol <= 0.0
							;//Do nothing!
						}
					}else{
						registNuclearReaction(a_pX->electrodePtr, REACTION_CANT_DETECT, a_pX->collideType, a_pX->nuclearReactionStr, a_pX->massDefectMeV, a_pX->targetAtomMol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
					}
				}else{
					registNuclearReaction(a_pX->electrodePtr, REACTION_COULOMB_BARRIER, a_pX->collideType, a_pX->nuclearReactionStr, a_pX->massDefectMeV, 0.0, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);					
				}
			}else{
				registNuclearReaction(a_pX->electrodePtr, REACTION_ENDOTHERMIC, a_pX->collideType, a_pX->nuclearReactionStr, a_pX->massDefectMeV, 0.0, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
			}
		}else{
			if(a_pX->nuclearReactionStr[0]){
				registNuclearReaction(a_pX->electrodePtr, REACTION_ERROR, a_pX->collideType, a_pX->nuclearReactionStr, a_pX->massDefectMeV, 0.0, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
				//fprintf(stderr, "%s\n", a_pX->nuclearReactionStr);
			}
		}
	}
	//fprintf(stderr, "DEBUG:%s:}END:collidedMol:%lp\n", __FUNCTION__, collidedMol);
	return collidedMol;
}

extern int iterateCollide(void * a_total, struct objectNodeConst * a_nodePtr)
{
	struct collide * pX = (struct collide *)a_total;
	//fprintf(stderr, "DEBUG:%s:{BEGIN:%lp %lp %lp\n", __FUNCTION__, a_total, a_keyPtr, a_valuePtr);
	pX->targetIsotopePropertyPtr = ((struct atomKey *)a_nodePtr->keyPtr)->isotopePropertyPtr;
	pX->targetAtomValuePtr = (struct atomValue *)a_nodePtr->valuePtr;
	pX->collidedBulletMol += collideParticleAtom(pX);
	//fprintf(stderr, "DEBUG:%s:}END:collidedBulletMol:%lg\n", __FUNCTION__, pX->collidedBulletMol);
	return KEEP_NODE;
}
//int e_inComptonEffect = 0;//DEBUG
//int e_cntCollideBulletToElectrode = 0;//DEBUG
extern void collideBulletToElectrode(struct electrode * a_electrodePtr, int a_collideType, double a_appliedVoltageMeV, double a_arriveBulletMol, double a_collideBulletRate, struct atomNodeConst * a_decreaseBulletPtr)
{
	struct collide col;
	//fprintf(stderr, "DEBUG:%s:{BEGIN a_collideType:%d a_appliedVoltageMeV:%lg a_arriveBulletMol:%lg a_collideBulletRate:%lg\n", __FUNCTION__, a_collideType, a_appliedVoltageMeV, a_arriveBulletMol, a_collideBulletRate);
	//e_cntCollideBulletToElectrode++;//DEBUG
	memset(&col, 0, sizeof(struct collide));
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	col.electrodePtr = a_electrodePtr;
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	col.collideType = a_collideType;
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	switch(col.collideType){
		case COLLIDE_ELECTRON:
		case COLLIDE_ELECTRON_C: col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON); break;
		case COLLIDE_PROTON:    col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN); break;
		case COLLIDE_DEUTERIUM: col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_DEUTERIUM, MASSNUMBER_DEUTERIUM);  break;
		case COLLIDE_TRITIUM:   col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_TRITIUM, MASSNUMBER_TRITIUM); break;
		case COLLIDE_NEUTRON:   col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON);
		if(a_appliedVoltageMeV != 0.0){
			fprintf(stderr, "FATAL ERROR:%s:%s a_appliedVoltageMeV:%lg!= 0.0\n", __FUNCTION__, getCollideName(col.collideType), a_appliedVoltageMeV);
			exit(1);
		}
		break;
	}
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	col.appliedVoltageMeV = a_appliedVoltageMeV;
	col.arriveBulletMol = a_arriveBulletMol;
	col.collideBulletMol = col.arriveBulletMol * a_collideBulletRate;
	col.imperfectCollideMol = col.arriveBulletMol - col.collideBulletMol;
	col.collidedBulletMol = 0.0;
	col.remainBulletMol = 0.0;
	col.totalCollideCrossSection = 0.0;
	col.totalTargetMol = 0.0;
	//fprintf(stderr, "DEBUG:%s(%d):col.collideBulletMol %lg col.electrodePtr->detectLimitMolForIsotope %lg\n", __FUNCTION__, __LINE__, col.collideBulletMol, col.electrodePtr->detectLimitMolForIsotope);
	if(col.collideBulletMol >= col.electrodePtr->detectLimitMolForIsotope){
		useElectronOrNeutron(&col);
#define KEEP_ATMIC_ORDER 1 //we should keep the atomic order to keep the same reaction order
#if  KEEP_ATMIC_ORDER == 1
		//fprintf(stderr, "DEBUG:%s:KEEP_ATMIC_ORDER == 1\n", __FUNCTION__);
		{
			unsigned int size, i;
			struct objectNodeConst ** oPtr;
			oPtr = getFlatTable(&col.electrodePtr->atomHashTable, SORT_ASCEND, &size);
			if(oPtr){
				//fprintf(stderr, "DEBUG:%s:iterateCrossSection usedCnt:%d\n", __FUNCTION__, col.electrodePtr->atomHashTable.usedCnt);
				for(i = 0; i < size; ++i){
					//fprintf(stderr, "DEBUG:%s:call iterateCrossSection i:%d\n", __FUNCTION__, i);
					iterateCrossSection(&col, oPtr[i]);
				}
				//fprintf(stderr, "DEBUG:%s:inCompton:%d\n", __FUNCTION__, e_inComptonEffect);
				//fprintf(stderr, "DEBUG:%s:totalCollideCrossSection:%lg\n", __FUNCTION__, col.totalCollideCrossSection);
				col.printNuclearReaction = 1;
				//fprintf(stderr, "DEBUG:%s:iterateCrossSection usedCnt:%d\n", __FUNCTION__, col.electrodePtr->atomHashTable.usedCnt);
				for(i = 0; i < size; ++i){
					//fprintf(stderr, "DEBUG:%s:call iterateCollide i:%d\n", __FUNCTION__, i);
					iterateCollide(&col, oPtr[i]);
				}
				free(oPtr);
			}
		}
#else
		//fprintf(stderr, "DEBUG:%s:KEEP_ATMIC_ORDER != 1\n", __FUNCTION__);
		//fprintf(stderr, "DEBUG:%s:call iterateCrossSection\n", __FUNCTION__);
		iterateInHashTable(&col.electrodePtr->atomHashTable, &col, iterateCrossSection);
		col.printNuclearReaction = 1;
		//fprintf(stderr, "DEBUG:%s:call iterateCollide\n", __FUNCTION__);
		iterateInHashTable(&col.electrodePtr->atomHashTable, &col, iterateCollide);
#endif
	}else if(col.collideBulletMol == 0.0){
		//fprintf(stderr, "DEBUG:%s:col.collideBulletMol == 0.0\n", __FUNCTION__);
	}
	if(col.collidedBulletMol > 0.0 && a_decreaseBulletPtr){
		setMolSub(a_decreaseBulletPtr->vPtr, col.collidedBulletMol);//decrease mol of atoms after processing collision.
	}
	col.remainBulletMol = col.collideBulletMol - col.collidedBulletMol;
	//fprintf(stderr, "DEBUG:%s:col.remainBulletMol:%lg\n", __FUNCTION__, col.remainBulletMol);
	if(col.remainBulletMol < 0.0){
		col.remainBulletMol = 0.0;//collect tolelance.
	}
	col.imperfectCollideMol += col.remainBulletMol;
	//fprintf(stderr, "DEBUG:%s:col.imperfectCollideMol:%lg\n", __FUNCTION__, col.imperfectCollideMol);
	if(col.imperfectCollideMol > 0.0){
		if(col.appliedVoltageMeV > 0.0){
			registOutput(col.electrodePtr, SCAT_IMPERFECT, col.appliedVoltageMeV, col.imperfectCollideMol);
			//fprintf(stderr, "DEBUG:%s:appliedVoltageMeV:%lg imperfectCollideMol:%lg remainBulletMol:%lg collideBulletMol:%lg collidedBulletMol:%lg\n", __FUNCTION__, col.appliedVoltageMeV, col.imperfectCollideMol, col.remainBulletMol, col.collideBulletMol, col.collidedBulletMol);
		}
		if(col.collideType == COLLIDE_ELECTRON || col.collideType == COLLIDE_ELECTRON_C){
			;//The electrons flew from the negative electrode to the positive electorode. But the electric circuit rotate electrons from the positive electrode to the negative electorode. So the amount of electrons is always same.
		}else if(col.collideType == COLLIDE_PROTON || col.collideType == COLLIDE_DEUTERIUM || col.collideType == COLLIDE_TRITIUM){
			//When the protons, deuteriums and tritiums flew from the positive electrode to the negative electorode with huge energy greater than the energy of beta decay, some of them collide imperfectly, so they will change from neurcuses to atoms.
			//We need to append the isotopes of hydrogen onto the negative electorode.
			increaseAtom(col.electrodePtr, col.bulletPropertyPtr->atomicNumber, col.bulletPropertyPtr->massNumber, col.imperfectCollideMol);
		}else  if(col.collideType == COLLIDE_NEUTRON){
			;//The neutrons do not fly from a electrode to another electrodes, because they do not have any electoric charge.
		}
	}
	//fprintf(stderr, "DEBUG:%s:}END\n", __FUNCTION__);
}

extern void generateNeutronInSpace(double a_pastSecond, int a_DeuteriumMassnumber, double a_protonRate)
{
	//We may be able to re-use the algorithum of the function "collideBulletToElectrode".
	
	//fprintf(stderr, "DEBUG:%s:{BEGIN a_pastSecond:%lg a_DeuteriumMassnumber:%d a_protonRate:%lg\n", __FUNCTION__, a_pastSecond, a_DeuteriumMassnumber, a_protonRate);
	if(a_protonRate > 0.0){
		int eNum;
		for(eNum = 0; eNum < 2; ++eNum){
			struct collide col;
			memset(&col, 0, sizeof(col));
			col.electrodePtr = (eNum == 0) ? &e_negativeElectrode : &e_positiveElectrode;
			col.collideType = COLLIDE_ELECTRON;
			col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON);
			col.appliedVoltageMeV = e_appliedVoltageMeV;
			col.arriveBulletMol = a_pastSecond * e_genedNeutronInSpaceMol * a_protonRate * ((eNum == 0) ? e_neutronGenInSpaceFlyRate : e_neutronGenInSpaceFallRate);
			col.collideBulletMol = col.arriveBulletMol;
			col.imperfectCollideMol = col.arriveBulletMol - col.collideBulletMol;//This is 0.0
			col.collidedBulletMol = 0.0;
			col.remainBulletMol = 0.0;
			col.totalCollideCrossSection = 0.0;
			col.totalTargetMol = 0.0;
			//fprintf(stderr, "DEBUG:%s:e_genedNeutronInSpaceMol:%lg eNum:%d electrodePtr:%lp col.arriveBulletMol:%lg\n", __FUNCTION__, e_genedNeutronInSpaceMol, eNum, col.electrodePtr, col.arriveBulletMol);
			if(col.collideBulletMol >= col.electrodePtr->detectLimitMolForIsotope){
				struct objectNodeConst nodeConst;
				struct atomKey dummyTargetAtomKey;
				struct atomValue dummyTargetAtomValue;
				useElectronOrNeutron(&col);
				dummyTargetAtomKey.isotopePropertyPtr = getIsotopePropertyPtr(MASSNUMBER_HYDROGEN, a_DeuteriumMassnumber);
				memset(&dummyTargetAtomValue, 0, sizeof(dummyTargetAtomValue));
				dummyTargetAtomValue.molIni = col.arriveBulletMol;
				nodeConst.keyPtr = &dummyTargetAtomKey;
				nodeConst.valuePtr = &dummyTargetAtomValue;

				//fprintf(stderr, "DEBUG:%s:eNum:%d call iterateCrossSection\n", __FUNCTION__, eNum);
				iterateCrossSection(&col, &nodeConst);
				col.printNuclearReaction = 1;
				//fprintf(stderr, "DEBUG:%s:eNum:%d iterateCollide\n", __FUNCTION__, eNum);
				iterateCollide(&col, &nodeConst);
			}else if(col.collideBulletMol == 0.0){
				//fprintf(stderr, "DEBUG:%s:eNum:col.collideBulletMol:%lg == 0.0\n", __FUNCTION__, col.collideBulletMol);
			}
			col.remainBulletMol = col.collideBulletMol - col.collidedBulletMol;
			if(col.remainBulletMol < 0.0){
				col.remainBulletMol = 0.0;//collect tolelance.
			}
			col.imperfectCollideMol += col.remainBulletMol;
			if(col.imperfectCollideMol > 0.0){
				if(col.appliedVoltageMeV > 0.0){
					registOutput(col.electrodePtr, SCAT_IMPERFECT_IN_SPACE, col.appliedVoltageMeV, col.imperfectCollideMol);
				}
			}
		}
	}
	//fprintf(stderr, "DEBUG:%s:}END\n", __FUNCTION__);
}
extern void checkMassInNuclearReaction(double a_formulaOfMassDefect, double a_massDefect, const char * a_nuclearReactionPtr)
{
	double diff;
	diff = a_formulaOfMassDefect - a_massDefect;
	if(diff < 0.0){
		diff = - diff;
	}
	if(a_massDefect > 0.0){
		if(diff / a_massDefect > 1e-8){
			fprintf(stderr, "ERROR:%s:Wrong %lg != a_massDefect:%lg\nNuclearReaction:%s\n", __FUNCTION__, a_formulaOfMassDefect, a_massDefect, a_nuclearReactionPtr);
			fprintf(e_logFp, "ERROR:%s:Wrong %lg != a_massDefect:%lg\nNuclearReaction:%s\n", __FUNCTION__, a_formulaOfMassDefect, a_massDefect, a_nuclearReactionPtr);
			exit(1);
		}
	}else{
		fprintf(stderr, "ERROR:%s:Zero or negative a_massDefect:%lg\nNuclearReaction:%s\n", __FUNCTION__, a_massDefect, a_nuclearReactionPtr);
		fprintf(e_logFp, "ERROR:%s:Zero or negative a_massDefect:%lg\nNuclearReaction:%s\n", __FUNCTION__, a_massDefect, a_nuclearReactionPtr);
		exit(1);
	}
}
//double e_debugDecayMassDiff = 0.0;
extern void decayNeuclay(struct electrode * a_electrodePtr, const struct isotopeProperty * a_isotopePropertyPtr, struct atomValue * a_atomValuePtr, double a_pastSecond)
{
	//[PHYSICS REVIEW] What is the half life?
	// When the half life of an isotope passed, the half amount of isotope decay.
	//
	// if(halflife == 1){ (afterTime):(survivalRate), 1:0.5^1, 2:0.5^2, 3:0.5^3 4:0.5^4, .... }
	// if(halflife == 2){ (afterTime):(survivalRate), 2:0.5^1, 4:0.5^2, 6:0.5^3 8:0.5^4, .... }
	// if(halflife == h){ (afterTime):(survivalRate),
	//  h:0.5^(h/h), 2h:0.5^(2h/h), 3h:0.5^(3h/h), 4h:0.5^(4h/h), x:0.5^(x/h)....
	// }
	//1[year]  = 31536000[sec] = 60*60*24*365
	//1[hour]  = 3600[sec] = 60*60
	//isotopes, halfLife, survival rate after 1 second, survival rate after 1 year
	//neutron,  886.7[sec], pow(0.5, 1/886.7)             = 0.9992185899073248, pow(0.5, 31536000/886.7) = 0
	//tritium, 12.32[year], pow(0.5, 1/(12.32*31536000))  = 0.9999999982159454, pow(0.5, 1/12.32)  = 0.9452914876844036
	//60Co,   5.2713[year], pow(0.5, 1/(5.2713*31536000)) = 0.9999999958303354, pow(0.5, 1/5.2713) = 0.8767840603933564
	//65Ni,   2.5172[hour], pow(0.5, 1/(2.5172*3600))     = 0.999923512823796,  pow(0.5, 31536000/(2.5172*3600)) = 0
	//The almost all isotopes that have shorter half life than 7 days will decay and servive olny 2.0r-16 after 1 year.
	//pow(0.5, 31536000/(3600*24)) = 1.3306124500025471e-110
	//pow(0.5, 31536000/(3600*24*7)) = 2.0111105320273042e-16
	//pow(0.5, 31536000/(3600*24*15)) = 4.7308237909323e-8
	//decayRate = 1.0 - survivalRate
	if(a_isotopePropertyPtr->halfLifeUnit != HLU_STABLE){
		double halfLifeSec;
		double decayRate;
		double mol, decayMol;
		double partialMol;
		int calcDecay, i;
		halfLifeSec = calcHalfLifeSec(a_isotopePropertyPtr);
		decayRate = 1.0 - pow(0.5, (a_atomValuePtr->pastSecond + a_pastSecond) / halfLifeSec);
		if(decayRate < 0.0){
			decayRate = 0.0;
		}
		if(decayRate > 1.0){
			decayRate = 1.0;
		}
		mol = getMol(a_atomValuePtr);
		decayMol = mol * decayRate;
		calcDecay = 0;
		if(decayMol >= a_electrodePtr->detectLimitMolForIsotope){
			if(a_isotopePropertyPtr->decayModeSize > 0){
				for(i = 0; i < a_isotopePropertyPtr->decayModeSize; ++i){
					partialMol = a_isotopePropertyPtr->decayModeRate[i] * decayMol;
					if(partialMol >= a_electrodePtr->detectLimitMolForIsotope){
						calcDecay = 1;
						break;
					}
				}
			}else{
				fprintf(stderr, "ERROR:%s:no decay mode:%s\n", __FUNCTION__, a_isotopePropertyPtr->symbol);
				fprintf(e_logFp, "ERROR:%s:no decay mode:%s\n", __FUNCTION__, a_isotopePropertyPtr->symbol);
			}
		}
		if(calcDecay){
			a_atomValuePtr->pastSecond = 0.0;
			for(i = 0; i < a_isotopePropertyPtr->decayModeSize; ++i){
				partialMol = a_isotopePropertyPtr->decayModeRate[i] * decayMol;
				if(partialMol >= a_electrodePtr->detectLimitMolForIsotope){
					int daughterAtomicNumber, daughterMassNumber;
					struct isotopeProperty * daughterIsotopePropertyPtr = NULL;
					int showErrorOfUndefinedDaughter = 0;
					int showErrorNegativeMassDefect = 0;
					char nuclearReactionStr[REACTION_LEN];
					double massDefectMeV;
					if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_ALPHA){
						//The Alpha Nuclear Reaction is : X -> Y + He,  He = α + 2 e-
						// atomic number - 2, mass number - 4
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 2;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - 4;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_helium4Ptr->massMeV);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + He(%lg) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										e_helium4Ptr->massMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_helium4Ptr->massMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM, partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + 						e_helium4Ptr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_ELECTRON_CAPTURE){
						// decrease atomic number, same mass number
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 1;
						daughterMassNumber = a_isotopePropertyPtr->massNumber;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massElectronMeV + e_massElectronMeV);
							if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS 
							&& massDefectMeV > 0.0){
								//The beta plus is as same as the positron(e+) emission
								//The Nuclear Reaction is : X           -> Y + e- + e+ + νe
								//                        : p+ + energy -> n      + e+ + νe [!CAUTION!]the proton must has energy.
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s(%lg) -> %s(%lg) + e-(%lg) + e+(%lg) + νe(0) + %lg[MeV/c^2]", 
									a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
									daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
									e_massElectronMeV, e_massElectronMeV,
									massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massElectronMeV + e_massElectronMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, - partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV, partialMol);
							}else if((massDefectMeV = a_isotopePropertyPtr->massMeV - daughterIsotopePropertyPtr->massMeV) > 0.0){
								//The electron capture do not emmit a positron(e+), only emmits a neutorino.
								//The Nuclear Reaction is : X+ + e-          -> Y + νe
								//                        : X                -> Y + νe
								//                        : p+ + e- + energy -> n + νe   [!CAUTION!]the proton and electron must has energy.
								int decayMode = a_isotopePropertyPtr->decayMode[i];
								if(decayMode == DECAY_MODE_BETA_PLUS){
									decayMode = DECAY_MODE_ELECTRON_CAPTURE_DEGRADE;
								}
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s+(%lg) + e-(%lg) -> %s(%lg) + νe(0) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV - e_massElectronMeV,
								e_massElectronMeV,
								daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
								massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - e_massElectronMeV + e_massElectronMeV - daughterIsotopePropertyPtr->massMeV, massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, decayMode, nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, - partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutputOfElectronCapture(a_electrodePtr, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_PROTON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_2PROTON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_3PROTON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_ALPHA
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_ELECTRON_CAPTURE_PROTON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_ELECTRON_CAPTURE_ALPHA){
						double zMeV;
						char * zSymbol;
						if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_PROTON
						|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_ELECTRON_CAPTURE_PROTON){
							daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 2;
							daughterMassNumber = a_isotopePropertyPtr->massNumber - 1;
							zMeV = e_hydrogenPtr->massMeV;
							zSymbol = "H";
						}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_2PROTON){
							daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 3;
							daughterMassNumber = a_isotopePropertyPtr->massNumber - 2;
							zMeV = e_hydrogenPtr->massMeV * 2;
							zSymbol = "H2";
						}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_3PROTON){
							daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 4;
							daughterMassNumber = a_isotopePropertyPtr->massNumber - 3;
							zMeV = e_hydrogenPtr->massMeV * 3;
							zSymbol = "H3";
						}else{//DECAY_MODE_BETA_PLUS_ALPHA, DECAY_MODE_ELECTRON_CAPTURE_ALPHA
							daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 3;
							daughterMassNumber = a_isotopePropertyPtr->massNumber - 4;
							zMeV = e_helium4Ptr->massMeV;
							zSymbol = "He";
						}
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massElectronMeV + e_massElectronMeV + zMeV);
							if(
							(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_PROTON
							|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_2PROTON
							|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_3PROTON
							|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_ALPHA
							)
							&& massDefectMeV > 0.0){
								//The beta plus and proton/2proton/alpha emition
								//The Nuclear Reaction is : X -> Y + e- + e+ + νe + (p+ + e-)|(p+ + e-)2|(He)
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s(%lg) -> %s(%lg) + e-(%lg) + e+(%lg) + νe(0) + %s(%lg) + %lg[MeV/c^2]", 
									a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
									daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
									e_massElectronMeV, e_massElectronMeV,
									zSymbol, zMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massElectronMeV + e_massElectronMeV + zMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, - partialMol);
								if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_2PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 2.0);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_3PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 3.0);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_ALPHA){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM, partialMol);
								}
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								//The mass defect is devided into Beta Plus and Protons, But the rate is fixed just half.
								registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * 0.5, partialMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV * 0.5, partialMol);
							}else if((massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + zMeV)) > 0.0){
								//The electron capture and proton emition
								//The Nuclear Reaction is : X+ + e- -> Y + νe + (p+ + e-)
								int decayMode = a_isotopePropertyPtr->decayMode[i];
								if(decayMode == DECAY_MODE_BETA_PLUS_PROTON){
									decayMode = DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_PROTON;
								}else if(decayMode == DECAY_MODE_BETA_PLUS_2PROTON){
									decayMode = DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON;
								}else if(decayMode == DECAY_MODE_BETA_PLUS_3PROTON){
									decayMode = DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON;
								}else if(decayMode == DECAY_MODE_BETA_PLUS_ALPHA){
									decayMode = DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA;
								}
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s+(%lg) + e-(%lg) -> %s(%lg) + νe(0) + %s(%lg) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV - e_massElectronMeV,
								e_massElectronMeV,
								daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
								zSymbol, zMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - e_massElectronMeV + e_massElectronMeV - (daughterIsotopePropertyPtr->massMeV + zMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, decayMode, nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, - partialMol);
								if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_PROTON
								|| decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol);
								}else if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 2.0);
								}else if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 3.0);
								}else if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_ALPHA
								|| decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM, partialMol);
								}
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								//The mass defect is devided into Electron Capture and Protons, But the rate is fixed just half.
								registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * 0.5, partialMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV * 0.5, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_DOUBLE_BETA_MINUS){
						//The beta minus emmits an electron and an anti-neutrino. The Nuclear Reaction is : n -> p + e- + ~νe
						// increase atomic number, same mass number
						if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS){
							daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber + 1;
						}else{//DECAY_MODE_DOUBLE_BETA_MINUS
							daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber + 2;
						}
						daughterMassNumber = a_isotopePropertyPtr->massNumber;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							char * daughterSymbolPtr; 
							if(daughterAtomicNumber == ATOMICNUMBER_HYDROGEN && daughterMassNumber == MASSNUMBER_HYDROGEN){
								daughterSymbolPtr = daughterIsotopePropertyPtr->nucleusName;
							}else{
								daughterSymbolPtr = daughterIsotopePropertyPtr->symbol;
							}
							double daughterIonMassMeV;
							char * fmt;
							if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS){
								daughterIonMassMeV = daughterIsotopePropertyPtr->massMeV - e_massElectronMeV;
							}else{//DECAY_MODE_DOUBLE_BETA_MINUS
								daughterIonMassMeV = daughterIsotopePropertyPtr->massMeV - e_massElectronMeV * 2.0;
							}
							if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS){
								massDefectMeV = a_isotopePropertyPtr->massMeV - daughterIsotopePropertyPtr->massMeV;
								fmt = "%s(%lg) -> %s+(%lg) + e-(%lg) + ~νe(0) + %lg[MeV/c^2]";
							}else{//DECAY_MODE_DOUBLE_BETA_MINUS
								massDefectMeV = a_isotopePropertyPtr->massMeV - daughterIsotopePropertyPtr->massMeV;
								fmt = "%s(%lg) -> %s2+(%lg) + e-2(%lg * 2) + ~νe2(0) + %lg[MeV/c^2]";
							}
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
								fmt, 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
								daughterSymbolPtr, daughterIonMassMeV,
								e_massElectronMeV, 
								massDefectMeV);
								if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS){
									checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIonMassMeV + e_massElectronMeV), massDefectMeV, nuclearReactionStr);
								}else{
									checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIonMassMeV + e_massElectronMeV * 2.0), massDefectMeV, nuclearReactionStr);
								}
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, partialMol);
									
								}else{//DECAY_MODE_DOUBLE_BETA_MINUS
									increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, partialMol * 2.0);									
								}
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutputOfBetaMinus(a_electrodePtr, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_NEWTRON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_2NEWTRON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_3NEWTRON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_4NEWTRON
					){
						int newtronCnt;
						char * newtronTxtPtr;
						char * scalePtr;
						switch(a_isotopePropertyPtr->decayMode[i]){
							case DECAY_MODE_BETA_MINUS_AND_NEWTRON: newtronCnt = 1; newtronTxtPtr = "n"; scalePtr = ""; break;
							case DECAY_MODE_BETA_MINUS_AND_2NEWTRON: newtronCnt = 2; newtronTxtPtr = "n2"; scalePtr = " * 2"; break;
							case DECAY_MODE_BETA_MINUS_AND_3NEWTRON: newtronCnt = 3; newtronTxtPtr = "n3"; scalePtr = " * 3"; break;
							case DECAY_MODE_BETA_MINUS_AND_4NEWTRON: newtronCnt = 4; newtronTxtPtr = "n4"; scalePtr = " * 4"; break;
						}
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber + 1;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - newtronCnt;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							char * daughterSymbolPtr; 
							if(daughterAtomicNumber == ATOMICNUMBER_HYDROGEN && daughterMassNumber == MASSNUMBER_HYDROGEN){
								daughterSymbolPtr = daughterIsotopePropertyPtr->nucleusName;
							}else{
								daughterSymbolPtr = daughterIsotopePropertyPtr->symbol;
							}
							double daughterIonMassMeV = daughterIsotopePropertyPtr->massMeV - e_massElectronMeV;
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * newtronCnt);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s(%lg) -> %s+(%lg) + e-(%lg) + ~νe(0) + %s(%lg%s) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
								daughterSymbolPtr, daughterIonMassMeV,
								e_massElectronMeV,
								newtronTxtPtr, e_massNeutronMeV, scalePtr,
								massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIonMassMeV + e_massElectronMeV + e_massNeutronMeV * newtronCnt), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, partialMol);
								//The mass defect is devided into Beta Minus and Neutrons, But the rate is fixed just half.
								registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * 0.5, partialMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV * 0.5, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol * newtronCnt);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * newtronCnt - a_isotopePropertyPtr->massMeV) * partialMol;
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_ALPHA){
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber + 1 - 2;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - 4;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							char * daughterSymbolPtr; 
							if(daughterAtomicNumber == ATOMICNUMBER_HYDROGEN && daughterMassNumber == MASSNUMBER_HYDROGEN){
								daughterSymbolPtr = daughterIsotopePropertyPtr->nucleusName;
							}else{
								daughterSymbolPtr = daughterIsotopePropertyPtr->symbol;
							}
							double daughterIonMassMeV = daughterIsotopePropertyPtr->massMeV - e_massElectronMeV;
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_helium4Ptr->massMeV);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s(%lg) -> %s+(%lg) + e-(%lg) + ~νe(0) + He(%lg) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
								daughterSymbolPtr, daughterIonMassMeV,
								e_massElectronMeV,
								e_helium4Ptr->massMeV, 
								massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIonMassMeV + e_massElectronMeV + e_helium4Ptr->massMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM, partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_helium4Ptr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								//The mass defect is devided into Beta Minus and ALPHA, But the rate is fixed just half.
								registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * 0.5, partialMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV * 0.5, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_PROTON_EMISSION){
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 1;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - 1;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_hydrogenPtr->massMeV);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + H(%lg) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										e_hydrogenPtr->massMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_hydrogenPtr->massMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_hydrogenPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_2PROTON_EMISSION){
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 2;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - 2;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_hydrogenPtr->massMeV * 2.0);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + H2(%lg * 2) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										e_hydrogenPtr->massMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_hydrogenPtr->massMeV * 2.0), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 2.0);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_hydrogenPtr->massMeV * 2.0 - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_NEWTRON_EMISSION){
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - 1;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + n(%lg) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										e_massNeutronMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_SELF_FISSION_80KR){
						fprintf(stderr, "ERROR:%s:not programmed decay mode:DECAY_MODE_SELF_FISSION_80KR\n", __FUNCTION__);
						fprintf(e_logFp, "ERROR:%s:not programmed decay mode:DECAY_MODE_SELF_FISSION_80KR\n", __FUNCTION__);
						exit(1);
					}else{
						fprintf(stderr, "ERROR:%s:unknown decay mode:%d\n", __FUNCTION__, a_isotopePropertyPtr->decayMode[i]);
						fprintf(e_logFp, "ERROR:%s:unknown decay mode:%d\n", __FUNCTION__, a_isotopePropertyPtr->decayMode[i]);
						exit(1);
					}
					if(showErrorOfUndefinedDaughter){
						struct atomProperty * atomPropertyPtr = getAtomPropertyPtr(daughterAtomicNumber);
						if(atomPropertyPtr){
							snprintf(nuclearReactionStr, REACTION_LEN, 
								"undefined daughter isotope:%d%s %s for %s[%d] on %s\n", 
								daughterMassNumber, atomPropertyPtr->symbol, atomPropertyPtr->name,
								getDecayModeText(a_isotopePropertyPtr->decayMode[i]), i,
								a_isotopePropertyPtr->symbol);
						}else{
							snprintf(nuclearReactionStr, REACTION_LEN, 
								"undefined daughter isotope:AtomicNumber:%d MassNumber:%d for %s[%d] on %s\n", 
								daughterAtomicNumber, daughterMassNumber,
								getDecayModeText(a_isotopePropertyPtr->decayMode[i]), i,
								a_isotopePropertyPtr->symbol);
						}
						registNuclearReaction(a_electrodePtr, REACTION_ERROR, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, 0.0, partialMol, 0.0, 0.0);
					}else if(showErrorNegativeMassDefect){
						snprintf(nuclearReactionStr, REACTION_LEN, "negative mass defect %lg[MeV/c^2] for %s[%d] on %s\n",
								massDefectMeV,
								getDecayModeText(a_isotopePropertyPtr->decayMode[i]), i,
								a_isotopePropertyPtr->symbol);
						registNuclearReaction(a_electrodePtr, REACTION_ERROR, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
					}else{
						setMolSub(a_atomValuePtr, partialMol);
					}
				}
			}
		}else{
			a_atomValuePtr->pastSecond += a_pastSecond;
		}
	}
}
extern void absorbeOrDecayNeutronInElectrode(struct electrode * a_electrodePtr, double a_pastSecond)
{
	struct atomNodeConst * neutronPtr;
	//fprintf(stderr, "DEBUG:%s:{BEGIN\n", __FUNCTION__);
	neutronPtr = findNeutronInElectrode(a_electrodePtr);
	if(neutronPtr){
		//There is another assumption about the behavior of neutrons. We know that the half life time of an alone neutron is 886.7 [sec] and the standard gravity of earth is 9.80665 [m/s2]. The assumption is that almost all neutrons is absorbed in nucleus while the free fall by the gravity during one second.
		double neutronPastSecond = 1.0;
		double neutronVoltageMeV = 0.0;
		double neutronMol;
		double collideNeutronRate = 1.0;
		if(a_pastSecond < 1.0){
			neutronPastSecond = a_pastSecond;
		}
		decayNeuclay(a_electrodePtr, neutronPtr->key->isotopePropertyPtr, neutronPtr->vPtr, neutronPastSecond);
		
		neutronMol = getMol(neutronPtr->vPtr);//get mol of neutrons after decay.
		collideBulletToElectrode(a_electrodePtr, COLLIDE_NEUTRON, neutronVoltageMeV, neutronMol, collideNeutronRate, neutronPtr);
	}
	//fprintf(stderr, "DEBUG:%s:}END\n", __FUNCTION__);
}
#if  KEEP_ATMIC_ORDER == 1
#else
struct decayBase {
	struct electrode * electrodePtr;
	double pastSecond;
};
extern int iterateDecayWithoutNewtron(void * a_total, struct objectNodeConst * a_nodePtr)
{
	struct decayBase * pDb = (struct decayBase *)a_total;
	const struct isotopeProperty * isotopePropertyPtr = ((const struct atomKey *)a_nodePtr->key)->isotopePropertyPtr;
	struct atomValue * atomValuePtr = (struct atomValue *)a_nodePtr->value;
	if(isotopePropertyPtr->atomicNumber != ATOMICNUMBER_NEUTRON
	|| isotopePropertyPtr->massNumber != MASSNUMBER_NEUTRON){
		decayNeuclay(pDb->electrodePtr, isotopePropertyPtr, atomValuePtr, pDb->pastSecond);
		
	}
	return KEEP_NODE;
}
#endif
extern void decayNeuclayWithoutNewtron(struct electrode * a_electrodePtr, double a_pastSecond)
{
#if  KEEP_ATMIC_ORDER == 1
	//fprintf(stderr, "DEBUG:%s:KEEP_ATMIC_ORDER == 1\n", __FUNCTION__);
	{
		unsigned int size, i;
		struct objectNodeConst ** oPtr;
		oPtr = getFlatTable(&a_electrodePtr->atomHashTable, SORT_ASCEND, &size);
		if(oPtr){
			for(i = 0; i < size; ++i){
				const struct isotopeProperty * isotopePropertyPtr = ((struct atomKey *)(oPtr[i]->keyPtr))->isotopePropertyPtr;
				struct atomValue * atomValuePtr = (struct atomValue *)oPtr[i]->valuePtr;
				//fprintf(stderr, "DEBUG:%s:decayNeuclay i:%d\n", __FUNCTION__, i);
				if(isotopePropertyPtr->atomicNumber != ATOMICNUMBER_NEUTRON
				|| isotopePropertyPtr->massNumber != MASSNUMBER_NEUTRON){
					decayNeuclay(a_electrodePtr, isotopePropertyPtr, atomValuePtr, a_pastSecond);
				}
			}
			free(oPtr);
		}
	}
#else
	//fprintf(stderr, "DEBUG:%s:KEEP_ATMIC_ORDER != 1\n", __FUNCTION__);
	{
		struct decayBase Db;
		Db.electrodePtr = a_electrodePtr;
		Db.pastSecond = a_pastSecond;
		//fprintf(stderr, "DEBUG:%s:iterateInHashTable\n", __FUNCTION__);
		iterateInHashTable(&a_electrodePtr->atomHashTable, &Db, iterateDecayWithoutNewtron);
	}
#endif	
}

//---------------------------------------------------------------------
#define SCATTER_CNT ((sizeof(e_cosTable) / sizeof(double)) - 1) //
double e_cosTable[7];//1.0 - cos(M_PI * 0 / 6), 1.0 - cos(M_PI * 1 / 6), 1.0 - cos(M_PI * 2 / 6), ... 
extern void initCosTable()
{
	int i;
	for(i = 1; i < (sizeof(e_cosTable) / sizeof(double)); ++i){
		e_cosTable[i] = 1.0 - cos((M_PI / 6.0) * i);
		//fprintf(stderr, "DEBUG:%s:e_cosTable[%d]%lg\n", __FUNCTION__, i, e_cosTable[i]);
	}
}
extern double calcScatteredPhotonMeV(double a_incidentPhotonMeV, int a_angleIndex)
{
	//[PHYSICS REVIEW] The Passed Energy of Compton scattering
	//https://docs.google.com/document/d/1ZmPn4N57MOAG2C02d_nFATG-t7OUhTOJm-7tQdR2cZA/edit?usp=sharing
	// We use the formula number (p).
	
	//Do not use a_angleIndex = 0, it means not collide.
	double scatteredPhotonMeV = (a_incidentPhotonMeV * e_massElectronMeV)
		/ (e_cosTable[a_angleIndex] * a_incidentPhotonMeV + e_massElectronMeV);
	if(scatteredPhotonMeV < 0.0){
		fprintf(stderr, "ERROR:%s:negetive:%lg\n", __FUNCTION__, scatteredPhotonMeV);
		scatteredPhotonMeV = 0.0;
	}
	if(scatteredPhotonMeV > a_incidentPhotonMeV){
		fprintf(stderr, "ERROR:%s:over:%lg > %lg\n", __FUNCTION__, scatteredPhotonMeV, a_incidentPhotonMeV);
		scatteredPhotonMeV = a_incidentPhotonMeV;
	}
	//fprintf(stderr, "DEBUG:%s:SCATTER_CNT:%d scattered Photon rate:%lg = %lg / %lg\n", __FUNCTION__, SCATTER_CNT, scatteredPhotonMeV / a_incidentPhotonMeV, scatteredPhotonMeV, a_incidentPhotonMeV);
	return scatteredPhotonMeV;
}
//int e_cntScatterPhoton = 0;//DEBUG
//int e_cntScatterPhotonUnnest = 0;//DEBUG
//int e_cntScatterPhoton_s = 0;//DEBUG
//int e_cntScatterPhotonUnnest_s = 0;//DEBUG
extern void scatterPhoton(struct electrode * a_electrodePtr, double a_MeV, double a_mol, int a_nest)
{
	//if(a_nest == 0){
	//	++e_cntScatterPhotonUnnest;//DEBUG
	//	++e_cntScatterPhotonUnnest_s;//DEBUG
	//}
	++a_nest;
	//++e_cntScatterPhoton;//DEBUG
	//++e_cntScatterPhoton_s;//DEBUG
	//fprintf(stderr, "DEBUG:%s:{BEGIN N:%d a_MeV:%lg a_mol:%lg\n", __FUNCTION__, a_nest, a_MeV, a_mol);
	{
		static double s_molLimit = 0.1;
		if(a_mol > s_molLimit){
			fprintf(stderr, "INFO:%s(%d):electrode %s, a_mol:%lg > %lg\n", __FUNCTION__, __LINE__, a_electrodePtr->atomHashTable.tableName, a_mol, s_molLimit);
			fprintf(e_logFp, "INFO:%s(%d):electrode %s, a_mol:%lg > %lg\n", __FUNCTION__, __LINE__, a_electrodePtr->atomHashTable.tableName, a_mol, s_molLimit);
			s_molLimit *= 10.0;
			//exit(1);
		}
	}
	//The next if-statement should be checked out of this function. 
	//if(a_MeV >= e_collideElectronMiniMeV && a_mol >= a_electrodePtr->detectLimitMolForIsotope * SCATTER_CNT && a_nest < 16)
	{
		//double debugCheck = 0.0, debugRate;
		struct atomNodeConst * electronPtr;
		double molDivideBy6 = a_mol / SCATTER_CNT;
		int i;
		electronPtr = findElectronInElectrode(a_electrodePtr);
		for(i = 1; i < (sizeof(e_cosTable) / sizeof(double)); ++i){
			double scatteredPhotonMeV = calcScatteredPhotonMeV(a_MeV, i);
			double scatteredElectronMeV = a_MeV - scatteredPhotonMeV;
			//fprintf(stderr, "DEBUG:%s:N:%d SA:%d scattered ElectronMeV:%lg(%lg) PhotonMeV:%lg(%lg) molDivideBy6:%lg\n", 
			//	__FUNCTION__, a_nest, i * 30, 
			//	scatteredElectronMeV, scatteredElectronMeV / a_MeV,
			//	scatteredPhotonMeV, scatteredPhotonMeV / a_MeV,
			//	molDivideBy6);
			if(scatteredElectronMeV >= e_appliedVoltageMeV){		
				//fprintf(stderr, " LINE:%d{ \n", __LINE__);
				collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_C, scatteredElectronMeV, molDivideBy6, e_collideElectronRateOnElectrode, electronPtr);
				//fprintf(stderr, " LINE:%d} \n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}else if(scatteredElectronMeV >= e_collideElectronMidiMeV){				
				//fprintf(stderr, " LINE:%d{ \n", __LINE__);
				collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_C, scatteredElectronMeV, molDivideBy6, e_collideElectronRateForMidiMeV, electronPtr);
				electronPtr = findElectronInElectrode(a_electrodePtr);
				//fprintf(stderr, " LINE:%d} \n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}else if(scatteredElectronMeV >= e_collideElectronMiniMeV){				
				//fprintf(stderr, " LINE:%d{ \n", __LINE__);
				collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_C, scatteredElectronMeV, molDivideBy6, e_collideElectronRateForMiniMeV, electronPtr);
				electronPtr = findElectronInElectrode(a_electrodePtr);
				//fprintf(stderr, " LINE:%d} \n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}else if(scatteredElectronMeV > 0.0){
				//fprintf(stderr, " LINE:%d{ registOutput SCAT_SMALL_BY_COMPTON_E call\n", __LINE__);
				registOutput(a_electrodePtr, SCAT_SMALL_BY_COMPTON_E, scatteredElectronMeV, molDivideBy6);
				//fprintf(stderr, " LINE:%d} registOutput SCAT_SMALL_BY_COMPTON_E return\n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}
			if(scatteredPhotonMeV >= e_collideElectronMiniMeV){
				registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, scatteredPhotonMeV, molDivideBy6);	
			}else if(scatteredPhotonMeV > 0.0){
				//fprintf(stderr, " LINE:%d{ registOutput SCAT_SMALL_BY_COMPTON call\n", __LINE__);
				registOutput(a_electrodePtr, SCAT_SMALL_BY_COMPTON, scatteredPhotonMeV, molDivideBy6);
				//fprintf(stderr, " LINE:%d} registOutput SCAT_SMALL_BY_COMPTON return\n", __LINE__);
				//debugCheck += (scatteredPhotonMeV * molDivideBy6);
			}
		}
		//debugRate = debugCheck / (a_MeV * a_mol);
		//if(debugRate < 0.9999 || 1.0001 < debugRate){
		//	fprintf(stderr, "FATAL ERROR:%s:debugRate:%lg = debugCheck:%lg /(a_MeV:%lg * a_mol:%lg)\n", __FUNCTION__, debugRate, debugCheck, a_MeV, a_mol);
		//	exit(1);
		//}
	}
	//fprintf(stderr, "DEBUG:%s:}END\n", __FUNCTION__);
}
//---------------------------------------------------------------------
#define DEGRADE_MAG 3
#define DEGRADE_DIV 16

struct cleanUp {
	//{ set by comptonEffect
	struct electrode * electrodePtr;
	double detectLimitMolForIsotopeSC;
	//}
	//{ set by iterateClenupMassDefect
	int cntOfBigMolNode;
	double maxChangeMeVOfBigMol;
	double minChangeMeVOfBigMol;
	int cntOfSmallMolNode;
	double totalMeVMol, degradeMeVMol;//DEBUG
	//}
	//{ set by initMaxMinOfClenup
	double maxMeV[DEGRADE_MAG], minMeV[DEGRADE_MAG], step[DEGRADE_MAG];
	struct massDefect md[DEGRADE_MAG][DEGRADE_DIV];
	//}
};
extern int iterateDegreadeMassDefect(void * a_total, struct objectNodeConst * a_nodePtr)
{
	int action = KEEP_NODE;
	struct cleanUp * cleanUpPtr = (struct cleanUp *)a_total;
	struct massDefect * valuePtr = (struct massDefect *)a_nodePtr->valuePtr;
	if(valuePtr->mol <= 0.0){
		action = FREE_NODE;
		//fprintf(stderr, "DEBUG:%s:FREE_NODE\n", __FUNCTION__);
	}else{
		int j, k;
		double val, rem;
		double debugV = 0.0, debugR;
		//fprintf(stderr, "DEBUG:%s: a_MeV:%lg valuePtr->mol:%lg\n", __FUNCTION__, valuePtr->MeV, valuePtr->mol);
		val = valuePtr->MeV;
		for(j = 0; j < DEGRADE_MAG; ++j){
			if(val >= cleanUpPtr->minMeV[j]){
				k = (val - cleanUpPtr->minMeV[j]) / cleanUpPtr->step[j];
				if(k >= DEGRADE_DIV && j == 0 && valuePtr->mol < cleanUpPtr->detectLimitMolForIsotopeSC){
					//regist the part of overflowing for the small amount of 'mol'
					int tm;
					k = DEGRADE_DIV - 1;
					tm = (int)(val / cleanUpPtr->md[j][k].MeV);
					cleanUpPtr->md[j][k].mol += (valuePtr->mol * tm);
					cleanUpPtr->degradeMeVMol += (cleanUpPtr->md[j][k].MeV * valuePtr->mol * tm);//DEBUG
					debugV += (cleanUpPtr->md[j][k].MeV * tm);//DEBUG
					rem = val - (cleanUpPtr->md[j][k].MeV * tm);
					val = rem;
					//re-calculation
					if(val >= cleanUpPtr->minMeV[j]){
						k = (val - cleanUpPtr->minMeV[j]) / cleanUpPtr->step[j];
					}else{
						continue;
					}
				}
				if(k < DEGRADE_DIV){
					cleanUpPtr->md[j][k].mol += valuePtr->mol;
					cleanUpPtr->degradeMeVMol += (cleanUpPtr->md[j][k].MeV * valuePtr->mol);//DEBUG
					debugV += cleanUpPtr->md[j][k].MeV;//DEBUG
					rem = val - (cleanUpPtr->step[j] * k + cleanUpPtr->minMeV[j]);
					val = rem;
				}else{
					fprintf(stderr, "ERROR:%s(%d): k:%d = (val:%lg - cleanUpPtr->minMeV[%d]:%lg) / cleanUpPtr->step[%d]:%lg, cleanUpPtr->maxMeV[%d]:%lg cleanUpPtr->md[%d][DEGRADE_DIV - 1].MeV:%lg\n", __FUNCTION__, __LINE__, k, val, j, cleanUpPtr->minMeV[j], j, cleanUpPtr->step[j], j, cleanUpPtr->maxMeV[j], j, cleanUpPtr->md[j][DEGRADE_DIV - 1].MeV);
					exit(1);
				}
			}
		}
		if(val > cleanUpPtr->minMeV[DEGRADE_MAG - 1]){//DEBUG
			fprintf(stderr, "ERROR:%s: val:%lg > cleanUpPtr->minMeV[DEGRADE_MAG - 1]:%lg\n", __FUNCTION__, val, cleanUpPtr->minMeV[DEGRADE_MAG - 1]);
			exit(1);
		}
		debugV += val;//DEBUG
		debugR = (debugV - valuePtr->MeV) / valuePtr->MeV;//DEBUG
		if(abs((int)(debugR * 10000)) > 0){//DEBUG
			fprintf(stderr, "ERROR:%s:valuePtr->MeV:%lg debugR:%lg\n", __FUNCTION__, valuePtr->MeV, debugR);
			exit(1);
		}
		if(val > 0.0){
			int scat_degrade;
			if(valuePtr->mol < cleanUpPtr->detectLimitMolForIsotopeSC){
				scat_degrade = SCAT_DEGRADE_IN_COMPTON_S;
			}else{
				scat_degrade = SCAT_DEGRADE_IN_COMPTON_B;
			}
			registOutput(cleanUpPtr->electrodePtr, scat_degrade, val, valuePtr->mol);
			cleanUpPtr->degradeMeVMol += (val * valuePtr->mol);//DEBUG
		}
		valuePtr->mol = 0.0;
		action = FREE_NODE;
	}
	return action;
}
extern void moveFromClenup(struct cleanUp * a_cleanUpPtr)
{
	int j, k;

	for(j = 0; j < DEGRADE_MAG; ++j){
		for(k = 0; k < DEGRADE_DIV; ++k){
			if(a_cleanUpPtr->md[j][k].mol > 0.0){
				registOutput(a_cleanUpPtr->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, a_cleanUpPtr->md[j][k].MeV, a_cleanUpPtr->md[j][k].mol);
			}
		}
	}
}
extern void initMaxMinOfClenup(struct cleanUp * a_cleanUpPtr)
{
	int j, k;

	a_cleanUpPtr->step[0] = (a_cleanUpPtr->maxMeV[0] - a_cleanUpPtr->minMeV[2]) / DEGRADE_DIV;
	a_cleanUpPtr->minMeV[0] = a_cleanUpPtr->minMeV[2] + a_cleanUpPtr->step[0];
	a_cleanUpPtr->step[1] = (a_cleanUpPtr->minMeV[0] - a_cleanUpPtr->minMeV[2]) / (DEGRADE_DIV + 1);
	a_cleanUpPtr->maxMeV[1] = a_cleanUpPtr->minMeV[0] - a_cleanUpPtr->step[1];
	a_cleanUpPtr->minMeV[1] = a_cleanUpPtr->minMeV[2] + a_cleanUpPtr->step[1];
	a_cleanUpPtr->step[2] = (a_cleanUpPtr->minMeV[1] - a_cleanUpPtr->minMeV[2]) / DEGRADE_DIV;
	a_cleanUpPtr->maxMeV[2] = a_cleanUpPtr->minMeV[1] - a_cleanUpPtr->step[2];
	
	for(j = 0; j < DEGRADE_MAG; ++j){
		for(k = 0; k < DEGRADE_DIV; ++k){
			a_cleanUpPtr->md[j][k].MeV = (a_cleanUpPtr->maxMeV[j] * k + a_cleanUpPtr->minMeV[j] * (DEGRADE_DIV - 1 - k)) / (DEGRADE_DIV - 1);
			a_cleanUpPtr->md[j][k].mol = 0.0;
		}
	}
}
extern int iterateClenupMassDefect(void * a_total, struct objectNodeConst * a_nodePtr)
{
	int action = KEEP_NODE;
	struct cleanUp * cleanUpPtr = (struct cleanUp *)a_total;
	struct massDefect * valuePtr = (struct massDefect *)a_nodePtr->valuePtr;
	if(valuePtr->mol <= 0.0){
		action = FREE_NODE;
		//fprintf(stderr, "DEBUG:%s:FREE_NODE\n", __FUNCTION__);
	}else{
		cleanUpPtr->totalMeVMol += (valuePtr->MeV * valuePtr->mol);//DEBUG
		if(valuePtr->mol < cleanUpPtr->detectLimitMolForIsotopeSC){
			cleanUpPtr->cntOfSmallMolNode++;
		}else{
			if(cleanUpPtr->cntOfBigMolNode == 0){
				cleanUpPtr->maxMeV[0] = valuePtr->MeV;
				cleanUpPtr->minMeV[2] = e_collideElectronMiniMeV;
			}else{
				if(cleanUpPtr->maxMeV[0] < valuePtr->MeV){
					cleanUpPtr->maxMeV[0] = valuePtr->MeV;
				}
			}
			cleanUpPtr->cntOfBigMolNode++;
		}
	}
	return action;
}

extern double calcPastTime(struct timeval * tv2, const struct timeval * tv1)
{
	gettimeofday(tv2, NULL);
	double t1 = (double)tv1->tv_sec + (double)tv1->tv_usec / 1000000.0;
	double t2 = (double)tv2->tv_sec + (double)tv2->tv_usec / 1000000.0;
	double diff = t2 - t1;
	return diff;
}

extern void comptonEffect(struct electrode * a_electrodePtr)
{
	int doProcess;
	struct cleanUp cUp;
	//struct timeval tv1, tv2;//DEBUG
    //gettimeofday(&tv1, NULL);//DEBUG
	//e_inComptonEffect = 1;//DEBUG
	//fprintf(stderr, "DEBUG:%s:{BEGIN %s usedCnt:%d\n", __FUNCTION__, a_electrodePtr->atomHashTable.tableName, a_electrodePtr->atomHashTable.usedCnt);	
	cUp.electrodePtr = a_electrodePtr;
	cUp.detectLimitMolForIsotopeSC = a_electrodePtr->detectLimitMolForIsotope * SCATTER_CNT;
	for(doProcess = 0; doProcess < 10; ++doProcess){
		unsigned int size, i, j = 0;
		struct objectNodeConst ** oPtr;
		
		oPtr = getFlatTable(&a_electrodePtr->massDefectHashTable, SORT_ASCEND, &size);
		//e_cntScatterPhotonUnnest = 0;//DEBUG
		//e_cntScatterPhoton = 0;//DEBUG
		//fprintf(stderr, "DEBUG:%s:doProcess:%d getFlatTable size:%d\n", __FUNCTION__, doProcess, size);
		//fprintf(stderr, "DEBUG:%s:calcPastTime-1:%lg[sec]\n", __FUNCTION__, calcPastTime(&tv2, &tv1));
		if(oPtr){
			struct massDefect * valuePtr;

			for(i = 0; i < size; ++i){
				valuePtr = (struct massDefect *)oPtr[i]->valuePtr;
				if(valuePtr->MeV >= e_collideElectronMiniMeV){
					if(valuePtr->mol >= cUp.detectLimitMolForIsotopeSC){
						int nest = 0;
						double saveMol = valuePtr->mol;
						valuePtr->mol = 0.0;
						if(e_useComptonEffect){
						//fprintf(stderr, "DEBUG:%s:loop:%lg[sec] i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, calcPastTime(&tv2, &tv1), i, valuePtr->MeV, saveMol);
						//fprintf(stderr, "DEBUG:%s:loop:i:%d tableName:%s usedCnt:%d\n", __FUNCTION__, i, a_electrodePtr->atomHashTable.tableName, a_electrodePtr->atomHashTable.usedCnt);
						//fprintf(stderr, "DEBUG:%s:e_cntScatterPhotonUnnest:%d e_cntScatterPhoton:%d\n", __FUNCTION__, e_cntScatterPhotonUnnest, e_cntScatterPhoton);
						//fprintf(stderr, "DEBUG:%s:massDefectHashTable.usedCnt:%u\n", __FUNCTION__, a_electrodePtr->massDefectHashTable.usedCnt);
						//e_cntScatterPhotonUnnest_s = 0;//DEBUG
						//e_cntScatterPhoton_s = 0;//DEBUG
						scatterPhoton(a_electrodePtr, valuePtr->MeV, saveMol, nest);
						//fprintf(stderr, "DEBUG:%s:e_cntScatterPhotonUnnes_s:%d e_cntScatterPhoton_s:%d\n", __FUNCTION__, e_cntScatterPhotonUnnest_s, e_cntScatterPhoton_s);
						//fprintf(stderr, "DEBUG:%s:massDefectHashTable.usedCnt:%u\n", __FUNCTION__, a_electrodePtr->massDefectHashTable.usedCnt);
						/*
						[!!CAUTION!!!] We can't use "scatterPhoton" as the argument of "iterateInHashTable", 
						because "scatterPhoton" will uses "insertObjectInHashTable" or "iterateInHashTable" through nesting calls.
						*/
						}else{
							registOutput(a_electrodePtr, SCAT_DEGRADE_IN_COMPTON_B, valuePtr->MeV, saveMol);
						}
						++j;
					}
				}else{
					fprintf(stderr, "ERROR:%s:%s valuePtr->MeV:%lg < e_collideElectronMiniMeV:%lg\n", __FUNCTION__, 
					a_electrodePtr->massDefectHashTable.tableName, valuePtr->MeV, e_collideElectronMiniMeV);
				}
			}
			free(oPtr);
			//fprintf(stderr, "DEBUG:%s:calcPastTime-2:%lg[sec] scatterPhoton * size:%d\n", __FUNCTION__, calcPastTime(&tv2, &tv1), size);
		}
		//fprintf(stderr, "DEBUG:%s:calcPastTime-3:%lg[sec] scatterPhoton * size:%d\n", __FUNCTION__, calcPastTime(&tv2, &tv1), size);
		//fprintf(stderr, "DEBUG:%s:e_cntScatterPhotonUnnest:%d e_cntScatterPhoton:%d\n", __FUNCTION__, e_cntScatterPhotonUnnest, e_cntScatterPhoton);
		//fprintf(stderr, "DEBUG:%s:doProcess:%d size:%u j:%u massDefectHashTable.usedCnt:%u\n", __FUNCTION__, doProcess, size, j, a_electrodePtr->massDefectHashTable.usedCnt);
		if(j > 0){
			//It's better to clean up "massDefectHashTable" for both small usage of memory and faster processing.
			cUp.cntOfBigMolNode = 0;
			cUp.maxChangeMeVOfBigMol = 0.0;
			cUp.minChangeMeVOfBigMol = 0.0;
			cUp.cntOfSmallMolNode = 0;
			cUp.totalMeVMol = 0.0;
			cUp.degradeMeVMol = 0.0;
			iterateInHashTable(&a_electrodePtr->massDefectHashTable, &cUp, iterateClenupMassDefect);
			initMaxMinOfClenup(&cUp);
			//I recommend the degrading process of SCAT_BIG_MASS_DEFECT_NOW with using the aproximatic calculation for faster calculation.
			//fprintf(stderr, "DEBUG:%s:doProcess:%d cUp.cntOfBigMolNode:%d cUp.maxChangeMeVOfBigMol:%lg cUp.minChangeMeVOfBigMol:%lg\n", __FUNCTION__, doProcess, cUp.cntOfBigMolNode, cUp.maxChangeMeVOfBigMol, cUp.minChangeMeVOfBigMol);
			//fprintf(stderr, "DEBUG:%s:doProcess:%d cUp.cntOfSmallMolNode:%d\n", __FUNCTION__, doProcess, cUp.cntOfSmallMolNode);
			if(cUp.cntOfBigMolNode > DEGRADE_MAG * DEGRADE_DIV || cUp.cntOfSmallMolNode > DEGRADE_MAG * DEGRADE_DIV * 8){
				iterateInHashTable(&a_electrodePtr->massDefectHashTable, &cUp, iterateDegreadeMassDefect);
				if(abs((int)(10000 * (cUp.degradeMeVMol - cUp.totalMeVMol) / cUp.totalMeVMol)) > 1){
					fprintf(stderr, "ERROR:%s(%d):cUp.degradeMeVMol:%lf != cUp.totalMeVMol:%lf\n", __FUNCTION__, __LINE__, cUp.degradeMeVMol, cUp.totalMeVMol);
					exit(1);
				}
				moveFromClenup(&cUp);
			}
		}
		//fprintf(stderr, "DEBUG:%s:massDefectHashTable.usedCnt:%u j:%d\n", __FUNCTION__, a_electrodePtr->massDefectHashTable.usedCnt, j);

		if(j == 0){
			break;
		}
	}
	//fprintf(stderr, "DEBUG:%s:tableName:%s usedCnt:%d\n", __FUNCTION__, a_electrodePtr->atomHashTable.tableName, a_electrodePtr->atomHashTable.usedCnt);
	//fprintf(stderr, "DEBUG:%s:}END calcPastTime-E:%lg[sec]\n", __FUNCTION__, calcPastTime(&tv2, &tv1));	
	//e_inComptonEffect = 0;//DEBUG
}
extern double calcDeuteriumRate(struct electrode * a_electrodePtr, double * a_tritiumRatePtr, double * a_hydrogenRatePtr)
{
	double hydrogenMol;
	double deuteriumMol;
	double tritiumMol;
	double totalMol;
	double deuteriumRate;
	struct atomNodeConst * h = findAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN);
	struct atomNodeConst * d = findAtom(a_electrodePtr, ATOMICNUMBER_DEUTERIUM, MASSNUMBER_DEUTERIUM);
	struct atomNodeConst * t = findAtom(a_electrodePtr, ATOMICNUMBER_TRITIUM, MASSNUMBER_TRITIUM);
	if(h){
		hydrogenMol = getMol(h->vPtr); 
	}else{
		hydrogenMol = 0.0;
	}
	if(d){
		deuteriumMol = getMol(d->vPtr); 
	}else{
		deuteriumMol = 0.0;
	}
	if(t){
		tritiumMol = getMol(t->vPtr);
	}else{
		tritiumMol = 0.0;
	}
	totalMol = hydrogenMol + deuteriumMol + tritiumMol;
	if(totalMol > 0.0){
		deuteriumRate = deuteriumMol / totalMol;
		if(deuteriumRate > 1.0){
			deuteriumRate = 1.0;
		}
		if(deuteriumRate < 0.0){
			deuteriumRate = 0.0;
		}
		*a_tritiumRatePtr = tritiumMol / totalMol;
		if(*a_tritiumRatePtr > 1.0){
			*a_tritiumRatePtr = 1.0;
		}
		if(*a_tritiumRatePtr < 0.0){
			*a_tritiumRatePtr = 0.0;
		}
		*a_hydrogenRatePtr = hydrogenMol / totalMol;
		if(*a_hydrogenRatePtr > 1.0){
			*a_hydrogenRatePtr = 1.0;
		}
		if(*a_hydrogenRatePtr < 0.0){
			*a_hydrogenRatePtr = 0.0;
		}
	}else{
		deuteriumRate = 0.0;
		*a_tritiumRatePtr = 0.0;
		*a_hydrogenRatePtr = 0.0;
	}
	return deuteriumRate;
}

extern void pulseCurrent(double a_pastSecond)
{
	//There is a negative electrode over the positive electrode.
	//upper[e_negative Electrode]==>(electron)==>|collide|<==(proton)<==[e_positive Electrode]lower
	int i;
	double protonRate, deuteriumRate, tritiumRate;
	struct atomNodeConst * positiveElectronPtr;
	struct atomNodeConst * negativeElectronPtr;
	positiveElectronPtr = findElectronInElectrode(&e_positiveElectrode);	
	negativeElectronPtr = findElectronInElectrode(&e_negativeElectrode);	
	
	//fprintf(stderr, "\nDEBUG:%s:{BEGINa_pastSecond:%lg\n", __FUNCTION__, a_pastSecond);
	//While prcessing, The rates of proton, deuterium and tritium will be changing step by step.
	deuteriumRate = calcDeuteriumRate(&e_positiveElectrode, &tritiumRate, &protonRate);
	//fprintf(stderr, "DEBUG:%s:positiveElectrode protonRate:%lg deuteriumRate:%lg tritiumRate:%lg\n", __FUNCTION__, protonRate, deuteriumRate, tritiumRate);
	
	//Both some of the protons emitted from the positive electrode and some of the electrons emitted from the negative electrode will collide in the intermediate spece of electrodes. They will tranform into neutrons with enough kinetic energy greater than the beta energy.
	//Then almost all of neutrons will fall on the positive electrode. The remain of neutrons will fly to the negative electrode. 
	generateNeutronInSpace(a_pastSecond, MASSNUMBER_HYDROGEN, protonRate);
	generateNeutronInSpace(a_pastSecond, MASSNUMBER_DEUTERIUM, deuteriumRate);
	generateNeutronInSpace(a_pastSecond, MASSNUMBER_TRITIUM, tritiumRate);

	//And more, both some of the protons emitted from the positive electrode and some of the electrons emitted from the negative electrode will collide in the intermediate spece of electrodes. They will be plasma, next they lost kinetic energy step by step with un-perfect collision, finally they will tranform into hydrogen atoms with kinetic energy smaller than the beta energy.
	//The total kinetic energy of hydrogen plasma will heat both the negative electrode and the positive electrode.
	//Here, we simulate that the kinetic energy losting step by step is the half of e_appliedVoltageMeV. It is that the steps is simulated only two.

	for(i = 0; i < e_stepsOfLostingEnergy; ++i){
		//fprintf(stderr, "\nDEBUG:%s:i:%d e_stepByStepLostingEnergyMeV:%lg e_genedHydrogenInSpaceMol:%lg\n", __FUNCTION__, i, e_stepByStepLostingEnergyMeV, e_genedHydrogenInSpaceMol);
		registOutput(&e_negativeElectrode, SCAT_NOT_COLLIDE, e_stepByStepLostingEnergyMeV, a_pastSecond * e_genedHydrogenInSpaceMol * 0.5);
		registOutput(&e_positiveElectrode, SCAT_NOT_COLLIDE, e_stepByStepLostingEnergyMeV, a_pastSecond * e_genedHydrogenInSpaceMol * 0.5);
	}

	collideBulletToElectrode(&e_positiveElectrode, COLLIDE_ELECTRON, e_appliedVoltageMeV, 
		a_pastSecond * e_arrivedElectronMol, e_collideElectronRateOnElectrode, positiveElectronPtr);

	collideBulletToElectrode(&e_negativeElectrode, COLLIDE_PROTON, e_appliedVoltageMeV, 
		a_pastSecond * e_arrivedProtonMol * protonRate, e_collideProtonRateOnElectrode, negativeElectronPtr);
	collideBulletToElectrode(&e_negativeElectrode, COLLIDE_DEUTERIUM, e_appliedVoltageMeV,
		a_pastSecond * e_arrivedProtonMol * deuteriumRate, e_collideProtonRateOnElectrode, negativeElectronPtr);
	collideBulletToElectrode(&e_negativeElectrode, COLLIDE_TRITIUM, e_appliedVoltageMeV,
		a_pastSecond * e_arrivedProtonMol * tritiumRate, e_collideProtonRateOnElectrode, negativeElectronPtr);

	//fprintf(stderr, "DEBUG:%s:}END\n\n", __FUNCTION__);
}


//---------------------------------------------------------------------

#define VERSION_OF_serializeRunningCondition 1
#define SERIALIZE_RUNNING_CONDITION(RW) \
extern int RW ## RunningCondition(FILE * a_fp) \
{\
	SERIALIZE_VERSION_CHECK(RW, VERSION_OF_serializeRunningCondition, a_fp) \
	SERIALIZE_VALUE(RW, e_rc.startTime, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.termTime, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.intervalTime, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.loopBeg, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.loopEnd, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.logStep, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.logNecleusReactions, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.dataStep, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.inputDataFname, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.outputDataFname, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.logFname, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.paramFname, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.inputTotalEnergy, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.simulatedLastTime, a_fp)\
	return 1;\
}
#define CALL_SERIALIZE_RUNNING_CONDITION(RW, FP) \
	RW ## RunningCondition(FP)
SERIALIZE_RUNNING_CONDITION(fread)//freadRunningCondition
SERIALIZE_RUNNING_CONDITION(fwrite)//fwriteRunningCondition

#define VERSION_OF_serializeALL 1

#define SERIALIZE_ALL(RW) \
extern int RW ## All(const char * a_fileName)\
{\
	int ret = 0;\
	FILE * fp = NULL;\
	fp = fopen(a_fileName, ((void *)RW == (void *)fread) ? "r" : "w");\
	if(fp){\
		ret = 1;\
		SERIALIZE_VERSION_CHECK(RW, VERSION_OF_serializeALL, fp)\
		if(ret){ret = CALL_SERIALIZE_RUNNING_CONDITION(RW, fp);}\
		if(ret){ret = CALL_SERIALIZE_USER_CONDITIONS(RW, fp);}\
		if(ret){ret = CALL_SERIALIZE_UNREGISTED_ATOM(RW, fp);}\
		if(ret){ret = CALL_SERIALIZE_UNREGISTED_ISOTOPE(RW, fp);}\
		if(ret){ret = CALL_SERIALIZE_ELECTRODE(RW, &e_negativeElectrode, fp);}\
		if(ret){ret = CALL_SERIALIZE_ELECTRODE(RW, &e_positiveElectrode, fp);}\
		if(ret){\
			if((void *)RW == (void *)fread){\
				ret = freadHashTable(&e_nuclearReaction, fp,\
						allocFreadNuclearReactionKey, allocFreadNuclearReactionValue, \
						calcHashSeedOfNuclearReactionKey,\
						allocCopyNuclearReactionKey,\
						freeCopyNuclearReactionKey,\
						compareNuclearReactionKey, \
						allocCopyNuclearReactionValue,\
						free,\
						foundActionForNuclearReactionValue);\
			}else if((void *)RW == (void *)fwrite){\
				ret = fwritedHashTable(&e_nuclearReaction, fp, \
						fwriteNuclearReactionKey, fwriteNuclearReactionValue);\
			}\
		}\
		if(fclose(fp) != 0){ret = 0;}\
	}\
	return ret;\
}
#define CALL_SERIALIZE_ALL(RW, FILE_NAME) \
	RW ## All(FILE_NAME)
SERIALIZE_ALL(fread)
SERIALIZE_ALL(fwrite)
//---------------------------------------------------------------------

extern double getNumber(const char * a_numberTxt)
{
	double number;
	size_t len;
	if(a_numberTxt[0] == '-' && a_numberTxt[1] != 0 && a_numberTxt[2] == '='){
		a_numberTxt += 3;//skip 3 character in "-d=1234"
	}
	len = strlen(a_numberTxt);
	if(len > 0){
		sscanf(a_numberTxt, "%lf", &number);
	}else{
		fprintf(stderr, "FATAL ERROR:unnumber in the option:%s\n", a_numberTxt);
		exit(1);
	}
	return number;
}
extern double getSecond(const char * a_timeTxt)
{
	char numberTxt[80];
	double sec = 0.0;
	size_t len;
	if(a_timeTxt[0] == '-' && a_timeTxt[1] != 0 && a_timeTxt[2] == '='){
		strncpy(numberTxt, a_timeTxt + 3, 80);//skip 3 character in "-t=1234n"
	}else{
		strncpy(numberTxt, a_timeTxt, 80);
	}
	len = strlen(numberTxt);
	if(len > 0){
		double scale = 1.0;
		if(numberTxt[len - 1] < '0' || '9' < numberTxt[len - 1]){
			if(numberTxt[len - 1] == 's'){
				scale = 1.0;
			}else if(numberTxt[len - 1] == 'm'){
				scale = MINUTE_SEC;
			}else if(numberTxt[len - 1] == 'h'){
				scale = HOUR_SEC;
			}else if(numberTxt[len - 1] == 'd'){
				scale = DAY_SEC;
			}else if(numberTxt[len - 1] == 'W'){
				scale = WEEK_SEC;
			}else if(numberTxt[len - 1] == 'M'){
				scale = MONTH_SEC;
			}else if(numberTxt[len - 1] == 'Y'){
				scale = YEAR_SEC;
			}else{
				fprintf(stderr, "FATAL ERROR:unkown time unit in the option:%s\n", a_timeTxt);
				exit(1);
			}
			numberTxt[len - 1] = 0;
		}
		sscanf(numberTxt, "%lf", &sec);
		sec *= scale;
	}else{
		fprintf(stderr, "FATAL ERROR:unnumber time in the option:%s\n", a_timeTxt);
		exit(1);
	}
	return sec;
}

extern void makeLogFileName(struct runningCondition * a_rc)
{
	char startTimeTxt[80];
	char endTimeTxt[80];
	int pretty = 0;
	formatSecond(startTimeTxt, 80, a_rc->startTime, pretty);
	formatSecond(endTimeTxt, 80, a_rc->startTime + a_rc->termTime, pretty);	
	snprintf(a_rc->logFname, FNAME_LEN, "%s_%s_%lg_%d_%d.log", startTimeTxt, endTimeTxt, a_rc->intervalTime, a_rc->logStep, a_rc->dataStep);
}
#define DEFAULT_TERMTIMESEC 60.0
#define DEFAULT_INTERVALSEC 10.0
#define DEFAULT_LOGSTEP 6
#define DEFAULT_DATASTEP 60
#define DEFAULT_DATA_FNAME "data.dat"
extern void initRunningCondition(struct runningCondition * a_rc)
{
	a_rc->startTime = 0.0;
	a_rc->termTime = DEFAULT_TERMTIMESEC;
	a_rc->intervalTime = DEFAULT_INTERVALSEC;
	a_rc->loopBeg = 0;
	a_rc->loopEnd = 6;
	a_rc->logStep = DEFAULT_LOGSTEP;
	a_rc->logNecleusReactions = 0;
	a_rc->dataStep = DEFAULT_DATASTEP;
	strncpy(a_rc->inputDataFname, DEFAULT_DATA_FNAME, FNAME_LEN);
	strncpy(a_rc->outputDataFname, DEFAULT_DATA_FNAME, FNAME_LEN);
	makeLogFileName(a_rc);
	strncpy(a_rc->paramFname, "", FNAME_LEN);
	a_rc->inputTotalEnergy = 0.0;
	a_rc->simulatedLastTime = 0.0;
};

#define BEGIN_PRINT 0
#define MID_PRINT 1
#define END_PRINT 2
extern void printRunningCondition(FILE * a_fp, int a_BeginOrEnd, const char * a_title, const time_t * a_timePtr, int argc, char * argv[])
{
	int i, j;
	char timeTxt[80];
	int pretty = 1;
	fprintf(a_fp, "[%s] %s\n", (a_BeginOrEnd == BEGIN_PRINT) ? "BEGIN" : "END", ctime(a_timePtr));
	//if(a_BeginOrEnd == BEGIN_PRINT){
		//time_t timeFin;
		//double takeTime;
		//takeTime = (e_rc.loopEnd - e_rc.loopBeg + 1) * (24.0 / 7.0);
		//formatSecond(timeTxt, 80, takeTime, pretty);
		//timeFin = *a_timePtr + (int)takeTime;
		//fprintf(a_fp, " This program will take about %s.\n It may finish at %s\n", timeTxt, ctime(&timeFin));
	//}
	fprintf(a_fp, "%s[Running Conditions]\n", a_title);
	i = 1;
	if(a_BeginOrEnd == BEGIN_PRINT || a_title[0] != 0){
		fprintf(a_fp, "%sR%d command ", a_title, i++);
		for(j = 0; j < argc; ++j){
			fprintf(a_fp, "%s", argv[j]);
			if(j < argc - 1){
				fprintf(a_fp, " ");
			}else{
				fprintf(a_fp, "\n");
			}
		}
		fprintf(a_fp, "%sR%d simulatedStartTime %s\n", a_title, i++, formatSecond(timeTxt, 80, e_rc.startTime, pretty));
		fprintf(a_fp, "%sR%d simulatedTermTime %s\n", a_title, i++, formatSecond(timeTxt, 80, e_rc.termTime, pretty));
		fprintf(a_fp, "%sR%d intervalTime %s\n", a_title, i++, formatSecond(timeTxt, 80, e_rc.intervalTime, pretty));
		fprintf(a_fp, "%sR%d loopBeg %d\n", a_title, i++, e_rc.loopBeg);
		fprintf(a_fp, "%sR%d loopEnd %d\n", a_title, i++, e_rc.loopEnd);
		fprintf(a_fp, "%sR%d logStep %d\n", a_title, i++, e_rc.logStep);
		fprintf(a_fp, "%sR%d logNecleusReactions %d\n", a_title, i++, e_rc.logNecleusReactions);
		fprintf(a_fp, "%sR%d dataStep %d\n", a_title, i++, e_rc.dataStep);
		fprintf(a_fp, "%sR%d inputDataFname %s\n", a_title, i++, e_rc.inputDataFname);
		fprintf(a_fp, "%sR%d outputDataFname %s\n", a_title, i++, e_rc.outputDataFname);
		fprintf(a_fp, "%sR%d logFname %s\n", a_title, i++, e_rc.logFname);
		fprintf(a_fp, "%sR%d paramFname %s\n", a_title, i++, e_rc.paramFname);
	}
	fprintf(a_fp, "%sR%d inputTotalEnergy %Lg [MeV]\n", a_title, i++, e_rc.inputTotalEnergy);
	if(a_BeginOrEnd == END_PRINT || a_title[0] != 0){
		fprintf(a_fp, "%sR%d simulatedLastTime %s\n", a_title, i++, formatSecond(timeTxt, 80, e_rc.simulatedLastTime, pretty));
	}
}
//---------------------------------------------------------------------

extern void printSumup(FILE * a_fp, char * a_Titile, int a_loop, char * a_Tail, int a_BeginOrEnd)
{
	static double beforeTime = -1.0;
	double diffTime = 0.0;
	char timeText[80];
	char * timeTitle;
	char timeMess[256];
	char totalText[256];
	int pretty = 1;
	struct sumMeVMol massDefectNowN;
	struct sumMeVMol massDefectNowP;
	long double diffInputHeat = 0.0, diffInputKW = 0.0, outputHeat = 0.0, diffOutputHeat = 0.0, diffOutputKW = 0.0, COP = 0.0, diffCOP = 0.0;//COP is the coefficient of performance.
	static long double oldInputHeat = 0.0, oldOutputHeat = 0.0, oldCOP = 0.0;
	double /* sumOfMassUMol_NE, */ sumOfMassUMolIni_NE, sumOfMassUMolAdd_NE, sumOfMassUMolSub_NE, heatCapacity_NE;
	double /* sumOfMassUMol_PE, */ sumOfMassUMolIni_PE, sumOfMassUMolAdd_PE, sumOfMassUMolSub_PE, heatCapacity_PE;
	double totalOfMassUMol, totalOfMassUMolIni, totalOfMassUMolAdd, totalOfMassUMolSub, totalheatCapacity, upTempertureByInput, totalDiff;
	formatSecond(timeText, 80, e_rc.simulatedLastTime, pretty);
	if(e_rc.simulatedLastTime > beforeTime){
		timeTitle = "[TIME]";
		if(beforeTime >= 0.0){
			diffTime = e_rc.simulatedLastTime - beforeTime;
		}
		beforeTime = e_rc.simulatedLastTime;
	}else{
		timeTitle = "[SAME TIME]";
	}
	snprintf(timeMess, 256, "%s %s Title %s %d %s", timeTitle, timeText, a_Titile, a_loop, a_Tail);
	fprintf(a_fp, "%s <<<<<\n\n", timeMess);
	/* sumOfMassUMol_NE = */ printElectrode(a_fp, &e_negativeElectrode, &massDefectNowN, &sumOfMassUMolIni_NE, &sumOfMassUMolAdd_NE, &sumOfMassUMolSub_NE, &heatCapacity_NE);
	fprintf(a_fp, "%s -----\n\n", timeMess);
	/*sumOfMassUMol_PE = */ printElectrode(a_fp, &e_positiveElectrode, &massDefectNowP, &sumOfMassUMolIni_PE, &sumOfMassUMolAdd_PE, &sumOfMassUMolSub_PE, &heatCapacity_PE);
	fprintf(a_fp, "%s -----\n\n", timeMess);
	totalOfMassUMolIni = sumOfMassUMolIni_NE + sumOfMassUMolIni_PE;
	totalOfMassUMolAdd = sumOfMassUMolAdd_NE + sumOfMassUMolAdd_PE;
	totalOfMassUMolSub = sumOfMassUMolSub_NE + sumOfMassUMolSub_PE;
	totalheatCapacity = heatCapacity_NE + heatCapacity_PE;
	upTempertureByInput = e_emittedPulseEnergyJ / totalheatCapacity;
	totalDiff = totalOfMassUMolAdd - totalOfMassUMolSub;//for precision
	totalOfMassUMol = totalOfMassUMolIni + totalDiff;
	fprintf(a_fp, "TOTAL of Mass * Mol : %lg (= %lg + %lg - %lg) [U]\n", totalOfMassUMol, totalOfMassUMolIni, totalOfMassUMolAdd, totalOfMassUMolSub);
	fprintf(a_fp, " = %lg [MeV]\n", totalOfMassUMol * NAvogadro * e_coefMassUToMeV);
	fprintf(a_fp, " totalheatCapacity = %lg [J/k] upTempertureByInput = %lg [k]\n\n", totalheatCapacity, upTempertureByInput);
	if(e_rc.inputTotalEnergy > 0.0L){
		struct sumMeVMol massDefectNow;
		struct sumMeVMol massDefectAll;
		struct sumMeVMol byNeutrinoAll;
		struct sumMeVMol scattered[SIZE_OF_SCATTERD];
		mergeSumMeVMol(&massDefectNow, &massDefectNowN, &massDefectNowP);
		mergeSumMeVMol(&massDefectAll, &e_negativeElectrode.massDefectAll, &e_positiveElectrode.massDefectAll);
		mergeSumMeVMol(&byNeutrinoAll, &e_negativeElectrode.byNeutrinoAll, &e_positiveElectrode.byNeutrinoAll);
		mergeScattered(scattered, e_negativeElectrode.scattered, e_positiveElectrode.scattered);
		outputHeat = printMassDefectHeat(a_fp, "TOTAL ENERGY IN THE SYSTEM", &massDefectNow, &massDefectAll, &byNeutrinoAll, scattered);
		diffInputHeat = e_rc.inputTotalEnergy - oldInputHeat;
		diffOutputHeat = outputHeat - oldOutputHeat;
		if(diffTime > 0.0){
			diffInputKW = (diffInputHeat * 3600 / diffTime) * 1.0E6 * e_energyJouleOfeV / e_energyJouleOfKWH;
			diffOutputKW = (diffOutputHeat * 3600 / diffTime) * 1.0E6 * e_energyJouleOfeV / e_energyJouleOfKWH;
		}
		COP = outputHeat / e_rc.inputTotalEnergy;
		diffCOP = COP - oldCOP;
	}
	snprintf(totalText, 256, "\nTOTAL ENERGY @ %s COP %Lg (= OUTPUT %Lg / INPUT %Lg [MeV]) diffCOP %Lg  diffOutput %Lg [MeV](= %Lg [KWH/h]) diffInput %Lg [MeV](= %Lg [KWH/h]) \n\n", timeText, COP, outputHeat, e_rc.inputTotalEnergy, diffCOP, diffOutputHeat, diffOutputKW, diffInputHeat, diffInputKW);
	oldInputHeat = e_rc.inputTotalEnergy;
	oldOutputHeat = outputHeat;
	oldCOP = COP;
	fputs(totalText, a_fp);
	fputs(totalText, stderr);
	fprintf(a_fp, "%s >>>>>\n\n", timeMess);

	if(a_BeginOrEnd == END_PRINT || e_rc.logNecleusReactions){
		//fprintf(stderr, "DEBUG:%s:e_nuclearReaction.usedCnt=%d\n", __FUNCTION__, e_nuclearReaction.usedCnt);
		printAllNuclearReaction(a_fp);
		printUnregistedAtomNumbers(a_fp);
		printUnregistedIsotopeNumbers(a_fp);
	}
	fflush(a_fp);
	fsync(fileno(a_fp));
}
//---------------------------------------------------------------------
extern int checkFileExists(const char * a_filename)
{
   FILE * fp = fopen(a_filename, "r");
   if(fp){
	   fclose(fp);
   }
   return (fp != NULL) ? 1 : 0;
}
//---------------------------------------------------------------------

extern void checkArgs(int argc, char * argv[], time_t * a_timeBeginPtr)
{
	struct runningCondition rc;
	int resumeMode;
	int i;
	char * runMode = NULL;
	char * termTimeOpt = NULL;
	char * intervalMode = NULL;
	double divideCount;
	char * logIntervalMode = NULL;
	double logInterval;
	char * logNecleusReactionsMode = NULL;
	char * dataIntervalMode = NULL;
	double dataInterval;
	char * inputDataFnameOpt = NULL;
	char * outputDataFnameOpt = NULL;
	char * logFnameOpt = NULL;
	char * paramFnameOpt = NULL;
	int adjustHydrogenHelium = 0;

	time(a_timeBeginPtr);
	initRunningCondition(&rc);
	for(i = 1; i < argc; ++i){
		//fprintf(stderr, "DEBUG:%s:argv[%d]:%s\n", __FUNCTION__, i, argv[i]);
		if(strcmp(argv[i], "-0") == 0 || strcmp(argv[i], "-f") == 0
		|| strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "-c") == 0){
			//fprintf(stderr, "DEBUG:%s:runMode\n", __FUNCTION__);
			runMode = argv[i];
		}else if(memcmp(argv[i], "-t=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:termTimeOpt\n", __FUNCTION__);
			termTimeOpt = argv[i];
			rc.termTime = getSecond(argv[i]);
			//fprintf(stderr, "DEBUG:%s:rc.termTime:%lg\n", __FUNCTION__, rc.termTime);
			if(rc.termTime < 1.0){
				fprintf(stderr, "FATAL ERROR:Too short termTime:%lg, it can't run.\n", rc.termTime);
				exit(1);
			}
		}else if(memcmp(argv[i], "-i=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:intervalMode -i\n", __FUNCTION__);
			intervalMode = argv[i];
			rc.intervalTime = getSecond(argv[i]);
			//fprintf(stderr, "DEBUG:%s:rc.intervalTime:%lg\n", __FUNCTION__, rc.intervalTime);
			if(rc.intervalTime < 1.0){
				fprintf(stderr, "FATAL ERROR:Too short intervalTime:%lg, it can't run.\n", rc.intervalTime);
				exit(1);
			}
		}else if(memcmp(argv[i], "-d=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:intervalMode -d\n", __FUNCTION__);
			intervalMode = argv[i];
			divideCount = getNumber(argv[i]);
			fprintf(stderr, "DEBUG:%s:divideCount:%lg\n", __FUNCTION__, divideCount);
			if(divideCount < 1.0){
				fprintf(stderr, "FATAL ERROR:Too small divideCount:%lg, it can't run.\n", divideCount);
				exit(1);
			}
		}else if(memcmp(argv[i], "-j=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:logIntervalMode -j\n", __FUNCTION__);
			logIntervalMode = argv[i];
			logInterval = getSecond(argv[i]);
			//fprintf(stderr, "DEBUG:%s:logInterval:%lg\n", __FUNCTION__, logInterval);
			if(logInterval < 1.0){
				fprintf(stderr, "FATAL ERROR:Too small logInterval:%lg, it can't run.\n", logInterval);
				exit(1);
			}
		}else if(memcmp(argv[i], "-s=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:logIntervalMode -s\n", __FUNCTION__);
			logIntervalMode = argv[i];
			rc.logStep = (int)getNumber(argv[i]);
			//fprintf(stderr, "DEBUG:%s:rc.logStep:%d\n", __FUNCTION__, rc.logStep);
			if(rc.logStep < 1){
				fprintf(stderr, "FATAL ERROR:Too small logStep:%d, it can't run.\n", rc.logStep);
				exit(1);
			}
		}else if(strcmp(argv[i], "-n0") == 0){
			//fprintf(stderr, "DEBUG:%s:logIntervalMode -n0\n", __FUNCTION__);
			logNecleusReactionsMode = argv[i];
			rc.logNecleusReactions = 0;
		}else if(strcmp(argv[i], "-n1") == 0){
			//fprintf(stderr, "DEBUG:%s:logIntervalMode -n1\n", __FUNCTION__);
			logNecleusReactionsMode = argv[i];
			rc.logNecleusReactions = 1;
		}else if(memcmp(argv[i], "-k=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:dataIntervalMode -k\n", __FUNCTION__);
			dataIntervalMode = argv[i];
			dataInterval = getSecond(argv[i]);
			//fprintf(stderr, "DEBUG:%s:dataInterval:%lg\n", __FUNCTION__, dataInterval);
			if(dataInterval < 1.0){
				fprintf(stderr, "FATAL ERROR:Too small dataInterval:%lg, it can't run.\n", dataInterval);
				exit(1);
			}
		}else if(memcmp(argv[i], "-p=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:dataIntervalMode -p\n", __FUNCTION__);
			dataIntervalMode = argv[i];
			rc.dataStep = getNumber(argv[i]);
			//fprintf(stderr, "DEBUG:%s:rc.dataStep:%d\n", __FUNCTION__, rc.dataStep);
			if(rc.dataStep < 1.0){
				fprintf(stderr, "FATAL ERROR:Too small dataStep:%d, it can't run.\n", rc.dataStep);
				exit(1);
			}
		}else if(memcmp(argv[i], "-I=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:inputDataFnameOpt\n", __FUNCTION__);
			inputDataFnameOpt = argv[i];
			strncpy(rc.inputDataFname, argv[i] + 3, FNAME_LEN);
		}else if(memcmp(argv[i], "-O=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:outputDataFnameOpt\n", __FUNCTION__);
			outputDataFnameOpt = argv[i];
			strncpy(rc.outputDataFname, argv[i] + 3, FNAME_LEN);
		}else if(memcmp(argv[i], "-L=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:logFnameOpt\n", __FUNCTION__);
			logFnameOpt = argv[i];
			strncpy(rc.logFname, argv[i] + 3, FNAME_LEN);
		}else if(memcmp(argv[i], "-P=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:paramFnameOpt\n", __FUNCTION__);
			paramFnameOpt = argv[i];
			strncpy(rc.paramFname, argv[i] + 3, FNAME_LEN);
		}else if(strcmp(argv[i], "-H-") == 0){
			adjustHydrogenHelium = -6;
		}else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0){
			printf(
"[help of usage]\n\
\n\
 This program simulates the nucleus reactions of the neutron generator.\n\
\n\
./simulate\n\
\n\
 You can run with the defaults, they are :\n\
 start=0[sec], term=%lg[sec], interval=%lg[sec], logStep=%d, dataStep=%d\n\
 You can specify options like following:\n\
\n\
./simulete [-0|-f|-r|-c] [-t=term] [-i=interval|-d=divideCount]\n\
 [-j=logInterval|-s=logStep] [-n0|-n1] [-k=dataInterval|-p=dataStep] \n\
 [-I=fname] [-O=fname] [-L=fname] [-P=fname] [-H-|-H+]\n\
\n\
-0 : It sets the start time is 0.\n\
   But if there is the output data file 'data.dat', It can't run.\n\
-f : It sets the start time is 0, too.\n\
   And it will ignore the input data file 'data.dat' at the start of simulation.\n\
-r : Resume mode. \n\
   It sets the start time is as same as the end time \n\
   in the input data file 'data.dat'.\n\
   And it also extracts 'interval', 'logInterval', 'dataInterval' \n\
   and '[-n0|-n1]' from the input data file 'data.dat'.\n\
   (If there is 'data.dat', '-r' is default. If not, '-0' is default.)\n\
-c : Confirm mode. \n\
   It only prints the simulated last time in the input data file 'data.dat'.\n\
   And the program does nothing anymore, it stops imediately. \n\
-t=term : It sets the long term from the start time to the end time of your\n\
   simulation. \n\
   The 'term' is a positive floating number and a character of time unit.\n\
   The time unit is 's', 'm', 'h', 'd', 'W', 'M' or 'Y'.\n\
   They are second, minute, hour, day, week, month or year.\n\
   If there is not a charcter of time unit, the default uint is second.\n\
   The default term is 60 seconds\n\
-i=interval : It sets the interval time of your simulation.\n\
   The process of simulation will occur after each proceeding of an interval.\n\
   The proceeding is not actual time waiting.\n\
   It's only the programming techniques of simulation.\n\
   It simply adds the interval time to the past time register.\n\
   The 'interval' is a positive floating number and a charcter of time unit.\n\
   Its time unit is as same as option '-t=term'.\n\
   If there is not '-i=interval', the default is 10 second with \n\
   both '-0' option and '-f' option. \n\
   Or the default is extracted from the file 'data.dat' with '-r' option.\n\
   You can override 'interval' in the resume mode.\n\
   The difference of interval time will make the accuracy of simulation.\n\
   The minimum interval is 1 second that takes many long time of\n\
   running simulation. We recommend 10 second.\n\
   The difference of accuracy of simulation between 1 second and 10 second is about 20[%%].\n\
-d=divideCount : It sets the interval time of your simulation by dividing\n\
   calculation.\n\
   The calculation is : (interval = termSecond / divideCount)\n\
   'divideCount' is not saved into the output data file 'data.dat'.\n\
-j=logInterval : It sets the step of logging by a calculation.\n\
   The 'logInterval' is the interval time of logging.\n\
   Its time format is as same as the option '-t=term'.\n\
   The calculation is : (logStep = logInterval / interval)\n\
   The default 'logInterval' is 6 times of interval.\n\
   The 'logInterval' is not saved into the output data file 'data.dat'.\n\
-s=logStep : It sets the step of logging. \n\
   If 'logStep' is 1, it prints log into the log file at each interval.\n\
   If 'logStep' is 2, it prints log into the log file at once every 2 intervals.\n\
   The 'logStep' must be positive integer.\n\
   The 'logStep' is saved into the output data file 'data.dat' instead of 'logInterval'.\n\
-n0 : It specifies not to print the necleus reactions at each logging time.\n\
   It prints the necleus reactions only at the end of simulation.\n\
   The default is '-n0', not '-n1' \n\
-n1 : It specifies to print the necleus reactions at each logging time.\n\
-k=dataInterval : It sets the step of writing the output data file\n\
   by the calculation.\n\
   The 'dataInterval' is the interval time of writing.\n\
   Its time format is as same as the option '-t=term'.\n\
   The calculation is : (dataStep = dataInterval / interval)\n\
   The default 'dataInterval' is 60 times of interval.\n\
   The 'dataInterval' is not saved into the output data file 'data.dat'.\n\
-p=dataStep : It sets the step of writing the output data file.\n\
   If 'dataStep' is 1, it writes the output data file at each interval.\n\
   If 'dataStep' is 2, it writes the output data file at once every 2 intervals.\n\
   The 'dataStep' must be positive integer.\n\
   The 'dataStep' is saved into the output data file 'data.dat' instead of 'dataInterval'.\n\
-I=fname : It specifies the name of the input data file.\n\
   The default name of data file is 'data.dat'.\n\
   The input data file is used in the resume mode with option '-r'.\n\
-O=fname : It specifies the name of a output data file.\n\
   The default name of data file is 'data.dat'.\n\
   The output data file is written at the end of the simulation and \n\
   at the each 'dataInterval'.\n\
   And the program write one more output data file at the end.\n\
   Its name is \"(simulated-last-time).dat\".\n\
   Its contents is as same as 'data.dat'.\n\
-L=fname : It specifies the name of a output logging file.\n\
   The default name of logging file is automatically generated.\n\
   Its foramt is (start)_(end)_(interval)_(logIntreval)_(dataInterval).log\n\
   (start) is the start time of the simulation.\n\
   (end) is the end time. \n\
   (interval) is the interval time.\n\
   (logInterval) is the interval time of logging.\n\
   (dataInterval)is the interval time of writing the output data file.\n\
   If there is a logging file, it can't run.\n\
-P=fname : It specifies the name of the input parameter file.\n\
    The input parameter file can overide the parameter in the data file\n\
	or the default vales of program.\n\
	Tipically, this option is used for stopping the high voltage current \n\
	by the zero of both 'e_emittedElectronMol' and 'e_emittedProtonMol'.\n\
	If fname is 'default', all input parameter are set to default values.\n\
	The format of input parameter file is like the programing language C.\n\
	The foramt of line is \"TYPE NAME = VALUE;//COMMENT\"\n\
	'TYPE' is 'int' or 'double'.\n\
	'NAME' is the name of a parameter.\n\
	'VALUE' is the value of a parameter.\n\
	'COMMENT' is the comment of a parameter.\n\
	You can know the purpose of parameters in the source code, '" __FILE__ "'.\n\
	You can write some part of parameters in the parameter file.\n\
	Here is the sample of a input parameter file with default values.\n\
		//[User conditions and default values]\n\
		int e_rangeDecimalDigitPrecision = %d;\n\
		double e_detectLimitRateForIsotope = %lg;\n\
		char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"%s\";//This is only used at zero time of the simulation. \n\
		char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"%s\";//This is only used at zero time of the simulation. \n\
		double e_appliedVoltageScale = %lg;\n\
		double e_emittedElectronMol = %lg;//The number of moles of electrons emitted per one(1) second from the negative electrode by high-voltage pulse like The Cockcroft-Walton voltage multiplier circuit.\n\
		double e_collideElectronRateOnElectrode = %lg;\n\
		double e_collideElectronRateForMidiMeV = %lg;\n\
		double e_collideElectronRateForMiniMeV = %lg;\n\
		double e_emittedProtonMol = %lg;//The number of moles of protons emitted per one(1) second from the by positive electrode by high-voltage pulse like The Cockcroft-Walton voltage multiplier circuit.\n\
		double e_neutronGenInSpaceProtonRate = %lg;\n\
		double e_neutronGenInSpaceFallRate = %lg;\n\
		double e_hydrogenGenInSpaceProtonRate = %lg;\n\
		double e_collideProtonRateOnElectrode = %lg;\n\
		double e_collideElectronMidiMeV = %lg;\n\
		double e_collideElectronMiniMeV = %lg;\n\
		int e_usebulletCrossSection = %d;\n\
		int e_useComptonEffect = %d;\n\
-H- : It specifies to reduce the amount of hydrogen and helium in the positive electrode\n\
    into one-1,000,000th before running the simulation.\n\
./simulate -h\n\
 (Show this help with '-h' or '-help'.)\n", 
			DEFAULT_TERMTIMESEC, DEFAULT_INTERVALSEC, DEFAULT_LOGSTEP, DEFAULT_DATASTEP,
			e_rangeDecimalDigitPrecision,
			e_detectLimitRateForIsotope,
			e_negativeElectrodeAtomicMols,
			e_positiveElectrodeAtomicMols,
			e_appliedVoltageScale,
			e_emittedElectronMol,
			e_collideElectronRateOnElectrode,
			e_collideElectronRateForMidiMeV,
			e_collideElectronRateForMiniMeV,
			e_emittedProtonMol,
			e_neutronGenInSpaceProtonRate,
			e_neutronGenInSpaceFallRate,
			e_hydrogenGenInSpaceProtonRate,
			e_collideProtonRateOnElectrode,
			e_collideElectronMidiMeV,
			e_collideElectronMiniMeV,
			e_usebulletCrossSection,
			e_useComptonEffect);
			exit(0);
		}
	}
	if(runMode == NULL){
		//fprintf(stderr, "DEBUG:%s:(runMode == NULL)\n", __FUNCTION__);
		if(checkFileExists(rc.inputDataFname)){
			//fprintf(stderr, "DEBUG:%s:FORCE (runMode = -r)\n", __FUNCTION__);
			runMode = "-r";
		}else{
			//fprintf(stderr, "DEBUG:%s:FORCE (runMode = -0)\n", __FUNCTION__);
			runMode = "-0";
		}
	}
	if(strcmp(runMode, "-0") == 0 && checkFileExists(rc.outputDataFname)){
		fprintf(stderr, "FATAL ERROR:When there is the option -0 or the start time is 0 and there is an old output file '%s', it can't run.\n", rc.outputDataFname);
		exit(1);
	}
	if(strcmp(runMode, "-r") == 0 && !checkFileExists(rc.inputDataFname)){
		fprintf(stderr, "FATAL ERROR:When there is not a input file '%s' in the resume mode, it can't run.\n", rc.inputDataFname);
		exit(1);
	}
	

	if(strcmp(runMode, "-0") == 0 || strcmp(runMode, "-f") == 0){
		//fprintf(stderr, "DEBUG:%s:resumeMode = 0;\n", __FUNCTION__);
		resumeMode = 0;
		e_rc = rc;
	}else if(strcmp(runMode, "-r") == 0 || strcmp(runMode, "-c") == 0){
		//fprintf(stderr, "DEBUG:%s(%d):resumeMode = 1 rc.inputDataFname:%s;\n", __FUNCTION__, __LINE__, rc.inputDataFname);
		resumeMode = 1;
		if(!CALL_SERIALIZE_ALL(fread, rc.inputDataFname)){
			fprintf(stderr, "FATAL ERROR:Can't read the data file '%s'..\n", rc.inputDataFname);
			exit(1);
		}
		if(strcmp(runMode, "-c") == 0){
			int pretty = 0;
			char simulatedLastTimeTxt[80];
			formatSecond(simulatedLastTimeTxt, 80, e_rc.simulatedLastTime, pretty);
			fprintf(stderr, "simulatedLastTime %s\n", simulatedLastTimeTxt);
			exit(0);
		}
		//fprintf(stderr, "DEBUG:%s(%d):resumeMode = 1;\n", __FUNCTION__, __LINE__);
		e_rc.startTime = e_rc.simulatedLastTime;//resume mode!
		//fprintf(stderr, "DEBUG:%s:NEW e_rc.startTime:%lg;\n", __FUNCTION__, e_rc.startTime);
		//override by user's options
		if(termTimeOpt){
			e_rc.termTime = rc.termTime;
			//fprintf(stderr, "DEBUG:%s(%d):termTimeOpt e_rc.termTime:%lg\n", __FUNCTION__, __LINE__, e_rc.termTime);
		}
		if(intervalMode){
			if(intervalMode[1] == 'i'){
				e_rc.intervalTime = rc.intervalTime;
				//fprintf(stderr, "DEBUG:%s(%d):e_rc.intervalTime:%lg\n", __FUNCTION__, __LINE__, e_rc.intervalTime);
			}
		}
		if(logIntervalMode){
			if(logIntervalMode[1] == 's'){
				e_rc.logStep = rc.logStep;
				//fprintf(stderr, "DEBUG:%s(%d):e_rc.logStep:%d\n", __FUNCTION__, __LINE__, e_rc.logStep);
			}
		}
		if(logNecleusReactionsMode){
			e_rc.logNecleusReactions = rc.logNecleusReactions;
			//fprintf(stderr, "DEBUG:%s(%d):e_rc.logNecleusReactions:%d\n", __FUNCTION__, __LINE__, e_rc.logNecleusReactions);
		}
		if(dataIntervalMode){
			if(dataIntervalMode[1] == 'p'){
				e_rc.dataStep = rc.dataStep;
				//fprintf(stderr, "DEBUG:%s(%d):e_rc.dataStep:%d\n", __FUNCTION__, __LINE__, e_rc.dataStep);
			}
		}
		if(inputDataFnameOpt){
			strncpy(e_rc.inputDataFname, rc.inputDataFname, FNAME_LEN);
			//fprintf(stderr, "DEBUG:%s(%d):e_rc.inputDataFname:%s\n", __FUNCTION__, __LINE__, e_rc.inputDataFname);
		}
		if(outputDataFnameOpt){
			strncpy(e_rc.outputDataFname, rc.outputDataFname, FNAME_LEN);
			//fprintf(stderr, "DEBUG:%s(%d):e_rc.outputDataFname:%s\n", __FUNCTION__, __LINE__, e_rc.outputDataFname);
		}
		if(logFnameOpt){
			strncpy(e_rc.logFname, rc.logFname, FNAME_LEN);
			//fprintf(stderr, "DEBUG:%s(%d):e_rc.logFname:%s\n", __FUNCTION__, __LINE__, e_rc.logFname);
		}
		if(paramFnameOpt){
			strncpy(e_rc.paramFname, rc.paramFname, FNAME_LEN);
			//fprintf(stderr, "DEBUG:%s(%d):e_rc.paramFname:%s\n", __FUNCTION__, __LINE__, e_rc.paramFname);
		}
	}
	// calculate options.
	if(intervalMode){
		if(intervalMode[1] == 'd'){
			e_rc.intervalTime = e_rc.termTime / divideCount;
		}
	}
	if(logIntervalMode){
		if(logIntervalMode[1] == 'j'){
			e_rc.logStep = logInterval / e_rc.intervalTime;
		}
	}
	if(dataIntervalMode){
		if(dataIntervalMode[1] == 'k'){
			e_rc.dataStep = dataInterval / e_rc.intervalTime;
		}
	}
	e_rc.loopBeg = 0;
	e_rc.loopEnd = (int)ceil(e_rc.termTime / e_rc.intervalTime);
	if(!logFnameOpt){
		makeLogFileName(&e_rc);
	}
	if(checkFileExists(e_rc.logFname)){
		fprintf(stderr, "FATAL ERROR:There is an old log file '%s'.\n", e_rc.logFname);
		exit(1);
	}
	printRunningCondition(stderr, BEGIN_PRINT, "", a_timeBeginPtr, argc, argv);
	e_logFp = fopen(e_rc.logFname, "w");
	if(e_logFp){
		printRunningCondition(e_logFp, BEGIN_PRINT, "", a_timeBeginPtr, argc, argv);
	}else{
		fprintf(stderr, "FATAL ERROR:can't open log file:%s\n", argv[i + 1]);
		exit(1);
	}
	if(paramFnameOpt){
		//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
		readUserConditions(rc.paramFname);
		//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	}	
	initCosTable();
	initNearZero();
	initiallyCalculateUserContants();
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	if(!resumeMode){
		//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
		initElectrodes(&e_negativeElectrode, e_negativeElectrodeAtomicMols, e_detectLimitRateForIsotope, "Negative");
		//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
		initElectrodes(&e_positiveElectrode, e_positiveElectrodeAtomicMols, e_detectLimitRateForIsotope, "Positive");
		//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
		initHashTable(&e_nuclearReaction, 
			calcHashSeedOfNuclearReactionKey,
			allocCopyNuclearReactionKey,
			freeCopyNuclearReactionKey,
			compareNuclearReactionKey, 
			allocCopyNuclearReactionValue,
			free,
			foundActionForNuclearReactionValue,
			"nuclearReaction", 100);
	}
	if(adjustHydrogenHelium < 0){
		double scale;
		scale = pow(10.0, adjustHydrogenHelium);
		struct atomNodeConst * atomPtr;
		int atomicNumberAry[] = {ATOMICNUMBER_HYDROGEN, ATOMICNUMBER_DEUTERIUM, ATOMICNUMBER_TRITIUM, ATOMICNUMBER_HELIUM, ATOMICNUMBER_HELIUM2, ATOMICNUMBER_HELIUM3};
		int massNumberAry[] = {MASSNUMBER_HYDROGEN, MASSNUMBER_DEUTERIUM, MASSNUMBER_TRITIUM, MASSNUMBER_HELIUM, MASSNUMBER_HELIUM2, MASSNUMBER_HELIUM3};
		int i, j, k;
		for(i = 0; i < sizeof(atomicNumberAry) / sizeof(int); ++i){
			for(j = 0; j < 2; ++j){
				atomPtr = findAtom((j == 0) ? &e_positiveElectrode : &e_negativeElectrode, atomicNumberAry[i], massNumberAry[i]);
				if(atomPtr){
					double mol, molSub, molRed;				
					mol = getMol(atomPtr->vPtr);
					molSub = mol * (1.0 - scale);
					molRed = mol - molSub;
					setMolSub(atomPtr->vPtr, molSub);
					for(k = 0; k < 2; ++k){
						fprintf((k == 0) ? e_logFp : stderr, "Reduce %s(%s) on %s Electrode. %lg - %lg = %lg [mol]\n", atomPtr->key->isotopePropertyPtr->name, atomPtr->key->isotopePropertyPtr->symbol, (j == 0) ? "Positive" : "Negative", mol, molSub, molRed);
					}
				}
			}
		}
	}
}
//---------------------------------------------------------------------

#define TITLE_LEN 80
int main(int argc, char * argv[])
{
	time_t timeBegin;
	time_t timeMid;
	time_t timeEnd;
	int loop;
	char tail[TITLE_LEN];
	//static char rotationMark[] = "0123456789";
	
	initiallyCalculatePhysicsContants();//This must be the first.
	initAtomProperties();//This must be the second.
	initUserConditionsByDefault();//This must be the third.
	checkArgs(argc, argv, &timeBegin);//This must be the forth.
	
	snprintf(tail, TITLE_LEN, "in [%d, %d] by %lg [sec]", e_rc.loopBeg, e_rc.loopEnd, e_rc.intervalTime);

	printContants(e_logFp);
	printAtomProperties(e_logFp);
	
	//fprintf(stderr,"DEBUG:stop\n"); exit(1);//debug stop

	for(loop = e_rc.loopBeg; loop < e_rc.loopEnd; ++loop){
		e_rc.simulatedLastTime = e_rc.startTime + e_rc.intervalTime * loop;
		fprintf(stderr, "%d/%d\n", loop, e_rc.loopEnd);
		if((loop - e_rc.loopBeg) % e_rc.logStep == 0){
			if(loop == e_rc.loopBeg){
				printSumup(e_logFp, "Begin loop", loop, tail, BEGIN_PRINT);
				//printSumup(stderr, "Begin loop", loop, tail, BEGIN_PRINT);//DEBUG
			}else{
				printSumup(e_logFp, "Processing", loop, tail, MID_PRINT);
				time(&timeMid);
				timeEnd = ((timeMid - timeBegin) * (e_rc.loopEnd - e_rc.loopBeg) / (loop - e_rc.loopBeg)) + timeBegin;
				fprintf(stderr, "This program will finish at %s\n", ctime(&timeEnd));
				//fprintf(stderr, "DEBUG:e_debugDecayMassDiff: %lg\n", e_debugDecayMassDiff);
				//fprintf(e_logFp, "DEBUG:e_debugDecayMassDiff: %lg\n", e_debugDecayMassDiff);
			}
		}
		if(loop > 0 && (loop - e_rc.loopBeg) % e_rc.dataStep == 0){
			if(!CALL_SERIALIZE_ALL(fwrite, e_rc.outputDataFname)){
				fprintf(stderr, "ERROR:CALL_SERIALIZE_ALL(fwrite, %s)\n", e_rc.outputDataFname);
			}
		}
		pulseCurrent(e_rc.intervalTime);
		e_rc.inputTotalEnergy += (e_emittedPulseEnergyMeV * e_rc.intervalTime);
		//printSumup(e_logFp, "after pulseCurrent", loop, tail, MID_PRINT);

		absorbeOrDecayNeutronInElectrode(&e_negativeElectrode, e_rc.intervalTime);
		absorbeOrDecayNeutronInElectrode(&e_positiveElectrode, e_rc.intervalTime);
		//printSumup(e_logFp, "after absorbeOrDecayNeutronInElectrode", loop, tail, MID_PRINT);
		
		decayNeuclayWithoutNewtron(&e_negativeElectrode, e_rc.intervalTime);
		decayNeuclayWithoutNewtron(&e_positiveElectrode, e_rc.intervalTime);
		//printSumup(e_logFp, "after decayNeuclayWithoutNewtron", loop, tail, MID_PRINT);
		
		comptonEffect(&e_negativeElectrode);
		comptonEffect(&e_positiveElectrode);
	}
	e_rc.simulatedLastTime = e_rc.startTime + e_rc.intervalTime * loop;
	fprintf(stderr, "%d/%d\n", loop, e_rc.loopEnd);
	printSumup(e_logFp, "End loop", loop, tail, END_PRINT);

	time(&timeEnd);
	printRunningCondition(stderr, END_PRINT, "", &timeEnd, argc, argv);
	printRunningCondition(e_logFp, END_PRINT, "", &timeEnd, argc, argv);
	
	if(fclose(e_logFp) != 0){
		fprintf(stderr, "ERROR:fclose %s\n", e_rc.logFname);
	}
	if(!CALL_SERIALIZE_ALL(fwrite, e_rc.outputDataFname)){
		fprintf(stderr, "ERROR:CALL_SERIALIZE_ALL(fwrite, %s)\n", e_rc.outputDataFname);
	}
	{
		int pretty = 0;
		char simulatedLastTimeTxt_Dat[80];
		formatSecond(simulatedLastTimeTxt_Dat, 80, e_rc.simulatedLastTime, pretty);
		strcat(simulatedLastTimeTxt_Dat, ".dat");
		if(!CALL_SERIALIZE_ALL(fwrite, simulatedLastTimeTxt_Dat)){
			fprintf(stderr, "ERROR:CALL_SERIALIZE_ALL(fwrite, %s)\n", simulatedLastTimeTxt_Dat);
		}
	}

	return 0;
}