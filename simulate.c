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
double e_raduisH = 5.29e-11;//[m], alomost radius of hydrogen atom
double e_raduisNi = 1.24e-10;//[m], alomost radius of Nickel atom
double e_r0 = 1.2E-15;//[m], radius zero in "R = r0 * pow(A, 1/3)", alomost radius of proton
//double e_radiusOfProtonByPion = 1.4138E-15;//[m], radius of proton by pion(pi-meson)
double e_radiusOfProtonByRmsCharge = 0.8751E-15;//[m], radius of proton by RMS(root mean square) charge 
double e_radiusOfDeuteronByRmsCharge = 2.1413E-15;//[m], radius of deuteron by RMS(root mean square) charge
double e_radiusOfElectronByWeakBoson = 2.4547e-18;//[m], radius of electron by weak boson, also radius of up-quark / down-quark by weak boson
//[calcuated Physics constants]
double e_coefElectroPotentialMeV;//(e_elementaryCharge * e_elementaryCharge) / (4.0 * M_PI * e_vacuumPermittivity * e_energyJouleOfeV * 1.0E6)
double e_relativeRadiusOfElectronByWeakBoson;//(e_radiusOfElectronByWeakBoson / e_r0)
double e_coefMassUToMeV;
double e_coefMeVtoMassU;
double e_massNeutronMassU;//[u], the mass of a neutron
double e_massProtonMassU;//[u], the mass of a proton
double e_massElectronMassU;//[u], the mass of an electron
double e_betaEnergyMassU;//
double e_betaEnergyMeV;// the kinetic energy that is released in beta decay, (It is random split into electrons and the anti-electron neutrino)
double e_coulombBarrierOf2Protons;// 0.822743 = e_coefElectroPotentialMeV / (e_radiusOfProtonByRmsCharge + e_radiusOfProtonByRmsCharge);
//0.822743 / 0.782333 = 1.0516531962731983 
double e_maxApparentRelativeRadiusOfCrossSectionForElectron;
double e_maxApparentRelativeRadiusOfCrossSectionForNecleus;
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

double e_collideElectronRateOnElectrode;//Probability to collide between (the electron emitted from the negative electrode) and (the particles in the positive electrode), the electron has large kinetic energy than the energy of beta decay, the electron emitted from the negative electrode, that reaches on the positive electrode
double e_collideElectronRateForMidiMeV;//This is paired with 'e_collideMidiMeV'
double e_collideElectronRateForMiniMeV;//This is paired with 'e_collideMiniMeV'
double e_emittedProtonMol;//The number of moles of protons emitted per one(1) second from the by positive electrode by high-voltage pulse like the Cockcroft-Walton voltage multiplier circuit.
double e_neutronGenInSpaceProtonRate;//The probability applied for the proton that is emitted from the positive electrode and will collide with the electron emitted from the negative electrode in the intermediate space of electrodes and will transform into a neutron with enough kinetic energy greater than the energy of beta decay.
double e_neutronGenInSpaceFallRate;//Almost all of the neutrons generated in the intermediate space fall on the positive electrode. The neutrons have smaller energy(= e_appliedDefectMeV).
double e_hydrogenGenInSpaceProtonRate;//The probability applied for the proton that is emitted from the positive electrode and will collide with the electron emitted from the negative electrode in the intermediate space of electrodes and will transform into a hydrogen with small kinetic energy less than the energy of beta decay.
double e_collideProtonRateOnElectrode;//Probability to collide between the proton and ((the nucleus or the electron) in the negative electrode), the proton has large kinetic energy than the energy of beta decay, the proton emitted from the positive electrode reach the negative electrode
double e_collideMidiMeV;//This is almost half of e_betaEnergyMeV = 0.78 [MeV]
double e_collideMiniMeV; //10[keV] is 1.16e+8 [K]. It is almost as same as D-T fusion reactor. It is about 100 times of the heat of the corona of the Sun, 1.0e+6 [K]. It is about 1,000 times of the energy to ionize a hydrogen, 13.6 [eV].
int e_usebulletCrossSection;// 1 (dafault) or 0
int e_useProtonScattering;// 1 (dafault) or 0, also for deuterium and tritium
int e_useNeutonScattering;// 1 (dafault) or 0
int e_useAlphaScattering;// 1 (dafault) or 0
int e_useElectronScattering;// 1 (dafault) or 0
int e_useComptonEffect;// 1 (dafault) or 0
double e_rateForAlphaParticle;// 0.5 (dafault) (from 0.0 to 1.0), The alpha particle will take the part of the enegy of the mass defect by this rate.
double e_rateForProtonAtBetaPlus;// 1.0/3.0 (dafault) (from 0.0 to 1.0),
double e_rateFor2ProtonAtBetaPlus;// 0.5 (dafault) (from 0.0 to 1.0),
double e_rateFor3ProtonAtBetaPlus;// 0.6 (dafault) (from 0.0 to 1.0),
double e_rateForAlphaParticleAtBetaPlus;// 1.0/3.0 (dafault) (from 0.0 to 1.0),
double e_rateForProtonAtEC;// 0.5 (dafault) (from 0.0 to 1.0),
double e_rateFor2ProtonAtEC;// 2.0/3.0 (dafault) (from 0.0 to 1.0),
double e_rateFor3ProtonAtEC;// 0.75 (dafault) (from 0.0 to 1.0),
double e_rateForAlphaParticleAtEC;// 0.5 (dafault) (from 0.0 to 1.0),
double e_rateForNeytonAtBetaMinus;// 1.0/3.0 (dafault) (from 0.0 to 1.0),
double e_rateFor2NeytonAtBetaMinus;// 1.0/4.0 (dafault) (from 0.0 to 1.0),
double e_rateFor3NeytonAtBetaMinus;// 1.0/5.0 (dafault) (from 0.0 to 1.0),
double e_rateFor4NeytonAtBetaMinus;// 1.0/6.0 (dafault) (from 0.0 to 1.0),
double e_rateForAlphaParticleAtBetaMinus;// 1.0/3.0 (dafault) (from 0.0 to 1.0),

extern void initUserConditionsByDefault()
{
	e_rangeDecimalDigitPrecision = 2;
	e_detectLimitRateForIsotope = 1.0E-21; //1.0E-15;
	strcpy(e_negativeElectrodeAtomicMols, "Ni=1.0");
	strcpy(e_positiveElectrodeAtomicMols, "H=0.1, Ni=1.0");
	e_appliedVoltageScale = 1.06;// > 1.0516 = e_coulombBarrierOf2Protons / e_betaEnergyMeV, 4.5 => 3.5[MeV]for aluminum
	e_emittedElectronMol = 0.04E-10;
	e_collideElectronRateOnElectrode = 1.0;
	e_collideElectronRateForMidiMeV = 1.0;
	e_collideElectronRateForMiniMeV = 1.0;
	e_emittedProtonMol = 0.4E-10;
	e_neutronGenInSpaceProtonRate = 2 * (e_radiusOfElectronByWeakBoson * e_radiusOfElectronByWeakBoson) / (e_r0 * e_r0);
	e_neutronGenInSpaceFallRate = 0.95;
	e_hydrogenGenInSpaceProtonRate = 0.001;
	e_collideProtonRateOnElectrode = 1.0;
	e_collideMidiMeV = 0.36;
	e_collideMiniMeV = 0.01;
	e_usebulletCrossSection = 1;
	e_useProtonScattering = 1;
	e_useNeutonScattering = 1;
	e_useAlphaScattering = 1;
	e_useElectronScattering = 1;
	e_useComptonEffect = 1;
	e_rateForAlphaParticle = 0.5;
	e_rateForProtonAtBetaPlus = 1.0 / 3.0;
	e_rateFor2ProtonAtBetaPlus = 0.5;
	e_rateFor3ProtonAtBetaPlus = 0.6;
	e_rateForAlphaParticleAtBetaPlus = 1.0 / 3.0;
	e_rateForProtonAtEC = 0.5;
	e_rateFor2ProtonAtEC = 2.0/3.0;
	e_rateFor3ProtonAtEC = 0.75;
	e_rateForAlphaParticleAtEC = 0.5;
	e_rateForNeytonAtBetaMinus = 1.0/3.0;
	e_rateFor2NeytonAtBetaMinus = 1.0/4.0;
	e_rateFor3NeytonAtBetaMinus = 1.0/5.0;
	e_rateFor4NeytonAtBetaMinus = 1.0/6.0;
	e_rateForAlphaParticleAtBetaMinus = 1.0/3.0;
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
	SERIALIZE_VALUE(RW, e_collideMidiMeV, a_fp)\
	SERIALIZE_VALUE(RW, e_collideMiniMeV, a_fp)\
	SERIALIZE_VALUE(RW, e_usebulletCrossSection, a_fp)\
	SERIALIZE_VALUE(RW, e_useProtonScattering, a_fp)\
	SERIALIZE_VALUE(RW, e_useNeutonScattering, a_fp)\
	SERIALIZE_VALUE(RW, e_useAlphaScattering, a_fp)\
	SERIALIZE_VALUE(RW, e_useElectronScattering, a_fp)\
	SERIALIZE_VALUE(RW, e_useComptonEffect, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForAlphaParticle, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForProtonAtBetaPlus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor2ProtonAtBetaPlus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor3ProtonAtBetaPlus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForAlphaParticleAtBetaPlus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForProtonAtEC, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor2ProtonAtEC, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor3ProtonAtEC, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForAlphaParticleAtEC, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForNeytonAtBetaMinus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor2NeytonAtBetaMinus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor3NeytonAtBetaMinus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateFor4NeytonAtBetaMinus, a_fp)\
	SERIALIZE_VALUE(RW, e_rateForAlphaParticleAtBetaMinus, a_fp)\
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
	e_relativeRadiusOfElectronByWeakBoson = e_radiusOfElectronByWeakBoson / e_r0;
	e_coefMassUToMeV = e_massNeutronMeV * ( e_unifiedAtomicMassUnitKg / e_massNeutronKg);
	e_coefMeVtoMassU = 1 / e_coefMassUToMeV;
	e_massNeutronMassU = e_massNeutronMeV * e_coefMeVtoMassU;
	e_massProtonMassU = e_massProtonMeV * e_coefMeVtoMassU;
	e_massElectronMassU = e_massElectronMeV * e_coefMeVtoMassU;
	e_betaEnergyMassU = e_massNeutronMassU - (e_massProtonMassU + e_massElectronMassU);
	e_betaEnergyMeV = e_massNeutronMeV - (e_massProtonMeV + e_massElectronMeV);	
	e_coulombBarrierOf2Protons = e_coefElectroPotentialMeV / (e_radiusOfProtonByRmsCharge + e_radiusOfProtonByRmsCharge);
	//fprintf(stderr, "e_betaEnergyMeV %lg e_coulombBarrierOf2Protons %lg\n", e_betaEnergyMeV, e_coulombBarrierOf2Protons);
	e_maxApparentRelativeRadiusOfCrossSectionForElectron = pow((pow(e_raduisNi, 3.0) / 28.0), 1.0 / 3.0);//pow((pow(1.24e-10, 3.0)/28), 1.0 / 3.0) = 4.083529287251773e-11
	if(e_maxApparentRelativeRadiusOfCrossSectionForElectron > e_raduisH){
		e_maxApparentRelativeRadiusOfCrossSectionForElectron = e_raduisH;
	}
	e_maxApparentRelativeRadiusOfCrossSectionForElectron = e_maxApparentRelativeRadiusOfCrossSectionForElectron / e_r0;
	e_maxApparentRelativeRadiusOfCrossSectionForNecleus = e_raduisNi / e_r0;
	//exit(1);//debug
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
	if(e_appliedVoltageScale < e_coulombBarrierOf2Protons / e_betaEnergyMeV){
		fprintf(stderr,"WARNING:%s:e_appliedVoltageScale:%lg < %lg = %lg / %lg\n", __FUNCTION__, e_appliedVoltageScale, e_coulombBarrierOf2Protons / e_betaEnergyMeV, e_coulombBarrierOf2Protons, e_betaEnergyMeV);
	}
	e_stepsOfLostingEnergy = (int)(e_appliedVoltageScale * 2.0);
	e_appliedVoltageMassU = e_betaEnergyMassU * e_appliedVoltageScale;
	e_appliedVoltageMeV = e_betaEnergyMeV * e_appliedVoltageScale;
	if(e_appliedVoltageMeV < e_collideMidiMeV){
		fprintf(stderr,"FATAL ERROR:%s:e_appliedVoltageMeV:%lg < e_collideMidiMeV:%lg\n", __FUNCTION__, e_appliedVoltageMeV, e_collideMidiMeV);
		exit(1);
	}
	if(e_collideMiniMeV > e_collideMidiMeV){
		fprintf(stderr,"FATAL ERROR:%s:e_collideMiniMeV:%lg < e_collideMidiMeV:%lg\n", __FUNCTION__, e_collideMiniMeV, e_collideMidiMeV);
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
	fprintf(a_fp, "P%d e_raduisH %lg [m]\n", i++, e_raduisH);
	fprintf(a_fp, "P%d e_raduisNi %lg [m]\n", i++, e_raduisNi);
	fprintf(a_fp, "P%d e_r0 %lg [m]\n", i++, e_r0);
	//fprintf(a_fp, "P%d e_radiusOfProtonByPion %lg [m]\n", i++, e_radiusOfProtonByPion);
	fprintf(a_fp, "P%d e_radiusOfProtonByRmsCharge %lg [m]\n", i++, e_radiusOfProtonByRmsCharge);
	fprintf(a_fp, "P%d e_radiusOfDeuteronByRmsCharge %lg [m]\n", i++, e_radiusOfDeuteronByRmsCharge);
	fprintf(a_fp, "P%d e_radiusOfElectronByWeakBoson %lg [m]\n", i++, e_radiusOfElectronByWeakBoson);
	
	fprintf(a_fp, "[calculated Physics constants]\n");
	i = 1;
	fprintf(a_fp, "Q%d coefMeVtoKg(by Neutron)  %lg [kg]/[MeV/C^2]\n", i++, e_massNeutronKg / e_massNeutronMeV);
	fprintf(a_fp, "Q%d coefMeVtoKg(by Proton)   %lg [kg]/[MeV/C^2]\n", i++, e_massProtonKg / e_massProtonMeV);
	fprintf(a_fp, "Q%d coefMeVtoKg(by Electron) %lg [kg]/[MeV/C^2]\n", i++, e_massElectronKg / e_massElectronMeV);
	
	fprintf(a_fp, "Q%d e_coefElectroPotentialMeV %lg\n", i++, e_coefElectroPotentialMeV);
	fprintf(a_fp, "Q%d e_relativeRadiusOfElectronByWeakBoson %lg\n", i++, e_relativeRadiusOfElectronByWeakBoson);
	
	fprintf(a_fp, "Q%d e_coefMassUToMeV %lg\n", i++, e_coefMassUToMeV);
	fprintf(a_fp, "Q%d e_coefMeVtoMassU %lg\n", i++, e_coefMeVtoMassU);
	fprintf(a_fp, "Q%d e_massNeutronMassU  %lg [u]\n", i++, e_massNeutronMassU);
	fprintf(a_fp, "Q%d e_massProtonMassU   %lg [u]\n", i++, e_massProtonMassU);
	fprintf(a_fp, "Q%d e_massElectronMassU %lg [u]\n", i++, e_massElectronMassU);
	fprintf(a_fp, "Q%d e_betaEnergyMassU %lg [u]\n", i++, e_betaEnergyMassU);
	fprintf(a_fp, "Q%d e_betaEnergyMeV   %lg [MeV/C^2]\n", i++, e_betaEnergyMeV);
	fprintf(a_fp, "Q%d e_coulombBarrierOf2Protons %lg [MeV/C^2]\n", i++, e_coulombBarrierOf2Protons);
	fprintf(a_fp, "Q%d e_maxApparentRelativeRadiusOfCrossSectionForElectron %lg \n", i++, e_maxApparentRelativeRadiusOfCrossSectionForElectron);
	fprintf(a_fp, "Q%d e_maxApparentRelativeRadiusOfCrossSectionForNecleus %lg \n", i++, e_maxApparentRelativeRadiusOfCrossSectionForNecleus);

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
	fprintf(a_fp, "U%d e_collideMidiMeV %lg\n", i++, e_collideMidiMeV);
	fprintf(a_fp, "U%d e_collideMiniMeV %lg\n", i++, e_collideMiniMeV);
	fprintf(a_fp, "U%d e_usebulletCrossSection %d\n", i++, e_usebulletCrossSection);
	fprintf(a_fp, "U%d e_useProtonScattering %d\n", i++, e_useProtonScattering);
	fprintf(a_fp, "U%d e_useNeutonScattering %d\n", i++, e_useNeutonScattering);
	fprintf(a_fp, "U%d e_useAlphaScattering %d\n", i++, e_useAlphaScattering);
	fprintf(a_fp, "U%d e_useElectronScattering %d\n", i++, e_useElectronScattering);
	fprintf(a_fp, "U%d e_useComptonEffect %d\n", i++, e_useComptonEffect);
	fprintf(a_fp, "U%d e_rateForAlphaParticle %lg\n", i++, e_rateForAlphaParticle);
	fprintf(a_fp, "U%d e_rateForProtonAtBetaPlus %lg\n", i++, e_rateForProtonAtBetaPlus);
	fprintf(a_fp, "U%d e_rateFor2ProtonAtBetaPlus %lg\n", i++, e_rateFor2ProtonAtBetaPlus);
	fprintf(a_fp, "U%d e_rateFor3ProtonAtBetaPlus %lg\n", i++, e_rateFor3ProtonAtBetaPlus);
	fprintf(a_fp, "U%d e_rateForAlphaParticleAtBetaPlus %lg\n", i++, e_rateForAlphaParticleAtBetaPlus);
	fprintf(a_fp, "U%d e_rateForProtonAtEC %lg\n", i++, e_rateForProtonAtEC);
	fprintf(a_fp, "U%d e_rateFor2ProtonAtEC %lg\n", i++, e_rateFor2ProtonAtEC);
	fprintf(a_fp, "U%d e_rateFor3ProtonAtEC %lg\n", i++, e_rateFor3ProtonAtEC);
	fprintf(a_fp, "U%d e_rateForAlphaParticleAtEC %lg\n", i++, e_rateForAlphaParticleAtEC);
	fprintf(a_fp, "U%d e_rateForNeytonAtBetaMinus %lg\n", i++, e_rateForNeytonAtBetaMinus);
	fprintf(a_fp, "U%d e_rateFor2NeytonAtBetaMinus %lg\n", i++, e_rateFor2NeytonAtBetaMinus);
	fprintf(a_fp, "U%d e_rateFor3NeytonAtBetaMinus %lg\n", i++, e_rateFor3NeytonAtBetaMinus);
	fprintf(a_fp, "U%d e_rateFor4NeytonAtBetaMinus %lg\n", i++, e_rateFor4NeytonAtBetaMinus);
	fprintf(a_fp, "U%d e_rateForAlphaParticleAtBetaMinus %lg\n", i++, e_rateForAlphaParticleAtBetaMinus);
	
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
							SET_VALUE("%lg", double, e_collideMidiMeV)
							SET_VALUE("%lg", double, e_collideMiniMeV)
							SET_VALUE("%d", int, e_usebulletCrossSection)
							SET_VALUE("%d", int, e_useProtonScattering)
							SET_VALUE("%d", int, e_useNeutonScattering)
							SET_VALUE("%d", int, e_useAlphaScattering)
							SET_VALUE("%d", int, e_useElectronScattering)
							SET_VALUE("%d", int, e_useComptonEffect)
							SET_VALUE("%lg", double, e_rateForAlphaParticle)
							SET_VALUE("%lg", double, e_rateForProtonAtBetaPlus)
							SET_VALUE("%lg", double, e_rateFor2ProtonAtBetaPlus)
							SET_VALUE("%lg", double, e_rateFor3ProtonAtBetaPlus)
							SET_VALUE("%lg", double, e_rateForAlphaParticleAtBetaPlus)
							SET_VALUE("%lg", double, e_rateForProtonAtEC)
							SET_VALUE("%lg", double, e_rateFor2ProtonAtEC)
							SET_VALUE("%lg", double, e_rateFor3ProtonAtEC)
							SET_VALUE("%lg", double, e_rateForAlphaParticleAtEC)
							SET_VALUE("%lg", double, e_rateForNeytonAtBetaMinus)
							SET_VALUE("%lg", double, e_rateFor2NeytonAtBetaMinus)
							SET_VALUE("%lg", double, e_rateFor3NeytonAtBetaMinus)
							SET_VALUE("%lg", double, e_rateFor4NeytonAtBetaMinus)
							SET_VALUE("%lg", double, e_rateForAlphaParticleAtBetaMinus)
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
#define MIN_DECAY_MODE   101
#define DECAY_MODE_ALPHA 101
#define DECAY_MODE_BETA_PLUS         111 
#define DECAY_MODE_DOUBLE_BETA_PLUS  112
#define DECAY_MODE_BETA_PLUS_ALPHA   113
#define DECAY_MODE_BETA_PLUS_PROTON  114
#define DECAY_MODE_BETA_PLUS_2PROTON 115
#define DECAY_MODE_BETA_PLUS_3PROTON 116
#define DECAY_MODE_ELECTRON_CAPTURE        121
#define DECAY_MODE_DOUBLE_ELECTRON_CAPTURE 122
#define DECAY_MODE_ELECTRON_CAPTURE_ALPHA  123
#define DECAY_MODE_ELECTRON_CAPTURE_PROTON 124
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE         131 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_DOUBLE_ELECTRON_CAPTURE_DEGRADE  132 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA   133 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_PROTON  134 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON 135 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON 136 // only use in the situation, can't use in isotopes table
#define DECAY_MODE_BETA_MINUS              141 
#define DECAY_MODE_DOUBLE_BETA_MINUS       142 
#define DECAY_MODE_BETA_MINUS_AND_ALPHA    143 
#define DECAY_MODE_BETA_MINUS_AND_NEUTRON  144
#define DECAY_MODE_BETA_MINUS_AND_2NEUTRON 145
#define DECAY_MODE_BETA_MINUS_AND_3NEUTRON 146
#define DECAY_MODE_BETA_MINUS_AND_4NEUTRON 147
#define DECAY_MODE_PROTON_EMISSION  151
#define DECAY_MODE_2PROTON_EMISSION 152
#define DECAY_MODE_NEUTRON_EMISSION 161 
#define DECAY_MODE_2NEUTRON_EMISSION 162 
#define DECAY_MODE_3NEUTRON_EMISSION 163 
#define DECAY_MODE_4NEUTRON_EMISSION 164 
#define DECAY_MODE_CLUSTER_DECAY_14C         106014
//#define DECAY_MODE_CLUSTER_DECAY_20O       108020
//#define DECAY_MODE_CLUSTER_DECAY_23F       109023
//#define DECAY_MODE_CLUSTER_DECAY_24Ne      110024
//#define DECAY_MODE_CLUSTER_DECAY_25Ne_24Ne 110225
//#define DECAY_MODE_CLUSTER_DECAY_26Ne_24Ne 110226
//#define DECAY_MODE_CLUSTER_DECAY_28Mg      112028
//#define DECAY_MODE_CLUSTER_DECAY_30Mg      112030
//#define DECAY_MODE_CLUSTER_DECAY_30Mg_28Mg 112230
//#define DECAY_MODE_CLUSTER_DECAY_32Si      114032
//#define DECAY_MODE_CLUSTER_DECAY_34Si      114034
#define DECAY_MODE_SELF_FISSION_80KR         136080 // 180Tl -> 100Ru +  80Kr

extern const char * getCollideName(int a_collideType);
extern const char * getDecayModeText(int a_decayMode)
{
	const char * modeText = "unkown";
	switch(a_decayMode){
		case DECAY_MODE_ALPHA: modeText = "ALPHA"; break;
		case DECAY_MODE_BETA_PLUS: modeText = "BETA_PLUS"; break;
		case DECAY_MODE_DOUBLE_BETA_PLUS: modeText = "DOUBLE_BETA_PLUS"; break;
		case DECAY_MODE_BETA_PLUS_ALPHA: modeText = "BETA_PLUS_ALPHA"; break;
		case DECAY_MODE_BETA_PLUS_PROTON: modeText = "BETA_PLUS_PROTON"; break;
		case DECAY_MODE_BETA_PLUS_2PROTON: modeText = "BETA_PLUS_2PROTON"; break;
		case DECAY_MODE_BETA_PLUS_3PROTON: modeText = "BETA_PLUS_3PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE: modeText = "ELECTRON_CAPTURE"; break;
		case DECAY_MODE_DOUBLE_ELECTRON_CAPTURE: modeText = "DOUBLE_ELECTRON_CAPTURE"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_ALPHA: modeText = "ELECTRON_CAPTURE_ALPHA"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_PROTON: modeText = "ELECTRON_CAPTURE_PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE: modeText = "ELECTRON_CAPTURE_DEGRADE"; break;
		case DECAY_MODE_DOUBLE_ELECTRON_CAPTURE_DEGRADE: modeText = "DOUBLE_ELECTRON_CAPTURE_DEGRADE"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA: modeText = "ELECTRON_CAPTURE_DEGRADE_ALPHA"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_PROTON: modeText = "ELECTRON_CAPTURE_DEGRADE_PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON: modeText = "ELECTRON_CAPTURE_DEGRADE_2PROTON"; break;
		case DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON: modeText = "ELECTRON_CAPTURE_DEGRADE_3PROTON"; break;
		case DECAY_MODE_BETA_MINUS: modeText = "BETA_MINUS"; break;
		case DECAY_MODE_DOUBLE_BETA_MINUS: modeText = "DOUBLE_BETA_MINUS"; break;
		case DECAY_MODE_BETA_MINUS_AND_ALPHA: modeText = "BETA_MINUS_AND_ALPHA"; break;
		case DECAY_MODE_BETA_MINUS_AND_NEUTRON: modeText = "BETA_MINUS_AND_NEUTRON"; break;
		case DECAY_MODE_BETA_MINUS_AND_2NEUTRON: modeText = "BETA_MINUS_AND_2NEUTRON"; break;
		case DECAY_MODE_BETA_MINUS_AND_3NEUTRON: modeText = "BETA_MINUS_AND_3NEUTRON"; break;
		case DECAY_MODE_BETA_MINUS_AND_4NEUTRON: modeText = "BETA_MINUS_AND_4NEUTRON"; break;
		case DECAY_MODE_PROTON_EMISSION: modeText = "PROTON_EMISSION"; break;
		case DECAY_MODE_2PROTON_EMISSION: modeText = "2PROTON_EMISSION"; break;
		case DECAY_MODE_NEUTRON_EMISSION: modeText = "NEUTRON_EMISSION"; break;
		case DECAY_MODE_2NEUTRON_EMISSION: modeText = "2NEUTRON_EMISSION"; break;
		case DECAY_MODE_3NEUTRON_EMISSION: modeText = "3NEUTRON_EMISSION"; break;
		case DECAY_MODE_4NEUTRON_EMISSION: modeText = "4NEUTRON_EMISSION"; break;
		case DECAY_MODE_CLUSTER_DECAY_14C: modeText = "CLUSTER_DECAY_14C"; break;
		case DECAY_MODE_SELF_FISSION_80KR: modeText = "SELF_FISSION_80KR"; break;
		default:
			if(a_decayMode < MIN_DECAY_MODE){
				modeText = getCollideName(a_decayMode);
			}else{
				fprintf(stderr, "FATAL:%s(%d):UNKNOWN DECAY %d\n", __FUNCTION__, __LINE__, a_decayMode);
				exit(1);
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
	double relativeNucleusRadius;// pow(massNumber, 1.0/3.0)
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
#define ATOMICNUMBER_CARBON 6
#define MASSNUMBER_CARBON14 14
#define ATOMICNUMBER_Fe 26
#define ATOMICNUMBER_Co 27
#define ATOMICNUMBER_Ni 28
#define ATOMICNUMBER_Cu 29
#define ATOMICNUMBER_Zn 30
#define ATOMICNUMBER_Kr 36
#define MASSNUMBER_Kr80 80
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
			if(a_massNumber == 1){
				isotopePropertyPtr->relativeNucleusRadius = e_radiusOfProtonByRmsCharge / e_r0;
			}else if(a_massNumber < pow(e_radiusOfDeuteronByRmsCharge / e_r0, 3.0)){
				//5.68 = pow(e_radiusOfDeuteronByRmsCharge / e_r0, 3.0)
				isotopePropertyPtr->relativeNucleusRadius = e_radiusOfDeuteronByRmsCharge / e_r0;
			}else{
				isotopePropertyPtr->relativeNucleusRadius = pow(a_massNumber, 1.0 / 3.0);
			}
		}else if(a_massNumber == MASSNUMBER_ELECTRON && a_atomicNumber == ATOMICNUMBER_ELECTRON){
			isotopePropertyPtr->relativeNucleusRadius = e_radiusOfElectronByWeakBoson / e_r0;
		}else{
			isotopePropertyPtr->relativeNucleusRadius = 0.0;
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
struct isotopeProperty * e_C14Ptr;
struct isotopeProperty * e_Kr80Ptr;

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
		
		initAtomProperty(1, "H", "hydrogen", "水素", 1, 7, 1, 3, 28.836);
		initIsotopeProperty(1, 1, 1.00782503207, 0.0, HLU_STABLE, 0.999885);
		initIsotopeProperty(1, 2, 2.0141017778, 0.0, HLU_STABLE, 0.000115); 
		initIsotopeProperty(1, 3, 3.0160492777, 12.32, HLU_YEAR, 0.0);
		initIsotopeProperty(1, 4, 4.02643, 1.39e-22, HLU_SECOND, 0.0);
		initIsotopeProperty(1, 5, 5.03531, 9.1e-22, HLU_SECOND, 0.0);
		initIsotopeProperty(1, 6, 6.04496, 2.9e-22, HLU_SECOND, 0.0);
		initIsotopeProperty(1, 7, 7.05275, 2.3e-23, HLU_SECOND, 0.0);
		initIsotopePropertySymbols(1, 1, "H", "hydrogen", "水素", "p");
		initIsotopePropertySymbols(1, 2, "D", "deuterium", "重水素", "d");
		initIsotopePropertySymbols(1, 3, "T", "tritium", "三重水素", "t");
		initDecayMode(1, 3, DECAY_MODE_BETA_MINUS);
		initDecayMode(1, 4, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode(1, 5, DECAY_MODE_2NEUTRON_EMISSION);
		initDecayMode2(1, 6, DECAY_MODE_3NEUTRON_EMISSION, 50.0, DECAY_MODE_4NEUTRON_EMISSION);
		initDecayMode(1, 7, DECAY_MODE_4NEUTRON_EMISSION);
		
		initAtomProperty(2, "He", "helium", "ヘリウム", 3, 10, 3, 4, 20.786);
		initIsotopeProperty(2, 2, 2.015894, 1.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 3, 3.0160293191, 0.0, HLU_STABLE, 0.00000134);
		initIsotopeProperty(2, 4, 4.00260325415, 0.0, HLU_STABLE, 0.99999866);
		initIsotopeProperty(2, 5, 5.01222, 700.0e-24, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 6, 6.0188891, 806.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 7, 7.028021, 2.9e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 8, 8.033922, 119.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 9, 9.04395, 7.0e-21, HLU_SECOND, 0.0);
		initIsotopeProperty(2, 10, 10.05240, 2.7e-21, HLU_SECOND, 0.0);
		initDecayMode2(2, 2, DECAY_MODE_PROTON_EMISSION, 99.99, DECAY_MODE_BETA_PLUS);
		initDecayMode(2, 5, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode2(2, 6, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode(2, 7, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode2(2, 8, DECAY_MODE_BETA_MINUS, 84.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);//ignore "β-, fission (0.9%, 5He + 3H)"
		initDecayMode(2, 9, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode(2, 10, DECAY_MODE_2NEUTRON_EMISSION);

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
		initDecayMode2(3, 9, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 50.8, DECAY_MODE_BETA_MINUS);
		initDecayMode(3, 10, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode5(3, 11, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 84.9, DECAY_MODE_BETA_MINUS, 8.07, DECAY_MODE_BETA_MINUS_AND_2NEUTRON, 4.1, DECAY_MODE_BETA_MINUS_AND_3NEUTRON, 1.9, DECAY_MODE_BETA_MINUS_AND_ALPHA); //ignore "β-, fission (.014%) 8Li, 3H", "β-, fission (.013%), 9Li 2H"

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
		initDecayMode2(4, 12, DECAY_MODE_BETA_MINUS, 99.48, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(4, 13, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode3(4, 14, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 81.0, DECAY_MODE_BETA_MINUS, 14.0, DECAY_MODE_BETA_MINUS_AND_2NEUTRON);

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
		initDecayMode2(5, 13, DECAY_MODE_BETA_MINUS, 99.72, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(5, 14, DECAY_MODE_BETA_MINUS, 93.96, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode3(5, 15, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 93.6, DECAY_MODE_BETA_MINUS, 6.0, 	DECAY_MODE_BETA_MINUS_AND_2NEUTRON);
		initDecayMode(5, 16, DECAY_MODE_NEUTRON_EMISSION);
		initDecayMode5(5, 17, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 63.0, DECAY_MODE_BETA_MINUS, 22.1, DECAY_MODE_BETA_MINUS_AND_2NEUTRON, 11.0, DECAY_MODE_BETA_MINUS_AND_3NEUTRON, 3.5, DECAY_MODE_BETA_MINUS_AND_4NEUTRON);
		initDecayMode(5, 18, DECAY_MODE_NEUTRON_EMISSION);
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
		initDecayMode2(6, 16, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 97.9, DECAY_MODE_BETA_MINUS);
		initDecayMode2(6, 17, DECAY_MODE_BETA_MINUS, 71.59, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(6, 18, DECAY_MODE_BETA_MINUS, 68.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode3(6, 19, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 47.0, DECAY_MODE_BETA_MINUS, 46.0, DECAY_MODE_BETA_MINUS_AND_2NEUTRON);
		initDecayMode2(6, 20, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 72.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(6, 21, DECAY_MODE_NEUTRON_EMISSION);
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
		initDecayMode3(7, 17, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 95.0, DECAY_MODE_BETA_MINUS, 4.99, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode3(7, 18, DECAY_MODE_BETA_MINUS, 76.9, DECAY_MODE_BETA_MINUS_AND_ALPHA, 12.2, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(7, 19, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 54.6, DECAY_MODE_BETA_MINUS);
		initDecayMode2(7, 20, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 56.99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(7, 21, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 80.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(7, 22, DECAY_MODE_BETA_MINUS, 65.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(7, 23, DECAY_MODE_BETA_MINUS);
		initDecayMode(7, 24, DECAY_MODE_NEUTRON_EMISSION);

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
		initDecayMode2(8, 22, DECAY_MODE_BETA_MINUS, 78.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(8, 23, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 57.99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(8, 24, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 57.99, DECAY_MODE_BETA_MINUS);

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
		initDecayMode2(9, 22, DECAY_MODE_BETA_MINUS, 89.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(9, 23, DECAY_MODE_BETA_MINUS, 86.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(9, 24, DECAY_MODE_BETA_MINUS, 94.1, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(9, 25, DECAY_MODE_BETA_MINUS, 76.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(9, 26, DECAY_MODE_BETA_MINUS, 68.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(9, 27, DECAY_MODE_BETA_MINUS);
		initDecayMode(9, 28, DECAY_MODE_NEUTRON_EMISSION);
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
		initDecayMode2(10, 26, DECAY_MODE_BETA_MINUS, 99.87, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(10, 27, DECAY_MODE_BETA_MINUS, 98.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(10, 28, DECAY_MODE_BETA_MINUS, 78.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(10, 29, DECAY_MODE_BETA_MINUS);
		initDecayMode(10, 30, DECAY_MODE_BETA_MINUS);
		initDecayMode2(10, 31, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(10, 32, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 50.0, DECAY_MODE_BETA_MINUS);

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
		initDecayMode2(11, 27, DECAY_MODE_BETA_MINUS, 99.87, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(11, 28, DECAY_MODE_BETA_MINUS, 99.421, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(11, 29, DECAY_MODE_BETA_MINUS, 74.09, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode3(11, 30, DECAY_MODE_BETA_MINUS, 68.83, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 30.0, 	DECAY_MODE_BETA_MINUS_AND_2NEUTRON); //	DECAY_MODE_BETA_MINUS_AND_ALPHA
		initDecayMode3(11, 31, DECAY_MODE_BETA_MINUS, 62.05, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 33.0, 	DECAY_MODE_BETA_MINUS_AND_2NEUTRON);// DECAY_MODE_BETA_MINUS_AND_3NEUTRON
		initDecayMode3(11, 32, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 35.0, 	DECAY_MODE_BETA_MINUS_AND_2NEUTRON);
		initDecayMode3(11, 33, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 52.0, DECAY_MODE_BETA_MINUS, 36.0, 	DECAY_MODE_BETA_MINUS_AND_2NEUTRON);
		initDecayMode3(11, 34, DECAY_MODE_BETA_MINUS_AND_2NEUTRON, 50.0, DECAY_MODE_BETA_MINUS, 35.0, 	DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(11, 35, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(12, 30, DECAY_MODE_BETA_MINUS, 99.94, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(12, 31, DECAY_MODE_BETA_MINUS, 98.3, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(12, 32, DECAY_MODE_BETA_MINUS, 97.6, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(12, 33, DECAY_MODE_BETA_MINUS, 83.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(12, 34, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(12, 35, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 52.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(12, 36, DECAY_MODE_BETA_MINUS);
		initDecayMode2(12, 37, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(13, 31, DECAY_MODE_BETA_MINUS, 98.4, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(13, 32, DECAY_MODE_BETA_MINUS, 99.3, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(13, 33, DECAY_MODE_BETA_MINUS, 91.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(13, 34, DECAY_MODE_BETA_MINUS, 87.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(13, 35, DECAY_MODE_BETA_MINUS, 74.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(13, 36, DECAY_MODE_BETA_MINUS, 69.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
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
		initDecayMode2(14, 35, DECAY_MODE_BETA_MINUS, 94.74, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(14, 36, DECAY_MODE_BETA_MINUS, 88.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(14, 37, DECAY_MODE_BETA_MINUS, 83.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(14, 38, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 50.0, DECAY_MODE_BETA_MINUS);
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
		initDecayMode2(15, 38, DECAY_MODE_BETA_MINUS, 88.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(15, 39, DECAY_MODE_BETA_MINUS, 74.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(15, 40, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(15, 41, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(15, 42, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(15, 43, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
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
		initDecayMode2(16, 41, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(16, 42, DECAY_MODE_BETA_MINUS, 96.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(16, 43, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(16, 44, DECAY_MODE_BETA_MINUS, 82.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(16, 45, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 54.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 46, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 47, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 48, DECAY_MODE_BETA_MINUS);
		initDecayMode(16, 49, DECAY_MODE_NEUTRON_EMISSION);

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
		initDecayMode2(17, 43, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(17, 44, DECAY_MODE_BETA_MINUS, 92.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(17, 45, DECAY_MODE_BETA_MINUS, 76.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(17, 46, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 60.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(17, 47, DECAY_MODE_BETA_MINUS, 97.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
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
		initDecayMode2(18, 47, DECAY_MODE_BETA_MINUS, 99.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(18, 48, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 49, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 50, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 51, DECAY_MODE_BETA_MINUS);
		initDecayMode(18, 52, DECAY_MODE_BETA_MINUS);
		initDecayMode2(18, 53, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(19, 48, DECAY_MODE_BETA_MINUS, 98.86, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(19, 49, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 86.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(19, 50, DECAY_MODE_BETA_MINUS, 71.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(19, 51, DECAY_MODE_BETA_MINUS, 53.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode3(19, 52, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 64.0, DECAY_MODE_BETA_MINUS_AND_2NEUTRON, 21.0, DECAY_MODE_BETA_MINUS);
		initDecayMode3(19, 53, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 67.0, DECAY_MODE_BETA_MINUS_AND_2NEUTRON, 17.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(19, 54, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(19, 55, DECAY_MODE_BETA_MINUS, 99.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(20, 51, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(20, 52, DECAY_MODE_BETA_MINUS, 98.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(20, 53, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(20, 54, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 50.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 55, DECAY_MODE_BETA_MINUS);
		initDecayMode(20, 56, DECAY_MODE_BETA_MINUS);
		initDecayMode2(20, 57, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(21, 53, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(21,54, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(21,55, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(21, 56, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 57, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 58, DECAY_MODE_BETA_MINUS);
		initDecayMode(21, 59, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(22, 56, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(22, 57, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(22, 58, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 59, DECAY_MODE_BETA_MINUS);
		initDecayMode(22, 60, DECAY_MODE_BETA_MINUS);
		initDecayMode2(22, 61, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(23, 56, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(23, 57, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(23, 58, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(23, 59, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(23, 60, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
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
		initDecayMode2(28, 72, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(29, 73, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(29, 74, DECAY_MODE_BETA_MINUS);
		initDecayMode2(29, 75, DECAY_MODE_BETA_MINUS, 96.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(29, 76, DECAY_MODE_BETA_MINUS, 97.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		
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
		initDecayMode2(31, 79, DECAY_MODE_BETA_MINUS, 99.911, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(31, 80, DECAY_MODE_BETA_MINUS, 99.11, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(31, 81, DECAY_MODE_BETA_MINUS, 88.11, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(31, 82, DECAY_MODE_BETA_MINUS, 78.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(31, 83, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(31, 84, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 70.0, DECAY_MODE_BETA_MINUS);

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
		initDecayMode2(32, 84, DECAY_MODE_BETA_MINUS, 89.2, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(32, 85, DECAY_MODE_BETA_MINUS, 86.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(32, 86, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 50.0, DECAY_MODE_BETA_MINUS);

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
		initDecayMode2(33, 84, DECAY_MODE_BETA_MINUS, 99.721, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(33, 85, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 59.4, DECAY_MODE_BETA_MINUS);
		initDecayMode2(33, 86, DECAY_MODE_BETA_MINUS, 67.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(33, 87, DECAY_MODE_BETA_MINUS, 84.6, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(33, 88, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
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
		initDecayMode2(34, 87, DECAY_MODE_BETA_MINUS, 99.64, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(34, 88, DECAY_MODE_BETA_MINUS, 99.01, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(34, 89, DECAY_MODE_BETA_MINUS, 92.2, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(34, 90, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 50.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(34, 91, DECAY_MODE_BETA_MINUS, 79.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
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
		initDecayMode2(35, 87, DECAY_MODE_BETA_MINUS, 97.48, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 88, DECAY_MODE_BETA_MINUS, 93.42, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 89, DECAY_MODE_BETA_MINUS, 86.2, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 90, DECAY_MODE_BETA_MINUS, 74.8, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 91, DECAY_MODE_BETA_MINUS, 80.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 92, DECAY_MODE_BETA_MINUS, 66.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 93, DECAY_MODE_BETA_MINUS, 89.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(35, 94, DECAY_MODE_BETA_MINUS, 70.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initDecayMode2(36, 92, DECAY_MODE_BETA_MINUS, 99.96, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(36, 93, DECAY_MODE_BETA_MINUS, 98.05, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(36, 94, DECAY_MODE_BETA_MINUS, 94.3, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(36, 95, DECAY_MODE_BETA_MINUS);
		initDecayMode(36, 96, DECAY_MODE_BETA_MINUS);
		initDecayMode2(36, 97, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(37, "Rb", "rubidium", "ルビジウム", 72, 102, 72, 102, 31.060);
		initIsotopeProperty(37, 72, 71.95908, 1.5e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 73, 72.95056, 30.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 74, 73.944265, 64.76e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 75, 74.938570, 19.0, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 76, 75.9350722, 36.5, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 77, 76.930408, 3.77, HLU_MINUTE, 0.0);
		initIsotopeProperty(37, 78, 77.928141, 17.66, HLU_MINUTE, 0.0);
		initIsotopeProperty(37, 79, 78.923989, 22.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(37, 80, 79.922519, 33.4, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 81, 80.918996, 4.570, HLU_HOUR, 0.0);
		initIsotopeProperty(37, 82, 81.9182086, 1.273, HLU_MINUTE, 0.0);
		initIsotopeProperty(37, 83, 82.915110, 86.2, HLU_DAY, 0.0);
		initIsotopeProperty(37, 84, 83.914385, 33.1, HLU_DAY, 0.0);
		initIsotopeProperty(37, 85, 84.911789738, 0.0, HLU_STABLE, 0.7217);
		initIsotopeProperty(37, 86, 85.91116742, 18.642, HLU_DAY, 0.0);
		initIsotopeProperty(37, 87, 86.909180527, 4.923e10, HLU_YEAR, 0.2783);
		initIsotopeProperty(37, 88, 87.91131559, 17.773, HLU_MINUTE, 0.0);
		initIsotopeProperty(37, 89, 88.912278, 15.15, HLU_MINUTE, 0.0);
		initIsotopeProperty(37, 90, 89.914802, 158.0, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 91, 90.916537, 58.4, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 92, 91.919729, 4.492, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 93, 92.922042, 5.84, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 94, 93.926405, 2.702, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 95, 94.929303, 377.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 96, 95.93427, 202.8e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 97, 96.93735, 169.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 98, 97.94179, 114e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 99, 98.94538, 50.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 100, 99.94987, 51.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 101, 100.95320, 32.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(37, 102, 101.95887, 37.0e-3, HLU_SECOND, 0.0);
		initDecayMode(37, 72, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(37, 73, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(37, 74, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 75, DECAY_MODE_BETA_PLUS);
		initDecayMode2(37, 76, DECAY_MODE_BETA_PLUS, 100.0 - 3.8e-7, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(37, 77, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 78, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 79, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 80, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 81, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 82, DECAY_MODE_BETA_PLUS);
		initDecayMode(37, 83, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(37, 84, DECAY_MODE_BETA_PLUS, 96.2, DECAY_MODE_BETA_MINUS);
		initDecayMode2(37, 86, DECAY_MODE_BETA_MINUS, 99.9948, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(37, 87, DECAY_MODE_BETA_MINUS);
		initDecayMode(37, 88, DECAY_MODE_BETA_MINUS);
		initDecayMode(37, 89, DECAY_MODE_BETA_MINUS);
		initDecayMode(37, 90, DECAY_MODE_BETA_MINUS);
		initDecayMode(37, 91, DECAY_MODE_BETA_MINUS);
		initDecayMode2(37, 92, DECAY_MODE_BETA_MINUS, 99.98, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(37, 93, DECAY_MODE_BETA_MINUS, 98.65, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(37, 94, DECAY_MODE_BETA_MINUS, 89.99, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(37, 95, DECAY_MODE_BETA_MINUS, 91.27, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(37, 96, DECAY_MODE_BETA_MINUS, 86.6, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(37, 97, DECAY_MODE_BETA_MINUS, 74.3, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode3(37, 98, DECAY_MODE_BETA_MINUS, 86.14, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 13.8, DECAY_MODE_BETA_MINUS_AND_2NEUTRON);
		initDecayMode2(37, 99, DECAY_MODE_BETA_MINUS, 84.1, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode3(37, 100, DECAY_MODE_BETA_MINUS, 94.25, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 5.6, DECAY_MODE_BETA_MINUS_AND_2NEUTRON);
		initDecayMode2(37, 101, DECAY_MODE_BETA_MINUS, 69.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(37, 102, DECAY_MODE_BETA_MINUS, 82.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(38, "Sr", "strontium", "ストロンチウム", 73, 104, 73, 104, 26.4);
		initIsotopeProperty(38, 73, 72.96597, 25.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 74, 73.95631, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 75, 74.94995, 88.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 76, 75.94177, 7.89, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 77, 76.937945, 9.0, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 78, 77.932180, 159, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 79, 78.929708, 2.25, HLU_MINUTE, 0.0);
		initIsotopeProperty(38, 80, 79.924521, 106.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(38, 81, 80.923212, 22.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(38, 82, 81.918402, 25.36, HLU_DAY, 0.0);
		initIsotopeProperty(38, 83, 82.917557, 32.41, HLU_HOUR, 0.0);
		initIsotopeProperty(38, 84, 83.913425, 0.0, HLU_STABLE, 0.0056);
		initIsotopeProperty(38, 85, 84.912933, 64.853, HLU_DAY, 0.0);
		initIsotopeProperty(38, 86, 85.9092607309, 0.0, HLU_STABLE, 0.0986);
		initIsotopeProperty(38, 87, 86.9088774970, 0.0, HLU_STABLE, 0.0700);
		initIsotopeProperty(38, 88, 87.9056122571, 0.0, HLU_STABLE, 0.8258);
		initIsotopeProperty(38, 89, 88.9074507, 50.57, HLU_DAY, 0.0);
		initIsotopeProperty(38, 90, 89.907738, 28.90, HLU_YEAR, 0.0);
		initIsotopeProperty(38, 91, 90.910203, 9.63, HLU_HOUR, 0.0);
		initIsotopeProperty(38, 92, 91.911038, 2.66, HLU_HOUR, 0.0);
		initIsotopeProperty(38, 93, 92.914026, 7.423, HLU_MINUTE, 0.0);
		initIsotopeProperty(38, 94, 93.915361, 75.3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 95, 94.919359, 23.90, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 96, 95.921697, 1.07, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 97, 96.926153, 429.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 98, 97.928453, 0.653, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 99, 98.93324, 0.269, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 100, 99.93535, 202.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 101, 100.94052, 118.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 102, 101.94302, 69.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 103, 102.94895, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(38, 104, 103.95233, 30.0e-3, HLU_SECOND, 0.0);
		//initIsotopeProperty(38, 105, 104.95858, 20.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(38, 73, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(38, 74, DECAY_MODE_BETA_PLUS);
		initDecayMode2(38, 75, DECAY_MODE_BETA_PLUS, 93.5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(38, 76, DECAY_MODE_BETA_PLUS);
		initDecayMode2(38, 77, DECAY_MODE_BETA_PLUS, 99.75, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(38, 78, DECAY_MODE_BETA_PLUS);
		initDecayMode(38, 79, DECAY_MODE_BETA_PLUS);
		initDecayMode(38, 80, DECAY_MODE_BETA_PLUS);
		initDecayMode(38, 81, DECAY_MODE_BETA_PLUS);
		initDecayMode(38, 82, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(38, 83, DECAY_MODE_BETA_PLUS);
		initDecayMode(38, 85, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(38, 89, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 90, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 91, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 92, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 93, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 94, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 95, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 96, DECAY_MODE_BETA_MINUS);
		initDecayMode2(38, 97, DECAY_MODE_BETA_MINUS, 99.95, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(38, 98, DECAY_MODE_BETA_MINUS, 99.75, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(38, 99, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(38, 100, DECAY_MODE_BETA_MINUS, 99.02, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(38, 101, DECAY_MODE_BETA_MINUS, 97.63, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(38, 102, DECAY_MODE_BETA_MINUS, 94.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(38, 103, DECAY_MODE_BETA_MINUS);
		initDecayMode(38, 104, DECAY_MODE_BETA_MINUS);

		initAtomProperty(39, "Y", "yttrium", "イットリウム", 77, 106, 77, 106, 26.53);
		initIsotopeProperty(39, 77, 76.94965, 63.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 78, 77.94361, 54.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 79, 78.93735, 14.8, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 80, 79.93428, 30.1, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 81, 80.92913, 70.4, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 82, 81.92679, 8.30, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 83, 82.92235, 7.08, HLU_MINUTE, 0.0);
		initIsotopeProperty(39, 84, 83.92039, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 85, 84.916433, 2.68, HLU_HOUR, 0.0);
		initIsotopeProperty(39, 86, 85.914886, 14.74, HLU_HOUR, 0.0);
		initIsotopeProperty(39, 87, 86.9108757, 79.8, HLU_HOUR, 0.0);
		initIsotopeProperty(39, 88, 87.9095011, 106.616, HLU_DAY, 0.0);
		initIsotopeProperty(39, 89, 88.9058483, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(39, 90, 89.9071519, 64.053, HLU_HOUR, 0.0);
		initIsotopeProperty(39, 91, 90.907305, 58.51, HLU_DAY, 0.0);
		initIsotopeProperty(39, 92, 91.908949, 3.54, HLU_HOUR, 0.0);
		initIsotopeProperty(39, 93, 92.909583, 10.18, HLU_HOUR, 0.0);
		initIsotopeProperty(39, 94, 93.911595, 18.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(39, 95, 94.912821, 10.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(39, 96, 95.915891, 5.34, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 97, 96.918134, 3.75, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 98, 97.922203, 0.548, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 99, 98.924636, 1.470, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 100, 99.92776, 735.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 101, 100.93031, 426.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 102, 101.93356, 0.30, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 103, 102.93673, 224.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 104, 103.94105, 180.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 105, 104.94487, 60.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(39, 106, 105.94979, 50.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(39, 77, DECAY_MODE_PROTON_EMISSION, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 78, DECAY_MODE_BETA_PLUS);
		initDecayMode2(39, 79, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(39, 80, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 81, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 82, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 83, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 84, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 85, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 86, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 87, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 88, DECAY_MODE_BETA_PLUS);
		initDecayMode(39, 90, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 91, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 92, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 93, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 94, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 95, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 96, DECAY_MODE_BETA_MINUS);
		initDecayMode2(39, 97, DECAY_MODE_BETA_MINUS, 99.942, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(39, 98, DECAY_MODE_BETA_MINUS, 99.669, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(39, 99, DECAY_MODE_BETA_MINUS, 98.1, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(39, 100, DECAY_MODE_BETA_MINUS, 98.98, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(39, 101, DECAY_MODE_BETA_MINUS, 98.06, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(39, 102, DECAY_MODE_BETA_MINUS, 95.1, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(39, 103, DECAY_MODE_BETA_MINUS, 91.7, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(39, 104, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 105, DECAY_MODE_BETA_MINUS);
		initDecayMode(39, 106, DECAY_MODE_BETA_MINUS);

		initAtomProperty(40, "Zr", "zirconium", "ジルコニウム", 79, 108, 79, 108, 25.36);
		initIsotopeProperty(40, 79, 78.94916, 56.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 80, 79.9404, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 81, 80.93721, 5.5, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 82, 81.93109, 32.0, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 83, 82.92865, 41.6, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 84, 83.92325, 25.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(40, 85, 84.92147, 7.86, HLU_MINUTE, 0.0);
		initIsotopeProperty(40, 86, 85.91647, 16.5, HLU_HOUR, 0.0);
		initIsotopeProperty(40, 87, 86.914816, 1.68, HLU_HOUR, 0.0);
		initIsotopeProperty(40, 88, 87.910227, 83.4, HLU_DAY, 0.0);
		initIsotopeProperty(40, 89, 88.908890, 78.41, HLU_HOUR, 0.0);
		initIsotopeProperty(40, 90, 89.9047044, 0.0, HLU_STABLE, 0.5145);
		initIsotopeProperty(40, 91, 90.9056458, 0.0, HLU_STABLE, 0.1122);
		initIsotopeProperty(40, 92, 91.9050408, 0.0, HLU_STABLE, 0.1715);
		initIsotopeProperty(40, 93, 92.9064760, 1.53e6, HLU_YEAR, 0.0);
		initIsotopeProperty(40, 94, 93.9063152, 0.0, HLU_STABLE, 0.1738);
		initIsotopeProperty(40, 95, 94.9080426, 64.032, HLU_DAY, 0.0);
		initIsotopeProperty(40, 96, 95.9082734, 20.0e18, HLU_YEAR, 0.0280);
		initIsotopeProperty(40, 97, 96.9109531, 16.744, HLU_HOUR, 0.0);
		initIsotopeProperty(40, 98, 97.912735, 30.7, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 99, 98.916512, 2.1, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 100, 99.91776, 7.1, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 101, 100.92114, 2.3, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 102, 101.92298, 2.9, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 103, 102.92660, 1.3, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 104, 103.92878, 1.2, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 105, 104.93305, 0.6, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 106, 105.93591, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 107, 106.94075, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(40, 108, 107.94396, 80.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(40, 79, DECAY_MODE_BETA_PLUS_PROTON, 50.0, DECAY_MODE_BETA_PLUS);
		initDecayMode(40, 80, DECAY_MODE_BETA_PLUS);
		initDecayMode2(40, 81, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(40, 82, DECAY_MODE_BETA_PLUS);
		initDecayMode2(40, 83, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(40, 84, DECAY_MODE_BETA_PLUS);
		initDecayMode(40, 85, DECAY_MODE_BETA_PLUS);
		initDecayMode(40, 86, DECAY_MODE_BETA_PLUS);
		initDecayMode(40, 87, DECAY_MODE_BETA_PLUS);
		initDecayMode(40, 88, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(40, 89, DECAY_MODE_BETA_PLUS);
		initDecayMode(40, 93, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 95, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 96, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(40, 97, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 98, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 99, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 100, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 101, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 102, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 103, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 104, DECAY_MODE_BETA_MINUS);
		initDecayMode2(40, 105, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(40, 106, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 107, DECAY_MODE_BETA_MINUS);
		initDecayMode(40, 108, DECAY_MODE_BETA_MINUS);

		initAtomProperty(41, "Nb", "niobium", "ニオブ", 82, 110, 82, 110, 24.6);
		initIsotopeProperty(41, 82, 81.94313, 51.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 83, 82.93671, 4.1, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 84, 83.93357, 9.8, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 85, 84.92791, 20.9, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 86, 85.92504, 88.0, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 87, 86.92036, 3.75, HLU_MINUTE, 0.0);
		initIsotopeProperty(41, 88, 87.91833, 14.55, HLU_MINUTE, 0.0);
		initIsotopeProperty(41, 89, 88.913418, 2.03, HLU_HOUR, 0.0);
		initIsotopeProperty(41, 90, 89.911265, 14.60, HLU_HOUR, 0.0);
		initIsotopeProperty(41, 91, 90.906996, 680.0, HLU_YEAR, 0.0);
		initIsotopeProperty(41, 92, 91.907194, 3.47e7, HLU_YEAR, 0.0);
		initIsotopeProperty(41, 93, 92.9063781, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(41, 94, 93.9072839, 2.03e4, HLU_YEAR, 0.0);
		initIsotopeProperty(41, 95, 94.9068358, 34.991, HLU_DAY, 0.0);
		initIsotopeProperty(41, 96, 95.908101, 23.35, HLU_HOUR, 0.0);
		initIsotopeProperty(41, 97, 96.9080986, 72.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(41, 98, 97.910328, 2.86, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 99, 98.911618, 15.0, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 100, 99.914182, 1.5, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 101, 100.915252, 7.1, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 102, 101.91804, 1.3, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 103, 102.91914, 1.5, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 104, 103.92246, 4.9, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 105, 104.92394, 2.95, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 106, 105.92797, 920.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 107, 106.93031, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 108, 107.93484, 0.193, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 109, 108.93763, 190.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(41, 110, 109.94244, 170.0e-3, HLU_SECOND, 0.0);
		initDecayMode(41, 82, DECAY_MODE_BETA_PLUS);
		initDecayMode(41, 83, DECAY_MODE_BETA_PLUS);
		initDecayMode2(41, 84, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(41, 85, DECAY_MODE_BETA_PLUS);
		initDecayMode(41, 86, DECAY_MODE_BETA_PLUS);
		initDecayMode(41, 87, DECAY_MODE_BETA_PLUS);
		initDecayMode(41, 88, DECAY_MODE_BETA_PLUS);
		initDecayMode(41, 89, DECAY_MODE_BETA_PLUS);
		initDecayMode(41, 90, DECAY_MODE_BETA_PLUS);
		initDecayMode2(41, 91, DECAY_MODE_ELECTRON_CAPTURE, 99.98, DECAY_MODE_BETA_PLUS);
		initDecayMode2(41, 92, DECAY_MODE_BETA_PLUS, 99.95, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 94, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 95, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 96, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 97, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 98, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 99, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 100, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 101, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 102, DECAY_MODE_BETA_MINUS);
		initDecayMode(41, 103, DECAY_MODE_BETA_MINUS);
		initDecayMode2(41, 104, DECAY_MODE_BETA_MINUS, 99.94, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(41, 105, DECAY_MODE_BETA_MINUS, 98.3, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(41, 106, DECAY_MODE_BETA_MINUS, 95.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(41, 107, DECAY_MODE_BETA_MINUS, 94.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(41, 108, DECAY_MODE_BETA_MINUS, 93.8, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(41, 109, DECAY_MODE_BETA_MINUS, 69.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(41, 110, DECAY_MODE_BETA_MINUS, 60.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(42, "Mo", "molybdenum", "モリブデン", 84, 113, 84, 113, 24.06);
		initIsotopeProperty(42, 84, 83.94009, 3.8, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 85, 84.93655, 3.2, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 86, 85.93070, 19.6, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 87, 86.92733, 14.05, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 88, 87.921953, 8.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(42, 89, 88.919480, 2.11, HLU_MINUTE, 0.0);
		initIsotopeProperty(42, 90, 89.913937, 5.56, HLU_HOUR, 0.0);
		initIsotopeProperty(42, 91, 90.911750, 15.49, HLU_MINUTE, 0.0);
		initIsotopeProperty(42, 92, 91.906811, 0.0, HLU_STABLE, 0.14649);
		initIsotopeProperty(42, 93, 92.906813, 4000.0, HLU_YEAR, 0.0);
		initIsotopeProperty(42, 94, 93.9050883, 0.0, HLU_STABLE, 0.09187);
		initIsotopeProperty(42, 95, 94.9058421, 0.0, HLU_STABLE, 0.15873);
		initIsotopeProperty(42, 96, 95.9046795, 0.0, HLU_STABLE, 0.16673);
		initIsotopeProperty(42, 97, 96.9060215, 0.0, HLU_STABLE, 0.09582);
		initIsotopeProperty(42, 98, 97.9054082, 0.0, HLU_STABLE, 0.24292);
		initIsotopeProperty(42, 99, 98.9077119, 2.7489, HLU_DAY, 0.0);
		initIsotopeProperty(42, 100, 99.907477, 8.5e18, HLU_YEAR, 0.09744);
		initIsotopeProperty(42, 101, 100.910347, 14.61, HLU_MINUTE, 0.0);
		initIsotopeProperty(42, 102, 101.910297, 11.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(42, 103, 102.91321, 67.5, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 104, 103.91376, 60.0, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 105, 104.91697, 35.6, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 106, 105.918137, 8.73, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 107, 106.92169, 3.5, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 108, 107.92345, 1.09, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 109, 108.92781, 0.53, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 110, 109.92973, 0.27, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 111, 110.93441, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 112, 111.93684, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(42, 113, 112.94188, 100.0e-3, HLU_SECOND, 0.0);
		initDecayMode(42, 84, DECAY_MODE_BETA_PLUS);
		initDecayMode(42, 85, DECAY_MODE_BETA_PLUS);
		initDecayMode(42, 86, DECAY_MODE_BETA_PLUS);
		initDecayMode2(42, 87, DECAY_MODE_BETA_PLUS, 85.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(42, 88, DECAY_MODE_BETA_PLUS);
		initDecayMode(42, 89, DECAY_MODE_BETA_PLUS);
		initDecayMode(42, 90, DECAY_MODE_BETA_PLUS);
		initDecayMode(42, 91, DECAY_MODE_BETA_PLUS);
		initDecayMode(42, 93, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(42, 99, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 100, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(42, 101, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 102, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 103, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 104, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 105, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 106, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 107, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 108, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 109, DECAY_MODE_BETA_MINUS);
		initDecayMode2(42, 110, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(42, 111, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 112, DECAY_MODE_BETA_MINUS);
		initDecayMode(42, 113, DECAY_MODE_BETA_MINUS);

		initAtomProperty(43, "Tc", "technetium", "テクネチウム", 86, 115, 86, 115, 24.27);
		initIsotopeProperty(43, 86, 85.94288, 55.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 87, 86.93653, 2.18, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 88, 87.93268, 5.8, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 89, 88.92717, 12.8, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 90, 89.92356, 8.7, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 91, 90.91843, 3.14, HLU_MINUTE, 0.0);
		initIsotopeProperty(43, 92, 91.915260, 4.25, HLU_MINUTE, 0.0);
		initIsotopeProperty(43, 93, 92.910249, 2.75, HLU_HOUR, 0.0);
		initIsotopeProperty(43, 94, 93.909657, 293.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(43, 95, 94.907657, 20.0, HLU_HOUR, 0.0);
		initIsotopeProperty(43, 96, 95.907871, 4.28, HLU_DAY, 0.0);
		initIsotopeProperty(43, 97, 96.906365, 2.6e6, HLU_YEAR, 0.0);
		initIsotopeProperty(43, 98, 97.907216, 4.2e6, HLU_YEAR, 0.0);
		initIsotopeProperty(43, 99, 98.9062547, 2.111e5, HLU_YEAR, 0.0);
		initIsotopeProperty(43, 100, 99.9076578, 15.8, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 101, 100.907315, 14.22, HLU_MINUTE, 0.0);
		initIsotopeProperty(43, 102, 101.909215, 5.28, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 103, 102.909181, 54.2, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 104, 103.91145, 18.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(43, 105, 104.91166, 7.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(43, 106, 105.914358, 35.6, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 107, 106.91508, 21.2, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 108, 107.91846, 5.17, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 109, 108.91998, 860.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 110, 109.92382, 0.92, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 111, 110.92569, 290.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 112, 111.92915, 290.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 113, 112.93159, 170.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 114, 113.93588, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(43, 115, 114.93869, 100.0e-3, HLU_SECOND, 0.0);
		initDecayMode(43, 86, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 87, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 88, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 89, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 90, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 91, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 92, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 93, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 94, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 95, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 96, DECAY_MODE_BETA_PLUS);
		initDecayMode(43, 97, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(43, 98, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(43, 100, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(43, 101, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 102, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 103, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 104, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 105, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 106, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 107, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 108, DECAY_MODE_BETA_MINUS);
		initDecayMode2(43, 109, DECAY_MODE_BETA_MINUS, 99.92, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(43, 110, DECAY_MODE_BETA_MINUS, 99.96, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(43, 111, DECAY_MODE_BETA_MINUS, 99.15, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(43, 112, DECAY_MODE_BETA_MINUS, 97.4, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(43, 113, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 114, DECAY_MODE_BETA_MINUS);
		initDecayMode(43, 115, DECAY_MODE_BETA_MINUS);

		initAtomProperty(44, "Ru", "ruthenium", "ルテニウム", 87, 118, 87, 118, 24.06);
		initIsotopeProperty(44, 87, 86.94918, 50.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 88, 87.94026, 1.3, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 89, 88.93611, 1.38, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 90, 89.92989, 11.7, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 91, 90.92629, 7.9, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 92, 91.92012, 3.65, HLU_MINUTE, 0.0);
		initIsotopeProperty(44, 93, 92.91705, 59.7, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 94, 93.911360, 51.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(44, 95, 94.910413, 1.643, HLU_HOUR, 0.0);
		initIsotopeProperty(44, 96, 95.907598, 0.0, HLU_STABLE, 0.0554);
		initIsotopeProperty(44, 97, 96.907555, 2.791, HLU_DAY, 0.0);
		initIsotopeProperty(44, 98, 97.905287, 0.0, HLU_STABLE, 0.0187);
		initIsotopeProperty(44, 99, 98.9059393, 0.0, HLU_STABLE, 0.1276);
		initIsotopeProperty(44, 100, 99.9042195, 0.0, HLU_STABLE, 0.1260);
		initIsotopeProperty(44, 101, 100.9055821, 0.0, HLU_STABLE, 0.1706);
		initIsotopeProperty(44, 102, 101.9043493, 0.0, HLU_STABLE, 0.3155);
		initIsotopeProperty(44, 103, 102.9063238, 39.26, HLU_DAY, 0.0);
		initIsotopeProperty(44, 104, 103.905433, 0.0, HLU_STABLE, 0.1862);
		initIsotopeProperty(44, 105, 104.907753, 4.44, HLU_HOUR, 0.0);
		initIsotopeProperty(44, 106, 105.907329, 373.59, HLU_DAY, 0.0);
		initIsotopeProperty(44, 107, 106.90991, 3.75, HLU_MINUTE, 0.0);
		initIsotopeProperty(44, 108, 107.91017, 4.55, HLU_MINUTE, 0.0);
		initIsotopeProperty(44, 109, 108.91320, 34.5, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 110, 109.91414, 11.6, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 111, 110.91770, 2.12, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 112, 111.91897, 1.75, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 113, 112.92249, 0.80, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 114, 113.92428, 0.53, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 115, 114.92869, 740.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 116, 115.93081, 400.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 117, 116.93558, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(44, 118, 117.93782, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode(44, 87, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 88, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 89, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 90, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 91, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 92, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 93, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 94, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 95, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 97, DECAY_MODE_BETA_PLUS);
		initDecayMode(44, 103, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 105, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 106, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 107, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 108, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 109, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 110, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 111, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 112, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 113, DECAY_MODE_BETA_MINUS);
		initDecayMode2(44, 114, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(44, 115, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(44, 116, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 117, DECAY_MODE_BETA_MINUS);
		initDecayMode(44, 118, DECAY_MODE_BETA_MINUS);

		initAtomProperty(45, "Rh", "rhodium", "ロジウム", 89, 121, 89, 121, 24.98);
		initIsotopeProperty(45, 89, 88.94884, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 90, 89.94287, 15.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 91, 90.93655, 1.74, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 92, 91.93198, 4.3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 93, 92.92574, 11.9, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 94, 93.92170, 70.6, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 95, 94.91590, 5.02, HLU_MINUTE, 0.0);
		initIsotopeProperty(45, 96, 95.914461, 9.90, HLU_MINUTE, 0.0);
		initIsotopeProperty(45, 97, 96.91134, 30.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(45, 98, 97.910708, 8.72, HLU_MINUTE, 0.0);
		initIsotopeProperty(45, 99, 98.908132, 16.1, HLU_DAY, 0.0);
		initIsotopeProperty(45, 100, 99.908122, 20.8, HLU_HOUR, 0.0);
		initIsotopeProperty(45, 101, 100.906164, 3.3, HLU_YEAR, 0.0);
		initIsotopeProperty(45, 102, 101.906843, 207.0, HLU_DAY, 0.0);
		initIsotopeProperty(45, 103, 102.905504, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(45, 104, 103.906656, 42.3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 105, 104.905694, 35.36, HLU_HOUR, 0.0);
		initIsotopeProperty(45, 106, 105.907287, 29.80, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 107, 106.906748, 21.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(45, 108, 107.90873, 16.8, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 109, 108.908737, 80.0, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 110, 109.91114, 28.5, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 111, 110.91159, 11.0, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 112, 111.91439, 3.45, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 113, 112.91553, 2.80, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 114, 113.91881, 1.85, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 115, 114.92033, 0.99, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 116, 115.92406, 0.68, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 117, 116.92598, 0.44, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 118, 117.93007, 310.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 119, 118.93211, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 120, 119.93641, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(45, 121, 120.93872, 100.0e-3, HLU_SECOND, 0.0);
		initDecayMode(45, 89, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 90, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 91, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 92, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 93, DECAY_MODE_BETA_PLUS);
		initDecayMode2(45, 94, DECAY_MODE_BETA_PLUS, 98.2, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(45, 95, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 96, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 97, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 98, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 99, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 100, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 101, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(45, 102, DECAY_MODE_BETA_PLUS, 80.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(45, 104, DECAY_MODE_BETA_MINUS, 99.55, DECAY_MODE_BETA_PLUS);
		initDecayMode(45, 105, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 106, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 107, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 108, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 109, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 110, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 111, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 112, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 113, DECAY_MODE_BETA_MINUS);
		initDecayMode2(45, 114, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(45, 115, DECAY_MODE_BETA_MINUS);
		initDecayMode2(45, 116, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(45, 117, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 118, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 119, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 120, DECAY_MODE_BETA_MINUS);
		initDecayMode(45, 121, DECAY_MODE_BETA_MINUS);

		initAtomProperty(46, "Pd", "palladium", "パラジウム", 91, 123, 91, 123, 25.98);
		initIsotopeProperty(46, 91, 90.94911, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 92, 91.94042, 1.1, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 93, 92.93591, 1.07, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 94, 93.92877, 9.0, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 95, 94.92469, 10.0, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 96, 95.91816, 122.0, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 97, 96.91648, 3.10, HLU_MINUTE, 0.0);
		initIsotopeProperty(46, 98, 97.912721, 17.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(46, 99, 98.911768, 21.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(46, 100, 99.908506, 3.63, HLU_DAY, 0.0);
		initIsotopeProperty(46, 101, 100.908289, 8.47, HLU_HOUR, 0.0);
		initIsotopeProperty(46, 102, 101.905609, 0.0, HLU_STABLE, 0.0102);
		initIsotopeProperty(46, 103, 102.906087, 16.991, HLU_DAY, 0.0);
		initIsotopeProperty(46, 104, 103.904036, 0.0, HLU_STABLE, 0.1114);
		initIsotopeProperty(46, 105, 104.905085, 0.0, HLU_STABLE, 0.2233);
		initIsotopeProperty(46, 106, 105.903486, 0.0, HLU_STABLE, 0.2733);
		initIsotopeProperty(46, 107, 106.905133, 6.e6, HLU_YEAR, 0.0);
		initIsotopeProperty(46, 108, 107.903892, 0.0, HLU_STABLE, 0.2646);
		initIsotopeProperty(46, 109, 108.905950, 13.7012, HLU_HOUR, 0.0);
		initIsotopeProperty(46, 110, 109.905153, 0.0, HLU_STABLE, 0.1172);
		initIsotopeProperty(46, 111, 110.907671, 23.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(46, 112, 111.907314, 21.03, HLU_HOUR, 0.0);
		initIsotopeProperty(46, 113, 112.91015, 93.0, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 114, 113.910363, 2.42, HLU_MINUTE, 0.0);
		initIsotopeProperty(46, 115, 114.91368, 25.0, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 116, 115.91416, 11.8, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 117, 116.91784, 4.3, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 118, 117.91898, 1.9, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 119, 118.92311, 0.92, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 120, 119.92469, 0.5, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 121, 120.92887, 400.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 122, 121.93055, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(46, 123, 122.93493, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode(46, 91, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 92, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 93, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 94, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 95, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 96, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 97, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 98, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 99, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 100, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(46, 101, DECAY_MODE_BETA_PLUS);
		initDecayMode(46, 103, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(46, 107, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 109, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 111, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 112, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 113, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 114, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 115, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 116, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 117, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 118, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 119, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 120, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 121, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 122, DECAY_MODE_BETA_MINUS);
		initDecayMode(46, 123, DECAY_MODE_BETA_MINUS);

		initAtomProperty(47, "Ag", "silver", "銀", 94, 127, 94, 127, 25.35);
		initIsotopeProperty(47, 94, 93.94278, 37.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 95, 94.93548, 1.74, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 96, 95.93068, 4.45, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 97, 96.92397, 25.3, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 98, 97.92157, 47.5, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 99, 98.91760, 124, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 100, 99.91610, 2.01, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 101, 100.91280, 11.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 102, 101.91169, 12.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 103, 102.908973, 65.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 104, 103.908629, 69.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 105, 104.906529, 41.29, HLU_DAY, 0.0);
		initIsotopeProperty(47, 106, 105.906669, 23.96, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 107, 106.905097, 0.0, HLU_STABLE, 0.51839);
		initIsotopeProperty(47, 108, 107.905956, 2.37, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 109, 108.904752, 0.0, HLU_STABLE, 0.48161);
		initIsotopeProperty(47, 110, 109.906107, 24.6, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 111, 110.905291, 7.45, HLU_DAY, 0.0);
		initIsotopeProperty(47, 112, 111.907005, 3.130, HLU_HOUR, 0.0);
		initIsotopeProperty(47, 113, 112.906567, 5.37, HLU_HOUR, 0.0);
		initIsotopeProperty(47, 114, 113.908804, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 115, 114.90876, 20.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 116, 115.91136, 2.68, HLU_MINUTE, 0.0);
		initIsotopeProperty(47, 117, 116.91168, 73.6, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 118, 117.91458, 3.76, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 119, 118.91567, 6.0, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 120, 119.91879, 1.23, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 121, 120.91985, 0.79, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 122, 121.92353, 0.529, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 123, 122.92490, .300, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 124, 123.92864, 172.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 125, 124.93043, 166.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 126, 125.93450, 107.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(47, 127, 126.93677, 79.0e-3, HLU_SECOND, 0.0);
		initDecayMode(47, 94, DECAY_MODE_BETA_PLUS);
		initDecayMode2(47, 95, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(47, 96, DECAY_MODE_BETA_PLUS, 96.3, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(47, 97, DECAY_MODE_BETA_PLUS);
		initDecayMode2(47, 98, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(47, 99, DECAY_MODE_BETA_PLUS);
		initDecayMode(47, 100, DECAY_MODE_BETA_PLUS);
		initDecayMode(47, 101, DECAY_MODE_BETA_PLUS);
		initDecayMode(47, 102, DECAY_MODE_BETA_PLUS);
		initDecayMode(47, 103, DECAY_MODE_BETA_PLUS);
		initDecayMode(47, 104, DECAY_MODE_BETA_PLUS);
		initDecayMode(47, 105, DECAY_MODE_BETA_PLUS);
		initDecayMode2(47, 106, DECAY_MODE_BETA_PLUS, 99.5, DECAY_MODE_BETA_MINUS);
		initDecayMode2(47, 108, DECAY_MODE_BETA_MINUS, 97.15, DECAY_MODE_BETA_PLUS);
		initDecayMode2(47, 110, DECAY_MODE_BETA_MINUS, 99.7, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(47, 111, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 112, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 113, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 114, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 115, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 116, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 117, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 118, DECAY_MODE_BETA_MINUS);
		initDecayMode(47, 119, DECAY_MODE_BETA_MINUS);
		initDecayMode2(47, 120, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 121, DECAY_MODE_BETA_MINUS, 99.92, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 122, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 123, DECAY_MODE_BETA_MINUS, 99.45, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 124, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 125, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 126, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(47, 127, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(48, "Cd", "cadmium", "カドミウム", 96, 130, 96, 130, 26.02);
		initIsotopeProperty(48, 96, 95.93977, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 97, 96.93494, 2.8, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 98, 97.92740, 9.2, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 99, 98.92501, 16.0, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 100, 99.92029, 49.1, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 101, 100.91868, 1.36, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 102, 101.91446, 5.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 103, 102.913419, 7.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 104, 103.909849, 57.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 105, 104.909468, 55.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 106, 105.906459, 0.0, HLU_STABLE, 0.0125);
		initIsotopeProperty(48, 107, 106.906618, 6.50, HLU_HOUR, 0.0);
		initIsotopeProperty(48, 108, 107.904184, 0.0, HLU_STABLE, 0.0089);
		initIsotopeProperty(48, 109, 108.904982, 461.4, HLU_DAY, 0.0);
		initIsotopeProperty(48, 110, 109.9030021, 0.0, HLU_STABLE, 0.1249);
		initIsotopeProperty(48, 111, 110.9041781, 0.0, HLU_STABLE, 0.128);
		initIsotopeProperty(48, 112, 111.9027578, 0.0, HLU_STABLE, 0.2413);
		initIsotopeProperty(48, 113, 112.9044017, 8.04e15, HLU_YEAR, 0.1222);
		initIsotopeProperty(48, 114, 113.9033585, 0.0, HLU_STABLE, 0.2873);
		initIsotopeProperty(48, 115, 114.9054310, 53.46, HLU_HOUR, 0.0);
		initIsotopeProperty(48, 116, 115.904756, 2.8e19, HLU_YEAR, 0.0749);
		initIsotopeProperty(48, 117, 116.907219, 2.49, HLU_HOUR, 0.0);
		initIsotopeProperty(48, 118, 117.906915, 50.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 119, 118.90992, 2.69, HLU_MINUTE, 0.0);
		initIsotopeProperty(48, 120, 119.90985, 50.80, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 121, 120.91298, 13.5, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 122, 121.91333, 5.24, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 123, 122.91700, 2.10, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 124, 123.91765, 1.25, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 125, 124.92125, 0.65, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 126, 125.92235, 0.515, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 127, 126.92644, 0.37, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 128, 127.92776, 0.28, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 129, 128.93215, 242.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(48, 130, 129.9339, 162.0e-3, HLU_SECOND, 0.0);	
		initDecayMode(48, 96, DECAY_MODE_BETA_PLUS);
		initDecayMode2(48, 97, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(48, 98, DECAY_MODE_BETA_PLUS, 99.975, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(48, 99, DECAY_MODE_BETA_PLUS, 99.78, DECAY_MODE_BETA_PLUS_PROTON, 0.21, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(48, 100, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 101, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 102, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 103, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 104, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 105, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 107, DECAY_MODE_BETA_PLUS);
		initDecayMode(48, 109, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(48, 113, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 115, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 116, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(48, 117, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 118, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 119, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 120, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 121, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 122, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 123, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 124, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 125, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 126, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 127, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 128, DECAY_MODE_BETA_MINUS);
		initDecayMode(48, 129, DECAY_MODE_BETA_MINUS);
		initDecayMode2(48, 130, DECAY_MODE_BETA_MINUS, 96.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(49, "In", "indium", "インジウム", 98, 134, 98, 134, 26.74);
		initIsotopeProperty(49, 98, 97.94214, 45.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 99, 98.93422, 3.1, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 100, 99.93111, 5.9, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 101, 100.92634, 15.1, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 102, 101.92409, 23.3, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 103, 102.919914, 60, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 104, 103.91830, 1.80, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 105, 104.914674, 5.07, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 106, 105.913465, 6.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 107, 106.910295, 32.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 108, 107.909698, 58.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 109, 108.907151, 4.2, HLU_HOUR, 0.0);
		initIsotopeProperty(49, 110, 109.907165, 4.9, HLU_HOUR, 0.0);
		initIsotopeProperty(49, 111, 110.905103, 2.8047, HLU_DAY, 0.0);
		initIsotopeProperty(49, 112, 111.905532, 14.97, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 113, 112.904058, 0.0, HLU_STABLE, 0.0429);
		initIsotopeProperty(49, 114, 113.904914, 71.9, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 115, 114.903878, 4.41e14, HLU_YEAR, 0.9571);
		initIsotopeProperty(49, 116, 115.905260, 14.10, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 117, 116.904514, 43.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 118, 117.906354, 5.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 119, 118.905845, 2.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(49, 120, 119.90796, 3.08, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 121, 120.907846, 23.1, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 122, 121.91028, 1.5, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 123, 122.910438, 6.17, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 124, 123.91318, 3.11, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 125, 124.91360, 2.36, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 126, 125.91646, 1.53, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 127, 126.91735, 1.09, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 128, 127.92017, 0.84, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 129, 128.92170, 611.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 130, 129.92497, 0.29, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 131, 130.92685, 0.28, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 132, 131.93299, 206.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 133, 132.93781, 165.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(49, 134, 133.94415, 140.0e-3, HLU_SECOND, 0.0);	
		initDecayMode(49, 98, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 99, DECAY_MODE_BETA_PLUS);
		initDecayMode2(49, 100, DECAY_MODE_BETA_PLUS, 96.1, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(49, 101, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(49, 102, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(49, 103, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 104, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 105, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 106, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 107, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 108, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 109, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 110, DECAY_MODE_BETA_PLUS);
		initDecayMode(49, 111, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(49, 112, DECAY_MODE_BETA_PLUS, 56.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(49, 114, DECAY_MODE_BETA_PLUS, 99.5, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 115, DECAY_MODE_BETA_MINUS);
		initDecayMode2(49, 116, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(49, 117, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 118, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 119, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 120, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 121, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 122, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 123, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 124, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 125, DECAY_MODE_BETA_MINUS);
		initDecayMode(49, 126, DECAY_MODE_BETA_MINUS);
		initDecayMode2(49, 127, DECAY_MODE_BETA_MINUS, 99.97, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(49, 128, DECAY_MODE_BETA_MINUS, 99.96, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(49, 129, DECAY_MODE_BETA_MINUS, 99.75, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(49, 130, DECAY_MODE_BETA_MINUS, 98.35, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(49, 131, DECAY_MODE_BETA_MINUS, 97.8, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(49, 132, DECAY_MODE_BETA_MINUS, 94.8, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(49, 133, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 85.0, DECAY_MODE_BETA_MINUS);
		initDecayMode3(49, 124, DECAY_MODE_BETA_MINUS, 79.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON, 17.0, DECAY_MODE_BETA_MINUS_AND_2NEUTRON);

		initAtomProperty(50, "Sn", "Tin", "スズ", 100, 137, 100, 137, 27.112);
		initIsotopeProperty(50, 100, 99.93904, 1.1, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 101, 100.93606, 3.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 102, 101.93030, 4.5, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 103, 102.92810, 7.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 104, 103.92314, 20.8, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 105, 104.92135, 34.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 106, 105.91688, 115.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 107, 106.91564, 2.90, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 108, 107.911925, 10.30, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 109, 108.911283, 18.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 110, 109.907843, 4.11, HLU_HOUR, 0.0);
		initIsotopeProperty(50, 111, 110.907734, 35.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 112, 111.904818, 0.0, HLU_STABLE, 0.0097);
		initIsotopeProperty(50, 113, 112.905171, 115.09, HLU_DAY, 0.0);
		initIsotopeProperty(50, 114, 113.902779, 0.0, HLU_STABLE, 0.0066);
		initIsotopeProperty(50, 115, 114.903342, 0.0, HLU_STABLE, 0.0034);
		initIsotopeProperty(50, 116, 115.901741, 0.0, HLU_STABLE, 0.1454);
		initIsotopeProperty(50, 117, 116.902952, 0.0, HLU_STABLE, 0.0768);
		initIsotopeProperty(50, 118, 117.901603, 0.0, HLU_STABLE, 0.2422);
		initIsotopeProperty(50, 119, 118.903308, 0.0, HLU_STABLE, 0.0859);
		initIsotopeProperty(50, 120, 119.9021947, 0.0, HLU_STABLE, 0.3258);
		initIsotopeProperty(50, 121, 120.9042355, 27.03, HLU_HOUR, 0.0);
		initIsotopeProperty(50, 122, 121.9034390, 0.0, HLU_STABLE, 0.0463);
		initIsotopeProperty(50, 123, 122.9057208, 129.2, HLU_DAY, 0.0);
		initIsotopeProperty(50, 124, 123.9052739, 0.0, HLU_STABLE, 0.0579);
		initIsotopeProperty(50, 125, 124.9077841, 9.64, HLU_DAY, 0.0);
		initIsotopeProperty(50, 126, 125.907653, 2.30e5, HLU_YEAR, 0.0);
		initIsotopeProperty(50, 127, 126.910360, 2.10, HLU_HOUR, 0.0);
		initIsotopeProperty(50, 128, 127.910537, 59.07, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 129, 128.91348, 2.23, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 130, 129.913967, 3.72, HLU_MINUTE, 0.0);
		initIsotopeProperty(50, 131, 130.917000, 56.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 132, 131.917816, 39.7, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 133, 132.92383, 1.45, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 134, 133.92829, 1.050, HLU_SECOND, 0.0);	
		initIsotopeProperty(50, 135, 134.93473, 530.0e-3, HLU_SECOND, 0.0);		
		initIsotopeProperty(50, 136, 135.93934, 0.25, HLU_SECOND, 0.0);		
		initIsotopeProperty(50, 137, 136.94599, 190.0e-3, HLU_SECOND, 0.0);	
		initDecayMode2(50, 100, DECAY_MODE_BETA_PLUS, 83.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(50, 101, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(50, 102, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(50, 103, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(50, 104, DECAY_MODE_BETA_PLUS);
		initDecayMode2(50, 105, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(50, 106, DECAY_MODE_BETA_PLUS);
		initDecayMode(50, 107, DECAY_MODE_BETA_PLUS);
		initDecayMode(50, 108, DECAY_MODE_BETA_PLUS);
		initDecayMode(50, 109, DECAY_MODE_BETA_PLUS);
		initDecayMode(50, 110, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(50, 111, DECAY_MODE_BETA_PLUS);
		initDecayMode(50, 113, DECAY_MODE_BETA_PLUS);
		initDecayMode(50, 121, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 123, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 125, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 126, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 127, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 128, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 129, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 130, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 131, DECAY_MODE_BETA_MINUS);
		initDecayMode(50, 132, DECAY_MODE_BETA_MINUS);
		initDecayMode2(50, 133, DECAY_MODE_BETA_MINUS, 99.97, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(50, 134, DECAY_MODE_BETA_MINUS, 83.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(50, 135, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(50, 136, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(50, 137, DECAY_MODE_BETA_MINUS);

		initAtomProperty(51, "Sb", "antimony", "アンチモン", 103, 139, 103, 139, 25.23);
		initIsotopeProperty(51, 103, 102.93969, 100.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 104, 103.93647, 0.47, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 105, 104.93149, 1.12, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 106, 105.92879, 0.6, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 107, 106.92415, 4.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 108, 107.92216, 7.4, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 109, 108.918132, 17.3, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 110, 109.91675, 23.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 111, 110.91316, 75.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 112, 111.912398, 51.4, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 113, 112.909372, 6.67, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 114, 113.90927, 3.49, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 115, 114.906598, 32.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 116, 115.906794, 15.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 117, 116.904836, 2.800, HLU_HOUR, 0.0);
		initIsotopeProperty(51, 118, 117.905529, 3.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 119, 118.903942, 38.190, HLU_HOUR, 0.0);
		initIsotopeProperty(51, 120, 119.905072, 15.89, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 121, 120.9038157, 0.0, HLU_STABLE, 0.5721);
		initIsotopeProperty(51, 122, 121.9051737, 2.7238, HLU_DAY, 0.0);
		initIsotopeProperty(51, 123, 122.9042140, 0.0, HLU_STABLE, 0.4279);
		initIsotopeProperty(51, 124, 123.9059357, 60.20, HLU_DAY, 0.0);
		initIsotopeProperty(51, 125, 124.9052538, 2.75856, HLU_YEAR, 0.0);
		initIsotopeProperty(51, 126, 125.90725, 12.35, HLU_DAY, 0.0);
		initIsotopeProperty(51, 127, 126.906924, 3.85, HLU_DAY, 0.0);
		initIsotopeProperty(51, 128, 127.909169, 9.010, HLU_HOUR, 0.0);
		initIsotopeProperty(51, 129, 128.909148, 4.400, HLU_HOUR, 0.0);
		initIsotopeProperty(51, 130, 129.911656, 39.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 131, 130.911982, 23.03, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 132, 131.914467, 2.79, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 133, 132.915252, 2.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(51, 134, 133.92038, 0.78, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 135, 134.92517, 1.68, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 136, 135.93035, 0.923, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 137, 136.93531, 450.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 138, 137.94079, 500.0e-3, HLU_SECOND, 0.0);	
		initIsotopeProperty(51, 139, 138.94598, 300.0e-3, HLU_SECOND, 0.0);	
		initDecayMode(51, 103, DECAY_MODE_BETA_PLUS);
		initDecayMode4(51, 104, DECAY_MODE_BETA_PLUS, 85.5, DECAY_MODE_PROTON_EMISSION, 7.0, DECAY_MODE_BETA_PLUS_PROTON, 7.0, DECAY_MODE_ALPHA);
		initDecayMode3(51, 105, DECAY_MODE_BETA_PLUS, 98.5, DECAY_MODE_PROTON_EMISSION, 1.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(51, 106, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 107, DECAY_MODE_BETA_PLUS);
		initDecayMode2(51, 108, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(51, 109, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 110, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 111, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 112, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 113, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 114, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 115, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 116, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 117, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 118, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 119, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(51, 120, DECAY_MODE_BETA_PLUS);
		initDecayMode2(51, 122, DECAY_MODE_BETA_MINUS, 97.59, DECAY_MODE_BETA_PLUS);
		initDecayMode(51, 124, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 125, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 126, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 127, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 128, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 129, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 130, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 131, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 132, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 133, DECAY_MODE_BETA_MINUS);
		initDecayMode(51, 134, DECAY_MODE_BETA_MINUS);
		initDecayMode2(51, 135, DECAY_MODE_BETA_MINUS, 82.4, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(51, 136, DECAY_MODE_BETA_MINUS, 83.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(51, 137, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(51, 138, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(51, 139, DECAY_MODE_BETA_MINUS);

		initAtomProperty(52, "Te", "tellurium", "テルル", 106, 142, 106, 142, 25.73);
		initIsotopeProperty(52, 106, 105.93750, 70.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 107, 106.93501, 3.1e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 108, 107.92944, 2.1, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 109, 108.92742, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 110, 109.92241, 18.6, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 111, 110.92111, 19.3, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 112, 111.91701, 2.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 113, 112.91589, 1.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 114, 113.91209, 15.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 115, 114.91190, 5.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 116, 115.90846, 2.49, HLU_HOUR, 0.0);
		initIsotopeProperty(52, 117, 116.908645, 62.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 118, 117.905828, 6.00, HLU_DAY, 0.0);
		initIsotopeProperty(52, 119, 118.906404, 16.05, HLU_HOUR, 0.0);
		initIsotopeProperty(52, 120, 119.90402, 0.0, HLU_STABLE, 9.00E-04);
		initIsotopeProperty(52, 121, 120.904936, 19.16, HLU_DAY, 0.0);
		initIsotopeProperty(52, 122, 121.9030439, 0.0, HLU_STABLE, 0.0255);
		initIsotopeProperty(52, 123, 122.9042700, 0.0, HLU_STABLE, 0.0089);
		initIsotopeProperty(52, 124, 123.9028179, 0.0, HLU_STABLE, 0.0474);
		initIsotopeProperty(52, 125, 124.9044307, 0.0, HLU_STABLE, 0.0707);
		initIsotopeProperty(52, 126, 125.9033117, 0.0, HLU_STABLE, 0.1884);
		initIsotopeProperty(52, 127, 126.9052263, 9.35, HLU_HOUR, 0.0);
		initIsotopeProperty(52, 128, 127.9044631, 2.2e24, HLU_YEAR, 0.3174);
		initIsotopeProperty(52, 129, 128.9065982, 69.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 130, 129.9062244, 790e18, HLU_YEAR, 0.3408);
		initIsotopeProperty(52, 131, 130.9085239, 25.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 132, 131.908553, 3.204, HLU_DAY, 0.0);
		initIsotopeProperty(52, 133, 132.910955, 12.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 134, 133.911369, 41.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(52, 135, 134.91645, 19.0, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 136, 135.92010, 17.63, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 137, 136.92532, 2.49, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 138, 137.92922, 1.4, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 139, 138.93473, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 140, 139.93885, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 141, 140.94465, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(52, 142, 141.94908, 50.0e-3, HLU_SECOND, 0.0);
		initDecayMode(52, 106, DECAY_MODE_ALPHA);
		initDecayMode2(52, 107, DECAY_MODE_ALPHA, 70.0, DECAY_MODE_BETA_PLUS);
		initDecayMode4(52, 108, DECAY_MODE_BETA_PLUS, 51.0 - 0.065, DECAY_MODE_ALPHA, 49.0 - 2.4, DECAY_MODE_BETA_PLUS_PROTON, 2.4, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode4(52, 109, DECAY_MODE_ELECTRON_CAPTURE, 82.695, DECAY_MODE_BETA_PLUS_PROTON, 9.4, DECAY_MODE_ALPHA, 7.9, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode2(52, 110, DECAY_MODE_ELECTRON_CAPTURE, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(52, 111, DECAY_MODE_ELECTRON_CAPTURE, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(52, 112, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 113, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 114, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 115, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 116, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 117, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 118, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 119, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 121, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(52, 127, DECAY_MODE_BETA_MINUS);
		initDecayMode(52, 128, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(52, 129, DECAY_MODE_BETA_MINUS);
		initDecayMode(52, 130, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(52, 131, DECAY_MODE_BETA_MINUS);
		initDecayMode(52, 132, DECAY_MODE_BETA_MINUS);
		initDecayMode(52, 133, DECAY_MODE_BETA_MINUS);
		initDecayMode(52, 134, DECAY_MODE_BETA_MINUS);
		initDecayMode(52, 135, DECAY_MODE_BETA_MINUS);
		initDecayMode2(52, 136, DECAY_MODE_BETA_MINUS, 98.7, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(52, 137, DECAY_MODE_BETA_MINUS, 97.01, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(52, 138, DECAY_MODE_BETA_MINUS, 93.7, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(52, 139, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(52, 140, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(52, 141, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(52, 142, DECAY_MODE_BETA_MINUS);

		initAtomProperty(53, "I", " iodine", "ヨウ素", 108, 144, 108, 144, 54.44/2.0);
		initIsotopeProperty(53, 108, 107.94348, 36.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 109, 108.93815, 103.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 110, 109.93524, 650.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 111, 110.93028, 2.5, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 112, 111.92797, 3.42, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 113, 112.92364, 6.6, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 114, 113.92185, 2.1, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 115, 114.91805, 1.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 116, 115.91681, 2.91, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 117, 116.91365, 2.22, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 118, 117.913074, 13.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 119, 118.91007, 19.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 120, 119.910048, 81.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 121, 120.907367, 2.12, HLU_HOUR, 0.0);
		initIsotopeProperty(53, 122, 121.907589, 3.63, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 123, 122.905589, 13.2235, HLU_HOUR, 0.0);
		initIsotopeProperty(53, 124, 123.9062099, 4.1760, HLU_DAY, 0.0);
		initIsotopeProperty(53, 125, 124.9046302, 59.400, HLU_DAY, 0.0);
		initIsotopeProperty(53, 126, 125.905624, 12.93, HLU_DAY, 0.0);
		initIsotopeProperty(53, 127, 126.904473, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(53, 128, 127.905809, 24.99, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 129, 128.904988, 1.57e7, HLU_YEAR, 0.0);
		initIsotopeProperty(53, 130, 129.906674, 12.36, HLU_HOUR, 0.0);
		initIsotopeProperty(53, 131, 130.9061246, 8.02070, HLU_DAY, 0.0);
		initIsotopeProperty(53, 132, 131.907997, 2.295, HLU_HOUR, 0.0);
		initIsotopeProperty(53, 133, 132.907797, 20.8, HLU_HOUR, 0.0);
		initIsotopeProperty(53, 134, 133.909744, 52.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(53, 135, 134.910048, 6.57, HLU_HOUR, 0.0);
		initIsotopeProperty(53, 136, 135.91465, 83.4, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 137, 136.917871, 24.13, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 138, 137.92235, 6.23, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 139, 138.92610, 2.282, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 140, 139.93100, 860.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 141, 140.93503, 430.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 142, 141.94018, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 143, 142.94456, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(53, 144, 143.94999, 50.0e-3, HLU_SECOND, 0.0);
		initDecayMode3(53, 108, DECAY_MODE_ALPHA, 90.0, DECAY_MODE_BETA_PLUS, 9.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(53, 109, DECAY_MODE_PROTON_EMISSION, 99.5, DECAY_MODE_ALPHA);
		initDecayMode4(53, 110, DECAY_MODE_BETA_PLUS, 83.0 - 1.1, DECAY_MODE_ALPHA, 17.0 - 1.09, DECAY_MODE_BETA_PLUS_PROTON, 1.1, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode2(53, 111, DECAY_MODE_BETA_PLUS, 99.91, DECAY_MODE_ALPHA);
		initDecayMode4(53, 112, DECAY_MODE_BETA_PLUS, 99.01, DECAY_MODE_BETA_PLUS_PROTON, 0.88, DECAY_MODE_BETA_PLUS_ALPHA, 0.104, DECAY_MODE_ALPHA);
		initDecayMode3(53, 113, DECAY_MODE_BETA_PLUS, 100.0 - 3.3e-7 * 1.3, DECAY_MODE_ALPHA, 3.3e-7, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode2(53, 114, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(53, 115, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 116, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 117, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 118, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 119, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 120, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 121, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 122, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 123, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(53, 124, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 125, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(53, 126, DECAY_MODE_BETA_PLUS, 56.3, DECAY_MODE_BETA_MINUS);
		initDecayMode2(53, 128, DECAY_MODE_BETA_MINUS, 93.1, DECAY_MODE_BETA_PLUS);
		initDecayMode(53, 129, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 130, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 131, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 132, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 133, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 134, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 135, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 136, DECAY_MODE_BETA_MINUS);
		initDecayMode2(53, 137, DECAY_MODE_BETA_MINUS, 92.86, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(53, 138, DECAY_MODE_BETA_MINUS, 94.54, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(53, 139, DECAY_MODE_BETA_MINUS, 90.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(53, 140, DECAY_MODE_BETA_MINUS, 90.7, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(53, 141, DECAY_MODE_BETA_MINUS, 78.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(53, 142, DECAY_MODE_BETA_MINUS, 75.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(53, 143, DECAY_MODE_BETA_MINUS);
		initDecayMode(53, 144, DECAY_MODE_BETA_MINUS);

		initAtomProperty(54, "Xe", "xenon", "キセノン", 110, 147, 110, 147, 20.786);
		initIsotopeProperty(54, 110, 109.94428, 310.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 111, 110.94160, 740.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 112, 111.93562, 2.7, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 113, 112.93334, 2.74, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 114, 113.927980, 10.0, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 115, 114.926294, 18.0, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 116, 115.921581, 59.0, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 117, 116.920359, 61.0, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 118, 117.916179, 3.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(54, 119, 118.915411, 5.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(54, 120, 119.911784, 40, HLU_MINUTE, 0.0);
		initIsotopeProperty(54, 121, 120.911462, 40.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(54, 122, 121.908368, 20.1, HLU_HOUR, 0.0);
		initIsotopeProperty(54, 123, 122.908482, 2.08, HLU_HOUR, 0.0);
		initIsotopeProperty(54, 124, 123.905893, 0.0, HLU_STABLE, 0.000952);
		initIsotopeProperty(54, 125, 124.9063955, 16.9, HLU_HOUR, 0.0);
		initIsotopeProperty(54, 126, 125.904274, 0.0, HLU_STABLE, 0.00089);
		initIsotopeProperty(54, 127, 126.905184, 36.345, HLU_DAY, 0.0);
		initIsotopeProperty(54, 128, 127.9035313, 0.0, HLU_STABLE, 0.019102);
		initIsotopeProperty(54, 129, 128.9047794, 0.0, HLU_STABLE, 0.264006);
		initIsotopeProperty(54, 130, 129.9035080, 0.0, HLU_STABLE, 0.04071);
		initIsotopeProperty(54, 131, 130.9050824, 0.0, HLU_STABLE, 0.212324);
		initIsotopeProperty(54, 132, 131.9041535, 0.0, HLU_STABLE, 0.269086);
		initIsotopeProperty(54, 133, 132.9059107, 5.2475, HLU_DAY, 0.0);
		initIsotopeProperty(54, 134, 133.9053945, 0.0, HLU_STABLE, 0.104357);
		initIsotopeProperty(54, 135, 134.907227, 9.14, HLU_HOUR, 0.0);
		initIsotopeProperty(54, 136, 135.907219, 2.165e21, HLU_YEAR, 0.088573);
		initIsotopeProperty(54, 137, 136.911562, 3.818, HLU_MINUTE, 0.0);
		initIsotopeProperty(54, 138, 137.91395, 14.08, HLU_MINUTE, 0.0);
		initIsotopeProperty(54, 139, 138.918793, 39.68, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 140, 139.92164, 13.60, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 141, 140.92665, 1.73, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 142, 141.92971, 1.22, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 143, 142.93511, 0.511, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 144, 143.93851, 0.388, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 145, 144.94407, 188.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 146, 145.94775, 146.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(54, 147, 146.95356, 130.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(54, 110, DECAY_MODE_BETA_PLUS, 80.0, DECAY_MODE_ALPHA);
		initDecayMode2(54, 111, DECAY_MODE_BETA_PLUS, 90.0, DECAY_MODE_ALPHA);
		initDecayMode2(54, 112, DECAY_MODE_BETA_PLUS, 99.1, DECAY_MODE_ALPHA);
		initDecayMode4(54, 113, DECAY_MODE_BETA_PLUS, 92.98, DECAY_MODE_BETA_PLUS_PROTON, 7.0, DECAY_MODE_ALPHA, 0.011, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(54, 114, DECAY_MODE_BETA_PLUS);
		initDecayMode3(54, 115, DECAY_MODE_BETA_PLUS, 99.65, DECAY_MODE_BETA_PLUS_PROTON, 0.34, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(54, 116, DECAY_MODE_BETA_PLUS);
		initDecayMode2(54, 117, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(54, 118, DECAY_MODE_BETA_PLUS);
		initDecayMode(54, 119, DECAY_MODE_BETA_PLUS);
		initDecayMode(54, 120, DECAY_MODE_BETA_PLUS);
		initDecayMode(54, 121, DECAY_MODE_BETA_PLUS);
		initDecayMode(54, 122, DECAY_MODE_BETA_PLUS);
		initDecayMode(54, 123, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(54, 125, DECAY_MODE_BETA_PLUS);
		initDecayMode(54, 127, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(54, 133, DECAY_MODE_BETA_MINUS);
		initDecayMode(54, 135, DECAY_MODE_BETA_MINUS);
		initDecayMode(54, 136, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(54, 137, DECAY_MODE_BETA_MINUS);
		initDecayMode(54, 138, DECAY_MODE_BETA_MINUS);
		initDecayMode(54, 139, DECAY_MODE_BETA_MINUS);
		initDecayMode(54, 140, DECAY_MODE_BETA_MINUS);
		initDecayMode2(54, 141, DECAY_MODE_BETA_MINUS, 99.45, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(54, 142, DECAY_MODE_BETA_MINUS, 99.59, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(54, 143, DECAY_MODE_BETA_MINUS);
		initDecayMode2(54, 144, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(54, 145, DECAY_MODE_BETA_MINUS);
		initDecayMode(54, 146, DECAY_MODE_BETA_MINUS);
		initDecayMode2(54, 147, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(55, "Cs", "caesium", "セシウム", 112, 151, 112, 151, 32.21);
		initIsotopeProperty(55, 112, 111.95030, 500.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 113, 112.94449, 16.7e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 114, 113.94145, 0.57, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 115, 114.93591, 1.4, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 116, 115.93337, 0.70, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 117, 116.92867, 8.4, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 118, 117.926559, 14.0, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 119, 118.922377, 43.0, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 120, 119.920677, 61.2, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 121, 120.917229, 155, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 122, 121.91611, 21.18, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 123, 122.912996, 5.88, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 124, 123.912258, 30.9, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 125, 124.909728, 46.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 126, 125.909452, 1.64, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 127, 126.907418, 6.25, HLU_HOUR, 0.0);
		initIsotopeProperty(55, 128, 127.907749, 3.640, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 129, 128.906064, 32.06, HLU_HOUR, 0.0);
		initIsotopeProperty(55, 130, 129.906709, 29.21, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 131, 130.905464, 9.689, HLU_DAY, 0.0);
		initIsotopeProperty(55, 132, 131.9064343, 6.480, HLU_DAY, 0.0);
		initIsotopeProperty(55, 133, 132.905451933, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(55, 134, 133.906718475, 2.0652, HLU_YEAR, 0.0);
		initIsotopeProperty(55, 135, 134.9059770, 2.3e6, HLU_YEAR, 0.0);
		initIsotopeProperty(55, 136, 135.9073116, 13.16, HLU_DAY, 0.0);
		initIsotopeProperty(55, 137, 136.9070895, 30.1671, HLU_YEAR, 0.0);
		initIsotopeProperty(55, 138, 137.911017, 33.41, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 139, 138.913364, 9.27, HLU_MINUTE, 0.0);
		initIsotopeProperty(55, 140, 139.917282, 63.7, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 141, 140.920046, 24.84, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 142, 141.924299, 1.689, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 143, 142.927352, 1.791, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 144, 143.932077, 994.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 145, 144.935526, 582.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 146, 145.94029, 0.321, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 147, 146.94416, 0.235, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 148, 147.94922, 146.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 149, 148.95293, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 150, 149.95817, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(55, 151, 150.96219, 60.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(55, 112, DECAY_MODE_PROTON_EMISSION, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(55, 113, DECAY_MODE_PROTON_EMISSION, 99.97, DECAY_MODE_BETA_PLUS);
		initDecayMode4(55, 114, DECAY_MODE_BETA_PLUS, 91.09, DECAY_MODE_BETA_PLUS_PROTON, 8.69, DECAY_MODE_BETA_PLUS_ALPHA, 0.19, DECAY_MODE_ALPHA);
		initDecayMode2(55, 115, DECAY_MODE_BETA_PLUS, 99.93, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(55, 116, DECAY_MODE_BETA_PLUS, 99.67, DECAY_MODE_BETA_PLUS_PROTON, 0.279, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(55, 117, DECAY_MODE_BETA_PLUS);
		initDecayMode3(55, 118, DECAY_MODE_BETA_PLUS, 99.95, DECAY_MODE_BETA_PLUS_PROTON, 0.042, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode2(55, 119, DECAY_MODE_BETA_PLUS, 100.0 - 2e-6, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode3(55, 120, DECAY_MODE_BETA_PLUS, 100.0 - 2e-5 - 7e-6, DECAY_MODE_BETA_PLUS_ALPHA, 2e-5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(55, 121, DECAY_MODE_BETA_PLUS);
		initDecayMode2(55, 122, DECAY_MODE_BETA_PLUS, 100.0 - 2e-7, DECAY_MODE_BETA_PLUS_ALPHA);
		initDecayMode(55, 123, DECAY_MODE_BETA_PLUS);
		initDecayMode(55, 124, DECAY_MODE_BETA_PLUS);
		initDecayMode(55, 125, DECAY_MODE_BETA_PLUS);
		initDecayMode(55, 126, DECAY_MODE_BETA_PLUS);
		initDecayMode(55, 127, DECAY_MODE_BETA_PLUS);
		initDecayMode(55, 128, DECAY_MODE_BETA_PLUS);
		initDecayMode(55, 129, DECAY_MODE_BETA_PLUS);
		initDecayMode2(55, 130, DECAY_MODE_BETA_PLUS, 98.4, DECAY_MODE_BETA_MINUS);
		initDecayMode(55, 131, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(55, 132, DECAY_MODE_BETA_PLUS, 98.13, DECAY_MODE_BETA_MINUS);
		initDecayMode2(55, 134, DECAY_MODE_BETA_MINUS, 100.0 - 3e-4, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(55, 135, DECAY_MODE_BETA_MINUS);
		initDecayMode(55, 136, DECAY_MODE_BETA_MINUS);
		initDecayMode(55, 137, DECAY_MODE_BETA_MINUS);
		initDecayMode(55, 138, DECAY_MODE_BETA_MINUS);
		initDecayMode(55, 139, DECAY_MODE_BETA_MINUS);
		initDecayMode(55, 140, DECAY_MODE_BETA_MINUS);
		initDecayMode2(55, 141, DECAY_MODE_BETA_MINUS, 99.96, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 142, DECAY_MODE_BETA_MINUS, 99.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 143, DECAY_MODE_BETA_MINUS, 98.38, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 144, DECAY_MODE_BETA_MINUS, 96.8, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 145, DECAY_MODE_BETA_MINUS, 85.7, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 146, DECAY_MODE_BETA_MINUS, 85.8, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 147, DECAY_MODE_BETA_MINUS, 71.5, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 148, DECAY_MODE_BETA_MINUS, 74.9, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 149, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 150, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(55, 151, DECAY_MODE_BETA_MINUS, 50.0, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

		initAtomProperty(56, "Ba", "barium", "バリウム", 114, 153, 114, 153, 28.07);
		initIsotopeProperty(56, 114, 113.95068, 530.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 115, 114.94737, 0.45, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 116, 115.94138, 1.3, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 117, 116.93850, 1.75, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 118, 117.93304, 5.2, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 119, 118.93066, 5.4, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 120, 119.92604, 24, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 121, 120.92405, 29.7, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 122, 121.91990, 1.95, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 123, 122.918781, 2.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 124, 123.915094, 11.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 125, 124.914473, 3.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 126, 125.911250, 100.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 127, 126.911094, 12.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 128, 127.908318, 2.43, HLU_DAY, 0.0);
		initIsotopeProperty(56, 129, 128.908679, 2.23, HLU_HOUR, 0.0);
		initIsotopeProperty(56, 130, 129.9063208, 1.6e21, HLU_YEAR, 0.00106);
		initIsotopeProperty(56, 131, 130.906941, 11.50, HLU_DAY, 0.0);
		initIsotopeProperty(56, 132, 131.9050613, 0.0, HLU_STABLE, 0.00101);
		initIsotopeProperty(56, 133, 132.9060075, 10.51, HLU_YEAR, 0.0);
		initIsotopeProperty(56, 134, 133.9045084, 0.0, HLU_STABLE, 0.02417);
		initIsotopeProperty(56, 135, 134.9056886, 0.0, HLU_STABLE, 0.06592);
		initIsotopeProperty(56, 136, 135.9045759, 0.0, HLU_STABLE, 0.07854);
		initIsotopeProperty(56, 137, 136.9058274, 0.0, HLU_STABLE, 0.11232);
		initIsotopeProperty(56, 138, 137.9052472, 0.0, HLU_STABLE, 0.71698);
		initIsotopeProperty(56, 139, 138.9088413, 83.06, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 140, 139.910605, 12.752, HLU_DAY, 0.0);
		initIsotopeProperty(56, 141, 140.914411, 18.27, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 142, 141.916453, 10.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(56, 143, 142.920627, 14.5, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 144, 143.922953, 11.5, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 145, 144.92763, 4.31, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 146, 145.93022, 2.22, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 147, 146.93495, 0.893, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 148, 147.93772, 0.612, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 149, 148.94258, 344.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 150, 149.94568, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 151, 150.95081, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 152, 151.95427, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(56, 153, 152.95961, 80.0e-3, HLU_SECOND, 0.0);
		initDecayMode3(56, 114, DECAY_MODE_BETA_PLUS_PROTON, 99.59, DECAY_MODE_ALPHA, 0.37, DECAY_MODE_BETA_PLUS);
		initDecayMode2(56, 115, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(56, 116, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(56, 117, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_ALPHA, 30.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(56, 118, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(56, 119, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(56, 120, DECAY_MODE_BETA_PLUS);
		initDecayMode2(56, 121, DECAY_MODE_BETA_PLUS, 99.98, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(56, 122, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 123, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 124, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 125, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 126, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 127, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 128, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 129, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 130, DECAY_MODE_DOUBLE_ELECTRON_CAPTURE);
		initDecayMode(56, 131, DECAY_MODE_BETA_PLUS);
		initDecayMode(56, 133, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(56, 139, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 140, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 141, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 142, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 143, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 144, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 145, DECAY_MODE_BETA_MINUS);
		initDecayMode2(56, 146, DECAY_MODE_BETA_MINUS, 99.98, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(56, 147, DECAY_MODE_BETA_MINUS, 99.94, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(56, 148, DECAY_MODE_BETA_MINUS, 99.6, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(56, 149, DECAY_MODE_BETA_MINUS, 99.57, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(56, 150, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(56, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 152, DECAY_MODE_BETA_MINUS);
		initDecayMode(56, 153, DECAY_MODE_BETA_MINUS);

		initAtomProperty(57, "La", "lanthanum", "ランタン", 117, 155, 117, 155, 24.11);
		initIsotopeProperty(57, 117, 116.95007, 23.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 118, 117.94673, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 119, 118.94099, 1.0, HLU_SECOND, 0.0);	
		initIsotopeProperty(57, 120, 119.93807, 2.8, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 121, 120.93301, 5.3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 122, 121.93071, 8.6, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 123, 122.92624, 17.0, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 124, 123.92457, 29.21, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 125, 124.920816, 64.8, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 126, 125.91951, 54.0, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 127, 126.916375, 5.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 128, 127.91559, 5.18, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 129, 128.912693, 11.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 130, 129.912369, 8.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 131, 130.91007, 59.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 132, 131.91010, 4.8, HLU_HOUR, 0.0);
		initIsotopeProperty(57, 133, 132.90822, 3.912, HLU_HOUR, 0.0);
		initIsotopeProperty(57, 134, 133.908514, 6.45, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 135, 134.906977, 19.5, HLU_HOUR, 0.0);
		initIsotopeProperty(57, 136, 135.90764, 9.87, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 137, 136.906494, 6.0e4, HLU_YEAR, 0.0);
		initIsotopeProperty(57, 138, 137.907112, 1.02e11, HLU_YEAR, 0.0009);
		initIsotopeProperty(57, 139, 138.9063533, 0.0, HLU_STABLE, 0.9991);
		initIsotopeProperty(57, 140, 139.9094776, 1.6781, HLU_DAY, 0.0);
		initIsotopeProperty(57, 141, 140.910962, 3.92, HLU_HOUR, 0.0);
		initIsotopeProperty(57, 142, 141.914079, 91.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 143, 142.916063, 14.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(57, 144, 143.91960, 40.8, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 145, 144.92165, 24.8, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 146, 145.92579, 6.27, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 147, 146.92824, 4.015, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 148, 147.93223, 1.26, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 149, 148.93473, 1.05, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 150, 149.93877, 510.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 151, 150.94172, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 152, 151.94625, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 153, 152.94962, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 154, 153.95450, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(57, 155, 154.95835, 60.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(57, 117, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(57, 118, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 119, DECAY_MODE_BETA_PLUS);
		initDecayMode2(57, 120, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(57, 121, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(57, 122, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(57, 123, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 124, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 125, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 126, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 127, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 128, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 129, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 130, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 131, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 132, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 133, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 134, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 135, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 136, DECAY_MODE_BETA_PLUS);
		initDecayMode(57, 137, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(57, 138, DECAY_MODE_BETA_PLUS, 66.4, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 140, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 141, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 142, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 143, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 144, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 145, DECAY_MODE_BETA_MINUS);
		initDecayMode2(57, 146, DECAY_MODE_BETA_MINUS, 99.99, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(57, 147, DECAY_MODE_BETA_MINUS, 99.96, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(57, 148, DECAY_MODE_BETA_MINUS, 99.85, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(57, 149, DECAY_MODE_BETA_MINUS, 98.6, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode2(57, 150, DECAY_MODE_BETA_MINUS, 97.3, DECAY_MODE_BETA_MINUS_AND_NEUTRON);
		initDecayMode(57, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 152, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 153, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 154, DECAY_MODE_BETA_MINUS);
		initDecayMode(57, 155, DECAY_MODE_BETA_MINUS);

		initAtomProperty(58, "Ce", "cerium", "セリウム", 119, 157, 119, 157, 26.94);
		initIsotopeProperty(58, 119, 118.95276, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 120, 119.94664, 250.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 121, 120.94342, 1.1, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 122, 121.93791, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 123, 122.93540, 3.8, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 124, 123.93041, 9.1, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 125, 124.92844, 9.3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 126, 125.92397, 51.0, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 127, 126.92273, 29.0, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 128, 127.91891, 3.93, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 129, 128.91810, 3.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 130, 129.91474, 22.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 131, 130.91442, 10.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 132, 131.911460, 3.51, HLU_HOUR, 0.0);
		initIsotopeProperty(58, 133, 132.911515, 97.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 134, 133.908925, 3.16, HLU_DAY, 0.0);
		initIsotopeProperty(58, 135, 134.909151, 17.7, HLU_DAY, 0.0);
		initIsotopeProperty(58, 136, 135.907172, 0.0, HLU_STABLE, 0.00185);
		initIsotopeProperty(58, 137, 136.907806, 9.0, HLU_DAY, 0.0);
		initIsotopeProperty(58, 138, 137.905991, 0.0, HLU_STABLE, 0.00251);
		initIsotopeProperty(58, 139, 138.906653, 137.641, HLU_DAY, 0.0);
		initIsotopeProperty(58, 140, 139.9054387, 0.0, HLU_STABLE, 0.8845);
		initIsotopeProperty(58, 141, 140.9082763, 32.508, HLU_DAY, 0.0);
		initIsotopeProperty(58, 142, 141.909244, 0.0, HLU_STABLE, 0.11114);
		initIsotopeProperty(58, 143, 142.912386, 33.039, HLU_DAY, 0.0);
		initIsotopeProperty(58, 144, 143.913647, 284.91, HLU_DAY, 0.0);
		initIsotopeProperty(58, 145, 144.91723, 3.01, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 146, 145.91876, 13.52, HLU_MINUTE, 0.0);
		initIsotopeProperty(58, 147, 146.92267, 56.4, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 148, 147.92443, 56.0, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 149, 148.9284, 5.3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 150, 149.93041, 4.0, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 151, 150.93398, 1.02, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 152, 151.93654, 1.4, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 153, 152.94058, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 154, 153.94342, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 155, 154.94804, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 156, 155.95126, 150.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(58, 157, 156.95634, 50.0e-3, HLU_SECOND, 0.0);
		initDecayMode(58, 120, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 121, DECAY_MODE_BETA_PLUS);
		initDecayMode2(58, 122, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(58, 123, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(58, 124, DECAY_MODE_BETA_PLUS);
		initDecayMode2(58, 125, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(58, 126, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 127, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 128, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 129, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 130, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 131, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 132, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 133, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 134, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(58, 135, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 137, DECAY_MODE_BETA_PLUS);
		initDecayMode(58, 139, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(58, 141, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 143, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 144, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 145, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 146, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 147, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 148, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 149, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 150, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 152, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 153, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 154, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 155, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 156, DECAY_MODE_BETA_MINUS);
		initDecayMode(58, 157, DECAY_MODE_BETA_MINUS);

		initAtomProperty(59, "Pr", "praseodymium", "プラセオジム", 121, 159, 121, 159, 27.2);
		initIsotopeProperty(59, 121, 120.95536, 600.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 122, 121.95181, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 123, 122.94596, 800.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 124, 123.94296, 1.2, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 125, 124.93783, 3.3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 126, 125.93531, 3.12, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 127, 126.93083, 4.2, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 128, 127.92879, 2.84, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 129, 128.92510, 32.0, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 130, 129.92359, 40.0, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 131, 130.92026, 1.50, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 132, 131.91926, 1.49, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 133, 132.916331, 6.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 134, 133.91571, 11.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 135, 134.913112, 24.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 136, 135.912692, 13.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 137, 136.910705, 1.28, HLU_HOUR, 0.0);
		initIsotopeProperty(59, 138, 137.910755, 1.45, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 139, 138.908938, 4.41, HLU_HOUR, 0.0);
		initIsotopeProperty(59, 140, 139.909076, 3.39, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 141, 140.9076528, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(59, 142, 141.9100448, 19.12, HLU_HOUR, 0.0);
		initIsotopeProperty(59, 143, 142.9108169, 13.57, HLU_DAY, 0.0);
		initIsotopeProperty(59, 144, 143.913305, 17.28, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 145, 144.914512, 5.984, HLU_HOUR, 0.0);
		initIsotopeProperty(59, 146, 145.91764, 24.15, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 147, 146.918996, 13.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 148, 147.922135, 2.29, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 149, 148.92372, 2.26, HLU_MINUTE, 0.0);
		initIsotopeProperty(59, 150, 149.926673, 6.19, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 151, 150.928319, 18.90, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 152, 151.93150, 3.63, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 153, 152.93384, 4.28, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 154, 153.93752, 2.3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 155, 154.94012, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 156, 155.94427, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 157, 156.94743, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 158, 157.95198, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(59, 159, 158.95550, 100.0e-3, HLU_SECOND, 0.0);
		initDecayMode3(59, 121, DECAY_MODE_PROTON_EMISSION, 99.99, DECAY_MODE_BETA_PLUS, 0.009, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(59, 122, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 123, DECAY_MODE_BETA_PLUS);
		initDecayMode2(59, 124, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(59, 125, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(59, 126, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(59, 127, DECAY_MODE_BETA_PLUS);
		initDecayMode2(59, 128, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(59, 129, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 130, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 131, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 132, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 133, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 134, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 135, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 136, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 137, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 138, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 139, DECAY_MODE_BETA_PLUS);
		initDecayMode(59, 140, DECAY_MODE_BETA_PLUS);
		initDecayMode2(59, 142, DECAY_MODE_BETA_MINUS, 99.98, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(59, 143, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 144, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 145, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 146, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 147, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 148, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 149, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 150, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 152, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 153, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 154, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 155, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 156, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 157, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 158, DECAY_MODE_BETA_MINUS);
		initDecayMode(59, 159, DECAY_MODE_BETA_MINUS);
	
		initAtomProperty(60, "Nd", "neodymium", "ネオジム", 126, 161, 126, 161, 27.45);
		initIsotopeProperty(60, 126, 125.94322, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 127, 126.94050, 1.8, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 128, 127.93539, 5.0, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 129, 128.93319, 4.9, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 130, 129.92851, 21.0, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 131, 130.92725, 33.0, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 132, 131.923321, 1.56, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 133, 132.92235, 70.0, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 134, 133.918790, 8.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 135, 134.918181, 12.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 136, 135.914976, 50.65, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 137, 136.914567, 38.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 138, 137.911950, 5.04, HLU_HOUR, 0.0);
		initIsotopeProperty(60, 139, 138.911978, 29.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 140, 139.90955, 3.37, HLU_DAY, 0.0);
		initIsotopeProperty(60, 141, 140.909610, 2.49, HLU_HOUR, 0.0);
		initIsotopeProperty(60, 142, 141.9077233, 0.0, HLU_STABLE, 0.272);
		initIsotopeProperty(60, 143, 142.9098143, 0.0, HLU_STABLE, 0.122);
		initIsotopeProperty(60, 144, 143.9100873, 2.29e15, HLU_YEAR, 0.238);
		initIsotopeProperty(60, 145, 144.9125736, 0.0, HLU_STABLE, 0.083);
		initIsotopeProperty(60, 146, 145.9131169, 0.0, HLU_STABLE, 0.172);
		initIsotopeProperty(60, 147, 146.9161004, 10.98, HLU_DAY, 0.0);
		initIsotopeProperty(60, 148, 147.916893, 0.0, HLU_STABLE, 0.057);
		initIsotopeProperty(60, 149, 148.920149, 1.728, HLU_HOUR, 0.0);
		initIsotopeProperty(60, 150, 149.920891, 6.7e18, HLU_YEAR, 0.056);
		initIsotopeProperty(60, 151, 150.923829, 12.44, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 152, 151.924682, 11.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(60, 153, 152.927698, 31.6, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 154, 153.92948, 25.9, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 155, 154.93293, 8.9, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 156, 155.93502, 5.49, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 157, 156.93903, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 158, 157.94160, 700.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 159, 158.94609, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 160, 159.94909, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(60, 161, 160.95388, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode(60, 126, DECAY_MODE_BETA_PLUS);
		initDecayMode2(60, 127, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(60, 128, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(60, 129, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(60, 130, DECAY_MODE_BETA_PLUS);
		initDecayMode2(60, 131, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(60, 132, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 133, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 134, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 135, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 136, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 137, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 138, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 139, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 140, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(60, 141, DECAY_MODE_BETA_PLUS);
		initDecayMode(60, 144, DECAY_MODE_ALPHA);
		initDecayMode(60, 147, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 149, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 150, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode(60, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 152, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 153, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 154, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 155, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 156, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 157, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 158, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 159, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 160, DECAY_MODE_BETA_MINUS);
		initDecayMode(60, 161, DECAY_MODE_BETA_MINUS);

		initAtomProperty(61, "Pm", "promethium", "プロメチウム", 128, 163, 128, 163, 27.45);//Heat capacity is unkown
		initIsotopeProperty(61, 128, 127.94842, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 129, 128.94316, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 130, 129.94045, 2.6, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 131, 130.93587, 6.3, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 132, 131.93375, 6.2, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 133, 132.92978, 15.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 134, 133.92835, 22.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 135, 134.92488, 49.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 136, 135.92357, 107.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 137, 136.920479, 2.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(61, 138, 137.919548, 10.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 139, 138.916804, 4.15, HLU_MINUTE, 0.0);
		initIsotopeProperty(61, 140, 139.91604, 9.2, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 141, 140.913555, 20.90, HLU_MINUTE, 0.0);
		initIsotopeProperty(61, 142, 141.912874, 40.5, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 143, 142.910933, 265.0, HLU_DAY, 0.0);
		initIsotopeProperty(61, 144, 143.912591, 363.0, HLU_DAY, 0.0);
		initIsotopeProperty(61, 145, 144.912749, 17.7, HLU_YEAR, 0.0);
		initIsotopeProperty(61, 146, 145.914696, 5.53, HLU_YEAR, 0.0);
		initIsotopeProperty(61, 147, 146.9151385, 2.6234, HLU_YEAR, 0.0);
		initIsotopeProperty(61, 148, 147.917475, 5.368, HLU_DAY, 0.0);
		initIsotopeProperty(61, 149, 148.918334, 53.08, HLU_HOUR, 0.0);
		initIsotopeProperty(61, 150, 149.920984, 2.68, HLU_HOUR, 0.0);
		initIsotopeProperty(61, 151, 150.921207, 28.40, HLU_HOUR, 0.0);
		initIsotopeProperty(61, 152, 151.923497, 4.12, HLU_MINUTE, 0.0);
		initIsotopeProperty(61, 153, 152.924117, 5.25, HLU_MINUTE, 0.0);
		initIsotopeProperty(61, 154, 153.92646, 1.73, HLU_MINUTE, 0.0);
		initIsotopeProperty(61, 155, 154.92810, 41.5, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 156, 155.93106, 26.70, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 157, 156.93304, 10.56, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 158, 157.93656, 4.8, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 159, 158.93897, 1.47, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 160, 159.94299, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 161, 160.94586, 700.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 162, 161.95029, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(61, 163, 162.95368, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(61, 128, DECAY_MODE_BETA_PLUS, 60.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(61, 129, DECAY_MODE_BETA_PLUS);
		initDecayMode2(61, 130, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(61, 131, DECAY_MODE_BETA_PLUS_PROTON, 60.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(61, 132, DECAY_MODE_BETA_PLUS, 100.0 - 5.0e-5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(61, 133, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 134, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 135, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 136, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 137, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 138, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 139, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 140, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 141, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 142, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 143, DECAY_MODE_BETA_PLUS);
		initDecayMode(61, 144, DECAY_MODE_BETA_PLUS);
		initDecayMode2(61, 145, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 2.8e-7, DECAY_MODE_ALPHA);
		initDecayMode2(61, 146, DECAY_MODE_ELECTRON_CAPTURE, 66.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 147, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 148, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 149, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 150, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 152, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 153, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 154, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 155, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 156, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 157, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 158, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 159, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 160, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 161, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 162, DECAY_MODE_BETA_MINUS);
		initDecayMode(61, 163, DECAY_MODE_BETA_MINUS);
	
		initAtomProperty(62, "Sm", "samarium", "サマリウム", 130, 165, 130, 165, 29.54);
		//initIsotopeProperty(62, 128, 127.95808, 0.5, HLU_SECOND, 0.0);
		//initIsotopeProperty(62, 129, 128.95464, 550.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 130, 129.94892, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 131, 130.94611, 1.2, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 132, 131.94069, 4.0, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 133, 132.93867, 2.90, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 134, 133.93397, 10.0, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 135, 134.93252, 10.3, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 136, 135.928276, 47.0, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 137, 136.92697, 45.0, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 138, 137.923244, 3.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 139, 138.922297, 2.57, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 140, 139.918995, 14.82, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 141, 140.918476, 10.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 142, 141.915198, 72.49, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 143, 142.914628, 8.75, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 144, 143.911999, 0.0, HLU_STABLE, 0.0307);
		initIsotopeProperty(62, 145, 144.913410, 340.0, HLU_DAY, 0.0);
		initIsotopeProperty(62, 146, 145.913041, 6.8e7, HLU_YEAR, 0.0);
		initIsotopeProperty(62, 147, 146.9148979, 1.06e11, HLU_YEAR, 0.0);
		initIsotopeProperty(62, 148, 147.9148227, 7.0e15, HLU_YEAR, 0.1124);
		initIsotopeProperty(62, 149, 148.9171847, 0.0, HLU_STABLE, 0.1382);
		initIsotopeProperty(62, 150, 149.9172755, 0.0, HLU_STABLE, 0.0738);
		initIsotopeProperty(62, 151, 150.9199324, 88.8, HLU_YEAR, 0.0);
		initIsotopeProperty(62, 152, 151.9197324, 0.0, HLU_STABLE, 0.2675);
		initIsotopeProperty(62, 153, 152.9220974, 46.284, HLU_HOUR, 0.0);
		initIsotopeProperty(62, 154, 153.9222093, 0.0, HLU_STABLE, 0.2275);
		initIsotopeProperty(62, 155, 154.9246402, 22.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 156, 155.925528, 9.4, HLU_HOUR, 0.0);
		initIsotopeProperty(62, 157, 156.92836, 8.03, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 158, 157.92999, 5.30, HLU_MINUTE, 0.0);
		initIsotopeProperty(62, 159, 158.93321, 11.37, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 160, 159.93514, 9.6, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 161, 160.93883, 4.8, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 162, 161.94122, 2.4, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 163, 162.94536, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 164, 163.94828, 500.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(62, 165, 164.95298, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode(62, 130, DECAY_MODE_BETA_PLUS);
		initDecayMode2(62, 131, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(62, 132, DECAY_MODE_BETA_PLUS, 99.995, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(62, 133, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(62, 134, DECAY_MODE_BETA_PLUS);
		initDecayMode2(62, 135, DECAY_MODE_BETA_PLUS, 99.98, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(62, 136, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 137, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 138, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 139, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 140, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 141, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 142, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 143, DECAY_MODE_BETA_PLUS);
		initDecayMode(62, 145, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(62, 146, DECAY_MODE_ALPHA);
		initDecayMode(62, 147, DECAY_MODE_ALPHA);
		initDecayMode(62, 148, DECAY_MODE_ALPHA);
		initDecayMode(62, 151, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 153, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 155, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 156, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 157, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 158, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 159, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 160, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 161, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 162, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 163, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 164, DECAY_MODE_BETA_MINUS);
		initDecayMode(62, 165, DECAY_MODE_BETA_MINUS);

		initAtomProperty(63, "Eu", "europium", "ユウロピウム", 132, 167, 132, 167, 27.66);
		initIsotopeProperty(63, 132, 131.95437, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 133, 132.94924, 200.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 134, 133.94651, 0.5, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 135, 134.94182, 1.5, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 136, 135.93960, 3.3, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 137, 136.93557, 8.4, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 138, 137.93371, 12.1, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 139, 138.929792, 17.9, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 140, 139.92809, 1.51, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 141, 140.924931, 40.7, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 142, 141.92343, 2.36, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 143, 142.920298, 2.59, HLU_MINUTE, 0.0);
		initIsotopeProperty(63, 144, 143.918817, 10.2, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 145, 144.916265, 5.93, HLU_DAY, 0.0);
		initIsotopeProperty(63, 146, 145.917206, 4.61, HLU_DAY, 0.0);
		initIsotopeProperty(63, 147, 146.916746, 24.1, HLU_DAY, 0.0);
		initIsotopeProperty(63, 148, 147.918086, 54.5, HLU_DAY, 0.0);
		initIsotopeProperty(63, 149, 148.917931, 93.1, HLU_DAY, 0.0);
		initIsotopeProperty(63, 150, 149.919702, 36.9, HLU_YEAR, 0.0);
		initIsotopeProperty(63, 151, 150.9198502, 4.62e18, HLU_YEAR, 0.4781);
		initIsotopeProperty(63, 152, 151.9217445, 13.537, HLU_YEAR, 0.0);
		initIsotopeProperty(63, 153, 152.9212303, 0.0, HLU_STABLE, 0.5219);
		initIsotopeProperty(63, 154, 153.9229792, 8.593, HLU_YEAR, 0.0);
		initIsotopeProperty(63, 155, 154.9228933, 4.7611, HLU_YEAR, 0.0);
		initIsotopeProperty(63, 156, 155.924752, 15.19, HLU_DAY, 0.0);
		initIsotopeProperty(63, 157, 156.925424, 15.18, HLU_HOUR, 0.0);
		initIsotopeProperty(63, 158, 157.92785, 45.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(63, 159, 158.929089, 18.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(63, 160, 159.93197, 38.0, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 161, 160.93368, 26.0, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 162, 161.93704, 10.6, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 163, 162.93921, 6.0, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 164, 163.94299, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 165, 164.94572, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 166, 165.94997, 400.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(63, 167, 166.95321, 200.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(63, 132, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(63, 133, DECAY_MODE_BETA_PLUS);
		initDecayMode2(63, 134, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(63, 135, DECAY_MODE_BETA_PLUS, 99.995, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(63, 136, DECAY_MODE_BETA_PLUS, 99.91, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(63, 137, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 138, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 139, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 140, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 141, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 142, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 143, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 144, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 145, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 146, DECAY_MODE_BETA_PLUS);
		initDecayMode2(63, 147, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_ALPHA);
		initDecayMode2(63, 148, DECAY_MODE_BETA_PLUS, 100.0 - 9.39e-7, DECAY_MODE_ALPHA);
		initDecayMode(63, 149, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(63, 150, DECAY_MODE_BETA_PLUS);
		initDecayMode(63, 151, DECAY_MODE_ALPHA);
		initDecayMode3(63, 152, DECAY_MODE_ELECTRON_CAPTURE, 72.09, DECAY_MODE_BETA_PLUS, 0.027, DECAY_MODE_BETA_MINUS);
		initDecayMode2(63, 154, DECAY_MODE_BETA_MINUS, 99.98, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(63, 155, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 156, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 157, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 158, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 159, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 160, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 161, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 162, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 163, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 164, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 165, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 166, DECAY_MODE_BETA_MINUS);
		initDecayMode(63, 167, DECAY_MODE_BETA_MINUS);

		initAtomProperty(64, "Gd", "gadolinium", "ガドリニウム", 136, 169, 136, 169, 37.03);
		initIsotopeProperty(64, 136, 135.94734, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 137, 136.94502, 2.2, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 138, 137.94012, 4.7, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 139, 138.93824, 5.7, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 140, 139.93367, 15.8, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 141, 140.932126, 14.0, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 142, 141.92812, 70.2, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 143, 142.92675, 39.0, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 144, 143.92296, 4.47, HLU_MINUTE, 0.0);
		initIsotopeProperty(64, 145, 144.921709, 23.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(64, 146, 145.918311, 48.27, HLU_DAY, 0.0);
		initIsotopeProperty(64, 147, 146.919094, 38.06, HLU_HOUR, 0.0);
		initIsotopeProperty(64, 148, 147.918115, 74.6, HLU_YEAR, 0.0);
		initIsotopeProperty(64, 149, 148.919341, 9.28, HLU_DAY, 0.0);
		initIsotopeProperty(64, 150, 149.918659, 1.79e6, HLU_YEAR, 0.0);
		initIsotopeProperty(64, 151, 150.920348, 124.0, HLU_DAY, 0.0);
		initIsotopeProperty(64, 152, 151.9197910, 1.08e14, HLU_YEAR, 0.002);
		initIsotopeProperty(64, 153, 152.9217495, 240.4, HLU_DAY, 0.0);
		initIsotopeProperty(64, 154, 153.9208656, 0.0, HLU_STABLE, 0.0218);
		initIsotopeProperty(64, 155, 154.9226220, 0.0, HLU_STABLE, 0.148);
		initIsotopeProperty(64, 156, 155.9221227, 0.0, HLU_STABLE, 0.2047);
		initIsotopeProperty(64, 157, 156.9239601, 0.0, HLU_STABLE, 0.1565);
		initIsotopeProperty(64, 158, 157.9241039, 0.0, HLU_STABLE, 0.2484);
		initIsotopeProperty(64, 159, 158.9263887, 18.479, HLU_HOUR, 0.0);
		initIsotopeProperty(64, 160, 159.9270541, 0.0, HLU_STABLE, 0.2186);
		initIsotopeProperty(64, 161, 160.9296692, 3.646, HLU_MINUTE, 0.0);
		initIsotopeProperty(64, 162, 161.930985, 8.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(64, 163, 162.93399, 68.0, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 164, 163.93586, 45.0, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 165, 164.93938, 10.3, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 166, 165.94160, 4.8, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 167, 166.94557, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 168, 167.94836, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(64, 169, 168.95287, 1.0, HLU_SECOND, 0.0);
		initDecayMode(64, 136, DECAY_MODE_BETA_PLUS);
		initDecayMode2(64, 137, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(64, 138, DECAY_MODE_BETA_PLUS);
		initDecayMode2(64, 139, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(64, 140, DECAY_MODE_BETA_PLUS);
		initDecayMode2(64, 141, DECAY_MODE_BETA_PLUS, 99.97, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(64, 142, DECAY_MODE_BETA_PLUS);
		initDecayMode3(64, 143, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_ALPHA, 0.0005, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(64, 144, DECAY_MODE_BETA_PLUS);
		initDecayMode(64, 145, DECAY_MODE_BETA_PLUS);
		initDecayMode(64, 146, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(64, 147, DECAY_MODE_BETA_PLUS);
		initDecayMode2(64, 148, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode2(64, 149, DECAY_MODE_BETA_PLUS, 100.0 - 4.34e-4, DECAY_MODE_ALPHA);
		initDecayMode2(64, 150, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode2(64, 151, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 1.0e-6, DECAY_MODE_ALPHA);
		initDecayMode(64, 152, DECAY_MODE_ALPHA);
		initDecayMode(64, 153, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(64, 159, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 161, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 162, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 163, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 164, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 165, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 166, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 167, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 168, DECAY_MODE_BETA_MINUS);
		initDecayMode(64, 169, DECAY_MODE_BETA_MINUS);

		initAtomProperty(65, "Tb", "terbium", "テルビウム", 138, 171, 138, 171, 28.91);
		initIsotopeProperty(65, 138, 137.95316, 800.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 139, 138.94829, 1.6, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 140, 139.94581, 2.4, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 141, 140.94145, 3.5, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 142, 141.93874, 597.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 143, 142.93512, 12.0, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 144, 143.93305, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 145, 144.92927, 20.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(65, 146, 145.92725, 8.0, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 147, 146.924045, 1.64, HLU_HOUR, 0.0);
		initIsotopeProperty(65, 148, 147.924272, 60.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(65, 149, 148.923246, 4.118, HLU_HOUR, 0.0);
		initIsotopeProperty(65, 150, 149.923660, 3.48, HLU_HOUR, 0.0);
		initIsotopeProperty(65, 151, 150.923103, 17.609, HLU_HOUR, 0.0);
		initIsotopeProperty(65, 152, 151.92407, 17.5, HLU_HOUR, 0.0);
		initIsotopeProperty(65, 153, 152.923435, 2.34, HLU_DAY, 0.0);
		initIsotopeProperty(65, 154, 153.92468, 21.5, HLU_HOUR, 0.0);
		initIsotopeProperty(65, 155, 154.923505, 5.32, HLU_DAY, 0.0);
		initIsotopeProperty(65, 156, 155.924747, 5.35, HLU_DAY, 0.0);
		initIsotopeProperty(65, 157, 156.9240246, 71.0, HLU_YEAR, 0.0);
		initIsotopeProperty(65, 158, 157.9254131, 180.0, HLU_YEAR, 0.0);
		initIsotopeProperty(65, 159, 158.9253468, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(65, 160, 159.9271676, 72.3, HLU_DAY, 0.0);
		initIsotopeProperty(65, 161, 160.9275699, 6.906, HLU_DAY, 0.0);
		initIsotopeProperty(65, 162, 161.92949, 7.60, HLU_MINUTE, 0.0);
		initIsotopeProperty(65, 163, 162.930648, 19.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(65, 164, 163.93335, 3.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(65, 165, 164.93488, 2.11, HLU_MINUTE, 0.0);
		initIsotopeProperty(65, 166, 165.93799, 25.6, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 167, 166.94005, 19.4, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 168, 167.94364, 8.2, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 169, 168.94622, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 170, 169.95025, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(65, 171, 170.95330, 500.0e-3, HLU_SECOND, 0.0);
		initDecayMode2(65, 138, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(65, 139, DECAY_MODE_BETA_PLUS);
		initDecayMode2(65, 140, DECAY_MODE_BETA_PLUS, 99.74, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(65, 141, DECAY_MODE_BETA_PLUS);
		initDecayMode2(65, 142, DECAY_MODE_BETA_PLUS, 99.5, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(65, 143, DECAY_MODE_BETA_PLUS);
		initDecayMode2(65, 144, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(65, 145, DECAY_MODE_BETA_PLUS);
		initDecayMode(65, 146, DECAY_MODE_BETA_PLUS);
		initDecayMode(65, 147, DECAY_MODE_BETA_PLUS);
		initDecayMode(65, 148, DECAY_MODE_BETA_PLUS);
		initDecayMode2(65, 149, DECAY_MODE_BETA_PLUS, 83.3, DECAY_MODE_ALPHA);
		initDecayMode2(65, 150, DECAY_MODE_BETA_PLUS, 99.95, DECAY_MODE_ALPHA);
		initDecayMode2(65, 151, DECAY_MODE_BETA_PLUS, 100.0 - 0.0095, DECAY_MODE_ALPHA);
		initDecayMode2(65, 152, DECAY_MODE_BETA_PLUS, 100.0 - 7.0e-7, DECAY_MODE_ALPHA);
		initDecayMode(65, 153, DECAY_MODE_BETA_PLUS);
		initDecayMode2(65, 154, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 155, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(65, 156, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 157, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(65, 158, DECAY_MODE_BETA_PLUS, 83.4, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 160, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 161, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 162, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 163, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 164, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 165, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 166, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 167, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 168, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 169, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 170, DECAY_MODE_BETA_MINUS);
		initDecayMode(65, 171, DECAY_MODE_BETA_MINUS);

		initAtomProperty(66, "Dy", "dysprosium", "ジスプロシウム", 140, 173, 140, 173, 27.7);
		initIsotopeProperty(66, 140, 139.95401, 700.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 141, 140.95135, 0.9, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 142, 141.94637, 2.3, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 143, 142.94383, 5.6, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 144, 143.93925, 9.1, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 145, 144.93743, 9.5, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 146, 145.932845, 33.2, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 147, 146.931092, 40.0, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 148, 147.927150, 3.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(66, 149, 148.927305, 4.20, HLU_MINUTE, 0.0);
		initIsotopeProperty(66, 150, 149.925585, 7.17, HLU_MINUTE, 0.0);
		initIsotopeProperty(66, 151, 150.926185, 17.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(66, 152, 151.924718, 2.38, HLU_HOUR, 0.0);
		initIsotopeProperty(66, 153, 152.925765, 6.4, HLU_HOUR, 0.0);
		initIsotopeProperty(66, 154, 153.924424, 3.0e6, HLU_YEAR, 0.0);
		initIsotopeProperty(66, 155, 154.925754, 9.9, HLU_HOUR, 0.0);
		initIsotopeProperty(66, 156, 155.924283, 0.0, HLU_STABLE, 5.60E-040);
		initIsotopeProperty(66, 157, 156.925466, 8.14, HLU_HOUR, 0.0);
		initIsotopeProperty(66, 158, 157.924409, 0.0, HLU_STABLE, 9.50E-040);
		initIsotopeProperty(66, 159, 158.9257392, 144.4, HLU_DAY, 0.0);
		initIsotopeProperty(66, 160, 159.9251975, 0.0, HLU_STABLE, 0.023290);
		initIsotopeProperty(66, 161, 160.9269334, 0.0, HLU_STABLE, 0.188890);
		initIsotopeProperty(66, 162, 161.9267984, 0.0, HLU_STABLE, 0.254750);
		initIsotopeProperty(66, 163, 162.9287312, 0.0, HLU_STABLE, 0.248960);
		initIsotopeProperty(66, 164, 163.9291748, 0.0, HLU_STABLE, 0.28260);
		initIsotopeProperty(66, 165, 164.9317033, 2.334, HLU_HOUR, 0.0);
		initIsotopeProperty(66, 166, 165.9328067, 81.6, HLU_HOUR, 0.0);
		initIsotopeProperty(66, 167, 166.93566, 6.20, HLU_MINUTE, 0.0);
		initIsotopeProperty(66, 168, 167.93713, 8.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(66, 169, 168.94031, 39.0, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 170, 169.94239, 30.0, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 171, 170.94620, 6.0, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 172, 171.94876, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(66, 173, 172.95300, 2.0, HLU_SECOND, 0.0);
		initDecayMode(66, 140, DECAY_MODE_BETA_PLUS);
		initDecayMode2(66, 141, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(66, 142, DECAY_MODE_BETA_PLUS, 99.94, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(66, 143, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(66, 144, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(66, 145, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(66, 146, DECAY_MODE_BETA_PLUS);
		initDecayMode2(66, 147, DECAY_MODE_BETA_PLUS, 99.95, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(66, 148, DECAY_MODE_BETA_PLUS);
		initDecayMode(66, 149, DECAY_MODE_BETA_PLUS);
		initDecayMode2(66, 150, DECAY_MODE_BETA_PLUS, 64.0, DECAY_MODE_ALPHA);
		initDecayMode2(66, 151, DECAY_MODE_BETA_PLUS, 94.4, DECAY_MODE_ALPHA);
		initDecayMode2(66, 152, DECAY_MODE_ELECTRON_CAPTURE, 99.9, DECAY_MODE_ALPHA);
		initDecayMode2(66, 153, DECAY_MODE_BETA_PLUS, 100.0 - 0.00939, DECAY_MODE_ALPHA);
		initDecayMode2(66, 154, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode(66, 155, DECAY_MODE_BETA_PLUS);
		initDecayMode(66, 157, DECAY_MODE_BETA_PLUS);
		initDecayMode(66, 159, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(66, 165, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 166, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 167, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 168, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 169, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 170, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 171, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 172, DECAY_MODE_BETA_MINUS);
		initDecayMode(66, 173, DECAY_MODE_BETA_MINUS);

		initAtomProperty(67, "Ho", "holmium", "ホルミウム", 142, 173, 142, 173, 27.15);
		initIsotopeProperty(67, 142, 141.95977, 400.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 143, 142.95461, 300.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 144, 143.95148, 0.7, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 145, 144.94720, 2.4, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 146, 145.94464, 3.6, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 147, 146.94006, 5.8, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 148, 147.93772, 2.2, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 149, 148.933775, 21.1, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 150, 149.933496, 76.8, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 151, 150.931688, 35.2, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 152, 151.931714, 161.8, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 153, 152.930199, 2.01, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 154, 153.930602, 11.76, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 155, 154.929103, 48.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 156, 155.92984, 56.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 157, 156.928256, 12.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 158, 157.928941, 11.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 159, 158.927712, 33.05, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 160, 159.928729, 25.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 161, 160.927855, 2.48, HLU_HOUR, 0.0);
		initIsotopeProperty(67, 162, 161.929096, 15.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 163, 162.9287339, 4570.0, HLU_YEAR, 0.0);
		initIsotopeProperty(67, 164, 163.9302335, 29.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 165, 164.9303221, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(67, 166, 165.9322842, 26.83, HLU_HOUR, 0.0);
		initIsotopeProperty(67, 167, 166.933133, 3.003, HLU_HOUR, 0.0);
		initIsotopeProperty(67, 168, 167.93552, 2.99, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 169, 168.936872, 4.72, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 170, 169.93962, 2.76, HLU_MINUTE, 0.0);
		initIsotopeProperty(67, 171, 170.94147, 53.0, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 172, 171.94482, 25.0, HLU_SECOND, 0.0);
		initIsotopeProperty(67, 173, 172.94729, 10.0, HLU_SECOND, 0.0);
		//initIsotopeProperty(67, 174, 173.95115, 8.0, HLU_SECOND, 0.0);
		//initIsotopeProperty(67, 175, 174.95405, 5.0, HLU_SECOND, 0.0);
		initDecayMode2(67, 142, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(67, 143, DECAY_MODE_BETA_PLUS);
		initDecayMode2(67, 144, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(67, 145, DECAY_MODE_BETA_PLUS);
		initDecayMode2(67, 146, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(67, 147, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(67, 148, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 149, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 150, DECAY_MODE_BETA_PLUS);
		initDecayMode2(67, 151, DECAY_MODE_BETA_PLUS, 78.0, DECAY_MODE_ALPHA);
		initDecayMode2(67, 152, DECAY_MODE_BETA_PLUS, 88.0, DECAY_MODE_ALPHA);
		initDecayMode2(67, 153, DECAY_MODE_BETA_PLUS, 99.94, DECAY_MODE_ALPHA);
		initDecayMode2(67, 154, DECAY_MODE_BETA_PLUS, 99.98, DECAY_MODE_ALPHA);
		initDecayMode(67, 155, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 156, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 157, DECAY_MODE_BETA_PLUS);
		initDecayMode2(67, 158, DECAY_MODE_BETA_PLUS, 93.0, DECAY_MODE_ALPHA);
		initDecayMode(67, 159, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 160, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 161, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(67, 162, DECAY_MODE_BETA_PLUS);
		initDecayMode(67, 163, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(67, 164, DECAY_MODE_ELECTRON_CAPTURE, 60.0, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 166, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 167, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 168, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 169, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 170, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 171, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 172, DECAY_MODE_BETA_MINUS);
		initDecayMode(67, 173, DECAY_MODE_BETA_MINUS);

		initAtomProperty(68, "Er", "erbium", "エルビウム", 144, 177, 144, 177, 28.12);
		initIsotopeProperty(68, 144, 143.96038, 400.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 145, 144.95739, 900.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 146, 145.95200, 1.7, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 147, 146.94949, 2.5, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 148, 147.94455, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 149, 148.94231, 4.0, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 150, 149.937914, 18.5, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 151, 150.937449, 23.5, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 152, 151.935050, 10.3, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 153, 152.935063, 37.1, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 154, 153.932783, 3.73, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 155, 154.933209, 5.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 156, 155.931065, 19.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 157, 156.93192, 18.65, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 158, 157.929893, 2.29, HLU_HOUR, 0.0);
		initIsotopeProperty(68, 159, 158.930684, 36.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 160, 159.929083, 28.58, HLU_HOUR, 0.0);
		initIsotopeProperty(68, 161, 160.929995, 3.21, HLU_HOUR, 0.0);
		initIsotopeProperty(68, 162, 161.928778, 0.0, HLU_STABLE, 0.00139);
		initIsotopeProperty(68, 163, 162.930033, 75.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 164, 163.929200, 0.0, HLU_STABLE, 0.01601);
		initIsotopeProperty(68, 165, 164.930726, 10.36, HLU_HOUR, 0.0);
		initIsotopeProperty(68, 166, 165.9302931, 0.0, HLU_STABLE, 0.33503);
		initIsotopeProperty(68, 167, 166.9320482, 0.0, HLU_STABLE, 0.22869);
		initIsotopeProperty(68, 168, 167.9323702, 0.0, HLU_STABLE, 0.26978);
		initIsotopeProperty(68, 169, 168.9345904, 9.392, HLU_DAY, 0.0);
		initIsotopeProperty(68, 170, 169.9354643, 0.0, HLU_STABLE, 0.1491);
		initIsotopeProperty(68, 171, 170.9380298, 7.516, HLU_HOUR, 0.0);
		initIsotopeProperty(68, 172, 171.939356, 49.3, HLU_HOUR, 0.0);
		initIsotopeProperty(68, 173, 172.94240, 0.434, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 174, 173.94423, 3.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 175, 174.94777, 1.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(68, 176, 175.95008, 20.0, HLU_SECOND, 0.0);
		initIsotopeProperty(68, 177, 176.95405, 3.0, HLU_SECOND, 0.0);
		initDecayMode(68, 144, DECAY_MODE_BETA_PLUS);
		initDecayMode2(68, 145, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(68, 146, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(68, 147, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(68, 148, DECAY_MODE_BETA_PLUS, 99.85, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(68, 149, DECAY_MODE_BETA_PLUS, 93.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(68, 150, DECAY_MODE_BETA_PLUS);
		initDecayMode(68, 151, DECAY_MODE_BETA_PLUS);
		initDecayMode2(68, 152, DECAY_MODE_ALPHA, 90.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(68, 153, DECAY_MODE_ALPHA, 53.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(68, 154, DECAY_MODE_BETA_PLUS, 99.53, DECAY_MODE_ALPHA);
		initDecayMode2(68, 155, DECAY_MODE_BETA_PLUS, 99.98, DECAY_MODE_ALPHA);
		initDecayMode(68, 156, DECAY_MODE_BETA_PLUS);
		initDecayMode(68, 157, DECAY_MODE_BETA_PLUS);
		initDecayMode(68, 158, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(68, 159, DECAY_MODE_BETA_PLUS);
		initDecayMode(68, 160, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(68, 161, DECAY_MODE_BETA_PLUS);
		initDecayMode(68, 163, DECAY_MODE_BETA_PLUS);
		initDecayMode(68, 165, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(68, 169, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 171, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 172, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 173, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 174, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 175, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 176, DECAY_MODE_BETA_MINUS);
		initDecayMode(68, 177, DECAY_MODE_BETA_MINUS);

		initAtomProperty(69, "Tm", "thulium", "ツリウム", 146, 179, 146, 179, 27.03);
		initIsotopeProperty(69, 146, 145.96643, 240.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 147, 146.96096, 0.58, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 148, 147.95784, 0.7, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 149, 148.95272, 0.9, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 150, 149.94996, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 151, 150.945483, 4.17, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 152, 151.94442, 8.0, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 153, 152.942012, 1.48, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 154, 153.941568, 8.1, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 155, 154.939199, 21.6, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 156, 155.938980, 83.8, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 157, 156.93697, 3.63, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 158, 157.936980, 3.98, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 159, 158.93498, 9.13, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 160, 159.93526, 9.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 161, 160.93355, 30.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 162, 161.933995, 21.70, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 163, 162.932651, 1.810, HLU_HOUR, 0.0);
		initIsotopeProperty(69, 164, 163.93356, 2.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 165, 164.932435, 30.06, HLU_HOUR, 0.0);
		initIsotopeProperty(69, 166, 165.933554, 7.70, HLU_HOUR, 0.0);
		initIsotopeProperty(69, 167, 166.9328516, 9.25, HLU_DAY, 0.0);
		initIsotopeProperty(69, 168, 167.934173, 93.1, HLU_DAY, 0.0);
		initIsotopeProperty(69, 169, 168.9342133, 0.0, HLU_STABLE, 1.0);
		initIsotopeProperty(69, 170, 169.9358014, 128.6, HLU_DAY, 0.0);
		initIsotopeProperty(69, 171, 170.9364294, 1.92, HLU_YEAR, 0.0);
		initIsotopeProperty(69, 172, 171.938400, 63.6, HLU_HOUR, 0.0);
		initIsotopeProperty(69, 173, 172.939604, 8.24, HLU_HOUR, 0.0);
		initIsotopeProperty(69, 174, 173.94217, 5.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 175, 174.94384, 15.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 176, 175.94699, 1.85, HLU_MINUTE, 0.0);
		initIsotopeProperty(69, 177, 176.94904, 90.0, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 178, 177.95264, 30.0, HLU_SECOND, 0.0);
		initIsotopeProperty(69, 179, 178.95534, 20.0, HLU_SECOND, 0.0);
		initDecayMode2(69, 146, DECAY_MODE_PROTON_EMISSION, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(69, 147, DECAY_MODE_BETA_PLUS, 85.0, DECAY_MODE_PROTON_EMISSION);
		initDecayMode(69, 148, DECAY_MODE_BETA_PLUS);
		initDecayMode2(69, 149, DECAY_MODE_BETA_PLUS, 99.74, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode(69, 150, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 151, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 152, DECAY_MODE_BETA_PLUS);
		initDecayMode2(69, 153, DECAY_MODE_ALPHA, 91.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(69, 154, DECAY_MODE_BETA_PLUS, 56.0, DECAY_MODE_ALPHA);
		initDecayMode2(69, 155, DECAY_MODE_BETA_PLUS, 98.1, DECAY_MODE_ALPHA);
		initDecayMode2(69, 156, DECAY_MODE_BETA_PLUS, 99.93, DECAY_MODE_ALPHA);
		initDecayMode(69, 157, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 158, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 159, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 160, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 161, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 162, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 163, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 164, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 165, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 166, DECAY_MODE_BETA_PLUS);
		initDecayMode(69, 167, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(69, 168, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_BETA_MINUS);
		initDecayMode2(69, 170, DECAY_MODE_BETA_MINUS, 99.86, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(69, 171, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 172, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 173, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 174, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 175, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 176, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 177, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 178, DECAY_MODE_BETA_MINUS);
		initDecayMode(69, 179, DECAY_MODE_BETA_MINUS);

		initAtomProperty(70, "Yb", "ytterbium", "イッテルビウム", 148, 181, 148, 181, 26.74);
		initIsotopeProperty(70, 148, 147.96742, 250.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 149, 148.96404, 0.7, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 150, 149.95842, 700.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 151, 150.95540, 1.6, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 152, 151.95029, 3.04, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 153, 152.94948, 4.2, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 154, 153.946394, 0.409, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 155, 154.945782, 1.793, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 156, 155.942818, 26.1, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 157, 156.942628, 38.6, HLU_SECOND, 0.0);
		initIsotopeProperty(70, 158, 157.939866, 1.49, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 159, 158.94005, 1.67, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 160, 159.937552, 4.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 161, 160.937902, 4.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 162, 161.935768, 18.87, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 163, 162.936334, 11.05, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 164, 163.934489, 75.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 165, 164.93528, 9.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 166, 165.933882, 56.7, HLU_HOUR, 0.0);
		initIsotopeProperty(70, 167, 166.934950, 17.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 168, 167.933897, 0.0, HLU_STABLE, 0.0013);
		initIsotopeProperty(70, 169, 168.935190, 32.026, HLU_DAY, 0.0);
		initIsotopeProperty(70, 170, 169.9347618, 0.0, HLU_STABLE, 0.0304);
		initIsotopeProperty(70, 171, 170.9363258, 0.0, HLU_STABLE, 0.1428);
		initIsotopeProperty(70, 172, 171.9363815, 0.0, HLU_STABLE, 0.2183);
		initIsotopeProperty(70, 173, 172.9382108, 0.0, HLU_STABLE, 0.1613);
		initIsotopeProperty(70, 174, 173.9388621, 0.0, HLU_STABLE, 0.3183);
		initIsotopeProperty(70, 175, 174.9412765, 4.185, HLU_DAY, 0.0);
		initIsotopeProperty(70, 176, 175.9425717, 0.0, HLU_STABLE, 0.1276);
		initIsotopeProperty(70, 177, 176.9452608, 1.911, HLU_HOUR, 0.0);
		initIsotopeProperty(70, 178, 177.946647, 74.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 179, 178.95017, 8.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 180, 179.95233, 2.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(70, 181, 180.95615, 1.0, HLU_MINUTE, 0.0);
		initDecayMode(70, 148, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 149, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 150, DECAY_MODE_BETA_PLUS);
		initDecayMode2(70, 151, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(70, 152, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode3(70, 153, DECAY_MODE_ALPHA, 50.0, DECAY_MODE_BETA_PLUS, 50.0 - 0.008, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(70, 154, DECAY_MODE_ALPHA, 100.0 - 7.119, DECAY_MODE_BETA_PLUS);
		initDecayMode2(70, 155, DECAY_MODE_ALPHA, 89.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(70, 156, DECAY_MODE_BETA_PLUS, 90.0, DECAY_MODE_ALPHA);
		initDecayMode2(70, 157, DECAY_MODE_BETA_PLUS, 99.5, DECAY_MODE_ALPHA);
		initDecayMode2(70, 158, DECAY_MODE_BETA_PLUS, 100.0 - 0.0021, DECAY_MODE_ALPHA);
		initDecayMode(70, 159, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 160, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 161, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 162, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 163, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 164, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(70, 165, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 166, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(70, 167, DECAY_MODE_BETA_PLUS);
		initDecayMode(70, 169, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(70, 175, DECAY_MODE_BETA_MINUS);
		initDecayMode(70, 177, DECAY_MODE_BETA_MINUS);
		initDecayMode(70, 178, DECAY_MODE_BETA_MINUS);
		initDecayMode(70, 179, DECAY_MODE_BETA_MINUS);
		initDecayMode(70, 180, DECAY_MODE_BETA_MINUS);
		initDecayMode(70, 181, DECAY_MODE_BETA_MINUS);

		initAtomProperty(71, "Lu", "lutetium", "ルテチウム", 152, 184, 152, 184, 26.86);
		//initIsotopeProperty(71, 150, 149.97323, 43.0e-3, HLU_SECOND, 0.0);
		//initIsotopeProperty(71, 151, 149.97323, 80.6e-3, HLU_SECOND, 0.0);//wrong mass 
		initIsotopeProperty(71, 152, 151.96412, 650.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 153, 152.95877, 0.9, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 154, 153.95752, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 155, 154.954316, 68.6e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 156, 155.95303, 494.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 157, 156.950098, 6.8, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 158, 157.949313, 10.6, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 159, 158.94663, 12.1, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 160, 159.94603, 36.1, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 161, 160.94357, 77.0, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 162, 161.94328, 1.37, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 163, 162.94118, 3.97, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 164, 163.94134, 3.14, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 165, 164.939407, 10.74, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 166, 165.93986, 2.65, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 167, 166.93827, 51.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 168, 167.93874, 5.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 169, 168.937651, 34.06, HLU_HOUR, 0.0);
		initIsotopeProperty(71, 170, 169.938475, 2.012, HLU_DAY, 0.0);
		initIsotopeProperty(71, 171, 170.9379131, 8.24, HLU_DAY, 0.0);
		initIsotopeProperty(71, 172, 171.939086, 6.70, HLU_DAY, 0.0);
		initIsotopeProperty(71, 173, 172.9389306, 1.37, HLU_YEAR, 0.0);
		initIsotopeProperty(71, 174, 173.9403375, 3.31, HLU_YEAR, 0.0);
		initIsotopeProperty(71, 175, 174.9407718, 0.0, HLU_STABLE, 0.9741);
		initIsotopeProperty(71, 176, 175.9426863, 38.5e9, HLU_YEAR, 0.0259);
		initIsotopeProperty(71, 177, 176.9437581, 6.6475, HLU_DAY, 0.0);
		initIsotopeProperty(71, 178, 177.945955, 28.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 179, 178.947327, 4.59, HLU_HOUR, 0.0);
		initIsotopeProperty(71, 180, 179.94988, 5.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 181, 180.95197, 3.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 182, 181.95504, 2.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(71, 183, 182.95757, 58.0, HLU_SECOND, 0.0);
		initIsotopeProperty(71, 184, 183.96091, 20.0, HLU_SECOND, 0.0);
		//initDecayMode2(71, 150, DECAY_MODE_PROTON_EMISSION, 80.0, DECAY_MODE_BETA_PLUS);
		//initDecayMode2(71, 151, DECAY_MODE_PROTON_EMISSION, 63.4, DECAY_MODE_BETA_PLUS);
		initDecayMode2(71, 152, DECAY_MODE_BETA_PLUS, 85.0, DECAY_MODE_BETA_PLUS_PROTON);
		initDecayMode2(71, 153, DECAY_MODE_ALPHA, 70.0, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 154, DECAY_MODE_BETA_PLUS);
		initDecayMode2(71, 155, DECAY_MODE_ALPHA, 76.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(71, 156, DECAY_MODE_ALPHA, 95.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(71, 157, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(71, 158, DECAY_MODE_BETA_PLUS, 99.09, DECAY_MODE_ALPHA);
		initDecayMode2(71, 159, DECAY_MODE_BETA_PLUS, 99.96, DECAY_MODE_ALPHA);
		initDecayMode2(71, 160, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-4, DECAY_MODE_ALPHA);
		initDecayMode(71, 161, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 162, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 163, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 164, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 165, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 166, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 167, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 168, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 169, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 170, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 171, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 172, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 173, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(71, 174, DECAY_MODE_BETA_PLUS);
		initDecayMode(71, 176, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 177, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 178, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 179, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 180, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 181, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 182, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 183, DECAY_MODE_BETA_MINUS);
		initDecayMode(71, 184, DECAY_MODE_BETA_MINUS);

		initAtomProperty(72, "Hf", "hafnium", "ハフニウム", 154, 186, 154, 186, 25.73);
		initIsotopeProperty(72, 154, 153.96486, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 155, 154.96339, 890.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 156, 155.95936, 23.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 157, 156.95840, 115.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 158, 157.954799, 2.84, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 159, 158.953995, 5.20, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 160, 159.950684, 13.6, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 161, 160.950275, 18.2, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 162, 161.94721, 39.4, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 163, 162.94709, 40.0, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 164, 163.944367, 111.0, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 165, 164.94457, 76.0, HLU_SECOND, 0.0);
		initIsotopeProperty(72, 166, 165.94218, 6.77, HLU_MINUTE, 0.0);
		initIsotopeProperty(72, 167, 166.94260, 2.05, HLU_MINUTE, 0.0);
		initIsotopeProperty(72, 168, 167.94057, 25.95, HLU_MINUTE, 0.0);
		initIsotopeProperty(72, 169, 168.94126, 3.24, HLU_MINUTE, 0.0);
		initIsotopeProperty(72, 170, 169.93961, 16.01, HLU_HOUR, 0.0);
		initIsotopeProperty(72, 171, 170.94049, 12.1, HLU_HOUR, 0.0);
		initIsotopeProperty(72, 172, 171.939448, 1.87, HLU_YEAR, 0.0);
		initIsotopeProperty(72, 173, 172.94051, 23.6, HLU_HOUR, 0.0);
		initIsotopeProperty(72, 174, 173.940046, 2.0e15, HLU_YEAR, 0.0016);
		initIsotopeProperty(72, 175, 174.941509, 70.0, HLU_DAY, 0.0);
		initIsotopeProperty(72, 176, 175.9414086, 0.0, HLU_STABLE, 0.0526);
		initIsotopeProperty(72, 177, 176.9432207, 0.0, HLU_STABLE, 0.186);
		initIsotopeProperty(72, 178, 177.9436988, 0.0, HLU_STABLE, 0.2728);
		initIsotopeProperty(72, 179, 178.9458161, 0.0, HLU_STABLE, 0.1362);
		initIsotopeProperty(72, 180, 179.9465500, 0.0, HLU_STABLE, 0.3508);
		initIsotopeProperty(72, 181, 180.9491012, 42.39, HLU_DAY, 0.0);
		initIsotopeProperty(72, 182, 181.950554, 8.90e6, HLU_YEAR, 0.0);
		initIsotopeProperty(72, 183, 182.95353, 1.067, HLU_HOUR, 0.0);
		initIsotopeProperty(72, 184, 183.95545, 4.12, HLU_HOUR, 0.0);
		initIsotopeProperty(72, 185, 184.95882, 3.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(72, 186, 185.96089, 2.6, HLU_MINUTE, 0.0);
		initDecayMode2(72, 154, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_ALPHA);
		initDecayMode2(72, 155, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_ALPHA);
		initDecayMode2(72, 156, DECAY_MODE_ALPHA, 97.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(72, 157, DECAY_MODE_ALPHA, 86.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(72, 158, DECAY_MODE_BETA_PLUS, 55.0, DECAY_MODE_ALPHA);
		initDecayMode2(72, 159, DECAY_MODE_BETA_PLUS, 59.0, DECAY_MODE_ALPHA);
		initDecayMode2(72, 160, DECAY_MODE_BETA_PLUS, 99.3, DECAY_MODE_ALPHA);
		initDecayMode2(72, 161, DECAY_MODE_BETA_PLUS, 99.7, DECAY_MODE_ALPHA);
		initDecayMode2(72, 162, DECAY_MODE_BETA_PLUS, 100.0 - 0.008, DECAY_MODE_ALPHA);
		initDecayMode2(72, 163, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-4, DECAY_MODE_ALPHA);
		initDecayMode(72, 164, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 165, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 166, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 167, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 168, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 169, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 170, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(72, 171, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 172, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(72, 173, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 174, DECAY_MODE_ALPHA);
		initDecayMode(72, 175, DECAY_MODE_BETA_PLUS);
		initDecayMode(72, 181, DECAY_MODE_BETA_MINUS);
		initDecayMode(72, 182, DECAY_MODE_BETA_MINUS);
		initDecayMode(72, 183, DECAY_MODE_BETA_MINUS);
		initDecayMode(72, 184, DECAY_MODE_BETA_MINUS);
		initDecayMode(72, 185, DECAY_MODE_BETA_MINUS);
		initDecayMode(72, 186, DECAY_MODE_BETA_MINUS);

		initAtomProperty(73, "Ta", "tantalum", "タンタル", 156, 188, 156, 188, 25.36);
		initIsotopeProperty(73, 156, 155.97230, 144.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 157, 156.96819, 10.1e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 158, 157.96670, 49.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 159, 158.963018, 1.04, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 160, 159.96149, 1.70, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 161, 160.95842, 3.0, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 162, 161.95729, 3.57, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 163, 162.95433, 10.6, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 164, 163.95353, 14.2, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 165, 164.950773, 31.0, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 166, 165.95051, 34.4, HLU_SECOND, 0.0);
		initIsotopeProperty(73, 167, 166.94809, 1.33, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 168, 167.94805, 2.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 169, 168.94601, 4.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 170, 169.94618, 6.76, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 171, 170.94448, 23.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 172, 171.94490, 36.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 173, 172.94375, 3.14, HLU_HOUR, 0.0);
		initIsotopeProperty(73, 174, 173.94445, 1.14, HLU_HOUR, 0.0);
		initIsotopeProperty(73, 175, 174.94374, 10.5, HLU_HOUR, 0.0);
		initIsotopeProperty(73, 176, 175.94486, 8.09, HLU_HOUR, 0.0);
		initIsotopeProperty(73, 177, 176.944472, 56.56, HLU_HOUR, 0.0);
		initIsotopeProperty(73, 178, 177.945778, 9.31, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 179, 178.9459295, 1.82, HLU_YEAR, 0.0);
		//initIsotopeProperty(73, 180, 179.9474648, 8.152, HLU_HOUR, 0.0);//180Ta
		//fprintf(stderr, "DEBUG:%s(%d) e_coefMeVtoMassU %lg\n", __FUNCTION__, __LINE__, e_coefMeVtoMassU);
		if(e_coefMeVtoMassU < 0.0){
			fprintf(stderr, "ERROR:%s(%d) wrong e_coefMeVtoMassU %lg\n", __FUNCTION__, __LINE__, e_coefMeVtoMassU);
			exit(0);
		}
		initIsotopeProperty(73, 180, 179.9474648 + 77.1e-3 * e_coefMeVtoMassU, 1.2e-4, HLU_STABLE, 0.0);//180m1Ta (77.1[keV])
		initIsotopeProperty(73, 181, 180.9479958, 0.0, HLU_STABLE, 0.99988);
		initIsotopeProperty(73, 182, 181.9501518, 114.43, HLU_DAY, 0.0);
		initIsotopeProperty(73, 183, 182.9513726, 5.1, HLU_DAY, 0.0);
		initIsotopeProperty(73, 184, 183.954008, 8.7, HLU_HOUR, 0.0);
		initIsotopeProperty(73, 185, 184.955559, 49.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 186, 185.95855, 10.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 187, 186.96053, 2.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(73, 188, 187.96370, 20.0, HLU_SECOND, 0.0);
		//initIsotopeProperty(73, 189, 188.96583, 3.0, HLU_SECOND, 0.0);
		//initIsotopeProperty(73, 190, 189.96923, 0.3, HLU_SECOND, 0.0);
		initDecayMode2(73, 156, DECAY_MODE_BETA_PLUS, 95.8, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(73, 157, DECAY_MODE_ALPHA, 91.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(73, 158, DECAY_MODE_ALPHA, 96.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(73, 159, DECAY_MODE_BETA_PLUS, 66.0, DECAY_MODE_ALPHA);
		initDecayMode2(73, 160, DECAY_MODE_ALPHA, 30.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(73, 161, DECAY_MODE_BETA_PLUS, 95.0, DECAY_MODE_ALPHA);
		initDecayMode2(73, 162, DECAY_MODE_BETA_PLUS, 99.92, DECAY_MODE_ALPHA);
		initDecayMode2(73, 163, DECAY_MODE_BETA_PLUS, 99.8, DECAY_MODE_ALPHA);
		initDecayMode(73, 164, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 165, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 166, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 167, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 168, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 169, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 170, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 171, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 172, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 173, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 174, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 175, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 176, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 177, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 178, DECAY_MODE_BETA_PLUS);
		initDecayMode(73, 179, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(73, 182, DECAY_MODE_BETA_MINUS);
		initDecayMode(73, 183, DECAY_MODE_BETA_MINUS);
		initDecayMode(73, 184, DECAY_MODE_BETA_MINUS);
		initDecayMode(73, 185, DECAY_MODE_BETA_MINUS);
		initDecayMode(73, 186, DECAY_MODE_BETA_MINUS);
		initDecayMode(73, 187, DECAY_MODE_BETA_MINUS);
		initDecayMode(73, 188, DECAY_MODE_BETA_MINUS);

		initAtomProperty(74, "W", "tungsten", "タングステン", 158, 190, 158, 190, 24.27);
		initIsotopeProperty(74, 158, 157.97456, 1.37e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 159, 158.97292, 8.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 160, 159.96848, 90.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 161, 160.96736, 409.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 162, 161.963497, 1.36, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 163, 162.96252, 2.8, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 164, 163.958954, 6.3, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 165, 164.958280, 5.1, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 166, 165.955027, 19.2, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 167, 166.954816, 19.9, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 168, 167.951808, 51.0, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 169, 168.951779, 76.0, HLU_SECOND, 0.0);
		initIsotopeProperty(74, 170, 169.949228, 2.42, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 171, 170.94945, 2.38, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 172, 171.94729, 6.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 173, 172.94769, 7.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 174, 173.94608, 33.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 175, 174.94672, 35.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 176, 175.94563, 2.5, HLU_HOUR, 0.0);
		initIsotopeProperty(74, 177, 176.94664, 132.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 178, 177.945876, 21.6, HLU_DAY, 0.0);
		initIsotopeProperty(74, 179, 178.947070, 37.05, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 180, 179.946704, 1.8e18, HLU_YEAR, 0.0012);
		initIsotopeProperty(74, 181, 180.948197, 121.2, HLU_DAY, 0.0);
		initIsotopeProperty(74, 182, 181.9482042, 0.0, HLU_STABLE, 0.265);
		initIsotopeProperty(74, 183, 182.9502230, 0.0, HLU_STABLE, 0.1431);
		initIsotopeProperty(74, 184, 183.9509312, 0.0, HLU_STABLE, 0.3064);
		initIsotopeProperty(74, 185, 184.9534193, 75.1, HLU_DAY, 0.0);
		initIsotopeProperty(74, 186, 185.9543641, 0.0, HLU_STABLE, 0.2843);
		initIsotopeProperty(74, 187, 186.9571605, 23.72, HLU_HOUR, 0.0);
		initIsotopeProperty(74, 188, 187.958489, 69.78, HLU_DAY, 0.0);
		initIsotopeProperty(74, 189, 188.96191, 11.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(74, 190, 189.96318, 30.0, HLU_MINUTE, 0.0);
		initDecayMode(74, 158, DECAY_MODE_ALPHA);
		initDecayMode2(74, 159, DECAY_MODE_ALPHA, 82.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(74, 160, DECAY_MODE_ALPHA, 87.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(74, 161, DECAY_MODE_ALPHA, 73.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(74, 162, DECAY_MODE_BETA_PLUS, 53.0, DECAY_MODE_ALPHA);
		initDecayMode2(74, 163, DECAY_MODE_BETA_PLUS, 59.0, DECAY_MODE_ALPHA);
		initDecayMode2(74, 164, DECAY_MODE_BETA_PLUS, 97.4, DECAY_MODE_ALPHA);
		initDecayMode2(74, 165, DECAY_MODE_BETA_PLUS, 99.8, DECAY_MODE_ALPHA);
		initDecayMode2(74, 166, DECAY_MODE_BETA_PLUS, 100.0 - 0.035, DECAY_MODE_ALPHA);
		initDecayMode2(74, 167, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_ALPHA);
		initDecayMode2(74, 168, DECAY_MODE_BETA_PLUS, 100.0 - 0.0319, DECAY_MODE_ALPHA);
		initDecayMode(74, 169, DECAY_MODE_BETA_PLUS);
		initDecayMode2(74, 170, DECAY_MODE_BETA_PLUS, 99.0, DECAY_MODE_ALPHA);
		initDecayMode(74, 171, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 172, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 173, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 174, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 175, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 176, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(74, 177, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 178, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(74, 179, DECAY_MODE_BETA_PLUS);
		initDecayMode(74, 180, DECAY_MODE_ALPHA);
		initDecayMode(74, 181, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(74, 185, DECAY_MODE_BETA_MINUS);
		initDecayMode(74, 187, DECAY_MODE_BETA_MINUS);
		initDecayMode(74, 188, DECAY_MODE_BETA_MINUS);
		initDecayMode(74, 189, DECAY_MODE_BETA_MINUS);
		initDecayMode(74, 190, DECAY_MODE_BETA_MINUS);

		initAtomProperty(75, "Re", "rhenium", "レニウム", 160, 192, 160, 192, 25.48);
		initIsotopeProperty(75, 160, 159.98212, 860.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 161, 160.97759, 0.37e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 162, 161.97600, 107.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 163, 162.972081, 390.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 164, 163.97032, 0.53, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 165, 164.967089, 1.0, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 166, 165.96581, 2.0, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 167, 166.96260, 3.4, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 168, 167.96157, 4.4, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 169, 168.95879, 8.1, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 170, 169.958220, 9.2, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 171, 170.95572, 15.2, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 172, 171.95542, 15.0, HLU_SECOND, 0.0);
		initIsotopeProperty(75, 173, 172.95324, 1.98, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 174, 173.95312, 2.40, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 175, 174.95138, 5.89, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 176, 175.95162, 5.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 177, 176.95033, 14.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 178, 177.95099, 13.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 179, 178.949988, 19.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 180, 179.950789, 2.44, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 181, 180.950068, 19.9, HLU_HOUR, 0.0);
		initIsotopeProperty(75, 182, 181.95121, 64.0, HLU_HOUR, 0.0);
		initIsotopeProperty(75, 183, 182.950820, 70.0, HLU_DAY, 0.0);
		initIsotopeProperty(75, 184, 183.952521, 38.0, HLU_DAY, 0.0);
		initIsotopeProperty(75, 185, 184.9529550, 0.0, HLU_STABLE, 0.374);
		initIsotopeProperty(75, 186, 185.9549861, 3.7186, HLU_DAY, 0.0);
		initIsotopeProperty(75, 187, 186.9557531, 41.2e9, HLU_YEAR, 0.626);
		initIsotopeProperty(75, 188, 187.9581144, 17.0040, HLU_HOUR, 0.0);
		initIsotopeProperty(75, 189, 188.959229, 24.3, HLU_HOUR, 0.0);
		initIsotopeProperty(75, 190, 189.96182, 3.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 191, 190.963125, 9.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(75, 192, 191.96596, 16.0, HLU_SECOND, 0.0);
		initDecayMode2(75, 160, DECAY_MODE_PROTON_EMISSION, 91.0, DECAY_MODE_ALPHA);
		initDecayMode(75, 161, DECAY_MODE_PROTON_EMISSION);
		initDecayMode2(75, 162, DECAY_MODE_ALPHA, 94.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(75, 163, DECAY_MODE_BETA_PLUS, 68.0, DECAY_MODE_ALPHA);
		initDecayMode2(75, 164, DECAY_MODE_ALPHA, 58.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(75, 165, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(75, 166, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(75, 167, DECAY_MODE_ALPHA, 50.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(75, 168, DECAY_MODE_BETA_PLUS, 100.0 - 0.005, DECAY_MODE_ALPHA);
		initDecayMode2(75, 169, DECAY_MODE_BETA_PLUS, 100.0 - 0.005, DECAY_MODE_ALPHA);
		initDecayMode2(75, 170, DECAY_MODE_BETA_PLUS, 99.99, DECAY_MODE_ALPHA);
		initDecayMode(75, 171, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 172, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 173, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 174, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 175, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 176, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 177, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 178, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 179, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 180, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 181, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 182, DECAY_MODE_BETA_PLUS);
		initDecayMode(75, 183, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(75, 184, DECAY_MODE_BETA_PLUS);
		initDecayMode2(75, 186, DECAY_MODE_BETA_MINUS, 93.1, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(75, 187, DECAY_MODE_BETA_MINUS, 100.0 - 1.0e-4, DECAY_MODE_ALPHA);
		initDecayMode(75, 188, DECAY_MODE_BETA_MINUS);
		initDecayMode(75, 189, DECAY_MODE_BETA_MINUS);
		initDecayMode(75, 190, DECAY_MODE_BETA_MINUS);
		initDecayMode(75, 191, DECAY_MODE_BETA_MINUS);
		initDecayMode(75, 192, DECAY_MODE_BETA_MINUS);

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
		initDecayMode3(81, 180, DECAY_MODE_ALPHA, 75.0 - 10.0e-4, DECAY_MODE_BETA_PLUS, 25.0, DECAY_MODE_SELF_FISSION_80KR);
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
		initDecayMode2(81,210, DECAY_MODE_BETA_MINUS, 99.991, DECAY_MODE_BETA_MINUS_AND_NEUTRON);

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
		initIsotopeProperty(82, 205, 204.9744818, 15.3e6, HLU_HOUR, 0.0);
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

		initAtomProperty(83, "Bi", "bismuth", "ビスマス", 185, 216, 185, 216, 25.52);
		initIsotopeProperty(83, 185, 184.99763, 2.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 186, 185.99660, 14.8e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 187, 186.993158, 32.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 188, 187.99227, 44.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 189, 188.98920, 674.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 190, 189.9883, 6.3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 191, 190.985786, 12.3, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 192, 191.98546, 34.6, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 193, 192.98296, 67.0, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 194, 193.98283, 95.0, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 195, 194.980651, 183.0, HLU_SECOND, 0.0);
		initIsotopeProperty(83, 196, 195.980667, 5.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 197, 196.978864, 9.33, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 198, 197.97921, 10.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 199, 198.977672, 27.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 200, 199.978132, 36.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 201, 200.977009, 108.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 202, 201.977742, 1.72, HLU_HOUR, 0.0);
		initIsotopeProperty(83, 203, 202.976876, 11.76, HLU_HOUR, 0.0);
		initIsotopeProperty(83, 204, 203.977813, 11.22, HLU_HOUR, 0.0);
		initIsotopeProperty(83, 205, 204.977389, 15.31, HLU_DAY, 0.0);
		initIsotopeProperty(83, 206, 205.978499, 6.243, HLU_DAY, 0.0);
		initIsotopeProperty(83, 207, 206.9784707, 32.9, HLU_YEAR, 0.0);
		initIsotopeProperty(83, 208, 207.9797422, 3.68e5, HLU_YEAR, 0.0);
		initIsotopeProperty(83, 209, 208.9803987, 1.9e19, HLU_YEAR, 1.0);
		initIsotopeProperty(83, 210, 209.9841204, 5.012, HLU_DAY, 0.0);
		initIsotopeProperty(83, 211, 210.987269, 2.14, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 212, 211.9912857, 60.55, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 213, 212.994385, 45.59, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 214, 213.998712, 19.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 215, 215.001770, 7.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(83, 216, 216.006306, 2.17, HLU_MINUTE, 0.0);
		initDecayMode2(83, 185, DECAY_MODE_PROTON_EMISSION, 99.999, DECAY_MODE_ALPHA);
		initDecayMode2(83, 186, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 187, DECAY_MODE_ALPHA, 50.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 188, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 189, DECAY_MODE_ALPHA, 51.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 190, DECAY_MODE_ALPHA, 77.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 191, DECAY_MODE_ALPHA, 60.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 192, DECAY_MODE_BETA_PLUS, 82.0, DECAY_MODE_ALPHA);
		initDecayMode2(83, 193, DECAY_MODE_BETA_PLUS, 95.0, DECAY_MODE_ALPHA);
		initDecayMode2(83, 194, DECAY_MODE_BETA_PLUS, 99.54, DECAY_MODE_ALPHA);
		initDecayMode2(83, 195, DECAY_MODE_BETA_PLUS, 99.97, DECAY_MODE_ALPHA);
		initDecayMode2(83, 196, DECAY_MODE_BETA_PLUS, 100.0 - 0.00115, DECAY_MODE_ALPHA);
		initDecayMode2(83, 197, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-4, DECAY_MODE_ALPHA);
		initDecayMode(83, 198, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 199, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 200, DECAY_MODE_BETA_PLUS);
		initDecayMode2(83, 201, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-4, DECAY_MODE_ALPHA);
		initDecayMode2(83, 202, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-5, DECAY_MODE_ALPHA);
		initDecayMode2(83, 203, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-5, DECAY_MODE_ALPHA);
		initDecayMode(83, 204, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 205, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 206, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 207, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 208, DECAY_MODE_BETA_PLUS);
		initDecayMode(83, 209, DECAY_MODE_ALPHA);
		initDecayMode2(83, 210, DECAY_MODE_BETA_MINUS, 100.0 - 1.32e-4, DECAY_MODE_ALPHA);
		initDecayMode2(83, 211, DECAY_MODE_ALPHA, 100.0 - 0.276, DECAY_MODE_BETA_MINUS);
		initDecayMode3(83, 212, DECAY_MODE_BETA_MINUS, 64.06 - 0.014, DECAY_MODE_ALPHA, 35.94, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode2(83, 213, DECAY_MODE_BETA_MINUS, 97.91, DECAY_MODE_ALPHA);
		initDecayMode3(83, 214, DECAY_MODE_BETA_MINUS, 100.0 - 0.021 - 0.003, DECAY_MODE_ALPHA, 0.021, DECAY_MODE_BETA_MINUS_AND_ALPHA);
		initDecayMode(83, 215, DECAY_MODE_BETA_MINUS);
		initDecayMode(83, 216, DECAY_MODE_BETA_MINUS);

		initAtomProperty(84, "Po", "polonium", "ポロニウム", 190, 218, 190, 218, 26.4);
		initIsotopeProperty(84, 190, 189.995101, 2.46e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 191, 190.994574, 22.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 192, 191.991335, 32.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 193, 192.99103, 420.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 194, 193.988186, 0.392, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 195, 194.98811, 4.64, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 196, 195.985535, 5.56, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 197, 196.98566, 53.6, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 198, 197.983389, 1.77, HLU_MINUTE, 0.0);
		initIsotopeProperty(84, 199, 198.983666, 5.48, HLU_MINUTE, 0.0);
		initIsotopeProperty(84, 200, 199.981799, 11.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(84, 201, 200.982260, 15.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(84, 202, 201.980758, 44.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(84, 203, 202.981420, 36.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(84, 204, 203.980318, 3.53, HLU_HOUR, 0.0);
		initIsotopeProperty(84, 205, 204.981203, 1.66, HLU_HOUR, 0.0);
		initIsotopeProperty(84, 206, 205.980481, 8.8, HLU_DAY, 0.0);
		initIsotopeProperty(84, 207, 206.981593, 5.80, HLU_HOUR, 0.0);
		initIsotopeProperty(84, 208, 207.9812457, 2.898, HLU_YEAR, 0.0);
		initIsotopeProperty(84, 209, 208.9824304, 125.2, HLU_YEAR, 0.0);
		initIsotopeProperty(84, 210, 209.9828737, 138.376, HLU_DAY, 0.0);
		initIsotopeProperty(84, 211, 210.9866532, 0.516, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 212, 211.9888680, 299e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 213, 212.992857, 3.65e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 214, 213.9952014, 164.3e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 215, 214.9994200, 1.781e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 216, 216.0019150, 0.145, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 217, 217.006335, 1.47, HLU_SECOND, 0.0);
		initIsotopeProperty(84, 218, 218.0089730, 3.10, HLU_MINUTE, 0.0);
		initDecayMode2(84, 190, DECAY_MODE_ALPHA, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 191, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 192, DECAY_MODE_ALPHA, 99.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 193, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 194, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 195, DECAY_MODE_ALPHA, 75.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 196, DECAY_MODE_ALPHA, 94.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 197, DECAY_MODE_BETA_PLUS, 54.0, DECAY_MODE_ALPHA);
		initDecayMode2(84, 198, DECAY_MODE_ALPHA, 57.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 199, DECAY_MODE_BETA_PLUS, 92.5, DECAY_MODE_ALPHA);
		initDecayMode2(84, 200, DECAY_MODE_BETA_PLUS, 88.8, DECAY_MODE_ALPHA);
		initDecayMode2(84, 201, DECAY_MODE_BETA_PLUS, 98.4, DECAY_MODE_ALPHA);
		initDecayMode2(84, 202, DECAY_MODE_BETA_PLUS, 98.0, DECAY_MODE_ALPHA);
		initDecayMode2(84, 203, DECAY_MODE_BETA_PLUS, 99.89, DECAY_MODE_ALPHA);
		initDecayMode2(84, 204, DECAY_MODE_BETA_PLUS, 99.33, DECAY_MODE_ALPHA);
		initDecayMode2(84, 205, DECAY_MODE_BETA_PLUS, 99.96, DECAY_MODE_ALPHA);
		initDecayMode2(84, 206, DECAY_MODE_BETA_PLUS, 94.55, DECAY_MODE_ALPHA);
		initDecayMode2(84, 207, DECAY_MODE_BETA_PLUS, 100.0 - 0.021, DECAY_MODE_ALPHA);
		initDecayMode2(84, 208, DECAY_MODE_ALPHA, 100.0 - 0.00277, DECAY_MODE_BETA_PLUS);
		initDecayMode2(84, 209, DECAY_MODE_ALPHA, 99.52, DECAY_MODE_BETA_PLUS);
		initDecayMode(84, 210, DECAY_MODE_ALPHA);
		initDecayMode(84, 211, DECAY_MODE_ALPHA);
		initDecayMode(84, 212, DECAY_MODE_ALPHA);
		initDecayMode(84, 213, DECAY_MODE_ALPHA);
		initDecayMode(84, 214, DECAY_MODE_ALPHA);
		initDecayMode2(84, 215, DECAY_MODE_ALPHA, 100.0 - 2.3e-4, DECAY_MODE_BETA_MINUS);
		initDecayMode2(84, 216, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode2(84, 217, DECAY_MODE_ALPHA, 95.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(84, 218, DECAY_MODE_ALPHA, 99.98, DECAY_MODE_BETA_MINUS);

		initAtomProperty(85, "At", "astatine", "アスタチン", 193, 222, 193, 222, 26.4);//unkown heat capacity
		initIsotopeProperty(85, 193, 192.99984, 28.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 194, 193.99873, 40.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 195, 194.996268, 328.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 196, 195.99579, 253.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 197, 196.99319, 0.390, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 198, 197.99284, 4.2, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 199, 198.99053, 6.92, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 200, 199.990351, 43.2, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 201, 200.988417, 85.0, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 202, 201.98863, 184.0, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 203, 202.986942, 7.37, HLU_MINUTE, 0.0);
		initIsotopeProperty(85, 204, 203.987251, 9.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(85, 205, 204.986074, 26.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(85, 206, 205.986667, 30.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(85, 207, 206.985784, 1.80, HLU_HOUR, 0.0);
		initIsotopeProperty(85, 208, 207.986590, 1.63, HLU_HOUR, 0.0);
		initIsotopeProperty(85, 209, 208.986173, 5.41, HLU_HOUR, 0.0);
		initIsotopeProperty(85, 210, 209.987148, 8.1, HLU_HOUR, 0.0);
		initIsotopeProperty(85, 211, 210.9874963, 7.214, HLU_HOUR, 0.0);
		initIsotopeProperty(85, 212, 211.990745, 0.314, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 213, 212.992937, 125.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 214, 213.996372, 558.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 215, 214.998653, 0.10e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 216, 216.002423, 0.30e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 217, 217.004719, 32.3e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 218, 218.008694, 1.5, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 219, 219.011162, 56.0, HLU_SECOND, 0.0);
		initIsotopeProperty(85, 220, 220.01541, 3.71, HLU_MINUTE, 0.0);
		initIsotopeProperty(85, 221, 221.01805, 2.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(85, 222, 222.02233, 54.0, HLU_SECOND, 0.0);
		initDecayMode(85, 193, DECAY_MODE_ALPHA);
		initDecayMode2(85, 194, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 195, DECAY_MODE_ALPHA, 75.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 196, DECAY_MODE_ALPHA, 96.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 197, DECAY_MODE_ALPHA, 96.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 198, DECAY_MODE_ALPHA, 94.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 199, DECAY_MODE_ALPHA, 89.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 200, DECAY_MODE_ALPHA, 57.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 201, DECAY_MODE_ALPHA, 71.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(85, 202, DECAY_MODE_BETA_PLUS, 88.0, DECAY_MODE_ALPHA);
		initDecayMode2(85, 203, DECAY_MODE_BETA_PLUS, 69.0, DECAY_MODE_ALPHA);
		initDecayMode2(85, 204, DECAY_MODE_BETA_PLUS, 96.2, DECAY_MODE_ALPHA);
		initDecayMode2(85, 205, DECAY_MODE_BETA_PLUS, 90.0, DECAY_MODE_ALPHA);
		initDecayMode2(85, 206, DECAY_MODE_BETA_PLUS, 99.11, DECAY_MODE_ALPHA);
		initDecayMode2(85, 207, DECAY_MODE_BETA_PLUS, 100.0 - 8.6, DECAY_MODE_ALPHA);
		initDecayMode2(85, 208, DECAY_MODE_BETA_PLUS, 100.0 - 0.55, DECAY_MODE_ALPHA);
		initDecayMode2(85, 209, DECAY_MODE_BETA_PLUS, 96.0, DECAY_MODE_ALPHA);
		initDecayMode2(85, 210, DECAY_MODE_BETA_PLUS, 100.0 - 0.18, DECAY_MODE_ALPHA);
		initDecayMode2(85, 211, DECAY_MODE_ELECTRON_CAPTURE, 58.2, DECAY_MODE_ALPHA);
		initDecayMode3(85, 212, DECAY_MODE_ALPHA, 99.95 - 2.0e-6, DECAY_MODE_BETA_PLUS, 0.05, DECAY_MODE_BETA_MINUS);
		initDecayMode(85, 213, DECAY_MODE_ALPHA);
		initDecayMode(85, 214, DECAY_MODE_ALPHA);
		initDecayMode(85, 215, DECAY_MODE_ALPHA);
		initDecayMode3(85, 216, DECAY_MODE_ALPHA, 100.0 - 0.006 - 3.0e-7, DECAY_MODE_BETA_MINUS, 0.006, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(85, 217, DECAY_MODE_ALPHA, 99.98, DECAY_MODE_BETA_MINUS);
		initDecayMode2(85, 218, DECAY_MODE_ALPHA, 99.9, DECAY_MODE_BETA_MINUS);
		initDecayMode2(85, 219, DECAY_MODE_ALPHA, 97.0, DECAY_MODE_BETA_MINUS);
		initDecayMode2(85, 220, DECAY_MODE_BETA_MINUS, 92.0, DECAY_MODE_ALPHA);
		initDecayMode(85, 221, DECAY_MODE_BETA_MINUS);
		initDecayMode(85, 222, DECAY_MODE_BETA_MINUS);

		initAtomProperty(86, "Rn", "radon", "ラドン", 196, 228, 196, 228, 20.786);
		initIsotopeProperty(86, 196, 196.002115, 4.7e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 197, 197.00158, 66.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 198, 197.998679, 65.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 199, 198.99837, 620.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 200, 199.995699, 0.96, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 201, 200.99563, 7.0, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 202, 201.993263, 9.94, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 203, 202.993387, 44.2, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 204, 203.991429, 1.17, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 205, 204.99172, 170.0, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 206, 205.990214, 5.67, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 207, 206.990734, 9.25, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 208, 207.989642, 24.35, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 209, 208.990415, 28.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 210, 209.989696, 2.4, HLU_HOUR, 0.0);
		initIsotopeProperty(86, 211, 210.990601, 14.6, HLU_HOUR, 0.0);
		initIsotopeProperty(86, 212, 211.990704, 23.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 213, 212.993883, 19.5e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 214, 213.995363, 0.27e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 215, 214.998745, 2.30e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 216, 216.000274, 45.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 217, 217.003928, 0.54e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 218, 218.0056013, 35.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 219, 219.0094802, 3.96, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 220, 220.0113940, 55.6, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 221, 221.015537, 25.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 222, 222.0175777, 3.8235, HLU_DAY, 0.0);
		initIsotopeProperty(86, 223, 223.02179, 24.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 224, 224.02409, 107.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 225, 225.02844, 4.66, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 226, 226.03089, 7.4, HLU_MINUTE, 0.0);
		initIsotopeProperty(86, 227, 227.03541, 20.8, HLU_SECOND, 0.0);
		initIsotopeProperty(86, 228, 228.03799, 65.0, HLU_SECOND, 0.0);
		initDecayMode2(86, 196, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 197, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 198, DECAY_MODE_ALPHA, 99.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 199, DECAY_MODE_ALPHA, 94.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 200, DECAY_MODE_ALPHA, 98.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 201, DECAY_MODE_ALPHA, 80.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 202, DECAY_MODE_ALPHA, 85.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 203, DECAY_MODE_ALPHA, 66.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 204, DECAY_MODE_ALPHA, 73.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 205, DECAY_MODE_BETA_PLUS, 77.0, DECAY_MODE_ALPHA);
		initDecayMode2(86, 206, DECAY_MODE_ALPHA, 62.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 207, DECAY_MODE_BETA_PLUS, 79.0, DECAY_MODE_ALPHA);
		initDecayMode2(86, 208, DECAY_MODE_ALPHA, 62.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 209, DECAY_MODE_BETA_PLUS, 83.0, DECAY_MODE_ALPHA);
		initDecayMode2(86, 210, DECAY_MODE_ALPHA, 96.0,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 211, DECAY_MODE_ALPHA, 72.6,  DECAY_MODE_BETA_PLUS);
		initDecayMode2(86, 212, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode(86, 213, DECAY_MODE_ALPHA);
		initDecayMode2(86, 214, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode(86, 215, DECAY_MODE_ALPHA);
		initDecayMode(86, 216, DECAY_MODE_ALPHA);
		initDecayMode(86, 217, DECAY_MODE_ALPHA);
		initDecayMode(86, 218, DECAY_MODE_ALPHA);
		initDecayMode(86, 219, DECAY_MODE_ALPHA);
		initDecayMode2(86, 220, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_MINUS);
		initDecayMode2(86, 221, DECAY_MODE_BETA_MINUS, 78.0, DECAY_MODE_ALPHA);
		initDecayMode(86, 222, DECAY_MODE_ALPHA);
		initDecayMode(86, 223, DECAY_MODE_BETA_MINUS);
		initDecayMode(86, 224, DECAY_MODE_BETA_MINUS);
		initDecayMode(86, 225, DECAY_MODE_BETA_MINUS);
		initDecayMode(86, 226, DECAY_MODE_BETA_MINUS);
		initDecayMode(86, 227, DECAY_MODE_BETA_MINUS);
		initDecayMode(86, 228, DECAY_MODE_BETA_MINUS);

		initAtomProperty(87, "Fr", "francium", "フランシウム", 200, 232, 200, 232, 20.786);//heat capacity is unknown
		initIsotopeProperty(87, 200, 200.00657, 24.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 201, 201.00386, 67.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 202, 202.00337, 290.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 203, 203.000925, 0.55, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 204, 204.000653, 1.7, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 205, 204.998594, 3.80, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 206, 205.99867, 16.0, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 207, 206.99695, 14.8, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 208, 207.99714, 59.1, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 209, 208.995954, 50.0, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 210, 209.996408, 3.18, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 211, 210.995537, 3.10, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 212, 211.996202, 20.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 213, 212.996189, 34.6, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 214, 213.998971, 5.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 215, 215.000341, 86.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 216, 216.003198, 0.70e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 217, 217.004632, 16.8e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 218, 218.007578, 1.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 219, 219.009252, 20.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 220, 220.012327, 27.4, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 221, 221.014255, 4.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 222, 222.017552, 14.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 223, 223.0197359, 22.00, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 224, 224.02325, 3.33, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 225, 225.02557, 4.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 226, 226.02939, 49.0, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 227, 227.03184, 2.47, HLU_MINUTE, 0.0);
		initIsotopeProperty(87, 228, 228.03573, 38.0, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 229, 229.03845, 50.2, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 230, 230.04251, 19.1, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 231, 231.04544, 17.6, HLU_SECOND, 0.0);
		initIsotopeProperty(87, 232, 232.04977, 5.0, HLU_SECOND, 0.0);
		initDecayMode(87, 200, DECAY_MODE_ALPHA);
		initDecayMode2(87, 201, DECAY_MODE_ALPHA, 99.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 202, DECAY_MODE_ALPHA, 97.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 203, DECAY_MODE_ALPHA, 95.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 204, DECAY_MODE_ALPHA, 96.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 205, DECAY_MODE_ALPHA, 99.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 206, DECAY_MODE_BETA_PLUS, 58.0, DECAY_MODE_ALPHA);
		initDecayMode2(87, 207, DECAY_MODE_ALPHA, 95.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 208, DECAY_MODE_ALPHA, 90.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 209, DECAY_MODE_ALPHA, 89.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 210, DECAY_MODE_ALPHA, 60.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 211, DECAY_MODE_ALPHA, 80.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(87, 212, DECAY_MODE_BETA_PLUS, 57.0, DECAY_MODE_ALPHA);
		initDecayMode2(87, 213, DECAY_MODE_ALPHA, 99.45, DECAY_MODE_BETA_PLUS);
		initDecayMode(87, 214, DECAY_MODE_ALPHA);
		initDecayMode(87, 215, DECAY_MODE_ALPHA);
		initDecayMode2(87, 216, DECAY_MODE_ALPHA, 100.0 - 2.0e-7, DECAY_MODE_BETA_PLUS);
		initDecayMode(87, 217, DECAY_MODE_ALPHA);
		initDecayMode(87, 218, DECAY_MODE_ALPHA);
		initDecayMode(87, 219, DECAY_MODE_ALPHA);
		initDecayMode2(87, 220, DECAY_MODE_ALPHA, 99.65, DECAY_MODE_BETA_MINUS);
		initDecayMode3(87, 221, DECAY_MODE_ALPHA, 99.9 - 8.79e-11, DECAY_MODE_BETA_MINUS, 0.1, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode(87, 222, DECAY_MODE_BETA_MINUS);
		initDecayMode2(87, 223, DECAY_MODE_BETA_MINUS, 100.0 - 0.006, DECAY_MODE_ALPHA);
		initDecayMode(87, 224, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 225, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 226, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 227, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 228, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 229, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 230, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 231, DECAY_MODE_BETA_MINUS);
		initDecayMode(87, 232, DECAY_MODE_BETA_MINUS);

		initAtomProperty(88, "Ra", "radium", "ラジウム", 203, 234, 203, 234, 20.786);//heat capacity is unknown
		initIsotopeProperty(88, 203, 203.00927, 4.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 204, 204.006500, 60.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 205, 205.00627, 220.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 206, 206.003827, 0.24, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 207, 207.00380, 1.3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 208, 208.001840, 1.3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 209, 209.00199, 4.6, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 210, 210.000495, 3.7, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 211, 211.000898, 13.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 212, 211.999794, 13.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 213, 213.000384, 2.74, HLU_MINUTE, 0.0);
		initIsotopeProperty(88, 214, 214.000108, 2.46, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 215, 215.002720, 1.55e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 216, 216.003533, 182.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 217, 217.006320, 1.63e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 218, 218.007140, 25.2e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 219, 219.010085, 10.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 220, 220.011028, 17.9e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 221, 221.013917, 28.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 222, 222.015375, 38.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 223, 223.0185022, 11.43, HLU_DAY, 0.0);
		initIsotopeProperty(88, 224, 224.0202118, 3.6319, HLU_DAY, 0.0);
		initIsotopeProperty(88, 225, 225.023612, 14.9, HLU_DAY, 0.0);
		initIsotopeProperty(88, 226, 226.0254098, 1600.0, HLU_YEAR, 0.0);
		initIsotopeProperty(88, 227, 227.0291778, 42.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(88, 228, 228.0310703, 5.75, HLU_YEAR, 0.0);
		initIsotopeProperty(88, 229, 229.034958, 4.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(88, 230, 230.037056, 93.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(88, 231, 231.04122, 103.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 232, 232.04364, 250.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 233, 233.04806, 30.0, HLU_SECOND, 0.0);
		initIsotopeProperty(88, 234, 234.05070, 30.0, HLU_SECOND, 0.0);
		initDecayMode2(88, 203, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 204, DECAY_MODE_ALPHA, 99.7, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 205, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode(88, 206, DECAY_MODE_ALPHA);
		initDecayMode2(88, 207, DECAY_MODE_ALPHA, 90.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 208, DECAY_MODE_ALPHA, 95.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 209, DECAY_MODE_ALPHA, 90.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 210, DECAY_MODE_ALPHA, 96.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 211, DECAY_MODE_ALPHA, 97.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 212, DECAY_MODE_ALPHA, 85.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 213, DECAY_MODE_ALPHA, 80.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(88, 214, DECAY_MODE_ALPHA, 99.94, DECAY_MODE_BETA_PLUS);
		initDecayMode(88, 215, DECAY_MODE_ALPHA);
		initDecayMode2(88, 216, DECAY_MODE_ALPHA, 100.0 - 1.0e-8, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(88, 217, DECAY_MODE_ALPHA);
		initDecayMode2(88, 218, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode(88, 219, DECAY_MODE_ALPHA);
		initDecayMode(88, 220, DECAY_MODE_ALPHA);
		initDecayMode2(88, 221, DECAY_MODE_ALPHA, 100.0 - 1.2e-10, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode2(88, 222, DECAY_MODE_ALPHA, 100.0 - 3.0e-8, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode2(88, 223, DECAY_MODE_ALPHA, 100.0 - 6.4e-8, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode2(88, 224, DECAY_MODE_ALPHA, 100.0 - 4.3e-9, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode(88, 225, DECAY_MODE_BETA_MINUS);
		initDecayMode3(88, 226, DECAY_MODE_ALPHA, 99.999 - 2.6e-9, DECAY_MODE_DOUBLE_BETA_MINUS, 0.001, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode(88, 227, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 228, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 229, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 230, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 231, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 232, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 233, DECAY_MODE_BETA_MINUS);
		initDecayMode(88, 234, DECAY_MODE_BETA_MINUS);

		initAtomProperty(89, "Ac", "actinium", "アクチニウム", 207, 236, 207, 236, 27.2);
		initIsotopeProperty(89, 207, 207.01195, 31.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 208, 208.01155, 97.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 209, 209.00949, 92.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 210, 210.00944, 350.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 211, 211.00773, 213.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 212, 212.00781, 920.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 213, 213.00661, 731.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 214, 214.006902, 8.2, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 215, 215.006454, 0.17, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 216, 216.008720, 0.440e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 217, 217.009347, 69.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 218, 218.01164, 1.08e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 219, 219.01242, 11.8e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 220, 220.014763, 26.36e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 221, 221.01559, 52.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 222, 222.017844, 5.0, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 223, 223.019137, 2.10, HLU_MINUTE, 0.0);
		initIsotopeProperty(89, 224, 224.021723, 2.78, HLU_HOUR, 0.0);
		initIsotopeProperty(89, 225, 225.023230, 10.0, HLU_DAY, 0.0);
		initIsotopeProperty(89, 226, 226.026098, 29.37, HLU_HOUR, 0.0);
		initIsotopeProperty(89, 227, 227.0277521, 21.772, HLU_YEAR, 0.0);
		initIsotopeProperty(89, 228, 228.0310211, 6.13, HLU_HOUR, 0.0);
		initIsotopeProperty(89, 229, 229.03302, 62.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(89, 230, 230.03629, 122.0, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 231, 231.03856, 7.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(89, 232, 232.04203, 119.0, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 233, 233.04455, 145.0, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 234, 234.04842, 44.0, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 235, 235.05123, 40.0, HLU_SECOND, 0.0);
		initIsotopeProperty(89, 236, 236.05530, 2.0, HLU_MINUTE, 0.0);
		initDecayMode(89, 207, DECAY_MODE_ALPHA);
		initDecayMode2(89, 208, DECAY_MODE_ALPHA, 99.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 209, DECAY_MODE_ALPHA, 99.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 210, DECAY_MODE_ALPHA, 96.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 211, DECAY_MODE_ALPHA, 99.8, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 212, DECAY_MODE_ALPHA, 97.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 213, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 214, DECAY_MODE_ALPHA, 89.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 215, DECAY_MODE_ALPHA, 99.91, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 216, DECAY_MODE_ALPHA, 100.0 - 7.0e-5, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 217, DECAY_MODE_ALPHA, 98.0, DECAY_MODE_BETA_PLUS);
		initDecayMode(89, 218, DECAY_MODE_ALPHA);
		initDecayMode2(89, 219, DECAY_MODE_ALPHA, 100.0 - 1.0e-6, DECAY_MODE_BETA_PLUS);
		initDecayMode2(89, 220, DECAY_MODE_ALPHA, 100.0 - 5.0e-4, DECAY_MODE_BETA_PLUS);
		initDecayMode(89, 221, DECAY_MODE_ALPHA);
		initDecayMode2(89, 222, DECAY_MODE_ALPHA, 99.0, DECAY_MODE_BETA_PLUS);
		initDecayMode3(89, 223, DECAY_MODE_ALPHA, 99.0 - 3.2e-9, DECAY_MODE_ELECTRON_CAPTURE, 1.0, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode3(89, 224, DECAY_MODE_BETA_PLUS, 100.0 - 9.1 - 1.6, DECAY_MODE_ALPHA, 9.1, DECAY_MODE_BETA_MINUS);
		initDecayMode2(89, 225, DECAY_MODE_ALPHA, 100.0 - 6.0e-10, DECAY_MODE_CLUSTER_DECAY_14C);
		initDecayMode3(89, 226, DECAY_MODE_BETA_MINUS, 83.0 - 0.006, DECAY_MODE_ELECTRON_CAPTURE, 17.0, DECAY_MODE_ALPHA);
		initDecayMode2(89, 227, DECAY_MODE_BETA_MINUS, 98.61, DECAY_MODE_ALPHA);
		initDecayMode2(89, 228, DECAY_MODE_BETA_MINUS, 100.0 - 5.5e-6, DECAY_MODE_ALPHA);
		initDecayMode(89, 229, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 230, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 231, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 232, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 233, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 234, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 235, DECAY_MODE_BETA_MINUS);
		initDecayMode(89, 236, DECAY_MODE_BETA_MINUS);

		initAtomProperty(90, "Th", "thorium", "トリウム", 210, 238, 210, 238, 26.230);
		initIsotopeProperty(90, 210, 210.015075, 17.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 211, 211.01493, 48.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 212, 212.01298, 36.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 213, 213.01301, 140.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 214, 214.011500, 100.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 215, 215.011730, 1.2, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 216, 216.011062, 26.8e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 217, 217.013114, 240.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 218, 218.013284, 109.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 219, 219.01554, 1.05e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 220, 220.015748, 9.7e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 221, 221.018184, 1.73e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 222, 222.018468, 2.237e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 223, 223.020811, 0.60, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 224, 224.021467, 1.05, HLU_SECOND, 0.0);
		initIsotopeProperty(90, 225, 225.023951, 8.72, HLU_MINUTE, 0.0);
		initIsotopeProperty(90, 226, 226.024903, 30.57, HLU_MINUTE, 0.0);
		initIsotopeProperty(90, 227, 227.0277041, 18.68, HLU_DAY, 0.0);
		initIsotopeProperty(90, 228, 228.0287411, 1.9116, HLU_YEAR, 0.0);
		initIsotopeProperty(90, 229, 229.031762, 7.34e3, HLU_YEAR, 0.0);
		initIsotopeProperty(90, 230, 230.0331338, 7.538e4, HLU_YEAR, 0.0);
		initIsotopeProperty(90, 231, 231.0363043, 25.52, HLU_HOUR, 0.0);
		initIsotopeProperty(90, 232, 232.0380553, 1.405e10, HLU_YEAR, 1.0);
		initIsotopeProperty(90, 233, 233.0415818, 21.83, HLU_MINUTE, 0.0);
		initIsotopeProperty(90, 234, 234.043601, 24.10, HLU_DAY, 0.0);
		initIsotopeProperty(90, 235, 235.04751, 7.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(90, 236, 236.04987, 37.5, HLU_MINUTE, 0.0);
		initIsotopeProperty(90, 237, 237.05389, 4.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(90, 238, 238.0565, 9.4, HLU_MINUTE, 0.0);
		initDecayMode2(90, 210, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(90, 211, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(90, 212, DECAY_MODE_ALPHA, 99.7, DECAY_MODE_BETA_PLUS);
		initDecayMode2(90, 213, DECAY_MODE_ALPHA,  99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode(90, 214, DECAY_MODE_ALPHA);
		initDecayMode(90, 215, DECAY_MODE_ALPHA);
		initDecayMode2(90, 216, DECAY_MODE_ALPHA, 100.0 - 0.006, DECAY_MODE_BETA_PLUS);
		initDecayMode(90, 217, DECAY_MODE_ALPHA);
		initDecayMode(90, 218, DECAY_MODE_ALPHA);
		initDecayMode2(90, 219, DECAY_MODE_ALPHA, 100.0 - 1.0e-7, DECAY_MODE_BETA_PLUS);
		initDecayMode2(90, 220, DECAY_MODE_ALPHA, 100.0 - 2.0e-7, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(90, 221, DECAY_MODE_ALPHA);
		initDecayMode2(90, 222, DECAY_MODE_ALPHA, 100.0 - 1.3e-8, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(90, 223, DECAY_MODE_ALPHA);
		initDecayMode2(90, 224, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);
		initDecayMode2(90, 225, DECAY_MODE_ALPHA, 90.0, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(90, 226, DECAY_MODE_ALPHA);
		initDecayMode(90, 227, DECAY_MODE_ALPHA);
		initDecayMode(90, 228, DECAY_MODE_ALPHA);//ignore, DECAY_MODE_CLUSTER_DECAY_20O(1.3e-11)
		initDecayMode(90, 229, DECAY_MODE_ALPHA);
		initDecayMode(90, 230, DECAY_MODE_ALPHA);//ignore, DECAY_MODE_CLUSTER_DECAY_24Ne(5.6e-11), SELF_FISSION(various 5.0e-11)
		initDecayMode2(90, 231, DECAY_MODE_BETA_MINUS, 100.0 - 1.0e-8, DECAY_MODE_ALPHA);
		initDecayMode2(90, 232, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_MINUS);//ignore, SELF_FISSION(various 1.1e-9), DECAY_MODE_CLUSTER_DECAY_26Ne_24Ne(2.78e-10)
		initDecayMode(90, 233, DECAY_MODE_BETA_MINUS);
		initDecayMode(90, 234, DECAY_MODE_BETA_MINUS);
		initDecayMode(90, 235, DECAY_MODE_BETA_MINUS);
		initDecayMode(90, 236, DECAY_MODE_BETA_MINUS);
		initDecayMode(90, 237, DECAY_MODE_BETA_MINUS);
		initDecayMode(90, 238, DECAY_MODE_BETA_MINUS);

		initAtomProperty(91, "Pa", "protactinium", "プロトアクチニウム", 213, 240, 213, 240, 26.230);//heat capacity is unknown
		initIsotopeProperty(91, 213, 213.02111, 7.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 214, 214.02092, 17.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 215, 215.01919, 14.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 216, 216.01911, 105.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 217, 217.01832, 3.48e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 218, 218.020042, 0.113e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 219, 219.01988, 53.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 220, 220.02188, 780.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 221, 221.02188, 4.9e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 222, 222.02374, 3.2e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 223, 223.02396, 5.1e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 224, 224.025626, 844.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 225, 225.02613, 1.7, HLU_SECOND, 0.0);
		initIsotopeProperty(91, 226, 226.027948, 1.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(91, 227, 227.028805, 38.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(91, 228, 228.031051, 22.0, HLU_HOUR, 0.0);
		initIsotopeProperty(91, 229, 229.0320968, 1.50, HLU_DAY, 0.0);
		initIsotopeProperty(91, 230, 230.034541, 17.4, HLU_DAY, 0.0);
		initIsotopeProperty(91, 231, 231.0358840, 3.276e4, HLU_YEAR, 1.0000);
		initIsotopeProperty(91, 232, 232.038592, 1.31, HLU_DAY, 0.0);
		initIsotopeProperty(91, 233, 233.0402473, 26.975, HLU_DAY, 0.0);
		initIsotopeProperty(91, 234, 234.043308, 6.70, HLU_HOUR, 0.0);
		initIsotopeProperty(91, 235, 235.04544, 24.44, HLU_MINUTE, 0.0);
		initIsotopeProperty(91, 236, 236.04868, 9.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(91, 237, 237.05115, 8.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(91, 238, 238.05450, 2.27, HLU_MINUTE, 0.0);
		initIsotopeProperty(91, 239, 239.05726, 1.8, HLU_HOUR, 0.0);
		initIsotopeProperty(91, 240, 240.06098, 2.0, HLU_MINUTE, 0.0);
		initDecayMode(91, 213, DECAY_MODE_ALPHA);
		initDecayMode(91, 214, DECAY_MODE_ALPHA);
		initDecayMode(91, 215, DECAY_MODE_ALPHA);
		initDecayMode2(91, 216, DECAY_MODE_ALPHA, 80.0, DECAY_MODE_BETA_PLUS);
		initDecayMode(91, 217, DECAY_MODE_ALPHA);
		initDecayMode(91, 218, DECAY_MODE_ALPHA);
		initDecayMode2(91, 219, DECAY_MODE_ALPHA, 100.0 - 5.0e-9, DECAY_MODE_BETA_PLUS);
		initDecayMode(91, 220, DECAY_MODE_ALPHA);
		initDecayMode(91, 221, DECAY_MODE_ALPHA);
		initDecayMode(91, 222, DECAY_MODE_ALPHA);
		initDecayMode2(91, 223, DECAY_MODE_ALPHA, 100.0 - 0.001, DECAY_MODE_BETA_PLUS);
		initDecayMode2(91, 224, DECAY_MODE_ALPHA, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode(91, 225, DECAY_MODE_ALPHA);
		initDecayMode2(91, 226, DECAY_MODE_ALPHA, 74.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(91, 227, DECAY_MODE_ALPHA, 85.0, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(91, 228, DECAY_MODE_BETA_PLUS, 98.15, DECAY_MODE_ALPHA);
		initDecayMode2(91, 229, DECAY_MODE_ELECTRON_CAPTURE, 99.52, DECAY_MODE_ALPHA);
		initDecayMode3(91, 230, DECAY_MODE_BETA_PLUS, 91.6 - 0.00319, DECAY_MODE_BETA_MINUS, 8.4, DECAY_MODE_ALPHA);
		initDecayMode(91, 231, DECAY_MODE_ALPHA);//ignore, DECAY_MODE_CLUSTER_DECAY_24Ne(1.34e-9), SELF_FISSION(various, 3.0e-10), DECAY_MODE_CLUSTER_DECAY_23F(1.0e-12)
		initDecayMode2(91, 232, DECAY_MODE_BETA_MINUS, 100.0 - 0.003, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(91, 233, DECAY_MODE_BETA_MINUS);
		initDecayMode(91, 234, DECAY_MODE_BETA_MINUS);//ignore, SELF_FISSION(various, 3.0e-10)
		initDecayMode(91, 235, DECAY_MODE_BETA_MINUS);
		initDecayMode(91, 236, DECAY_MODE_BETA_MINUS);//ignore, BETA_MINUS_SELF_FISSION(various, 6.0e-8)
		initDecayMode(91, 237, DECAY_MODE_BETA_MINUS);
		initDecayMode(91, 238, DECAY_MODE_BETA_MINUS);//ignore, BETA_MINUS_SELF_FISSION(various, 2.6e-6)
		initDecayMode(91, 239, DECAY_MODE_BETA_MINUS);
		initDecayMode(91, 240, DECAY_MODE_BETA_MINUS);

		initAtomProperty(92, "U", "uranium", "ウラン", 217, 242, 217, 242, 27.665);
		initIsotopeProperty(92, 217, 217.02437, 26.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 218, 218.02354, 6.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 219, 219.02492, 55.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 220, 220.02472, 60.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 221, 221.02640, 700.0e-9, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 222, 222.02609, 1.4e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 223, 223.02774, 21.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 224, 224.027605, 940.0e-6, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 225, 225.02939, 61.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 226, 226.029339, 269.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(92, 227, 227.031156, 1.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(92, 228, 228.031374, 9.1, HLU_MINUTE, 0.0);
		initIsotopeProperty(92, 229, 229.033506, 58.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(92, 230, 230.033940, 20.8, HLU_DAY, 0.0);
		initIsotopeProperty(92, 231, 231.036294, 4.2, HLU_DAY, 0.0);
		initIsotopeProperty(92, 232, 232.0371562, 68.9, HLU_YEAR, 0.0);
		initIsotopeProperty(92, 233, 233.0396352, 1.59e5, HLU_YEAR, 0.0);
		initIsotopeProperty(92, 234, 234.0409521, 2.455e5, HLU_YEAR, 0.000054);
		initIsotopeProperty(92, 235, 235.0439299, 7.04e8, HLU_YEAR, 0.007204);
		initIsotopeProperty(92, 236, 236.045568, 2.342e7, HLU_YEAR, 0.0);
		initIsotopeProperty(92, 237, 237.0487302, 6.75, HLU_DAY, 0.0);
		initIsotopeProperty(92, 238, 238.0507882, 4.468e9, HLU_YEAR, 0.992742);
		initIsotopeProperty(92, 239, 239.0542933, 23.45, HLU_MINUTE, 0.0);
		initIsotopeProperty(92, 240, 240.056592, 14.1, HLU_HOUR, 0.0);
		initIsotopeProperty(92, 241, 241.06033, 5.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(92, 242, 242.06293, 16.8, HLU_MINUTE, 0.0);
		initDecayMode(92, 217, DECAY_MODE_ALPHA);
		initDecayMode(92, 218, DECAY_MODE_ALPHA);
		initDecayMode(92, 219, DECAY_MODE_ALPHA);
		initDecayMode2(92, 220, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(92, 221, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(92, 222, DECAY_MODE_ALPHA, 100.0 - 1.0e-6, DECAY_MODE_BETA_PLUS);
		initDecayMode(92, 223, DECAY_MODE_ALPHA);
		initDecayMode(92, 224, DECAY_MODE_ALPHA);
		initDecayMode(92, 225, DECAY_MODE_ALPHA);
		initDecayMode(92, 226, DECAY_MODE_ALPHA);
		initDecayMode2(92, 227, DECAY_MODE_ALPHA, 100.0 - 0.001, DECAY_MODE_BETA_PLUS);
		initDecayMode2(92, 228, DECAY_MODE_ALPHA, 95.0, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode2(92, 229, DECAY_MODE_BETA_PLUS, 80.0, DECAY_MODE_ALPHA);
		initDecayMode2(92, 230, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_DOUBLE_BETA_PLUS);//ignore, SELF_FISSION(various, 1.4e-10) 
		initDecayMode2(92, 231, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 0.004, DECAY_MODE_ALPHA);
		initDecayMode(92, 232, DECAY_MODE_ALPHA);//ignore, DECAY_MODE_CLUSTER_DECAY_24Ne(8.9e-10), DECAY_MODE_CLUSTER_DECAY_28Mg(5.0e-12), SELF_FISSION(various, 1.0e-11) 
		initDecayMode(92, 233, DECAY_MODE_ALPHA);//ignore SELF_FISSION(various, 6.0e-9), DECAY_MODE_CLUSTER_DECAY_24Ne(7.2e-11), DECAY_MODE_CLUSTER_DECAY_28Mg(1.3e-13)
		initDecayMode(92, 234, DECAY_MODE_ALPHA);//ignore SELF_FISSION(various, 1.73e-9), DECAY_MODE_CLUSTER_DECAY_28Mg(1.4e-11), DECAY_MODE_CLUSTER_DECAY_26Ne_24Ne(9.0e-12)
		initDecayMode(92, 235, DECAY_MODE_ALPHA);//ignore SELF_FISSION(various, 7.0e-9), DECAY_MODE_CLUSTER_DECAY_25Ne_24Ne(8.0e-10)
		initDecayMode(92, 236, DECAY_MODE_ALPHA);//ignore SELF_FISSION(various, 9.6e-8)
		initDecayMode(92, 237, DECAY_MODE_BETA_MINUS);
		initDecayMode2(92, 238, DECAY_MODE_ALPHA, 100.0 - 2.19e-10, DECAY_MODE_DOUBLE_BETA_MINUS);//ignore SELF_FISSION(various, 5.45e-5)
		initDecayMode(92, 239, DECAY_MODE_BETA_MINUS);
		initDecayMode2(92, 240, DECAY_MODE_BETA_MINUS, 100.0 - 1.0e-10, DECAY_MODE_ALPHA);
		initDecayMode(92, 241, DECAY_MODE_BETA_MINUS);
		initDecayMode(92, 242, DECAY_MODE_BETA_MINUS);

		initAtomProperty(93, "Np", "neptunium", "ネプツニウム", 225, 244, 225, 244, 29.46);
		initIsotopeProperty(93, 225, 225.03391, 3.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(93, 226, 226.03515, 35.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(93, 227, 227.03496, 510.0e-3, HLU_SECOND, 0.0);
		initIsotopeProperty(93, 228, 228.03618, 61.4, HLU_SECOND, 0.0);
		initIsotopeProperty(93, 229, 229.03626, 4.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 230, 230.03783, 4.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 231, 231.03825, 48.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 232, 232.04011, 14.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 233, 233.04074, 36.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 234, 234.042895, 4.4, HLU_DAY, 0.0);
		initIsotopeProperty(93, 235, 235.0440633, 396.1, HLU_DAY, 0.0);
		initIsotopeProperty(93, 236, 236.04657, 1.54e5, HLU_YEAR, 0.0);
		initIsotopeProperty(93, 237, 237.0481734, 2.144e6, HLU_YEAR, 0.0);
		initIsotopeProperty(93, 238, 238.0509464, 2.117, HLU_DAY, 0.0);
		initIsotopeProperty(93, 239, 239.0529390, 2.356, HLU_DAY, 0.0);
		initIsotopeProperty(93, 240, 240.056162, 61.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 241, 241.05825, 13.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 242, 242.06164, 2.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 243, 243.06428, 1.85, HLU_MINUTE, 0.0);
		initIsotopeProperty(93, 244, 244.06785, 2.29, HLU_MINUTE, 0.0);
		initDecayMode(93, 225, DECAY_MODE_ALPHA);
		initDecayMode(93, 226, DECAY_MODE_ALPHA);
		initDecayMode2(93, 227, DECAY_MODE_ALPHA, 99.95, DECAY_MODE_BETA_PLUS);
		initDecayMode2(93, 228, DECAY_MODE_BETA_PLUS, 59.0, DECAY_MODE_ALPHA);//ignore,  BETA_PLUS_SELF_FISSION(various, 0.012)
		initDecayMode2(93, 229, DECAY_MODE_ALPHA, 51.0, DECAY_MODE_BETA_PLUS);
		initDecayMode2(93, 230, DECAY_MODE_BETA_PLUS, 97.0, DECAY_MODE_ALPHA);
		initDecayMode2(93, 231, DECAY_MODE_BETA_PLUS, 98.0, DECAY_MODE_ALPHA);
		initDecayMode2(93, 232, DECAY_MODE_BETA_PLUS, 99.997, DECAY_MODE_ALPHA);
		initDecayMode2(93, 233, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_ALPHA);
		initDecayMode(93, 234, DECAY_MODE_BETA_PLUS);
		initDecayMode2(93, 235, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 0.0026, DECAY_MODE_ALPHA);
		initDecayMode3(93, 236, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 12.5 - 0.16, DECAY_MODE_BETA_MINUS, 12.5, DECAY_MODE_ALPHA);
		initDecayMode(93, 237, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 2.0e-10), DECAY_MODE_CLUSTER_DECAY_30Mg(4.0e-12)
		initDecayMode(93, 238, DECAY_MODE_BETA_MINUS);
		initDecayMode(93, 239, DECAY_MODE_BETA_MINUS);
		initDecayMode(93, 240, DECAY_MODE_BETA_MINUS);
		initDecayMode(93, 241, DECAY_MODE_BETA_MINUS);
		initDecayMode(93, 242, DECAY_MODE_BETA_MINUS);
		initDecayMode(93, 243, DECAY_MODE_BETA_MINUS);
		initDecayMode(93, 244, DECAY_MODE_BETA_MINUS);

		initAtomProperty(94, "Pu", "plutonium", "プルトニウム", 228, 247, 228, 247, 35.5);
		initIsotopeProperty(94, 228, 228.03874, 1.1, HLU_SECOND, 0.0);
		initIsotopeProperty(94, 229, 229.04015, 120.0, HLU_SECOND, 0.0);
		initIsotopeProperty(94, 230, 230.039650, 1.70, HLU_MINUTE, 0.0);
		initIsotopeProperty(94, 231, 231.041101, 8.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(94, 232, 232.041187, 33.7, HLU_MINUTE, 0.0);
		initIsotopeProperty(94, 233, 233.04300, 20.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(94, 234, 234.043317, 8.8, HLU_HOUR, 0.0);
		initIsotopeProperty(94, 235, 235.045286, 25.3, HLU_MINUTE, 0.0);
		initIsotopeProperty(94, 236, 236.0460580, 2.858, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 237, 237.0484097, 45.2, HLU_DAY, 0.0);
		initIsotopeProperty(94, 238, 238.0495599, 87.7, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 239, 239.0521634, 2.411e4, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 240, 240.0538135, 6.561e3, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 241, 241.0568515, 14.290, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 242, 242.0587426, 3.75e5, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 243, 243.062003, 4.956, HLU_HOUR, 0.0);
		initIsotopeProperty(94, 244, 244.064204, 8.00e7, HLU_YEAR, 0.0);
		initIsotopeProperty(94, 245, 245.067747, 10.5, HLU_HOUR, 0.0);
		initIsotopeProperty(94, 246, 246.070205, 10.84, HLU_DAY, 0.0);
		initIsotopeProperty(94, 247, 247.07407, 2.27, HLU_DAY, 0.0);
		initDecayMode2(94, 228, DECAY_MODE_ALPHA, 99.9, DECAY_MODE_BETA_PLUS);
		initDecayMode(94, 229, DECAY_MODE_ALPHA);
		initDecayMode2(94, 230, DECAY_MODE_ALPHA, 99.999, DECAY_MODE_BETA_PLUS);
		initDecayMode2(94, 231, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_ALPHA);
		initDecayMode2(94, 232, DECAY_MODE_ELECTRON_CAPTURE, 89.0, DECAY_MODE_ALPHA);
		initDecayMode2(94, 233, DECAY_MODE_BETA_PLUS, 99.88, DECAY_MODE_ALPHA);
		initDecayMode2(94, 234, DECAY_MODE_ELECTRON_CAPTURE, 94.0, DECAY_MODE_ALPHA);
		initDecayMode2(94, 235, DECAY_MODE_BETA_PLUS, 100.00 - 0.0027, DECAY_MODE_ALPHA);
		initDecayMode2(94, 236, DECAY_MODE_ALPHA, 99.999999, DECAY_MODE_DOUBLE_BETA_PLUS);//ignore, SELF_FISSION(various, 1.37e-7), DECAY_MODE_CLUSTER_DECAY_28Mg(2.0e-12)
		initDecayMode2(94, 237, DECAY_MODE_ELECTRON_CAPTURE, 100.0 - 0.0042, DECAY_MODE_ALPHA);
		initDecayMode(94, 238, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 1.9e-7), DECAY_MODE_CLUSTER_DECAY_32Si(1.4e-14), DECAY_MODE_CLUSTER_DECAY_30Mg_28Mg(6.0e-15)
		initDecayMode(94, 239, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 3.1e-10)
		initDecayMode(94, 240, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 5.7e-6), DECAY_MODE_CLUSTER_DECAY_34Si(1.3e-13)
		initDecayMode2(94, 241, DECAY_MODE_BETA_MINUS, 100.0 - 0.00245, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 2.4e-14)
		initDecayMode(94, 242, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 5.5e-4)
		initDecayMode(94, 243, DECAY_MODE_BETA_MINUS);
		initDecayMode2(94, 244, DECAY_MODE_ALPHA, 100.0 - 7.3e-9, DECAY_MODE_DOUBLE_BETA_MINUS);//ignore, SELF_FISSION(various, 0.123)
		initDecayMode(94, 245, DECAY_MODE_BETA_MINUS);
		initDecayMode(94, 246, DECAY_MODE_BETA_MINUS);
		initDecayMode(94, 247, DECAY_MODE_BETA_MINUS);

		initAtomProperty(95, "Am", "americium", "アメリシウム", 231, 249, 231, 249, 62.7);
		initIsotopeProperty(95, 231, 231.04556, 30.0, HLU_SECOND, 0.0);
		initIsotopeProperty(95, 232, 232.04659, 79.0, HLU_SECOND, 0.0);
		initIsotopeProperty(95, 233, 233.04635, 3.2, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 234, 234.04781, 2.32, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 235, 235.04795, 9.9, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 236, 236.04958, 3.6, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 237, 237.05000, 73.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 238, 238.05198, 98.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 239, 239.0530245, 11.9, HLU_HOUR, 0.0);
		initIsotopeProperty(95, 240, 240.055300, 50.8, HLU_HOUR, 0.0);
		initIsotopeProperty(95, 241, 241.0568291, 432.2, HLU_YEAR, 0.0);
		initIsotopeProperty(95, 242, 242.0595492, 16.02, HLU_HOUR, 0.0);
		initIsotopeProperty(95, 243, 243.0613811, 7370.0, HLU_YEAR, 0.0);
		initIsotopeProperty(95, 244, 244.0642848, 10.1, HLU_HOUR, 0.0);
		initIsotopeProperty(95, 245, 245.066452, 2.05, HLU_HOUR, 0.0);
		initIsotopeProperty(95, 246, 246.069775, 39.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 247, 247.07209, 23.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 248, 248.07575, 3.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(95, 249, 249.07848, 1.0, HLU_MINUTE, 0.0);
		initDecayMode2(95, 231, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_ALPHA);
		initDecayMode2(95, 232, DECAY_MODE_BETA_PLUS, 98.0, DECAY_MODE_ALPHA);//ignore, BETA_PLUS_SELF_FISSION(various, 0.069)
		initDecayMode2(95, 233, DECAY_MODE_BETA_PLUS, 98.0, DECAY_MODE_ALPHA);
		initDecayMode2(95, 234, DECAY_MODE_BETA_PLUS, 100.0 - 0.04, DECAY_MODE_ALPHA);//ignore, BETA_PLUS_SELF_FISSION(various, 0.0066)
		initDecayMode2(95, 235, DECAY_MODE_BETA_PLUS, 99.999, DECAY_MODE_ALPHA);
		initDecayMode2(95, 236, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_ALPHA);
		initDecayMode2(95, 237, DECAY_MODE_BETA_PLUS, 100.0 - 0.025, DECAY_MODE_ALPHA);
		initDecayMode2(95, 238, DECAY_MODE_BETA_PLUS, 100.0 - 1.0e-4, DECAY_MODE_ALPHA);
		initDecayMode2(95, 239, DECAY_MODE_ELECTRON_CAPTURE, 99.99, DECAY_MODE_ALPHA);
		initDecayMode2(95, 240, DECAY_MODE_BETA_PLUS, 100.0 - 1.9e-4, DECAY_MODE_ALPHA);
		initDecayMode(95, 241, DECAY_MODE_ALPHA);//ignore, DECAY_MODE_CLUSTER_DECAY_34Si(7.4e-10), SELF_FISSION(various, 4.3e-10)
		initDecayMode2(95, 242, DECAY_MODE_BETA_MINUS, 82.7, DECAY_MODE_ELECTRON_CAPTURE);
		initDecayMode(95, 243, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 3.7e-9)
		initDecayMode(95, 244, DECAY_MODE_BETA_MINUS);
		initDecayMode(95, 245, DECAY_MODE_BETA_MINUS);
		initDecayMode(95, 246, DECAY_MODE_BETA_MINUS);
		initDecayMode(95, 247, DECAY_MODE_BETA_MINUS);
		initDecayMode(95, 248, DECAY_MODE_BETA_MINUS);
		initDecayMode(95, 249, DECAY_MODE_BETA_MINUS);

		initAtomProperty(96, "Cm", "curium", "キュリウム", 233, 252, 233, 252, 62.7);//heat apacity is unknown
		initIsotopeProperty(96, 233, 233.05077, 1.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(96, 234, 234.05016, 51.0, HLU_SECOND, 0.0);
		initIsotopeProperty(96, 235, 235.05143, 5.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(96, 236, 236.05141, 10.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(96, 237, 237.05290, 20.0, HLU_MINUTE, 0.0);
		initIsotopeProperty(96, 238, 238.05303, 2.4, HLU_HOUR, 0.0);
		initIsotopeProperty(96, 239, 239.05496, 2.9, HLU_HOUR, 0.0);
		initIsotopeProperty(96, 240, 240.0555295, 27.0, HLU_DAY, 0.0);
		initIsotopeProperty(96, 241, 241.0576530, 32.8, HLU_DAY, 0.0);
		initIsotopeProperty(96, 242, 242.0588358, 162.8, HLU_DAY, 0.0);
		initIsotopeProperty(96, 243, 243.0613891, 29.1, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 244, 244.0627526, 18.10, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 245, 245.0654912, 8.5e3, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 246, 246.0672237, 4.76e3, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 247, 247.070354, 1.56e7, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 248, 248.072349, 3.48e5, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 249, 249.075953, 64.15, HLU_MINUTE, 0.0);
		initIsotopeProperty(96, 250, 250.078357, 8300.0, HLU_YEAR, 0.0);
		initIsotopeProperty(96, 251, 251.082285, 16.8, HLU_MINUTE, 0.0);
		initIsotopeProperty(96, 252, 252.08487, 1.0, HLU_DAY, 0.0);
		initDecayMode2(96, 233, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 234, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 235, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 236, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 237, DECAY_MODE_BETA_PLUS, 50.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 238, DECAY_MODE_ELECTRON_CAPTURE, 90.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 239, DECAY_MODE_BETA_PLUS, 99.9, DECAY_MODE_ALPHA);
		initDecayMode2(96, 240, DECAY_MODE_ALPHA, 99.5, DECAY_MODE_ELECTRON_CAPTURE);//ignore, SELF_FISSION(various, 3.9e-6)
		initDecayMode2(96, 241, DECAY_MODE_ELECTRON_CAPTURE, 99.0, DECAY_MODE_ALPHA);
		initDecayMode2(96, 242, DECAY_MODE_ALPHA, 99.999999, DECAY_MODE_DOUBLE_BETA_PLUS);//ignore, SELF_FISSION(various, 6.33e-6), DECAY_MODE_CLUSTER_DECAY_34Si(1.0e-14)
		initDecayMode2(96, 243, DECAY_MODE_ALPHA, 99.71, DECAY_MODE_ELECTRON_CAPTURE);//ignore, SELF_FISSION(various, 5.3e-9)
		initDecayMode(96, 244, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 1.34e-4)
		initDecayMode(96, 245, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 6.1e-7)
		initDecayMode(96, 246, DECAY_MODE_ALPHA);//ignore, SELF_FISSION(various, 0.0261)
		initDecayMode(96, 247, DECAY_MODE_ALPHA);
		initDecayMode2(96, 248, DECAY_MODE_ALPHA, 99.9999, DECAY_MODE_DOUBLE_BETA_MINUS);//ignore, SELF_FISSION(various, 8.26)
		initDecayMode(96, 249, DECAY_MODE_BETA_MINUS);
		initDecayMode2(96, 250, DECAY_MODE_ALPHA, 55.0, DECAY_MODE_BETA_MINUS);//ignore, SELF_FISSION(various, 80.0)
		initDecayMode(96, 251, DECAY_MODE_BETA_MINUS);
		initDecayMode(96, 252, DECAY_MODE_BETA_MINUS);

		e_hydrogenPtr = getIsotopePropertyPtr(ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN);
		if(!e_hydrogenPtr){
			fprintf(stderr, "FATAL ERROR:%s:NO hydrogen\n", __FUNCTION__);
		}
		e_helium4Ptr = getIsotopePropertyPtr(ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM);
		if(!e_helium4Ptr){
			fprintf(stderr, "FATAL ERROR:%s:NO helium4\n", __FUNCTION__);
		}
		e_C14Ptr = getIsotopePropertyPtr(ATOMICNUMBER_CARBON, MASSNUMBER_CARBON14);
		if(!e_C14Ptr){
			fprintf(stderr, "FATAL ERROR:%s:NO C14\n", __FUNCTION__);
		}
		e_Kr80Ptr = getIsotopePropertyPtr(ATOMICNUMBER_Kr, MASSNUMBER_Kr80);
		if(!e_Kr80Ptr){
			fprintf(stderr, "FATAL ERROR:%s:NO Kr80\n", __FUNCTION__);
		}
	}
}
extern void printAllElectronCapture(FILE * a_fp)
{
	int atomicNumber;
	fputs("[Electron Capture Properties]\n", a_fp);
	fputs("\n", a_fp);
	fputs("AN : \"atomic number\"\n", a_fp);
	fputs("MN : \"mass number\"\n", a_fp);
	fputs("HLs : \"Half Lift time in second\"\n", a_fp);
	fputs("dRt : \"decay Rate of Electron Capture\"\n", a_fp);
	fputs("masU : \"mass [U]\"\n", a_fp);
	fputs("mas : \"mass [MeV]\"\n", a_fp);
	fputs("mDf : \"mass Defect in [MeV]\"\n", a_fp);
	fputs("r1 : \"Electron orbit radius of k shell [m]\"\n", a_fp);
	fputs("CP1 : \"Electrostatic potential energy of k shell [MeV]\"\n", a_fp);
	fputs("E1 : \"Ionization Energy of k shell [MeV]\"\n", a_fp);
	fputs("r0 : \"the radius of a nucleus [m]\"\n", a_fp);
	fputs("CP0 : \"Electrostatic potential energy of the radius of a nucleus [MeV]\"\n", a_fp);
	fputs("\n", a_fp);
	fputs("AN symbol name MN HLs dRt masU mas mDf r1 CP1 E1 CP1-E1 r0 CP0\n", a_fp);
	for(atomicNumber = MIN_ATOM_NUMBER; atomicNumber <= MAX_ATOM_NUMBER; ++atomicNumber){
		struct atomProperty * atomPropertyPtr;
		atomPropertyPtr = getAtomPropertyPtr(atomicNumber);
		if(atomPropertyPtr){
			if(atomPropertyPtr->atomicNumber != UNDEF_ATOMIC_NUMBER){
				int massNumber;
				for(massNumber = atomPropertyPtr->minMassNumber; massNumber <= atomPropertyPtr->maxMassNumber; ++massNumber){
					struct isotopeProperty * isotopePropertyPtr
					= getIsotopePropertyPtr(atomPropertyPtr->atomicNumber, massNumber);
					if(atomPropertyPtr->atomicNumber == isotopePropertyPtr->atomicNumber){
						if(isotopePropertyPtr->massU > 0.0){
							if(isotopePropertyPtr->decayModeSize > 0){
								int i;
								for(i = 0; i < isotopePropertyPtr->decayModeSize; ++i){
									if(isotopePropertyPtr->decayMode[i] == DECAY_MODE_ELECTRON_CAPTURE){
										double massDefectMeV, r1, CP1, E1, r0, CP0;
										int daughterAtomicNumber, daughterMassNumber;
										struct isotopeProperty * daughterIsotopePropertyPtr;
										daughterAtomicNumber = isotopePropertyPtr->atomicNumber - 1;
										daughterMassNumber = isotopePropertyPtr->massNumber;
										daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
										massDefectMeV = isotopePropertyPtr->massMeV - daughterIsotopePropertyPtr->massMeV;

										fprintf(a_fp, "%d", atomPropertyPtr->atomicNumber);
										fprintf(a_fp, " \"%s\"", atomPropertyPtr->symbol);
										fprintf(a_fp, " \"%s\"", atomPropertyPtr->name);
										//fprintf(a_fp, " \"%s\"", atomPropertyPtr->japaneseName);
										fprintf(a_fp, " %d", isotopePropertyPtr->massNumber);
										//fprintf(a_fp, " \"%lg [%s]\"", isotopePropertyPtr->halfLife, getHalfLifeUnitString(isotopePropertyPtr->halfLifeUnit));
										fprintf(a_fp, " %lg", calcHalfLifeSec(isotopePropertyPtr));
									    //fprintf(a_fp, " %s", getDecayModeText(isotopePropertyPtr->decayMode[i]));
									    fprintf(a_fp, " %lg", isotopePropertyPtr->decayModeRate[i]);
										fprintf(a_fp, " %lg", isotopePropertyPtr->massU);
										fprintf(a_fp, " %lg", isotopePropertyPtr->massMeV);
									    fprintf(a_fp, " %lg", massDefectMeV);

										r1 = 5.288e-11 / atomPropertyPtr->atomicNumber;
										CP1 = (e_elementaryCharge * e_elementaryCharge * atomPropertyPtr->atomicNumber)/(4.0* M_PI *e_vacuumPermittivity * r1) / e_energyJouleOfeV / 1.0E6;
										E1 = 13.6 * atomPropertyPtr->atomicNumber * atomPropertyPtr->atomicNumber * 1e-6;
										r0 = isotopePropertyPtr->relativeNucleusRadius * e_r0;
										CP0 = (e_elementaryCharge * e_elementaryCharge * atomPropertyPtr->atomicNumber)/(4.0* M_PI *e_vacuumPermittivity * r0) / e_energyJouleOfeV / 1.0E6;
										fprintf(a_fp, " %lg %lg %lg %lg %lg %lg\n",  r1, CP1, E1, CP1 - E1, r0, CP0);
									}
								}
							}
						}
					}
				}
			}
		}else{
		}
	}
	fputs("AN symbol name MN HLs dRt masU mas mDf r1 CP1 E1 CP1-E1 r0 CP0\n", a_fp);
	fputs("\n", a_fp);
	exit(0);
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
		fprintf(a_fp, " NucleusRadius %lg\n", a_isotopePropertyPtr->relativeNucleusRadius * e_r0);
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
#define SCAT_BIG_MASS_DEFECT_ALL -2//The big mass defect type is the part of energy by mass defection from either cillide of nucleus particles or decay of nuecuses and the amont of energy is eaual or bigger than "e_collideMiniMeV".
#define SCAT_LOST_BY_NEUTRINO -1 //the enegey type carryed out by neutrino in beta decays.

#define SCAT_NOT_COLLIDE 0 // The not collide type is a part of enegy of input high voltage pulse current. The part is not used to collide of nucleus perticles by high voltage pulse current.
#define SCAT_IMPERFECT 1  // The IMPERFECT type is a part of enegy of input high voltage pulse current.
#define SCAT_IMPERFECT_IN_SPACE 2  // The IMPERFECT type is a part of enegy of input high voltage pulse current.
#define SCAT_SMALL_MASS_DEFECT 3 // The small mass defect type is the output energy from neclear actions. But they are less than "e_collideMiniMeV". Therefore they can't cause the next neclear action anymore.
#define SCAT_SMALL_BY_COMPTON 4 // This Comton effect type is the reflection photon with less energy than "e_collideMiniMeV". Therefore they can't cause the next neclear action anymore.
#define SCAT_SMALL_BY_COMPTON_E 5 // This Comton effect type is the enegy passed to an electron from a photon. But the amount of energy is less than "e_collideMiniMeV". Therefore they can't cause the next neclear action anymore.


#define SCAT_DEGRADE_IN_COMPTON_B 6 //If there are too many SCAT_BIG_MASS_DEFECT_NOW nodes in the hash table, "massDefectHashTable[MASS_DEFECT_BY_GANMMA]", they waste the memory and power of a computer. So, we do degreade them by shopping the energy into 32 greades. The SCAT_DEGRADE_IN_COMPTON_B type contains shopped small remain part of BIG energy. The grading big energy is still in the "massDefectHashTable[MASS_DEFECT_BY_GANMMA]".
#define SCAT_DEGRADE_IN_COMPTON_S 7
#define SIZE_OF_SCATTERD 8 //The array size of scatterd photon types. Those types of photon are reason of heat.
extern const char * getMassDefectByName(int a_massDefectBy);

extern const char * getScatName(int a_scatType)
{
	char * retPtr = "";
	switch(a_scatType){
		case SCAT_BIG_MASS_DEFECT_NOW: retPtr = "BIG_MASS_DEFECT_NOW"; break;
		case SCAT_BIG_MASS_DEFECT_ALL: retPtr = "BIG_MASS_DEFECT_ALL"; break;
		case SCAT_LOST_BY_NEUTRINO: retPtr = "LOST_BY_NEUTRINO"; break;
		//
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
extern long double printSumMeVMol(FILE * a_fp, const char * a_operatorStr, int a_scatType, int a_massDefectBy, const struct sumMeVMol * a_sumMeVMol, long double a_grandTotal, const char * a_rateName)
{
	long double MeVMolNAvogadro, rate = 0.0;
	MeVMolNAvogadro = a_sumMeVMol->MeVMol * NAvogadro;
	if(a_scatType == SCAT_BIG_MASS_DEFECT_NOW){
		fprintf(a_fp, "%s %s %Lg [MeV] %lld[cnt]",
			a_operatorStr, getMassDefectByName(a_massDefectBy), MeVMolNAvogadro, a_sumMeVMol->cnt);
	}else{
		fprintf(a_fp, "%s %s %Lg [MeV] %lld[cnt]",
			a_operatorStr, getScatName(a_scatType), MeVMolNAvogadro, a_sumMeVMol->cnt);
	}
	if(a_grandTotal > 0.0){
		rate = 100.0 * MeVMolNAvogadro / a_grandTotal;
		fprintf(a_fp, " %Lg %s", rate, a_rateName);
	}
	fprintf(a_fp, "\n (range[ %lg , %lg ][MeV] * [ %lg , %lg ][mol])\n", 
			a_sumMeVMol->minMeV, a_sumMeVMol->maxMeV,
			a_sumMeVMol->minMol, a_sumMeVMol->maxMol);
	return MeVMolNAvogadro;
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
		printSumMeVMol(a_fp, (i == 0) ? " " : "+", i, 0, &a_scatteredAry[i], grandTotal, "[%]");
	}
	fprintf(a_fp, "= %Lg [MeV]\n",  grandTotal);
	return grandTotal;
}

//---------------------------------------------------------------------
#define MASS_DEFECT_BY_NEUTRINO (-1) //No count
#define MASS_DEFECT_BY_PROTON 0
#define MASS_DEFECT_BY_DEUTERIUM 1
#define MASS_DEFECT_BY_TRITIUM 2
#define MASS_DEFECT_BY_NEUTRON 3
#define MASS_DEFECT_BY_ALPHA 4
#define MASS_DEFECT_BY_BETA 5
#define MASS_DEFECT_BY_GANMMA 6
#define COUNT_OF_MASS_DEFECT_HASH_TABLE 7

extern const char * getMassDefectByName(int a_massDefectBy)
{
	const char * pRet = "";
	switch(a_massDefectBy){
	case MASS_DEFECT_BY_NEUTRINO: pRet = "MASS_DEFECT_BY_NEUTRINO"; break;
	case MASS_DEFECT_BY_PROTON: pRet = "MASS_DEFECT_BY_PROTON"; break;
	case MASS_DEFECT_BY_DEUTERIUM: pRet = "MASS_DEFECT_BY_DEUTERIUM"; break;
	case MASS_DEFECT_BY_TRITIUM: pRet = "MASS_DEFECT_BY_TRITIUM"; break;
	case MASS_DEFECT_BY_NEUTRON: pRet = "MASS_DEFECT_BY_NEUTRON"; break;
	case MASS_DEFECT_BY_ALPHA: pRet = "MASS_DEFECT_BY_ALPHA"; break;
	case MASS_DEFECT_BY_BETA: pRet = "MASS_DEFECT_BY_BETA"; break;
	case MASS_DEFECT_BY_GANMMA: pRet = "MASS_DEFECT_BY_GANMMA"; break;
	}
	return pRet;
}
struct electrode {
	double detectLimitMolForIsotope;
	struct hashTable atomHashTable;
	struct hashTable massDefectHashTable[COUNT_OF_MASS_DEFECT_HASH_TABLE];
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
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_PROTON], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_DEUTERIUM], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_TRITIUM], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_NEUTRON], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_ALPHA], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_BETA], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
		if(!fwritedHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_GANMMA], a_fp, fwriteRange, fwriteMassDefect)){return 0;}\
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
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_PROTON], a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_DEUTERIUM], a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_TRITIUM], a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_NEUTRON], a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_ALPHA], a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_BETA], a_fp,\
			allocFreadRange, allocFreadMassDefect, \
			calcHashSeedOfRange,\
			allocCopyRange,\
			free,\
			comparRange,\
			allocCopyMassDefect,\
			free,\
			foundActionForMassDefect)){return 0;}\
		if(!freadHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_GANMMA], a_fp,\
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
extern struct atomNodeConst * findProtonInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN);
}
extern struct atomNodeConst * findDeuteriumInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_DEUTERIUM, MASSNUMBER_DEUTERIUM);
}
extern struct atomNodeConst * findTritiumInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_TRITIUM, MASSNUMBER_TRITIUM);
}
extern struct atomNodeConst * findHeliumInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM);
}
extern struct atomNodeConst * findNeutronInElectrode(struct electrode * a_electrodePtr)
{
	return findAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON);
}

extern void increaseAtom(struct electrode * a_electrodePtr, int a_atomicNumber, int a_massNumber, double a_increaseMol)
{
	//fprintf(stderr, "DEBUG:%s:BEGIN{a_atomicNumber:%d massNumber:%d a_increaseMol:%lg \n", __FUNCTION__, a_atomicNumber, a_massNumber, a_increaseMol);
	if(a_increaseMol >= 0.0){
		struct atomKey akey;
		akey.isotopePropertyPtr = getIsotopePropertyPtr(a_atomicNumber, a_massNumber);
		if(akey.isotopePropertyPtr){
			struct atomNodeConst * atomPtr;
			double molNow;
			struct atomValue aValue;
			struct objectNodeConst nodeConst;
			atomPtr = findAtom(a_electrodePtr, a_atomicNumber, a_massNumber);
			if(atomPtr){
				molNow = getMol(atomPtr->vPtr);
			}else{
				molNow = 0.0;
			}
			memset(&aValue, 0, sizeof(aValue));
			aValue.molIni = a_increaseMol;
			nodeConst.keyPtr = &akey;
			nodeConst.valuePtr = &aValue;
			//if(molNow > 0.0){
				//double warnScale = 128.0;
				if(a_increaseMol > 4.0 /* molNow * warnScale*/){//DEBUG
					//double exitScale = 1.0e8;
					char * levelPtr = "WARN";
					//int doExit = 0;
					//if(a_increaseMol > molNow * exitScale){
					//	levelPtr = "FATAL";
					//	doExit = 1;
					//}
					fprintf(stderr, "%s:%s(%d):%s %s a_atomicNumber:%d massNumber:%d a_increaseMol:%lg molNow:%lg\n", levelPtr, __FUNCTION__, __LINE__, a_electrodePtr->atomHashTable.tableName, akey.isotopePropertyPtr->symbol, a_atomicNumber, a_massNumber, a_increaseMol, molNow);
					fprintf(e_logFp, "%s:%s(%d):%s %s a_atomicNumber:%d massNumber:%d a_increaseMol:%lg molNow:%lg\n", levelPtr, __FUNCTION__, __LINE__, a_electrodePtr->atomHashTable.tableName, akey.isotopePropertyPtr->symbol, a_atomicNumber, a_massNumber, a_increaseMol, molNow);
					//if(doExit){
					//	exit(1);
					//}
				}
			//}
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
	
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_PROTON], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_PROTON)),
		100);
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_DEUTERIUM], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_DEUTERIUM)),
		100);
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_TRITIUM], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_TRITIUM)),
		100);
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_NEUTRON], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_NEUTRON)),
		100);
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_ALPHA], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_ALPHA)),
		100);
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_BETA], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_BETA)),
		100);
	initHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_GANMMA], 
		calcHashSeedOfRange,
		allocCopyRange,
		free,
		comparRange,
		allocCopyMassDefect,
		free,
		foundActionForMassDefect,
		strcat(strcat(strcpy(name, a_electrodeName), " "), getMassDefectByName(MASS_DEFECT_BY_GANMMA)),
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

extern void registOutput(struct electrode * a_electrodePtr, int a_scatType, int a_massDefectBy, double a_MeV, double a_mol)
{
	//fprintf(stderr, "DEBUG:%s:{BEGIN %lp, %s a_MeV:%lg * a_mol:%lg = %lg\n\n", __FUNCTION__, a_electrodePtr, getScatName(a_scatType), a_MeV, a_mol, a_MeV * a_mol);
	if(a_scatType == SCAT_BIG_MASS_DEFECT_ALL){
		fprintf(stderr, "FATAL:%s:%p, %s %s a_MeV:%lg * a_mol:%lg = %lg\n", __FUNCTION__, a_electrodePtr, getScatName(a_scatType), getMassDefectByName(a_massDefectBy), a_MeV, a_mol, a_MeV * a_mol);
		exit(1);
	}
	if(a_mol > 0.0){
		//Check data!
		if(a_scatType == SCAT_BIG_MASS_DEFECT_NOW){
			if(a_MeV < e_collideMiniMeV){
				a_scatType = SCAT_SMALL_MASS_DEFECT;
			}
		}
		if(a_scatType == SCAT_BIG_MASS_DEFECT_NOW){
			if(MASS_DEFECT_BY_PROTON <= a_massDefectBy && a_massDefectBy <= MASS_DEFECT_BY_GANMMA){
				struct range key;
				struct massDefect value;
				struct objectNodeConst nodeConst;
				key = calcPrecisionRange(a_MeV);
				value.MeV = a_MeV;
				value.mol = a_mol;
				nodeConst.keyPtr = &key;
				nodeConst.valuePtr = &value;
				//if(a_massDefectBy != MASS_DEFECT_BY_GANMMA){
				//	fprintf(stderr, "DEBUG:%s(%d):%s %s a_MeV:%lg * a_mol:%lg = %lg\n\n", __FUNCTION__, __LINE__, getScatName(a_scatType), getMassDefectByName(a_massDefectBy), a_MeV, a_mol, a_MeV * a_mol);
				//}
				insertObjectInHashTable(&a_electrodePtr->massDefectHashTable[a_massDefectBy], &nodeConst);
				//insertObjectInHashTable(&a_electrodePtr->massDefectHashTable[MASS_DEFECT_BY_GANMMA], &nodeConst);//DEBUG
				addSumMeVMol(&a_electrodePtr->massDefectAll, a_MeV, a_mol);
			}else{
				fprintf(stderr, "FATAL:%s:%p, %s %s a_MeV:%lg * a_mol:%lg = %lg\n", __FUNCTION__, a_electrodePtr, getScatName(a_scatType), getMassDefectByName(a_massDefectBy), a_MeV, a_mol, a_MeV * a_mol);
				exit(1);				
			}
		}else if(a_scatType == SCAT_LOST_BY_NEUTRINO){
			if(a_massDefectBy == MASS_DEFECT_BY_NEUTRINO){
				addSumMeVMol(&a_electrodePtr->byNeutrinoAll, a_MeV, a_mol);
			}else{
				fprintf(stderr, "FATAL:%s:%p, %s %s a_MeV:%lg * a_mol:%lg = %lg\n", __FUNCTION__, a_electrodePtr, getScatName(a_scatType), getMassDefectByName(a_massDefectBy), a_MeV, a_mol, a_MeV * a_mol);
				exit(1);
			}
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
		registOutput(a_electrodePtr, SCAT_LOST_BY_NEUTRINO, MASS_DEFECT_BY_NEUTRINO, energyCarryedByNeutrino, molOfNeutrino);
		//debugSum += energyCarryedByNeutrino * molOfNeutrino;//DEBUG
		registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_BETA, remainMassDefect, molOfElectoron);
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
		registOutput(a_electrodePtr, SCAT_LOST_BY_NEUTRINO, MASS_DEFECT_BY_NEUTRINO, energyCarryedByNeutrino, molOfNeutrino);
		//debugSum += energyCarryedByNeutrino * molOfNeutrino;//DEBUG
		registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, remainMassDefect, molOfDaughter);	
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
		registOutput(a_electrodePtr, SCAT_LOST_BY_NEUTRINO, MASS_DEFECT_BY_NEUTRINO, energyCarryedByNeutrino, molOfNeutrino);
		debugSum += energyCarryedByNeutrino * molOfNeutrino;//DEBUG
		registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, energyOfGamma, molOfGamma);
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

extern void sumupBigMassDefect(struct electrode * a_electrodePtr, int a_i, struct sumMeVMol a_massDefectNow[COUNT_OF_MASS_DEFECT_HASH_TABLE])
{
	initSumMeVMol(&a_massDefectNow[a_i]);
	iterateInHashTable(&a_electrodePtr->massDefectHashTable[a_i], &a_massDefectNow[a_i], iterateSumupMassDefect);
}
extern long double printMassDefectHeat(FILE * a_fp, const char * a_titlePtr, const struct sumMeVMol a_massDefectNow[COUNT_OF_MASS_DEFECT_HASH_TABLE], const struct sumMeVMol * a_massDefectAllPtr, const struct sumMeVMol * a_byNeutrinoAllPtr, struct sumMeVMol a_scatteredAry[SIZE_OF_SCATTERD])
{
	long double outputHeat = 0.0;
	int i;
	fprintf(a_fp, "%s MASS DEFECT\n(gamma ray >= %lg [MeV], beta energy = %lg [MeV])\n", a_titlePtr, e_collideMiniMeV, e_betaEnergyMeV);
	printSumMeVMol(a_fp, " ", SCAT_BIG_MASS_DEFECT_ALL, 0, a_massDefectAllPtr, 0.0, "");
	for(i = 0; i < COUNT_OF_MASS_DEFECT_HASH_TABLE; ++i){
		outputHeat += printSumMeVMol(a_fp, " ", SCAT_BIG_MASS_DEFECT_NOW, i, &a_massDefectNow[i], 0.0, "");
	}
	fprintf(a_fp, "(BIG_MASS_DEFECT_NOW toatl %Lg [MeV])\n", outputHeat);
	outputHeat += printScattered(a_fp, a_titlePtr, a_scatteredAry);
	fprintf(a_fp, "%s LOST HEAT =\n", a_titlePtr);
	printSumMeVMol(a_fp, " ", SCAT_LOST_BY_NEUTRINO, 0, a_byNeutrinoAllPtr, outputHeat, "[%(/Output heat)]");
	return outputHeat;
}
extern double printElectrode(FILE * a_fp, struct electrode * a_electrodePtr, struct sumMeVMol a_massDefectNow[COUNT_OF_MASS_DEFECT_HASH_TABLE], double * a_sumOfMassUMolIniPtr, double * a_sumOfMassUMolAddPtr, double * a_sumOfMassUMolSubPtr, double * a_heatCapacityPtr)
{
	//long double grandTotal;
	double sumOfMassUMol;
	int i;

	fprintf(a_fp, "%s detectLimitMolForIsotope %lg\n", a_electrodePtr->atomHashTable.tableName, a_electrodePtr->detectLimitMolForIsotope);
	sumOfMassUMol = printAtomList(a_fp, &a_electrodePtr->atomHashTable, a_sumOfMassUMolIniPtr, a_sumOfMassUMolAddPtr, a_sumOfMassUMolSubPtr, a_heatCapacityPtr);
	for(i = 0; i < COUNT_OF_MASS_DEFECT_HASH_TABLE; ++i){
		sumupBigMassDefect(a_electrodePtr, i, a_massDefectNow);
	}
	/*grandTotal = */printMassDefectHeat(a_fp, a_electrodePtr->atomHashTable.tableName, a_massDefectNow, &a_electrodePtr->massDefectAll, &a_electrodePtr->byNeutrinoAll, a_electrodePtr->scattered);
	fputs("\n", a_fp);
	return sumOfMassUMol;
}
//---------------------------------------------------------------------

#define COLLIDE_PROTON      1 //by high volatage
#define COLLIDE_DEUTERIUM   2 //by high volatage
#define COLLIDE_TRITIUM     3 //by high volatage
#define COLLIDE_ELECTRON    4 //by high volatage
#define COLLIDE_NEUTRON     5 //absorbing cold neutrons
#define COLLIDE_PROTON_S    6 //by scattering
#define COLLIDE_DEUTERIUM_S 7 //by scattering
#define COLLIDE_TRITIUM_S   8 //by scattering
#define COLLIDE_ALPHA_S     9 //by scattering
#define COLLIDE_ELECTRON_S 10 //by scattering, the compton Effect

extern const char * getCollideName(int a_collideType)
{
	const char * collideNamePtr = "unkown";
	switch(a_collideType){
	case COLLIDE_PROTON:      collideNamePtr = "COLLIDE_PROTON_BY_HIGH_VOLTAGE"; break;
	case COLLIDE_DEUTERIUM:   collideNamePtr = "COLLIDE_DEUTERIUM_BY_HIGH_VOLTAGE"; break;
	case COLLIDE_TRITIUM:     collideNamePtr = "COLLIDE_TRITIUM_BY_HIGH_VOLTAGE";  break;
	case COLLIDE_ELECTRON:	  collideNamePtr = "COLLIDE_ELECTRON_BY_HIGH_VOLTAGE"; break;
	case COLLIDE_NEUTRON:     collideNamePtr = "COLLIDE_NEUTRON";  break;
	case COLLIDE_PROTON_S:    collideNamePtr = "COLLIDE_PROTON_BY_SCATTERING"; break;
	case COLLIDE_DEUTERIUM_S: collideNamePtr = "COLLIDE_DEUTERIUM_BY_SCATTERING"; break;
	case COLLIDE_TRITIUM_S:   collideNamePtr = "COLLIDE_TRITIUM_BY_SCATTERING"; break;
	case COLLIDE_ALPHA_S:     collideNamePtr = "COLLIDE_ALPHA_BY_SCATTERING"; break;
	case COLLIDE_ELECTRON_S:  collideNamePtr = "COLLIDE_ELECTRON_BY_SCATTERING"; break;
	default:
		if(a_collideType > MIN_DECAY_MODE){
			collideNamePtr = getDecayModeText(a_collideType);
		}else{
			fprintf(stderr, "FATAL:%s(%d):UNKNOWN COLLIDE %d\n", __FUNCTION__, __LINE__, a_collideType);
			exit(1);
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
#define MAX_NEWCOUNT 4

#define EMISSION_IS_GAMMA 0 //default
#define EMISSION_IS_ELECTRO_NEUTRINO 1
#define EMISSION_IS_PROTON 4
#define EMISSION_IS_NEUTRON 5
#define EMISSION_IS_2_ELECTRONS 7 //collision between two electrons

	int newCount;
	int emissionType[MAX_NEWCOUNT];
	int newAtomicNumber[MAX_NEWCOUNT];
	int newMassNumber[MAX_NEWCOUNT];
	const struct isotopeProperty * newIsotopePropertyPtr[MAX_NEWCOUNT];
	double minAppliedVoltageMeV[MAX_NEWCOUNT];
	double massDefectMeV[MAX_NEWCOUNT];
	double minMassDefectMeV[MAX_NEWCOUNT];
	char * plusAlpha[MAX_NEWCOUNT];
	char nuclearReactionStr[MAX_NEWCOUNT][REACTION_LEN + 1];
	double rate[MAX_NEWCOUNT];
	int subAtomicNumber[MAX_NEWCOUNT];
	int subMassNumber[MAX_NEWCOUNT];
	const struct isotopeProperty * subIsotopePropertyPtr[MAX_NEWCOUNT];
	//}

	double newIsotopeMol;// set by collideParticleAtom
};
extern void useElectronOrNeutron(struct collide * a_pX)
{
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_S){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 0;
	}else if(a_pX->collideType == COLLIDE_PROTON || a_pX->collideType == COLLIDE_PROTON_S){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 1;
	}else if(a_pX->collideType == COLLIDE_DEUTERIUM || a_pX->collideType == COLLIDE_DEUTERIUM_S){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 1;
	}else if(a_pX->collideType == COLLIDE_TRITIUM || a_pX->collideType == COLLIDE_TRITIUM_S){
		a_pX->useElectron = 1;
		a_pX->useNeutron = 1;
	}else if(a_pX->collideType == COLLIDE_ALPHA_S){
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
#define VERSION_OF_serializeNuclearReactionKey 2 //I changed the order of decayOrCollideType. (Oct 20, 2016)
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
	int numberOfBreakingIntervals;
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
extern void printAllNuclearReaction(FILE * a_fp, const char * a_timeMessPtr)
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
					fprintf(a_fp, "%s\n[ %s ReactionType %s ]>>>\n", a_timeMessPtr, electrodeTypeStr, reactionTypeStr);
				}
				oldElectrodeType = keyPtr->electrodeType;
				oldReactionType = keyPtr->reactionType;
				electrodeTypeStr = getElectrodeName(keyPtr->electrodeType);
				electrodeChar = getElectrodeChar(keyPtr->electrodeType);
				reactionTypeStr = getNuclearReactionName(keyPtr->reactionType);
				reactionChar = getNuclearReactionChar(keyPtr->reactionType);
				fprintf(a_fp, "\n[ %s ReactionType %s ]<<<\n%s\n", electrodeTypeStr, reactionTypeStr, a_timeMessPtr);
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
			fprintf(a_fp, "%s\n[ %s ReactionType %s ]>>>\n", a_timeMessPtr, electrodeTypeStr, reactionTypeStr);
		}
		free(nodePtrs);
	}
}


//---------------------------------------------------------------------
extern void sprintNoDataOfAtom(struct collide * a_pX, struct atomProperty * a_atomPropertyPtr, int a_index)
{
	if(a_pX->printNuclearReaction){
		if(a_atomPropertyPtr){
			snprintf(a_pX->nuclearReactionStr[a_index], REACTION_LEN,
				"There is no data of atom, %s: %s + %s -> %d%s\n", 
				getCollideName(a_pX->collideType), 
				a_pX->targetIsotopePropertyPtr->symbol, a_pX->bulletPropertyPtr->nucleusName,
				a_pX->newMassNumber[a_index], a_atomPropertyPtr->symbol);
		}else{
			snprintf(a_pX->nuclearReactionStr[a_index], REACTION_LEN,
				"There is no data of atom, %s: %s + %s -> %d~%d(mass~atom)\n", 
				getCollideName(a_pX->collideType), 
				a_pX->targetIsotopePropertyPtr->symbol, a_pX->bulletPropertyPtr->nucleusName,
				a_pX->newMassNumber[a_index], a_pX->newAtomicNumber[a_index]);		
		}
	}
}
extern void checkNoDataOfAtom(struct collide * a_pX, int a_index)
{
	struct atomProperty * atomPropertyPtr = getAtomPropertyPtr(a_pX->newAtomicNumber[a_index]);
	if(atomPropertyPtr){
		if(atomPropertyPtr->min1secMassNumber <= a_pX->newMassNumber[a_index] 
		&& a_pX->newMassNumber[a_index] <= atomPropertyPtr->max1secMassNumber){
			sprintNoDataOfAtom(a_pX, atomPropertyPtr, a_index);
		}else{
			;//The new isotope has very short half-time. So, we ignore them.
		}
	}else{
		sprintNoDataOfAtom(a_pX, atomPropertyPtr, a_index);
	}
}
int e_debugElectroPotential = 0;
int e_logElectroPotential = 0;
extern void calcElectroPotentialMeVByCollide(struct collide * a_pX)
{
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_S){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			if(a_pX->appliedVoltageMeV <= 0.0){
				fprintf(stderr, "FATAL:%s(%d):invalid appliedVoltageMeV:%lg\n", __FUNCTION__, __LINE__, a_pX->appliedVoltageMeV);
				exit(1);
			}
			//The coulomb barrier between an electron and an electron
			a_pX->electroPotentialMeV = a_pX->appliedVoltageMeV * 2.0;//Coulomb potential to allow a moving electron and a static electron to approach
		}else if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_NEUTRON 
		&& (a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_NEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_DINEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_TRINEUTRON
			)
		){
			a_pX->electroPotentialMeV = 0.0;
		}else{
			a_pX->electroPotentialMeV = -1.0 * e_coefElectroPotentialMeV * (a_pX->targetIsotopePropertyPtr->atomicNumber) / ((a_pX->targetIsotopePropertyPtr->relativeNucleusRadius + a_pX->bulletPropertyPtr->relativeNucleusRadius) * e_r0);
		}
	}else if(a_pX->collideType == COLLIDE_PROTON || a_pX->collideType == COLLIDE_PROTON_S
	|| a_pX->collideType == COLLIDE_DEUTERIUM || a_pX->collideType == COLLIDE_DEUTERIUM_S
	|| a_pX->collideType == COLLIDE_TRITIUM || a_pX->collideType == COLLIDE_TRITIUM_S
	|| a_pX->collideType == COLLIDE_ALPHA_S){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			a_pX->electroPotentialMeV = -1.0 * e_coefElectroPotentialMeV * (a_pX->bulletPropertyPtr->atomicNumber) / ((a_pX->bulletPropertyPtr->relativeNucleusRadius + a_pX->targetIsotopePropertyPtr->relativeNucleusRadius) * e_r0);
		}else if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_NEUTRON 
		&& (a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_NEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_DINEUTRON
			|| a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_TRINEUTRON
			)
		){
			a_pX->electroPotentialMeV = 0.0;
		}else{
			//The coulomb barrier between a proton/deuterium/tritium and a nucleus
			a_pX->electroPotentialMeV = e_coefElectroPotentialMeV * (a_pX->targetIsotopePropertyPtr->atomicNumber * a_pX->bulletPropertyPtr->atomicNumber) / ((a_pX->targetIsotopePropertyPtr->relativeNucleusRadius + a_pX->bulletPropertyPtr->relativeNucleusRadius) * e_r0);
		}
	}else{//COLLIDE_NEUTRON
		a_pX->electroPotentialMeV = 0.0;
	}
	if(!a_pX->printNuclearReaction){//DEBUG
		if(e_debugElectroPotential){
			fprintf(stderr, "DEBUG:%d:%s + %s, electroPotential %lg %s appliedVoltage %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->electroPotentialMeV, ((a_pX->electroPotentialMeV > a_pX->appliedVoltageMeV) ? ">" : ((a_pX->electroPotentialMeV == a_pX->appliedVoltageMeV) ? "==" : "<")), a_pX->appliedVoltageMeV);
		}
		if(e_logElectroPotential){
			fprintf(e_logFp, "[INFO]:%d:%s + %s, electroPotential %lg %s appliedVoltage %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->electroPotentialMeV, ((a_pX->electroPotentialMeV > a_pX->appliedVoltageMeV) ? ">" : ((a_pX->electroPotentialMeV == a_pX->appliedVoltageMeV) ? "==" : "<")), a_pX->appliedVoltageMeV);
		}
	}
	//fprintf(stderr, "DEBUG:%s:%s target:%s electroPotentialMeV=%lg\n", __FUNCTION__, getCollideName(a_pX->collideType), a_pX->targetIsotopePropertyPtr->symbol, a_pX->electroPotentialMeV);
}
extern void set_minMassDefectMeV_plusAlpha(struct collide * a_pX, int a_index)
{
	//fprintf(stderr, "DEBUG:%s:BEGIN{a_pX %p\n", __FUNCTION__, a_pX);	
	if(a_pX->minAppliedVoltageMeV[a_index] >= 0.0){
		a_pX->minMassDefectMeV[a_index] = 0.0;
	}else{
		a_pX->minMassDefectMeV[a_index] = - a_pX->minAppliedVoltageMeV[a_index];
		a_pX->minAppliedVoltageMeV[a_index] = 0.0;
	}
	if(a_pX->massDefectMeV[a_index] >= 0.0){
		a_pX->plusAlpha[a_index] = " + x";
	}else{
		a_pX->plusAlpha[a_index] = " - x";
	}
	//fprintf(stderr, "DEBUG:%s:END}\n", __FUNCTION__);
}
int e_debugNewIsotope = 0;
int e_logNewIsotope = 0;
extern void getCollidedNewIsotopePropertyPtr(struct collide * a_pX)
{
	//fprintf(stderr, "DEBUG:%s:BEGIN{a_pX %p\n", __FUNCTION__, a_pX);
	calcElectroPotentialMeVByCollide(a_pX);
	a_pX->newCount = 1;
	a_pX->newAtomicNumber[0] = UNDEF_ATOMIC_NUMBER;//for safety
	a_pX->newMassNumber[0] = UNDEF_MASS_NUMBER;//for safety
	a_pX->newIsotopePropertyPtr[0] = NULL;//for safety
	a_pX->minAppliedVoltageMeV[0] = 0.0;//for safety
	a_pX->massDefectMeV[0] = 0.0;//for safety
	a_pX->minMassDefectMeV[0] = 0.0;//for safety
	a_pX->plusAlpha[0] = "";//for safety
	a_pX->nuclearReactionStr[0][0] = 0;//for safety
	a_pX->rate[0] = 1.0;
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_S){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			a_pX->emissionType[0] = EMISSION_IS_2_ELECTRONS;
			a_pX->newAtomicNumber[0] = a_pX->targetIsotopePropertyPtr->atomicNumber;
			a_pX->newMassNumber[0] = a_pX->targetIsotopePropertyPtr->massNumber;
			a_pX->newIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->newAtomicNumber[0], a_pX->newMassNumber[0]);
			if(a_pX->newIsotopePropertyPtr[0]){
				a_pX->minAppliedVoltageMeV[0] = 0.0;
				a_pX->massDefectMeV[0] = 0.0;//a_pX->appliedVoltageMeV;
				set_minMassDefectMeV_plusAlpha(a_pX, 0);
				if(a_pX->printNuclearReaction){
					snprintf(a_pX->nuclearReactionStr[0], REACTION_LEN,
						"e-(%lg) + e-(%lg) + (X+0)[MeV] -> e-(%lg) + (Y)[MeV] + e-(%lg) + (Z)[MeV]", 
						e_massElectronMeV, e_massElectronMeV, 
						e_massElectronMeV, e_massElectronMeV);
					//fprintf(stderr, "DEBUG:%s:snprintf %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
				}
			}else{
				checkNoDataOfAtom(a_pX, 0);
			}
		}else{
			a_pX->emissionType[0] = EMISSION_IS_ELECTRO_NEUTRINO;
			a_pX->newAtomicNumber[0] = a_pX->targetIsotopePropertyPtr->atomicNumber - 1;
			a_pX->newMassNumber[0] = a_pX->targetIsotopePropertyPtr->massNumber;
			if(a_pX->newAtomicNumber[0] != ATOMICNUMBER_ELECTRON || a_pX->newMassNumber[0] != MASSNUMBER_ELECTRON){
				a_pX->newIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->newAtomicNumber[0], a_pX->newMassNumber[0]);
				if(a_pX->newIsotopePropertyPtr[0]){
					a_pX->minAppliedVoltageMeV[0] = a_pX->newIsotopePropertyPtr[0]->massMeV - a_pX->targetIsotopePropertyPtr->massMeV;
					a_pX->massDefectMeV[0] = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV[0];
					set_minMassDefectMeV_plusAlpha(a_pX, 0);
					if(a_pX->printNuclearReaction){
						double targetIonMassMeV = a_pX->targetIsotopePropertyPtr->massMeV - e_massElectronMeV;
						char * targetSymbolPtr;
						if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_HYDROGEN){
							targetSymbolPtr = a_pX->targetIsotopePropertyPtr->nucleusName;
						}else{
							targetSymbolPtr = a_pX->targetIsotopePropertyPtr->symbol;
						}
						snprintf(a_pX->nuclearReactionStr[0], REACTION_LEN,
							"%s+(%lg) + e-(%lg) + (%lg%s)[MeV] -> %s(%lg) + νe + (%lg%s)[MeV]", 
							targetSymbolPtr, targetIonMassMeV,
							e_massElectronMeV,
							a_pX->minAppliedVoltageMeV[0], a_pX->plusAlpha[0],
							a_pX->newIsotopePropertyPtr[0]->symbol, a_pX->newIsotopePropertyPtr[0]->massMeV,
							a_pX->minMassDefectMeV[0], a_pX->plusAlpha[0]);
						//fprintf(stderr, "DEBUG:%s:snprintf %s\n", __FUNCTION__, a_pX->nuclearReactionStr[0]);
					}
				}else{
					checkNoDataOfAtom(a_pX, 0);
				}
			}else{
				fprintf(stderr, "FATAL:%s(%d):COLLIDE_ELECTRON:%d\n", __FUNCTION__, __LINE__, a_pX->collideType);
			}
		}
	}else if(a_pX->collideType == COLLIDE_PROTON || a_pX->collideType == COLLIDE_PROTON_S
	|| a_pX->collideType == COLLIDE_DEUTERIUM || a_pX->collideType == COLLIDE_DEUTERIUM_S
	|| a_pX->collideType == COLLIDE_TRITIUM || a_pX->collideType == COLLIDE_TRITIUM_S
	|| a_pX->collideType == COLLIDE_ALPHA_S){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			//fprintf(stderr, "DEBUG:%s:COLLIDE_PROTON... ATOMICNUMBER_ELECTRON\n", __FUNCTION__);
			a_pX->emissionType[0] = EMISSION_IS_ELECTRO_NEUTRINO;
			switch(a_pX->collideType){
				case COLLIDE_PROTON:
				case COLLIDE_PROTON_S:
					a_pX->newAtomicNumber[0] = ATOMICNUMBER_NEUTRON;
					a_pX->newMassNumber[0] = MASSNUMBER_NEUTRON;
					break;
				case COLLIDE_DEUTERIUM:
				case COLLIDE_DEUTERIUM_S:
					a_pX->newAtomicNumber[0] = ATOMICNUMBER_NEUTRON;
					a_pX->newMassNumber[0] = MASSNUMBER_DINEUTRON;
					break;
				case COLLIDE_TRITIUM:
				case COLLIDE_TRITIUM_S:
					a_pX->newAtomicNumber[0] = ATOMICNUMBER_NEUTRON;
					a_pX->newMassNumber[0] = MASSNUMBER_TRINEUTRON;
					break;
				case COLLIDE_ALPHA_S:
					a_pX->newAtomicNumber[0] = MASSNUMBER_HYDROGEN;
					a_pX->newMassNumber[0] = 4;
					break;
			}
			a_pX->newIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->newAtomicNumber[0], a_pX->newMassNumber[0]);
			if(a_pX->newIsotopePropertyPtr[0]){
				a_pX->minAppliedVoltageMeV[0] = a_pX->newIsotopePropertyPtr[0]->massMeV - a_pX->bulletPropertyPtr->massMeV;
				a_pX->massDefectMeV[0] = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV[0];
				set_minMassDefectMeV_plusAlpha(a_pX, 0);
				if(a_pX->printNuclearReaction){
					double bulletIonMassMeV = a_pX->bulletPropertyPtr->massMeV - e_massElectronMeV;
					//fprintf(stderr, "DEBUG:%s:snprintf COLLIDE_PROTON\n", __FUNCTION__);
					snprintf(a_pX->nuclearReactionStr[0], REACTION_LEN,
						"e-(%lg) + %s+(%lg) + (%lg%s)[MeV] -> %s(%lg) + νe + (%lg%s)[MeV]", 
						e_massElectronMeV,
						a_pX->bulletPropertyPtr->nucleusName, bulletIonMassMeV,
						a_pX->minAppliedVoltageMeV[0], a_pX->plusAlpha[0],
						a_pX->newIsotopePropertyPtr[0]->symbol, a_pX->newIsotopePropertyPtr[0]->massMeV,
						a_pX->minMassDefectMeV[0], a_pX->plusAlpha[0]);
					//fprintf(stderr, "DEBUG:%s:snprintf COLLIDE_PROTON %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
				}
			}else{
				checkNoDataOfAtom(a_pX, 0);
			}
		}else{
			//fprintf(stderr, "DEBUG:%s:COLLIDE_PROTON... other\n", __FUNCTION__);
			if(a_pX->collideType == COLLIDE_DEUTERIUM
			&& a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_DEUTERIUM 
			&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_DEUTERIUM){
				//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
				a_pX->newCount = 3;
				a_pX->emissionType[0] = EMISSION_IS_PROTON;
				a_pX->newAtomicNumber[0] = ATOMICNUMBER_TRITIUM;
				a_pX->newMassNumber[0] = MASSNUMBER_TRITIUM;
				a_pX->newIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->newAtomicNumber[0], a_pX->newMassNumber[0]);
				a_pX->rate[0] = 0.5 - 0.5e-6;
				a_pX->subAtomicNumber[0] = ATOMICNUMBER_HYDROGEN;
				a_pX->subMassNumber[0] = MASSNUMBER_HYDROGEN;
				a_pX->subIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->subAtomicNumber[0], a_pX->subMassNumber[0]);
				a_pX->minAppliedVoltageMeV[0] = (a_pX->newIsotopePropertyPtr[0]->massMeV + a_pX->subIsotopePropertyPtr[0]->massMeV) - (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV);
				a_pX->massDefectMeV[0] = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV[0];
				set_minMassDefectMeV_plusAlpha(a_pX, 0);
				if(a_pX->printNuclearReaction){
					//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
					snprintf(a_pX->nuclearReactionStr[0], REACTION_LEN,
						"D(%lg) + D(%lg) + (%lg%s)[MeV] -> T(%lg) + H(%lg) + (%lg%s)[MeV] (rate:%lg)", 
						a_pX->targetIsotopePropertyPtr->massMeV,
						a_pX->bulletPropertyPtr->massMeV,
						a_pX->minAppliedVoltageMeV[0], a_pX->plusAlpha[0],
						a_pX->newIsotopePropertyPtr[0]->massMeV,
						a_pX->subIsotopePropertyPtr[0]->massMeV,
						a_pX->minMassDefectMeV[0], a_pX->plusAlpha[0], a_pX->rate[0]);
					//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, a_pX->nuclearReactionStr[0]);
				}

				//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
				a_pX->emissionType[1] = EMISSION_IS_NEUTRON;
				a_pX->newAtomicNumber[1] = ATOMICNUMBER_HELIUM;
				a_pX->newMassNumber[1] = MASSNUMBER_HELIUM3;
				a_pX->newIsotopePropertyPtr[1] = getIsotopePropertyPtr(a_pX->newAtomicNumber[1], a_pX->newMassNumber[1]);
				a_pX->rate[1] = 0.5 - 0.5e-6;
				a_pX->subAtomicNumber[1] = ATOMICNUMBER_NEUTRON;
				a_pX->subMassNumber[1] = MASSNUMBER_NEUTRON;
				a_pX->subIsotopePropertyPtr[1] = getIsotopePropertyPtr(a_pX->subAtomicNumber[1], a_pX->subMassNumber[1]);
				a_pX->minAppliedVoltageMeV[1] = (a_pX->newIsotopePropertyPtr[1]->massMeV + a_pX->subIsotopePropertyPtr[1]->massMeV) - (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV);
				a_pX->massDefectMeV[1] = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV[1];
				set_minMassDefectMeV_plusAlpha(a_pX, 1);
				//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
				if(a_pX->printNuclearReaction){
					//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
					snprintf(a_pX->nuclearReactionStr[1], REACTION_LEN,
						"D(%lg) + D(%lg) + (%lg%s)[MeV] -> 3He(%lg) + n(%lg) + (%lg%s)[MeV] (rate:%lg)", 
						a_pX->targetIsotopePropertyPtr->massMeV,
						a_pX->bulletPropertyPtr->massMeV,
						a_pX->minAppliedVoltageMeV[1], a_pX->plusAlpha[1],
						a_pX->newIsotopePropertyPtr[1]->massMeV,
						a_pX->subIsotopePropertyPtr[1]->massMeV,
						a_pX->minMassDefectMeV[1], a_pX->plusAlpha[1], a_pX->rate[1]);
					//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, a_pX->nuclearReactionStr[1]);
				}

				//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
				a_pX->emissionType[2] = EMISSION_IS_GAMMA;
				a_pX->newAtomicNumber[2] = ATOMICNUMBER_HELIUM;
				a_pX->newMassNumber[2] = MASSNUMBER_HELIUM;
				a_pX->newIsotopePropertyPtr[2] = getIsotopePropertyPtr(a_pX->newAtomicNumber[2], a_pX->newMassNumber[2]);
				a_pX->rate[2] = 1.0 - a_pX->rate[0] - a_pX->rate[1];
				a_pX->subAtomicNumber[2] = UNDEF_ATOMIC_NUMBER;
				a_pX->subMassNumber[2] = UNDEF_MASS_NUMBER;
				a_pX->subIsotopePropertyPtr[2] = NULL;
				a_pX->minAppliedVoltageMeV[2] = a_pX->newIsotopePropertyPtr[2]->massMeV - (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV);
				a_pX->massDefectMeV[2] = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV[2];
				set_minMassDefectMeV_plusAlpha(a_pX, 2);
				if(a_pX->printNuclearReaction){
					//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
					snprintf(a_pX->nuclearReactionStr[2], REACTION_LEN,
						"D(%lg) + D(%lg) + (%lg%s)[MeV] -> 4He(%lg) + (%lg%s)[MeV] (rate:%lg)", 
						a_pX->targetIsotopePropertyPtr->massMeV,
						a_pX->bulletPropertyPtr->massMeV,
						a_pX->minAppliedVoltageMeV[2], a_pX->plusAlpha[2],
						a_pX->newIsotopePropertyPtr[2]->massMeV,
						a_pX->minMassDefectMeV[2], a_pX->plusAlpha[2], a_pX->rate[2]);
					//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, a_pX->nuclearReactionStr[2]);
				}
				//fprintf(stderr, "DEBUG:%s(%d):D + D\n", __FUNCTION__, __LINE__);
			}else{
				a_pX->emissionType[0] = EMISSION_IS_GAMMA;
				a_pX->newAtomicNumber[0] = a_pX->targetIsotopePropertyPtr->atomicNumber + a_pX->bulletPropertyPtr->atomicNumber;
				a_pX->newMassNumber[0] = a_pX->targetIsotopePropertyPtr->massNumber + a_pX->bulletPropertyPtr->massNumber;
				a_pX->newIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->newAtomicNumber[0], a_pX->newMassNumber[0]);
				if(a_pX->newIsotopePropertyPtr[0]){
					a_pX->minAppliedVoltageMeV[0] = a_pX->newIsotopePropertyPtr[0]->massMeV - (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV);
					a_pX->massDefectMeV[0] = a_pX->appliedVoltageMeV - a_pX->minAppliedVoltageMeV[0];
					set_minMassDefectMeV_plusAlpha(a_pX, 0);
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
							newAtomMassMeV = a_pX->newIsotopePropertyPtr[0]->massMeV - e_massElectronMeV;
							newSymbol = a_pX->newIsotopePropertyPtr[0]->nucleusName;
							newIonSign = "+";
						}else{
							targetSymbol = a_pX->targetIsotopePropertyPtr->symbol;
							targetIonSign = "-";
							targetIonMassMeV = a_pX->targetIsotopePropertyPtr->massMeV + e_massElectronMeV;
							bulletIonMassMeV = a_pX->bulletPropertyPtr->massMeV - e_massElectronMeV;
							newAtomMassMeV = a_pX->newIsotopePropertyPtr[0]->massMeV;
							newSymbol = a_pX->newIsotopePropertyPtr[0]->symbol;
							newIonSign = "";
						}
						//fprintf(stderr, "DEBUG:%s:snprintf other\n", __FUNCTION__);
						snprintf(a_pX->nuclearReactionStr[0], REACTION_LEN,
							"%s%s(%lg) + %s+(%lg) + (%lg%s)[MeV] -> %s%s(%lg) + (%lg%s)[MeV]", 
							targetSymbol, targetIonSign, targetIonMassMeV,
							a_pX->bulletPropertyPtr->nucleusName, bulletIonMassMeV,
							a_pX->minAppliedVoltageMeV[0], a_pX->plusAlpha[0],
							newSymbol, newIonSign, newAtomMassMeV,
							a_pX->minMassDefectMeV[0], a_pX->plusAlpha[0]);
						//fprintf(stderr, "DEBUG:%s:snprintf other %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
					}
				}else{
					checkNoDataOfAtom(a_pX, 0);
				}
			}
		}
	}else if(a_pX->collideType == COLLIDE_NEUTRON){ 
		a_pX->emissionType[0] = EMISSION_IS_GAMMA;
		a_pX->newAtomicNumber[0] = a_pX->targetIsotopePropertyPtr->atomicNumber + a_pX->bulletPropertyPtr->atomicNumber;
		a_pX->newMassNumber[0] = a_pX->targetIsotopePropertyPtr->massNumber + a_pX->bulletPropertyPtr->massNumber;
		a_pX->newIsotopePropertyPtr[0] = getIsotopePropertyPtr(a_pX->newAtomicNumber[0], a_pX->newMassNumber[0]);
		if(a_pX->newIsotopePropertyPtr[0]){
			a_pX->massDefectMeV[0] = (a_pX->targetIsotopePropertyPtr->massMeV + a_pX->bulletPropertyPtr->massMeV) - a_pX->newIsotopePropertyPtr[0]->massMeV;
			if(a_pX->printNuclearReaction){
				//fprintf(stderr, "DEBUG:%s:snprintf other\n", __FUNCTION__);
				snprintf(a_pX->nuclearReactionStr[0], REACTION_LEN,
					"%s(%lg) + %s(%lg) -> %s(%lg) + %lg[MeV]", 
					a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->massMeV,
					a_pX->bulletPropertyPtr->nucleusName, a_pX->bulletPropertyPtr->massMeV,
					a_pX->newIsotopePropertyPtr[0]->symbol, a_pX->newIsotopePropertyPtr[0]->massMeV,
					a_pX->massDefectMeV[0]);
				//fprintf(stderr, "DEBUG:%s:snprintf other %s\n", __FUNCTION__, a_pX->nuclearReactionStr);
			}
		}else{
			checkNoDataOfAtom(a_pX, 0);
		}
	}
	if(!a_pX->printNuclearReaction){//DEBUG
		int i;
		for(i = 0; i < a_pX->newCount; ++i){
			if(a_pX->newIsotopePropertyPtr[i]){
				if(e_debugNewIsotope){
					fprintf(stderr, "DEBUG:%d:%s + %s -> %s massDefect %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->newIsotopePropertyPtr[i]->symbol, a_pX->massDefectMeV[i]);
				}
				if(e_logNewIsotope){
					fprintf(e_logFp, "[INFO]:%d:%s + %s -> %s massDefect %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->newIsotopePropertyPtr[i]->symbol, a_pX->massDefectMeV[i]);
				}
			}else{
				if(e_debugNewIsotope){
					fprintf(stderr, "DEBUG:%d:%s + %s -> NO_NEW_ISOTOPES\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol);
				}
				if(e_logNewIsotope){
					fprintf(e_logFp, "[INFO]:%d:%s + %s -> NO_NEW_ISOTOPES\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol);
				}
			}
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
	double relativeR;
	double CrossSection;
	if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_S){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			// U = - (Z * e^2) / (4 * π * ε0 * r)
			// r = - (Z * e^2) / (4 * π * ε0 * U)
			// Z = 1, U [J] = U' [MeV] * 1.0e6 * 1.60217662E-19 [J/eV]
			
			//I think that the radius of electrons is treated larger up to the Coulomb Barrier in the case of collision between an electron and an electron.
			relativeR = e_coefElectroPotentialMeV / (a_pX->appliedVoltageMeV * e_r0);
			if(relativeR > e_maxApparentRelativeRadiusOfCrossSectionForElectron){
				fprintf(stderr, "WARN:%s(%d):relativeR:%lg > e_maxApparentRelativeRadiusOfCrossSectionForElectron:%lg, R:%lg > %lg \n", __FUNCTION__, __LINE__, relativeR, e_maxApparentRelativeRadiusOfCrossSectionForElectron, relativeR* e_r0, e_maxApparentRelativeRadiusOfCrossSectionForElectron * e_r0);
				relativeR = e_maxApparentRelativeRadiusOfCrossSectionForElectron;
			}
			CrossSection = relativeR * relativeR * a_pX->targetAtomMol;
		}else{
			//I think that the radius of nucleus is treated smaller down to the runing distunce of a Weak Boson in the case of collision between an electron and a nucleus.
			//The electron will collide two up-quark only in protons of a nucleus, not any up-quark in a nutron.
			
			//Although the lifetime of particles that move at a speed close to the speed of light increases due to the effect of the theory of relativity, this effect may not change the collision cross section that has significance in the direction perpendicular to the motion.
			
			//double relativeTheoryExpand = (e_massElectronMeV + a_pX->appliedVoltageMeV) / e_massElectronMeV;
			//static double s_maxRelativeTheoryExpand = 1.0;
			//if(relativeTheoryExpand > s_maxRelativeTheoryExpand){
			//	s_maxRelativeTheoryExpand = relativeTheoryExpand;
			//	fprintf(stderr, "WARN:%s(%d):s_maxRelativeTheoryExpand:%lg \n", __FUNCTION__, __LINE__, s_maxRelativeTheoryExpand);
			//}
			if(e_usebulletCrossSection){
				relativeR = e_relativeRadiusOfElectronByWeakBoson + a_pX->bulletPropertyPtr->relativeNucleusRadius /* * relativeTheoryExpand */;
			}else{
				relativeR = e_relativeRadiusOfElectronByWeakBoson /* * relativeTheoryExpand */;
			}
			CrossSection = relativeR * relativeR * a_pX->targetIsotopePropertyPtr->atomicNumber * 2 * a_pX->targetAtomMol;
		}
	}else if(a_pX->collideType == COLLIDE_PROTON || a_pX->collideType == COLLIDE_PROTON_S
	|| a_pX->collideType == COLLIDE_DEUTERIUM || a_pX->collideType == COLLIDE_DEUTERIUM_S
	|| a_pX->collideType == COLLIDE_TRITIUM || a_pX->collideType == COLLIDE_TRITIUM_S
	|| a_pX->collideType == COLLIDE_ALPHA_S){
		if(a_pX->targetIsotopePropertyPtr->atomicNumber == ATOMICNUMBER_ELECTRON 
		&& a_pX->targetIsotopePropertyPtr->massNumber == MASSNUMBER_ELECTRON){
			//Although the lifetime of particles that move at a speed close to the speed of light increases due to the effect of the theory of relativity, this effect may not change the collision cross section that has significance in the direction perpendicular to the motion.
			
			//double relativeTheoryExpand = (a_pX->bulletPropertyPtr->massMeV + a_pX->appliedVoltageMeV) / a_pX->bulletPropertyPtr->massMeV; 
			//static double s_maxRelativeTheoryExpand = 1.0;
			//if(relativeTheoryExpand > s_maxRelativeTheoryExpand){
			//	s_maxRelativeTheoryExpand = relativeTheoryExpand;
			//	fprintf(stderr, "WARN:%s(%d):s_maxRelativeTheoryExpand:%lg \n", __FUNCTION__, __LINE__, s_maxRelativeTheoryExpand);
			//}
			if(e_usebulletCrossSection){
				relativeR = a_pX->targetIsotopePropertyPtr->relativeNucleusRadius + e_relativeRadiusOfElectronByWeakBoson /* * relativeTheoryExpand */;
			}else{
				relativeR = a_pX->targetIsotopePropertyPtr->relativeNucleusRadius /* * relativeTheoryExpand */;
			}
			CrossSection = relativeR * relativeR * a_pX->bulletPropertyPtr->atomicNumber * 2 * a_pX->targetAtomMol;
		}else{
			if(a_pX->appliedVoltageMeV >= a_pX->electroPotentialMeV){
				relativeR = a_pX->targetIsotopePropertyPtr->relativeNucleusRadius;
				if(e_usebulletCrossSection){
					relativeR += a_pX->bulletPropertyPtr->relativeNucleusRadius;
				}
				CrossSection = relativeR * relativeR * a_pX->targetAtomMol;
			}else{
				// U = - (Z * e^2) / (4 * π * ε0 * r)
				// r = - (Z * e^2) / (4 * π * ε0 * U)
				// Z = Zbullet * Ztarget, U [J] = U' [MeV] * 1.0e6 * 1.60217662E-19 [J/eV]
				
				//I think that the radius of electrons is treated larger up to the Coulomb Barrier in the case of collision between an electron and an electron.

				relativeR = a_pX->bulletPropertyPtr->atomicNumber * a_pX->targetIsotopePropertyPtr->atomicNumber * e_coefElectroPotentialMeV / (a_pX->appliedVoltageMeV * e_r0);
				if(relativeR > e_maxApparentRelativeRadiusOfCrossSectionForNecleus){
					fprintf(stderr, "WARN:%s(%d):relativeR:%lg > e_maxApparentRelativeRadiusOfCrossSectionForNecleus:%lg, R:%lg > %lg \n", __FUNCTION__, __LINE__, relativeR, e_maxApparentRelativeRadiusOfCrossSectionForNecleus, relativeR* e_r0, e_maxApparentRelativeRadiusOfCrossSectionForNecleus * e_r0);
					relativeR = e_maxApparentRelativeRadiusOfCrossSectionForNecleus;
				}
				if(e_usebulletCrossSection){
					relativeR += a_pX->bulletPropertyPtr->relativeNucleusRadius;
				}
				CrossSection = relativeR * relativeR * a_pX->targetAtomMol;
				{//DEBUG
					double relativeR2;
					double CrossSection2;
					relativeR2 = a_pX->targetIsotopePropertyPtr->relativeNucleusRadius;
					if(e_usebulletCrossSection){
						relativeR2 += a_pX->bulletPropertyPtr->relativeNucleusRadius;
					}
					CrossSection2 = relativeR * relativeR * a_pX->targetAtomMol;
					if(CrossSection < CrossSection2){
						fprintf(stderr, "FATAL:DEBUG:%s(%d):CrossSection:%lg < CrossSection2:%lg relativeR:%lg relativeR2:%lg\n", __FUNCTION__, __LINE__, CrossSection, CrossSection2, relativeR, relativeR2);
						exit(1);
					}
				}//DEBUG
			}
		}
	}else{//COLLIDE_NEUTRON
		if(a_pX->massDefectMeV[0] > 0.0){
			relativeR = a_pX->targetIsotopePropertyPtr->relativeNucleusRadius;
			if(e_usebulletCrossSection){
				relativeR +=a_pX->bulletPropertyPtr->relativeNucleusRadius;
			}else{
			}
			CrossSection = relativeR * relativeR * a_pX->targetAtomMol;
		}else{
			CrossSection = 0.0;
		}
	}
	return CrossSection;
}
int e_debugIgnoreCrossSection = 0;
int e_logIgnoreCrossSection = 0;
extern void sumTotalCollideCrossSection(struct collide * a_pX)
{
	a_pX->targetAtomMol = getMol(a_pX->targetAtomValuePtr);
	if(a_pX->targetAtomMol >= a_pX->electrodePtr->detectLimitMolForIsotope){
		a_pX->totalCollideCrossSection += getCrossSection(a_pX);
		a_pX->totalTargetMol += a_pX->targetAtomMol;
		//fprintf(stderr, "DEBUG:%s:totalCollideCrossSection:%lg\n", __FUNCTION__, a_pX->totalCollideCrossSection);
	}else{
		if(e_debugIgnoreCrossSection){
			fprintf(stderr, "DEBUG:%d:%s + %s IgnoreCrossSection by target %lg < detectLimit %lg [mol]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetAtomMol, a_pX->electrodePtr->detectLimitMolForIsotope);
		}
		if(e_logIgnoreCrossSection){
			fprintf(e_logFp, "[INFO]:%d:%s + %s IgnoreCrossSection by target %lg < detectLimit %lg [mol]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetAtomMol, a_pX->electrodePtr->detectLimitMolForIsotope);
		}
	}
}
extern void debugIgnoreCrossSection(struct collide * a_pX)
{
	if(e_debugIgnoreCrossSection){
		fprintf(stderr, "DEBUG:%d:%s + %s IgnoreCrossSection by massDefect %lg <= 0.0 [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->massDefectMeV[0]);
	}
	if(e_logIgnoreCrossSection){
		fprintf(e_logFp, "[INFO]:%d:%s + %s IgnoreCrossSection by massDefect %lg <= 0.0 [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->massDefectMeV[0]);
	}
}
extern int iterateCrossSection(void * a_total, struct objectNodeConst * a_nodePtr)
{
	struct collide * pX = (struct collide *)a_total;
	pX->targetIsotopePropertyPtr = ((struct atomKey *)a_nodePtr->keyPtr)->isotopePropertyPtr;
	pX->targetAtomValuePtr = (struct atomValue *)a_nodePtr->valuePtr;
	//if(pX->collideType == COLLIDE_ELECTRON){
	//	fprintf(stderr, "DEBUG:%s:{BEGIN:%p %s\n", __FUNCTION__, pX, pX->targetIsotopePropertyPtr->symbol);
	//}
	if(checkCollidingUsage(pX)){
		getCollidedNewIsotopePropertyPtr(pX);
		sumTotalCollideCrossSection(pX);
		/* OLD
		if(pX->newIsotopePropertyPtr[0]){
			if(pX->collideType == COLLIDE_ELECTRON || pX->collideType == COLLIDE_ELECTRON_S){
				//if(pX->massDefectMeV[0] > 0.0){
					sumTotalCollideCrossSection(pX);
				//}else{
				//	debugIgnoreCrossSection(pX);
				//}
				//if(pX->collideType == COLLIDE_ELECTRON){
				//	fprintf(stderr, "DEBUG:%s(%d):%p ELECTRON+target:%s new:%s massDefectMeV[0]:%lg totalCollideCrossSection:%lg\n", __FUNCTION__, __LINE__, pX, pX->targetIsotopePropertyPtr->symbol, pX->newIsotopePropertyPtr[0]->symbol, pX->massDefectMeV[0], pX->totalCollideCrossSection);
				//}
			}else if(pX->collideType == COLLIDE_PROTON || pX->collideType == COLLIDE_PROTON_S
			|| pX->collideType == COLLIDE_DEUTERIUM || pX->collideType == COLLIDE_DEUTERIUM_S
			|| pX->collideType == COLLIDE_TRITIUM || pX->collideType == COLLIDE_TRITIUM_S
			|| pX->collideType == COLLIDE_ALPHA_S){
				sumTotalCollideCrossSection(pX);
			}else{//COLLIDE_NEUTRON
				//if(pX->massDefectMeV[0] > 0.0){
					sumTotalCollideCrossSection(pX);
				//}else{
				//	debugIgnoreCrossSection(pX);
				//}
			}
		}else{
			sumTotalCollideCrossSection(pX);
		}
		*/
	}
	//if(pX->collideType == COLLIDE_ELECTRON){
	//	fprintf(stderr, "DEBUG:%s:}END:totalCollideCrossSection:%lg\n", __FUNCTION__, pX->totalCollideCrossSection);
	//}
	return KEEP_NODE;
}
int e_debugCollideCrossSectionRate = 0;
int e_logCollideCrossSectionRate = 0;
extern double calcCollideCrossSectionRate(struct collide * a_pX)
{
	double collideCrossSectionRate;
	if(a_pX->totalCollideCrossSection > 0.0){
		collideCrossSectionRate = getCrossSection(a_pX) / a_pX->totalCollideCrossSection;
		if(e_debugCollideCrossSectionRate){
			fprintf(stderr, "DEBUG:%d:%s + %s, CrossSectionRate %lg\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, collideCrossSectionRate);
		}
		if(e_logCollideCrossSectionRate){
			fprintf(e_logFp, "[INFO]:%d:%s + %s, CrossSectionRate %lg\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, collideCrossSectionRate);
		}
	}else{
		collideCrossSectionRate = 0.0;
	}
	return collideCrossSectionRate;
}
int e_debugCollide = 0;
int e_LogCollide = 0;
extern int calcTargetAtomMolAndNewIsotopeMol(struct collide * a_pX)
{
	int iRet = 0;
	a_pX->targetAtomMol = getMol(a_pX->targetAtomValuePtr);
	if(a_pX->targetAtomMol >= a_pX->electrodePtr->detectLimitMolForIsotope){
		double collideCrossSectionRate = calcCollideCrossSectionRate(a_pX);
		a_pX->newIsotopeMol = a_pX->collideBulletMol * collideCrossSectionRate;
		if(a_pX->newIsotopeMol < 0.0){
			a_pX->newIsotopeMol = 0.0;//collect tolelance.
		}
		if(a_pX->newIsotopeMol > a_pX->targetAtomMol){
			//if(a_pX->newIsotopeMol > a_pX->targetAtomMol * 8.0){
			//	fprintf(e_logFp, "WARN:%s(%d):%s newIsotope %s Mol:%lg > targetAtom %s mol:%lg\n",
			//		__FUNCTION__, __LINE__, getCollideName(a_pX->collideType), 
			//		a_pX->newIsotopePropertyPtr[0]->symbol, a_pX->newIsotopeMol, 
			//		a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetAtomMol);
			//}
			//if(a_pX->newIsotopeMol > a_pX->targetAtomMol * 4.0){
			//	fprintf(stderr, "WARN:%s(%d):%s newIsotope %s Mol:%lg > targetAtom %s mol:%lg\n",
			//		__FUNCTION__, __LINE__, getCollideName(a_pX->collideType), 
			//		a_pX->newIsotopePropertyPtr[0]->symbol, a_pX->newIsotopeMol, 
			//		a_pX->targetIsotopePropertyPtr->symbol, a_pX->targetAtomMol);
			//}
			a_pX->newIsotopeMol = a_pX->targetAtomMol;
		}
		iRet = 1;
	}
	return 	iRet;
}
extern void debugScatVoltage(struct collide * a_pX, double a_scatVoltage)
{
	if(e_debugCollide){
		fprintf(stderr, "DEBUG:%d:%s + %s, SCAT NEXT appliedVoltage %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_scatVoltage);
	}
	if(e_LogCollide){
		fprintf(e_logFp, "[INFO]:%d:%s + %s, SCAT NEXT appliedVoltage %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_scatVoltage);
	}
}
extern void scatBullet(struct collide * a_pX)
{
	int sumMassNumber = a_pX->bulletPropertyPtr->massNumber + a_pX->targetIsotopePropertyPtr->massNumber;
	if(sumMassNumber > 0){
		int scatType;
		double v1Scale, scatVoltage, gammaVoltage;
		if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_S){
			fprintf(stderr, "FATAL:%s(%d):a_pX->collideType:%s\n", __FUNCTION__, __LINE__, getCollideName(a_pX->collideType));
			exit(1);
		}else if(a_pX->collideType == COLLIDE_PROTON || a_pX->collideType == COLLIDE_PROTON_S){
			scatType = MASS_DEFECT_BY_PROTON;
		}else if(a_pX->collideType == COLLIDE_DEUTERIUM || a_pX->collideType == COLLIDE_DEUTERIUM_S){
			scatType = MASS_DEFECT_BY_DEUTERIUM;
		}else if(a_pX->collideType == COLLIDE_TRITIUM || a_pX->collideType == COLLIDE_TRITIUM_S){
			scatType = MASS_DEFECT_BY_TRITIUM;
		}else if(a_pX->collideType == COLLIDE_ALPHA_S){
			scatType = MASS_DEFECT_BY_ALPHA;
		}else if(a_pX->collideType == COLLIDE_NEUTRON){
			fprintf(stderr, "FATAL:%s(%d):a_pX->collideType:%s\n", __FUNCTION__, __LINE__, getCollideName(a_pX->collideType));
			exit(1);
		}
		/* Perfect elastic collision of phisics.
		There are two mass points.
		The first mass point has mass, m1, and velocity v1.
		The second mass point has has mass, m2, and velocity v2 = 0.
		We think that the first mass point collides the second mass point.
		The velocity after collision are v1’, v2’ for each mass point.
		The momentum is preserved.
			m1 v1  = m1 v1’ + m2 v2’
		The kinetic energy of two mass points is preserved.
			(1/2)m1 (v1)^2 = (1/2)m1 (v1’)^2  + (1/2)m2 (v2’)^2 
		We can calculate v1’, v2’:
			v1’ = v1 and v2’ = 0 (Trivial solution)
			or
			v1’ = ((m1 - m2)/(m1 + m2))v1 and v2’ =  (2 m1/(m1 + m2))v1
		We can know:
		if m2 << m1 <=> m2 = almost 0, v1’ = almost v1 and v2’ = almost 2 v1.
		if m1 almost = m2, v1’ = almost 0 and v2’ = almost v1.
		if m1 << m2 <=> m1 = almost 0, v1’ = almost - v1 and v2’ = almost 0.
		We can also know:
			((1/2)m1 (v1’)^2) / ((1/2)m1 (v1)^2) 
			= (v1’/v1)^2
			= ((m1 - m2)/(m1 + m2))^2
		*/
		if(a_pX->targetIsotopePropertyPtr->massNumber > 0){
			v1Scale = (double)(a_pX->bulletPropertyPtr->massNumber - a_pX->targetIsotopePropertyPtr->massNumber) / (double)sumMassNumber;
		}else{
			v1Scale = (sumMassNumber - (e_massElectronMeV / e_massProtonMeV)) / sumMassNumber;
		}
		scatVoltage = a_pX->appliedVoltageMeV * v1Scale * v1Scale;
		gammaVoltage = a_pX->appliedVoltageMeV - scatVoltage;
		registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, scatType, scatVoltage, a_pX->newIsotopeMol);
		registOutput(a_pX->electrodePtr, SCAT_NOT_COLLIDE, MASS_DEFECT_BY_GANMMA, gammaVoltage, a_pX->newIsotopeMol);
		debugScatVoltage(a_pX, scatVoltage);
	}else{
		fprintf(stderr, "FATAL:%s(%d):sumMassNumber:%d\n", __FUNCTION__, __LINE__, sumMassNumber);
		exit(1);
	}
}

extern double collideParticleAtom(struct collide * a_pX)
{
	//The collision of a particle like an electron, a proton, a deterium or a tritium, it may have kinetic energy,  it will be absorbed by neuclaus.
	int i;
	double collidedMol = 0.0;
	//if(a_pX->collideType == COLLIDE_ELECTRON){
	//	fprintf(stderr, "DEBUG:%s:{BEGIN %s\n", __FUNCTION__, a_pX->targetIsotopePropertyPtr->symbol);
	//}
	if(checkCollidingUsage(a_pX)){
		getCollidedNewIsotopePropertyPtr(a_pX);
		if(calcTargetAtomMolAndNewIsotopeMol(a_pX)){
			if(a_pX->newIsotopePropertyPtr[0]){
				if(a_pX->appliedVoltageMeV >= a_pX->electroPotentialMeV){
					if(e_debugCollide){
						fprintf(stderr, "DEBUG:%d:%s + %s, appliedVoltage %lg >= electroPotential %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
					}
					if(e_LogCollide){
						fprintf(e_logFp, "[INFO]:%d:%s + %s, appliedVoltage %lg >= electroPotential %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
					}
					if(a_pX->massDefectMeV[0] > 0.0){
						if(e_debugCollide){
							fprintf(stderr, "DEBUG:%d:%s + %s, massDefect %lg > 0.0 [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->massDefectMeV[0]);
						}
						if(e_LogCollide){
							fprintf(e_logFp, "[INFO]:%d:%s + %s, massDefect %lg > 0.0 [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->massDefectMeV[0]);
						}
						for(i = 0; i < a_pX->newCount; ++i){
							double newIsotopeMolRate = a_pX->newIsotopeMol * a_pX->rate[i];
							//if(a_pX->collideType == COLLIDE_ELECTRON){//DEBUG
								//fprintf(stderr, "DEBUG:%s(%d):%s (i:%d/newCount:%d) newIsotopeMol:%lg, rate:%lg, newIsotopeMolRate:%lg detectLimitMolForIsotope:%lg\n", __FUNCTION__, __LINE__, a_pX->newIsotopePropertyPtr[0]->symbol, i, a_pX->newCount, a_pX->newIsotopeMol, a_pX->rate[i], newIsotopeMolRate, a_pX->electrodePtr->detectLimitMolForIsotope);
							//}
							if(newIsotopeMolRate >= a_pX->electrodePtr->detectLimitMolForIsotope){
								if(a_pX->emissionType[i] == EMISSION_IS_2_ELECTRONS){
									fprintf(stderr, "FATAL:%s(%d):a_pX->collideType:%s EMISSION_IS_2_ELECTRONS i:%d\n", __FUNCTION__, __LINE__, getCollideName(a_pX->collideType), i);
									exit(1);
								}else{
									registNuclearReaction(a_pX->electrodePtr, REACTION_DETECT, a_pX->collideType, a_pX->nuclearReactionStr[i], a_pX->massDefectMeV[i], newIsotopeMolRate, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
									if(a_pX->emissionType[i] == EMISSION_IS_GAMMA){
										registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, a_pX->massDefectMeV[i], newIsotopeMolRate);
									}else if(a_pX->emissionType[i] == EMISSION_IS_ELECTRO_NEUTRINO){
										registOutputOfElectronCapture(a_pX->electrodePtr, a_pX->massDefectMeV[i], newIsotopeMolRate);
									}else if(a_pX->emissionType[i] == EMISSION_IS_PROTON){
										registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, a_pX->massDefectMeV[i], newIsotopeMolRate);
									}else if(a_pX->emissionType[i] == EMISSION_IS_NEUTRON){
										registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, a_pX->massDefectMeV[i], newIsotopeMolRate);
									}else{
										fprintf(stderr, "FATAL:%s(%d):invalid emissionType:%d\n", __FUNCTION__, __LINE__, a_pX->emissionType[i]);
										exit(1);
									}
									collidedMol += newIsotopeMolRate;
									setMolSub(a_pX->targetAtomValuePtr, newIsotopeMolRate);
									if(a_pX->newAtomicNumber[i] == ATOMICNUMBER_NEUTRON 
									&& (a_pX->newMassNumber[i] == MASSNUMBER_DINEUTRON 
									|| a_pX->newMassNumber[i] == MASSNUMBER_TRINEUTRON)){
										// The dineutron and trineutron will imediately separate into neutrons.
										newIsotopeMolRate *= a_pX->newMassNumber[i];
										a_pX->newMassNumber[i] = MASSNUMBER_NEUTRON;
										a_pX->newIsotopePropertyPtr[i] = getIsotopePropertyPtr(ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON);
									}
									increaseAtom(a_pX->electrodePtr, a_pX->newAtomicNumber[i], a_pX->newMassNumber[i], newIsotopeMolRate);
									if(a_pX->subIsotopePropertyPtr[i]){
										increaseAtom(a_pX->electrodePtr, a_pX->subAtomicNumber[i], a_pX->subMassNumber[i], newIsotopeMolRate);
									}
								}
							}else if(newIsotopeMolRate > 0.0){
								registNuclearReaction(a_pX->electrodePtr, REACTION_CANT_DETECT, a_pX->collideType, a_pX->nuclearReactionStr[i], a_pX->massDefectMeV[i], newIsotopeMolRate, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
								//fprintf(stderr, "DEBUG:%s:%s newIsotopeMolRate:%lg < detectLimitMolForIsotope:%lg\n %s\n",
								//__FUNCTION__, a_pX->newIsotopePropertyPtr->symbol, 
								//newIsotopeMolRate, a_pX->electrodePtr->detectLimitMolForIsotope, a_pX->nuclearReactionStr);
							}else{ //newIsotopeMolRate <= 0.0
								;//Do nothing!
							}
						}
					}else{//a_pX->massDefectMeV[0] <= 0.0
						registNuclearReaction(a_pX->electrodePtr, REACTION_ENDOTHERMIC, a_pX->collideType, a_pX->nuclearReactionStr[0], a_pX->massDefectMeV[0], 0.0, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
						if(e_debugCollide){
							fprintf(stderr, "DEBUG:%d:%s + %s, massDefect %lg <= 0.0 [MeV] REACTION_ENDOTHERMIC\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->massDefectMeV[0]);
						}
						if(e_LogCollide){
							fprintf(e_logFp, "[INFO]:%d:%s + %s, massDefect %lg <= 0.0 [MeV] REACTION_ENDOTHERMIC\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->massDefectMeV[0]);
						}
						if(a_pX->collideType == COLLIDE_ELECTRON || a_pX->collideType == COLLIDE_ELECTRON_S){
							if(a_pX->emissionType[0] == EMISSION_IS_2_ELECTRONS){
								fprintf(stderr, "FATAL:%s(%d):a_pX->collideType:%s EMISSION_IS_2_ELECTRONS\n", __FUNCTION__, __LINE__, getCollideName(a_pX->collideType));
								exit(1);
							}else{
								double scatVoltage, gammaVoltage;
								scatVoltage = a_pX->appliedVoltageMeV * (a_pX->targetIsotopePropertyPtr->massNumber - (e_massElectronMeV / e_massProtonMeV)) / a_pX->targetIsotopePropertyPtr->massNumber;
								gammaVoltage = a_pX->appliedVoltageMeV - scatVoltage;
								registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_BETA, scatVoltage, a_pX->newIsotopeMol);
								registOutput(a_pX->electrodePtr, SCAT_NOT_COLLIDE, MASS_DEFECT_BY_GANMMA, gammaVoltage, a_pX->newIsotopeMol);
								debugScatVoltage(a_pX, scatVoltage);
								collidedMol += a_pX->newIsotopeMol;
							}
						}else if(a_pX->collideType == COLLIDE_NEUTRON){
							if(a_pX->newIsotopeMol > 0.0){
								fprintf(stderr, "FATAL:%s(%d):a_pX->collideType:%s a_pX->newIsotopeMol:%lg > 0.0\n", __FUNCTION__, __LINE__, getCollideName(a_pX->collideType), a_pX->newIsotopeMol);
								exit(1);
							}else{
								;//DO NOTHING!
							}
						}else{
							scatBullet(a_pX);
							collidedMol += a_pX->newIsotopeMol;
						}
					}
				}else{//a_pX->appliedVoltageMeV < a_pX->electroPotentialMeV
					if(e_debugCollide){
						fprintf(stderr, "DEBUG:%d:%s + %s, appliedVoltage %lg < electroPotential %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
					}
					if(e_LogCollide){
						fprintf(e_logFp, "[INFO]:%d:%s + %s, appliedVoltage %lg < electroPotential %lg [MeV]\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
					}
					if(a_pX->emissionType[0] == EMISSION_IS_2_ELECTRONS){
						//The bullet electrons lost the half of energy and pass it to target electrons.
						registOutput(a_pX->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_BETA, a_pX->appliedVoltageMeV * 0.5, a_pX->newIsotopeMol * 2.0);
						debugScatVoltage(a_pX, a_pX->appliedVoltageMeV * 0.5);
						collidedMol += a_pX->newIsotopeMol;
					}else{
						registNuclearReaction(a_pX->electrodePtr, REACTION_COULOMB_BARRIER, a_pX->collideType, a_pX->nuclearReactionStr[0], a_pX->massDefectMeV[0], 0.0, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
						scatBullet(a_pX);
						collidedMol += a_pX->newIsotopeMol;
					}
				}
			}else{//!a_pX->newIsotopePropertyPtr[0]
				if(a_pX->nuclearReactionStr[0][0]){
					registNuclearReaction(a_pX->electrodePtr, REACTION_ERROR, a_pX->collideType, a_pX->nuclearReactionStr[0], a_pX->massDefectMeV[0], a_pX->targetAtomMol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);
					//fprintf(stderr, "%s\n", a_pX->nuclearReactionStr);
				}
				if(e_debugCollide){
					fprintf(stderr, "DEBUG:%d:%s + %s -> REACTION_ERROR\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol);
				}
				if(e_LogCollide){
					fprintf(e_logFp, "[INFO]:%d:%s + %s -> REACTION_ERROR\n", __LINE__, a_pX->bulletPropertyPtr->symbol, a_pX->targetIsotopePropertyPtr->symbol);
				}
			}
		}else{
			if(a_pX->nuclearReactionStr[0][0]){
				registNuclearReaction(a_pX->electrodePtr, REACTION_CANT_DETECT, a_pX->collideType, a_pX->nuclearReactionStr[0], a_pX->massDefectMeV[0], a_pX->targetAtomMol, a_pX->appliedVoltageMeV, a_pX->electroPotentialMeV);		
			}
		}
	}
	//if(a_pX->collideType == COLLIDE_ELECTRON){
	//	fprintf(stderr, "DEBUG:%s:}END:collidedMol:%lg\n", __FUNCTION__, collidedMol);
	//}
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
int e_timesOfCollide[2][COLLIDE_ELECTRON_S + 1];

extern void collideBulletToElectrode(struct electrode * a_electrodePtr, int a_collideType, double a_appliedVoltageMeV, double a_arriveBulletMol, double a_collideBulletRate, struct atomNodeConst * a_decreaseBulletPtr, int a_nest)
{
	struct collide col;
	int logCrossSectionTimes = 100, remTimes;//DEBUG
	int electrodeId;//DEBUG
	//fprintf(stderr, "DEBUG:%s:{BEGIN a_collideType:%s a_appliedVoltageMeV:%lg a_arriveBulletMol:%lg a_collideBulletRate:%lg\n", __FUNCTION__, getCollideName(a_collideType), a_appliedVoltageMeV, a_arriveBulletMol, a_collideBulletRate);
	//e_cntCollideBulletToElectrode++;//DEBUG
	memset(&col, 0, sizeof(struct collide));
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	col.electrodePtr = a_electrodePtr;
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	col.collideType = a_collideType;
	//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
	switch(col.collideType){//DEBUG
	case COLLIDE_ELECTRON: logCrossSectionTimes = 100; break;
	case COLLIDE_PROTON: logCrossSectionTimes = 100; break;
	case COLLIDE_DEUTERIUM: logCrossSectionTimes = 100; break;
	case COLLIDE_TRITIUM: logCrossSectionTimes = 100; break;
	case COLLIDE_NEUTRON: logCrossSectionTimes = 100; break;
	case COLLIDE_ELECTRON_S: logCrossSectionTimes = 100000; break;
	case COLLIDE_PROTON_S: logCrossSectionTimes = 10000; break;
	case COLLIDE_DEUTERIUM_S: logCrossSectionTimes = 10000; break;
	case COLLIDE_TRITIUM_S: logCrossSectionTimes = 10000; break;
	case COLLIDE_ALPHA_S: logCrossSectionTimes = 10000; break;
	}
	electrodeId = (a_electrodePtr == &e_negativeElectrode) ? 0 : 1;//DEBUG
	remTimes = e_timesOfCollide[electrodeId][a_collideType] % logCrossSectionTimes;//DEBUG
	if(remTimes < 0){//DEBUG
		e_debugElectroPotential = 1;
		e_logElectroPotential = 1;
		e_debugNewIsotope = 1;
		e_logNewIsotope = 1;
		e_debugIgnoreCrossSection = 1;
		e_logIgnoreCrossSection = 1;
		e_debugCollideCrossSectionRate = 1;
		e_logCollideCrossSectionRate = 1;
		e_debugCollide = 1;
		e_LogCollide = 1;
		fprintf(stderr, "DEBUG:%d:{%s %s appliedVoltage %lg [MeV] %d times \n", __LINE__, a_electrodePtr->atomHashTable.tableName, getCollideName(a_collideType), a_appliedVoltageMeV, e_timesOfCollide[electrodeId][a_collideType]);
		fprintf(e_logFp, "[INFO]:%d:{%s %s appliedVoltage %lg [MeV] %d times \n", __LINE__, a_electrodePtr->atomHashTable.tableName, getCollideName(a_collideType), a_appliedVoltageMeV, e_timesOfCollide[electrodeId][a_collideType]);
	}else{//DEBUG
		e_debugElectroPotential = 0;
		e_logElectroPotential = 0;
		e_debugNewIsotope = 0;
		e_logNewIsotope = 0;
		e_debugIgnoreCrossSection = 0;
		e_logIgnoreCrossSection = 0;
		e_debugCollideCrossSectionRate = 0;
		e_logCollideCrossSectionRate = 0;
		e_debugCollide = 0;
		e_LogCollide = 0;
	}
	switch(col.collideType){
		case COLLIDE_ELECTRON:
		case COLLIDE_ELECTRON_S:  col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON); break;
		case COLLIDE_PROTON:
		case COLLIDE_PROTON_S:    col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN); break;
		case COLLIDE_DEUTERIUM:
		case COLLIDE_DEUTERIUM_S: col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_DEUTERIUM, MASSNUMBER_DEUTERIUM);  break;
		case COLLIDE_TRITIUM:
		case COLLIDE_TRITIUM_S:   col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_TRITIUM, MASSNUMBER_TRITIUM); break;
		case COLLIDE_ALPHA_S:     col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM); break;
		case COLLIDE_NEUTRON:	  col.bulletPropertyPtr = getIsotopePropertyPtr(ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON); break;
	}
	if(col.collideType == COLLIDE_NEUTRON){
		if(a_appliedVoltageMeV != 0.0){
			fprintf(stderr, "FATAL ERROR:%s:%s a_appliedVoltageMeV:%lg!= 0.0\n", __FUNCTION__, getCollideName(col.collideType), a_appliedVoltageMeV);
			exit(1);
		}
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
					//fprintf(stderr, "DEBUG:%s:call iterateCrossSection i/size:%d/%d\n", __FUNCTION__, i, size);
					iterateCrossSection(&col, oPtr[i]);
				}
				//fprintf(stderr, "DEBUG:%s:inCompton:%d\n", __FUNCTION__, e_inComptonEffect);
				//fprintf(stderr, "DEBUG:%s:totalCollideCrossSection:%lg\n", __FUNCTION__, col.totalCollideCrossSection);
				col.printNuclearReaction = 1;
				if(e_debugCollide){
					fprintf(stderr, "DEBUG:%d:\n", __LINE__);
				}
				if(e_LogCollide){
					fprintf(e_logFp, "[INFO]:%d:\n", __LINE__);
				}
				//fprintf(stderr, "DEBUG:%s:iterateCrossSection usedCnt:%d\n", __FUNCTION__, col.electrodePtr->atomHashTable.usedCnt);
				for(i = 0; i < size; ++i){
					//fprintf(stderr, "DEBUG:%s:call iterateCollide i/size:%d/%d\n", __FUNCTION__, i, size);
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
	/* OLD
	if(col.remainBulletMol > col.electrodePtr->detectLimitMolForIsotope){
		if(a_nest < 8){
			collideBulletToElectrode(a_electrodePtr, a_collideType, a_appliedVoltageMeV, col.remainBulletMol, a_collideBulletRate, a_decreaseBulletPtr, a_nest + 1);
			col.remainBulletMol = 0.0;
		}else{
			fprintf(stderr, "WARN:%s(%d):a_nest:%d reached to the limit at remainBulletMol:%lg a_collideType:%s\n", __FUNCTION__, __LINE__, a_nest, col.remainBulletMol, getCollideName(a_collideType));
		}
	}
	*/
	col.imperfectCollideMol += col.remainBulletMol;
	//fprintf(stderr, "DEBUG:%s:col.imperfectCollideMol:%lg\n", __FUNCTION__, col.imperfectCollideMol);
	if(col.imperfectCollideMol > 0.0){
		if(col.appliedVoltageMeV > 0.0){
			registOutput(col.electrodePtr, SCAT_IMPERFECT, MASS_DEFECT_BY_GANMMA, col.appliedVoltageMeV, col.imperfectCollideMol);
		}
		if(col.collideType == COLLIDE_PROTON || col.collideType == COLLIDE_DEUTERIUM || col.collideType == COLLIDE_TRITIUM){
			//When the protons, deuteriums and tritiums flew from the positive electrode to the negative electorode with huge energy greater than the energy of beta decay, some of them collide imperfectly, so they will change from neurcuses to atoms.
			//We need to append the isotopes of hydrogen onto the negative electorode.
			increaseAtom(col.electrodePtr, col.bulletPropertyPtr->atomicNumber, col.bulletPropertyPtr->massNumber, col.imperfectCollideMol);
		}else if(col.collideType == COLLIDE_ELECTRON){
			;//The electrons flew from the negative electrode to the positive electorode. But the electric circuit rotate electrons from the positive electrode to the negative electorode. So the amount of electrons is always same.
		}else  if(col.collideType == COLLIDE_NEUTRON){
			;//The neutrons do not fly from a electrode to another electrodes, because they do not have any electoric charge.
		}else if(col.collideType == COLLIDE_PROTON_S || col.collideType == COLLIDE_DEUTERIUM_S || col.collideType == COLLIDE_TRITIUM_S || col.collideType == COLLIDE_ALPHA_S || col.collideType == COLLIDE_ELECTRON_S){
			;//The sacttering particles are in the electrode, they do not fly from a electrode to another electrodes.
		}
	}
	if(remTimes < 0){//DEBUG
		fprintf(stderr, "DEBUG:%d:}%s %s %d times appliedVoltage %lg [MeV]\n", __LINE__, a_electrodePtr->atomHashTable.tableName, getCollideName(a_collideType), e_timesOfCollide[electrodeId][a_collideType], a_appliedVoltageMeV);
		fprintf(e_logFp, "[INFO]:%d:}%s %s %d times appliedVoltage %lg [MeV]\n", __LINE__, a_electrodePtr->atomHashTable.tableName, getCollideName(a_collideType), e_timesOfCollide[electrodeId][a_collideType], a_appliedVoltageMeV);
	}
	e_timesOfCollide[electrodeId][a_collideType]++;//DEBUG
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
					registOutput(col.electrodePtr, SCAT_IMPERFECT_IN_SPACE, MASS_DEFECT_BY_GANMMA, col.appliedVoltageMeV, col.imperfectCollideMol);
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
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, massDefectMeV * (1.0 - e_rateForAlphaParticle), partialMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_ALPHA, massDefectMeV * e_rateForAlphaParticle, partialMol);
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
									
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * (1.0 - e_rateForProtonAtBetaPlus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * e_rateForProtonAtBetaPlus, partialMol);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_2PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 2.0);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor2ProtonAtBetaPlus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * e_rateFor2ProtonAtBetaPlus * 0.5 , partialMol * 2.0);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_3PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 3.0);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor3ProtonAtBetaPlus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * e_rateFor3ProtonAtBetaPlus / 3.0, partialMol * 3.0);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_PLUS_ALPHA){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM, partialMol);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * (1.0 - e_rateForAlphaParticleAtBetaPlus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_ALPHA, massDefectMeV * e_rateForAlphaParticleAtBetaPlus, partialMol);
								}
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
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * (1.0 - e_rateForProtonAtEC), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * e_rateForProtonAtEC, partialMol);
								}else if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_2PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 2.0);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor2ProtonAtEC), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * e_rateFor2ProtonAtEC * 0.5, partialMol * 2.0);
								}else if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_3PROTON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HYDROGEN, MASSNUMBER_HYDROGEN, partialMol * 3.0);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor3ProtonAtEC), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * e_rateFor3ProtonAtEC / 3.0, partialMol * 3.0);
								}else if(decayMode == DECAY_MODE_ELECTRON_CAPTURE_ALPHA
								|| decayMode == DECAY_MODE_ELECTRON_CAPTURE_DEGRADE_ALPHA){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_HELIUM, MASSNUMBER_HELIUM, partialMol);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + zMeV - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * (1.0 - e_rateForAlphaParticleAtEC), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_ALPHA, massDefectMeV * e_rateForAlphaParticleAtEC, partialMol);
								}
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_DOUBLE_BETA_PLUS
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_DOUBLE_ELECTRON_CAPTURE){
						// decrease atomic number, same mass number
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 2;
						daughterMassNumber = a_isotopePropertyPtr->massNumber;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massElectronMeV + e_massElectronMeV + e_massElectronMeV + e_massElectronMeV);
							if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_DOUBLE_BETA_PLUS 
							&& massDefectMeV > 0.0){
								//The beta plus is as same as the positron(e+) emission
								//The Nuclear Reaction is : X           -> Y + e- + e- + e+ + e+ + νe + νe
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s(%lg) -> %s(%lg) + e-(%lg) + e+(%lg) + νe(0) + %lg[MeV/c^2] + e-(%lg) + e+(%lg) + νe(0) + %lg[MeV/c^2]", 
									a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
									daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
									e_massElectronMeV, e_massElectronMeV, massDefectMeV * 0.5,
									e_massElectronMeV, e_massElectronMeV, massDefectMeV * 0.5);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massElectronMeV + e_massElectronMeV + e_massElectronMeV + e_massElectronMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, - partialMol * 2.0);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * 0.5, partialMol);
								registOutputOfBetaPlusWithPositron(a_electrodePtr, massDefectMeV * 0.5, partialMol);
							}else if((massDefectMeV = a_isotopePropertyPtr->massMeV - daughterIsotopePropertyPtr->massMeV) > 0.0){
								//The electron capture do not emmit a positron(e+), only emmits a neutorino.
								//The Nuclear Reaction is : X+2 + e- + e-  -> Y + νe
								int decayMode = a_isotopePropertyPtr->decayMode[i];
								if(decayMode == DECAY_MODE_DOUBLE_BETA_PLUS){
									decayMode = DECAY_MODE_DOUBLE_ELECTRON_CAPTURE_DEGRADE;
								}
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s2+(%lg) + e-(%lg) + e-(%lg) -> %s(%lg) + νe(0) + %lg[MeV/c^2] + νe(0) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV - e_massElectronMeV - e_massElectronMeV,
								e_massElectronMeV, e_massElectronMeV,
								daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
								massDefectMeV * 0.5, massDefectMeV * 0.5);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - e_massElectronMeV - e_massElectronMeV + e_massElectronMeV + e_massElectronMeV - daughterIsotopePropertyPtr->massMeV, massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, decayMode, nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, - partialMol * 2.0);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * 0.5, partialMol);
								registOutputOfElectronCapture(a_electrodePtr, massDefectMeV * 0.5, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					/* 
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_DOUBLE_ELECTRON_CAPTURE){
						// decrease atomic number, same mass number
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - 2;
						daughterMassNumber = a_isotopePropertyPtr->massNumber;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = (a_isotopePropertyPtr->massMeV - daughterIsotopePropertyPtr->massMeV) * 0.5;
							if(massDefectMeV > 0.0){
								//The DOUBLE electron capture do not emmit a positron(e+), only emmits two neutorinos.
								//The Nuclear Reaction is : X2+ + e- + e-  -> Y + νe + νe
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s2+(%lg) + e-(%lg) + e-(%lg) -> %s(%lg) + νe(0) + %lg[MeV/c^2] + νe(0) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, 
								a_isotopePropertyPtr->massMeV - e_massElectronMeV - e_massElectronMeV,
								e_massElectronMeV, e_massElectronMeV,
								daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
								massDefectMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - e_massElectronMeV  - e_massElectronMeV + e_massElectronMeV + e_massElectronMeV - daughterIsotopePropertyPtr->massMeV, massDefectMeV + massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV + massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, -partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, -partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutputOfElectronCapture(a_electrodePtr, massDefectMeV, partialMol);
								registOutputOfElectronCapture(a_electrodePtr, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
						*/
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
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_NEUTRON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_2NEUTRON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_3NEUTRON
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_4NEUTRON
					){
						int neutronCnt;
						char * neutronTxtPtr;
						char * scalePtr;
						switch(a_isotopePropertyPtr->decayMode[i]){
							case DECAY_MODE_BETA_MINUS_AND_NEUTRON: neutronCnt = 1; neutronTxtPtr = "n"; scalePtr = ""; break;
							case DECAY_MODE_BETA_MINUS_AND_2NEUTRON: neutronCnt = 2; neutronTxtPtr = "n2"; scalePtr = " * 2"; break;
							case DECAY_MODE_BETA_MINUS_AND_3NEUTRON: neutronCnt = 3; neutronTxtPtr = "n3"; scalePtr = " * 3"; break;
							case DECAY_MODE_BETA_MINUS_AND_4NEUTRON: neutronCnt = 4; neutronTxtPtr = "n4"; scalePtr = " * 4"; break;
						}
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber + 1;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - neutronCnt;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							char * daughterSymbolPtr; 
							if(daughterAtomicNumber == ATOMICNUMBER_HYDROGEN && daughterMassNumber == MASSNUMBER_HYDROGEN){
								daughterSymbolPtr = daughterIsotopePropertyPtr->nucleusName;
							}else{
								daughterSymbolPtr = daughterIsotopePropertyPtr->symbol;
							}
							double daughterIonMassMeV = daughterIsotopePropertyPtr->massMeV - e_massElectronMeV;
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * neutronCnt);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
								"%s(%lg) -> %s+(%lg) + e-(%lg) + ~νe(0) + %s(%lg%s) + %lg[MeV/c^2]", 
								a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
								daughterSymbolPtr, daughterIonMassMeV,
								e_massElectronMeV,
								neutronTxtPtr, e_massNeutronMeV, scalePtr,
								massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIonMassMeV + e_massElectronMeV + e_massNeutronMeV * neutronCnt), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_ELECTRON, MASSNUMBER_ELECTRON, partialMol);
								//The mass defect is devided into Beta Minus and Neutrons, But the rate is fixed just half.
								if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_NEUTRON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * neutronCnt - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * (1.0 - e_rateForNeytonAtBetaMinus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, massDefectMeV * e_rateForNeytonAtBetaMinus, partialMol);
									
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_2NEUTRON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol * neutronCnt);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * neutronCnt - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor2NeytonAtBetaMinus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, massDefectMeV * e_rateFor2NeytonAtBetaMinus / neutronCnt, partialMol * neutronCnt);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_3NEUTRON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol * neutronCnt);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * neutronCnt - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor3NeytonAtBetaMinus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, massDefectMeV * e_rateFor3NeytonAtBetaMinus / neutronCnt, partialMol * neutronCnt);
								}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_BETA_MINUS_AND_4NEUTRON){
									increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol * neutronCnt);
									//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * neutronCnt - a_isotopePropertyPtr->massMeV) * partialMol;
									registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * (1.0 - e_rateFor4NeytonAtBetaMinus), partialMol);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, massDefectMeV * e_rateFor4NeytonAtBetaMinus / neutronCnt, partialMol * neutronCnt);
								}
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
								registOutputOfBetaMinus(a_electrodePtr, massDefectMeV * (1.0 - e_rateForAlphaParticleAtBetaMinus), partialMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_ALPHA, massDefectMeV * e_rateForAlphaParticleAtBetaMinus, partialMol);
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
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV, partialMol);
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
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_PROTON, massDefectMeV * 0.5, partialMol * 2.0);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_NEUTRON_EMISSION){
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
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_2NEUTRON_EMISSION
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_3NEUTRON_EMISSION
					|| a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_4NEUTRON_EMISSION){
						int ncnt = 2;
						switch(a_isotopePropertyPtr->decayMode[i]){
						case DECAY_MODE_2NEUTRON_EMISSION: ncnt = 2; break;
						case DECAY_MODE_3NEUTRON_EMISSION: ncnt = 3; break;
						case DECAY_MODE_4NEUTRON_EMISSION: ncnt = 4; break;
						}
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - ncnt;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * ncnt);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + %d n(%d * %lg) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										ncnt, ncnt,
										e_massNeutronMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV * ncnt), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_NEUTRON, MASSNUMBER_NEUTRON, partialMol * ncnt);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_massNeutronMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_NEUTRON, massDefectMeV / ncnt, partialMol * ncnt);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_CLUSTER_DECAY_14C){
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - ATOMICNUMBER_CARBON;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - MASSNUMBER_CARBON14;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_C14Ptr->massMeV);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + 14C(%lg) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										e_C14Ptr->massMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_C14Ptr->massMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_CARBON, MASSNUMBER_CARBON14, partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_C14Ptr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
					}else if(a_isotopePropertyPtr->decayMode[i] == DECAY_MODE_SELF_FISSION_80KR){
						daughterAtomicNumber = a_isotopePropertyPtr->atomicNumber - ATOMICNUMBER_Kr;
						daughterMassNumber = a_isotopePropertyPtr->massNumber - MASSNUMBER_Kr80;
						daughterIsotopePropertyPtr = getIsotopePropertyPtr(daughterAtomicNumber, daughterMassNumber);
						if(daughterIsotopePropertyPtr){
							massDefectMeV = a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_Kr80Ptr->massMeV);
							if(massDefectMeV > 0.0){
								snprintf(nuclearReactionStr, REACTION_LEN,
									"%s(%lg) -> %s(%lg) + 80Kr(%lg) + %lg[MeV/c^2]", 
										a_isotopePropertyPtr->symbol, a_isotopePropertyPtr->massMeV,
										daughterIsotopePropertyPtr->symbol, daughterIsotopePropertyPtr->massMeV,
										e_Kr80Ptr->massMeV, massDefectMeV);
								checkMassInNuclearReaction(a_isotopePropertyPtr->massMeV - (daughterIsotopePropertyPtr->massMeV + e_Kr80Ptr->massMeV), massDefectMeV, nuclearReactionStr);
								registNuclearReaction(a_electrodePtr, REACTION_DETECT, a_isotopePropertyPtr->decayMode[i], nuclearReactionStr, massDefectMeV, partialMol, 0.0, 0.0);
								increaseAtom(a_electrodePtr, daughterAtomicNumber, daughterMassNumber, partialMol);
								increaseAtom(a_electrodePtr, ATOMICNUMBER_Kr, MASSNUMBER_Kr80, partialMol);
								//e_debugDecayMassDiff += (daughterIsotopePropertyPtr->massMeV + e_Kr80Ptr->massMeV - a_isotopePropertyPtr->massMeV) * partialMol;
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, massDefectMeV, partialMol);
							}else{
								showErrorNegativeMassDefect = 1;
							}
						}else{
							showErrorOfUndefinedDaughter = 1;
						}
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
		collideBulletToElectrode(a_electrodePtr, COLLIDE_NEUTRON, neutronVoltageMeV, neutronMol, collideNeutronRate, neutronPtr, 1);
	}
	//fprintf(stderr, "DEBUG:%s:}END\n", __FUNCTION__);
}
#if  KEEP_ATMIC_ORDER == 1
#else
struct decayBase {
	struct electrode * electrodePtr;
	double pastSecond;
};
extern int iterateDecayWithoutNeutron(void * a_total, struct objectNodeConst * a_nodePtr)
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
extern void decayNeuclayWithoutNeutron(struct electrode * a_electrodePtr, double a_pastSecond)
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
		iterateInHashTable(&a_electrodePtr->atomHashTable, &Db, iterateDecayWithoutNeutron);
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
	//if(a_MeV >= e_collideMiniMeV && a_mol >= a_electrodePtr->detectLimitMolForIsotope * SCATTER_CNT && a_nest < 16)
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
				collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_S, scatteredElectronMeV, molDivideBy6, e_collideElectronRateOnElectrode, electronPtr, 1);
				//fprintf(stderr, " LINE:%d} \n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}else if(scatteredElectronMeV >= e_collideMidiMeV){				
				//fprintf(stderr, " LINE:%d{ \n", __LINE__);
				collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_S, scatteredElectronMeV, molDivideBy6, e_collideElectronRateForMidiMeV, electronPtr, 1);
				electronPtr = findElectronInElectrode(a_electrodePtr);
				//fprintf(stderr, " LINE:%d} \n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}else if(scatteredElectronMeV >= e_collideMiniMeV){				
				//fprintf(stderr, " LINE:%d{ \n", __LINE__);
				collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_S, scatteredElectronMeV, molDivideBy6, e_collideElectronRateForMiniMeV, electronPtr, 1);
				electronPtr = findElectronInElectrode(a_electrodePtr);
				//fprintf(stderr, " LINE:%d} \n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}else if(scatteredElectronMeV > 0.0){
				//fprintf(stderr, " LINE:%d{ registOutput SCAT_SMALL_BY_COMPTON_E call\n", __LINE__);
				registOutput(a_electrodePtr, SCAT_SMALL_BY_COMPTON_E, MASS_DEFECT_BY_GANMMA, scatteredElectronMeV, molDivideBy6);
				//fprintf(stderr, " LINE:%d} registOutput SCAT_SMALL_BY_COMPTON_E return\n", __LINE__);
				//debugCheck += (scatteredElectronMeV * molDivideBy6);
			}
			if(scatteredPhotonMeV >= e_collideMiniMeV){
				registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, scatteredPhotonMeV, molDivideBy6);	
			}else if(scatteredPhotonMeV > 0.0){
				//fprintf(stderr, " LINE:%d{ registOutput SCAT_SMALL_BY_COMPTON call\n", __LINE__);
				registOutput(a_electrodePtr, SCAT_SMALL_BY_COMPTON, MASS_DEFECT_BY_GANMMA, scatteredPhotonMeV, molDivideBy6);
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
	//{ set by scatterIn
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
			registOutput(cleanUpPtr->electrodePtr, scat_degrade, MASS_DEFECT_BY_GANMMA, val, valuePtr->mol);
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
				registOutput(a_cleanUpPtr->electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, a_cleanUpPtr->md[j][k].MeV, a_cleanUpPtr->md[j][k].mol);
			}
		}
	}
}
extern int initMaxMinOfClenup(struct cleanUp * a_cleanUpPtr)
{
	int ret = 0;
	int j, k;

	a_cleanUpPtr->step[0] = (a_cleanUpPtr->maxMeV[0] - a_cleanUpPtr->minMeV[2]) / DEGRADE_DIV;
	if(a_cleanUpPtr->step[0] > 0.0){
		a_cleanUpPtr->minMeV[0] = a_cleanUpPtr->minMeV[2] + a_cleanUpPtr->step[0];
		a_cleanUpPtr->step[1] = (a_cleanUpPtr->minMeV[0] - a_cleanUpPtr->minMeV[2]) / (DEGRADE_DIV + 1);
		if(a_cleanUpPtr->step[1] > 0.0){
			a_cleanUpPtr->maxMeV[1] = a_cleanUpPtr->minMeV[0] - a_cleanUpPtr->step[1];
			a_cleanUpPtr->minMeV[1] = a_cleanUpPtr->minMeV[2] + a_cleanUpPtr->step[1];
			a_cleanUpPtr->step[2] = (a_cleanUpPtr->minMeV[1] - a_cleanUpPtr->minMeV[2]) / DEGRADE_DIV;
			if(a_cleanUpPtr->step[2] > 0.0){
				a_cleanUpPtr->maxMeV[2] = a_cleanUpPtr->minMeV[1] - a_cleanUpPtr->step[2];
				ret = 1;
				
				for(j = 0; j < DEGRADE_MAG; ++j){
					for(k = 0; k < DEGRADE_DIV; ++k){
						a_cleanUpPtr->md[j][k].MeV = (a_cleanUpPtr->maxMeV[j] * k + a_cleanUpPtr->minMeV[j] * (DEGRADE_DIV - 1 - k)) / (DEGRADE_DIV - 1);
						a_cleanUpPtr->md[j][k].mol = 0.0;
					}
				}
			}
		}
	}
	return ret;
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
				cleanUpPtr->minMeV[2] = e_collideMiniMeV;
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
extern void scatterProtonByNeutron(struct electrode * a_electrodePtr, double a_MeV, double a_scale, int a_collideType, double a_Mol, struct atomNodeConst * a_protonPtr)
{
	//The neutrons with high energy scatter same number of protons step-half by step-half.
	double MeV, op;
	double debugMeVMol = 0.0;//DEBUG
	//fprintf(stderr, "DEBUG:%s(%d):{a_MeV:%lg a_scale:%lg a_collideType:%d a_Mol:%lg, a_MeV * a_Mol=%lg\n", __FUNCTION__, __LINE__, a_MeV, a_scale, a_collideType, a_Mol, a_MeV * a_Mol);
	op = 1.0 - a_scale;
	for(MeV = a_MeV * a_scale; MeV >= e_collideMiniMeV; MeV *= op){
		collideBulletToElectrode(a_electrodePtr, a_collideType, MeV, a_Mol, e_collideProtonRateOnElectrode, a_protonPtr, 1);
		debugMeVMol += (MeV * a_Mol);//DEBUG
		//fprintf(stderr, "DEBUG:%s(%d):MeV:%lg * a_Mol:%lg = %lg, debugMeVMol:%lg\n", __FUNCTION__, __LINE__, MeV, a_Mol, MeV * a_Mol, debugMeVMol);
	}
	registOutput(a_electrodePtr, SCAT_SMALL_MASS_DEFECT, MASS_DEFECT_BY_GANMMA, MeV , a_Mol / a_scale);
	debugMeVMol += (MeV * a_Mol / a_scale);//DEBUG
	//fprintf(stderr, "DEBUG:%s(%d):(MeV:%lg * a_Mol:%lg / a_scale:%lg) = %lg, debugMeVMol:%lg\n", __FUNCTION__, __LINE__, MeV, a_Mol, a_scale, MeV * a_Mol / a_scale, debugMeVMol);
	
	//DEBUG, check code
	if(debugMeVMol > 0.0){
		double debugRate = ((a_MeV * a_Mol) - debugMeVMol) / debugMeVMol;
		if(debugRate < 0.0){
			debugRate = - debugRate;
		}
		if(debugRate > 0.00001){
			fprintf(stderr, "ERROR:%s(%d):a_MeV:%lg a_scale:%lg a_collideType:%d a_Mol:%lg debugRate:%lg\n", __FUNCTION__, __LINE__, a_MeV, a_scale, a_collideType, a_Mol, debugRate);
			exit(1);
		}
	}
	//fprintf(stderr, "DEBUG:%s(%d):}a_MeV:%lg a_scale:%lg a_collideType:%d a_Mol:%lg\n", __FUNCTION__, __LINE__, a_MeV, a_scale, a_collideType, a_Mol);
}
extern void scatterIn(struct electrode * a_electrodePtr, int a_massDefectBy)
{
	int doProcess;
	struct cleanUp cUp;
	//struct timeval tv1, tv2;//DEBUG
    //gettimeofday(&tv1, NULL);//DEBUG
	//e_inComptonEffect = 1;//DEBUG
	//fprintf(stderr, "DEBUG:%s(%d):{BEGIN %s usedCnt:%d %s\n", __FUNCTION__, __LINE__, a_electrodePtr->atomHashTable.tableName, a_electrodePtr->atomHashTable.usedCnt, getMassDefectByName(a_massDefectBy));	
	cUp.electrodePtr = a_electrodePtr;
	if(a_massDefectBy == MASS_DEFECT_BY_GANMMA){
		cUp.detectLimitMolForIsotopeSC = a_electrodePtr->detectLimitMolForIsotope * SCATTER_CNT;
	}else{
		cUp.detectLimitMolForIsotopeSC = a_electrodePtr->detectLimitMolForIsotope;	
	}
	for(doProcess = 0; doProcess < 10; ++doProcess){
		unsigned int size, i, j = 0;
		struct objectNodeConst ** oPtr;
		
		oPtr = getFlatTable(&a_electrodePtr->massDefectHashTable[a_massDefectBy], SORT_ASCEND, &size);
		//e_cntScatterPhotonUnnest = 0;//DEBUG
		//e_cntScatterPhoton = 0;//DEBUG
		if(a_massDefectBy != MASS_DEFECT_BY_GANMMA){
			//fprintf(stderr, "DEBUG:%s(%d):doProcess:%d getFlatTable size:%d\n", __FUNCTION__, __LINE__, doProcess, size);
		}
		//fprintf(stderr, "DEBUG:%s:calcPastTime-1:%lg[sec]\n", __FUNCTION__, calcPastTime(&tv2, &tv1));
		if(oPtr){
			struct massDefect * valuePtr;
			//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
			for(i = 0; i < size; ++i){
				valuePtr = (struct massDefect *)oPtr[i]->valuePtr;
				//if(a_massDefectBy != MASS_DEFECT_BY_GANMMA){//DEBUG
					//fprintf(stderr, "DEBUG:%s(%d):%s i:%d/size:%d \n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy), i, size);
				//}
				if(valuePtr->mol >= cUp.detectLimitMolForIsotopeSC){
					int nest = 0;
					double saveMol = valuePtr->mol;
					//if(a_massDefectBy != MASS_DEFECT_BY_GANMMA){//DEBUG
						//fprintf(stderr, "DEBUG:%s(%d):valuePtr->mol:%lg >= cUp.detectLimitMolForIsotopeSC:%lg\n", __FUNCTION__, __LINE__, valuePtr->mol, cUp.detectLimitMolForIsotopeSC);
					//}
					if(a_massDefectBy == MASS_DEFECT_BY_PROTON){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							struct atomNodeConst * protonPtr;
							protonPtr = findProtonInElectrode(a_electrodePtr);
							valuePtr->mol = 0.0;
							if(e_useProtonScattering && protonPtr){
								//fprintf(stderr, "DEBUG:%s(%d):e_useProtonScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useProtonScattering, i, valuePtr->MeV, saveMol);
								collideBulletToElectrode(a_electrodePtr, COLLIDE_PROTON_S, valuePtr->MeV, saveMol, e_collideProtonRateOnElectrode, protonPtr, 1);
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):e_useProtonScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useProtonScattering, i, valuePtr->MeV, saveMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}else if(a_massDefectBy == MASS_DEFECT_BY_DEUTERIUM){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							struct atomNodeConst * deuteriumPtr;
							deuteriumPtr = findDeuteriumInElectrode(a_electrodePtr);
							valuePtr->mol = 0.0;
							if(e_useProtonScattering && deuteriumPtr){
								//fprintf(stderr, "DEBUG:%s(%d):e_useProtonScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useProtonScattering, i, valuePtr->MeV, saveMol);
								collideBulletToElectrode(a_electrodePtr, COLLIDE_DEUTERIUM_S, valuePtr->MeV, saveMol, e_collideProtonRateOnElectrode, deuteriumPtr, 1);
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):e_useProtonScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useProtonScattering, i, valuePtr->MeV, saveMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}else if(a_massDefectBy == MASS_DEFECT_BY_TRITIUM){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							struct atomNodeConst * tritiumPtr;
							tritiumPtr = findTritiumInElectrode(a_electrodePtr);
							valuePtr->mol = 0.0;
							if(e_useProtonScattering && tritiumPtr){
								//fprintf(stderr, "DEBUG:%s(%d):e_useProtonScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useProtonScattering, i, valuePtr->MeV, saveMol);
								collideBulletToElectrode(a_electrodePtr, COLLIDE_TRITIUM_S, valuePtr->MeV, saveMol, e_collideProtonRateOnElectrode, tritiumPtr, 1);
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):e_useProtonScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useProtonScattering, i, valuePtr->MeV, saveMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}else if(a_massDefectBy == MASS_DEFECT_BY_NEUTRON){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							valuePtr->mol = 0.0;
							if(e_useNeutonScattering){
								struct atomNodeConst * protonPtr;
								struct atomNodeConst * deuteriumPtr;
								struct atomNodeConst * tritiumPtr;
								double protonMol, deuteriumMol, tritiumMol, totalMol, actMol, remMol;
								//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
								protonPtr = findProtonInElectrode(a_electrodePtr);
								deuteriumPtr = findDeuteriumInElectrode(a_electrodePtr);
								tritiumPtr = findTritiumInElectrode(a_electrodePtr);
								protonMol = (protonPtr) ? getMol(protonPtr->vPtr) : 0.0;
								deuteriumMol = (deuteriumPtr) ? getMol(deuteriumPtr->vPtr) : 0.0;
								tritiumMol = (tritiumPtr) ? getMol(tritiumPtr->vPtr) : 0.0;
								totalMol = protonMol + deuteriumMol + tritiumMol;
								if(totalMol > saveMol){
									//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
									actMol = saveMol;
									remMol = 0.0;
								}else{
									//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
									actMol = totalMol;
									remMol = saveMol - totalMol;
								}
								//fprintf(stderr, "DEBUG:%s(%d):e_useNeutonScattering:%d i:%d MeV:%lg saveMol:%lg actMol:%lg remMol:%lg\n", __FUNCTION__, __LINE__, e_useNeutonScattering, i, valuePtr->MeV, saveMol, actMol, remMol);
								
								if(actMol > 0.0){
									//The neutrons with high energy scatter protons, deutriums and tritiums.
									//We ignore that the neutrons with high energy scatter other atoms.
									if(protonPtr){
										//fprintf(stderr, "DEBUG:%s(%d):scatterProtonByNeutron(actMol * protonMol:%lg / totalMol:%lg = %lg)\n", __FUNCTION__, __LINE__, protonMol, totalMol, actMol * protonMol / totalMol);
										scatterProtonByNeutron(a_electrodePtr, valuePtr->MeV, 0.5, COLLIDE_PROTON_S, actMol * protonMol / totalMol, protonPtr);
									}
									if(deuteriumPtr){
										//fprintf(stderr, "DEBUG:%s(%d):scatterProtonByNeutron(actMol * deuteriumMol:%lg / totalMol:%lg = %lg)\n", __FUNCTION__, __LINE__, deuteriumMol, totalMol, actMol * deuteriumMol / totalMol);
										scatterProtonByNeutron(a_electrodePtr, valuePtr->MeV, 1.0 / 3.0, COLLIDE_DEUTERIUM_S, actMol * deuteriumMol / totalMol, deuteriumPtr);
									}
									if(tritiumPtr){
										//fprintf(stderr, "DEBUG:%s(%d):scatterProtonByNeutron(actMol * tritiumMol:%lg / totalMol:%lg = %lg)\n", __FUNCTION__, __LINE__, tritiumMol, totalMol, actMol * tritiumMol / totalMol);
										scatterProtonByNeutron(a_electrodePtr, valuePtr->MeV, 0.25, COLLIDE_TRITIUM_S, actMol * tritiumMol / totalMol, tritiumPtr);
									}
								}
								if(remMol > 0.0){
									//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
									registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, remMol);
								}
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}else if(a_massDefectBy == MASS_DEFECT_BY_ALPHA){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							struct atomNodeConst * heliumPtr;
							heliumPtr = findHeliumInElectrode(a_electrodePtr);
							valuePtr->mol = 0.0;
							if(e_useAlphaScattering && heliumPtr){
								//fprintf(stderr, "DEBUG:%s(%d):e_useAlphaScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useAlphaScattering, i, valuePtr->MeV, saveMol);
								collideBulletToElectrode(a_electrodePtr, COLLIDE_ALPHA_S, valuePtr->MeV, saveMol, e_collideProtonRateOnElectrode, heliumPtr, 1);
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):e_useAlphaScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useAlphaScattering, i, valuePtr->MeV, saveMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}else if(a_massDefectBy == MASS_DEFECT_BY_BETA){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							struct atomNodeConst * electronPtr;
							electronPtr = findElectronInElectrode(a_electrodePtr);
							valuePtr->mol = 0.0;
							if(e_useElectronScattering && electronPtr){
								//fprintf(stderr, "DEBUG:%s(%d):e_useElectronScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useElectronScattering, i, valuePtr->MeV, saveMol);
								collideBulletToElectrode(a_electrodePtr, COLLIDE_ELECTRON_S, valuePtr->MeV, saveMol, e_collideProtonRateOnElectrode, electronPtr, 1);								
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):e_useElectronScattering:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useElectronScattering, i, valuePtr->MeV, saveMol);
								registOutput(a_electrodePtr, SCAT_BIG_MASS_DEFECT_NOW, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}else if(a_massDefectBy == MASS_DEFECT_BY_GANMMA){
						//fprintf(stderr, "DEBUG:%s(%d):%s\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy));
						if(valuePtr->MeV >= e_collideMiniMeV){
							//fprintf(stderr, "DEBUG:%s(%d):\n", __FUNCTION__, __LINE__);
							valuePtr->mol = 0.0;
							if(e_useComptonEffect){
								//fprintf(stderr, "DEBUG:%s(%d):e_useComptonEffect:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useComptonEffect, i, valuePtr->MeV, saveMol);
								//fprintf(stderr, "DEBUG:%s:loop:%lg[sec] i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, calcPastTime(&tv2, &tv1), i, valuePtr->MeV, saveMol);
								//fprintf(stderr, "DEBUG:%s:loop:i:%d tableName:%s usedCnt:%d\n", __FUNCTION__, i, a_electrodePtr->atomHashTable.tableName, a_electrodePtr->atomHashTable.usedCnt);
								//fprintf(stderr, "DEBUG:%s:e_cntScatterPhotonUnnest:%d e_cntScatterPhoton:%d\n", __FUNCTION__, e_cntScatterPhotonUnnest, e_cntScatterPhoton);
								//fprintf(stderr, "DEBUG:%s:massDefectHashTable[a_massDefectBy].usedCnt:%u\n", __FUNCTION__, a_electrodePtr->massDefectHashTable[a_massDefectBy].usedCnt);
								//e_cntScatterPhotonUnnest_s = 0;//DEBUG
								//e_cntScatterPhoton_s = 0;//DEBUG
								scatterPhoton(a_electrodePtr, valuePtr->MeV, saveMol, nest);
								//fprintf(stderr, "DEBUG:%s:e_cntScatterPhotonUnnes_s:%d e_cntScatterPhoton_s:%d\n", __FUNCTION__, e_cntScatterPhotonUnnest_s, e_cntScatterPhoton_s);
								//fprintf(stderr, "DEBUG:%s:massDefectHashTable[a_massDefectBy].usedCnt:%u\n", __FUNCTION__, a_electrodePtr->massDefectHashTable[a_massDefectBy].usedCnt);
								/*
								[!!CAUTION!!!] We can't use "scatterPhoton" as the argument of "iterateInHashTable", 
								because "scatterPhoton" will uses "insertObjectInHashTable" or "iterateInHashTable" through nesting calls.
								*/
							}else{
								//fprintf(stderr, "DEBUG:%s(%d):e_useComptonEffect:%d i:%d MeV:%lg saveMol:%lg\n", __FUNCTION__, __LINE__, e_useComptonEffect, i, valuePtr->MeV, saveMol);
								registOutput(a_electrodePtr, SCAT_DEGRADE_IN_COMPTON_B, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
							}
							++j;
						}else{
							fprintf(stderr, "FATAL:%s(%d):%d %s valuePtr->MeV:%lg < e_collideMiniMeV:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, a_electrodePtr->massDefectHashTable[a_massDefectBy].tableName, valuePtr->MeV, e_collideMiniMeV);
							exit(1);
						}
					}
				}else{
					double saveMol = valuePtr->mol;
					//fprintf(stderr, "DEBUG:%s(%d):%s valuePtr->mol:%lg < cUp.detectLimitMolForIsotopeSC:%lg\n", __FUNCTION__, __LINE__, getMassDefectByName(a_massDefectBy), valuePtr->mol, cUp.detectLimitMolForIsotopeSC);
					valuePtr->mol = 0.0;
					++j;
					registOutput(a_electrodePtr, SCAT_DEGRADE_IN_COMPTON_B, MASS_DEFECT_BY_GANMMA, valuePtr->MeV, saveMol);
				}
			}
			free(oPtr);
			//fprintf(stderr, "DEBUG:%s:calcPastTime-2:%lg[sec] scatterPhoton * size:%d\n", __FUNCTION__, calcPastTime(&tv2, &tv1), size);
		}
		//fprintf(stderr, "DEBUG:%s:calcPastTime-3:%lg[sec] scatterPhoton * size:%d\n", __FUNCTION__, calcPastTime(&tv2, &tv1), size);
		//fprintf(stderr, "DEBUG:%s:e_cntScatterPhotonUnnest:%d e_cntScatterPhoton:%d\n", __FUNCTION__, e_cntScatterPhotonUnnest, e_cntScatterPhoton);
		//fprintf(stderr, "DEBUG:%s:doProcess:%d size:%u j:%u massDefectHashTable[a_massDefectBy].usedCnt:%u\n", __FUNCTION__, doProcess, size, j, a_electrodePtr->massDefectHashTable[a_massDefectBy].usedCnt);
		//fprintf(stderr, "DEBUG:%s(%d):j:%d\n", __FUNCTION__, __LINE__, j);
		if(j > 0){
			//It's better to clean up "massDefectHashTable[a_massDefectBy]" for both small usage of memory and faster processing.
			cUp.cntOfBigMolNode = 0;
			cUp.maxChangeMeVOfBigMol = 0.0;
			cUp.minChangeMeVOfBigMol = 0.0;
			cUp.cntOfSmallMolNode = 0;
			cUp.totalMeVMol = 0.0;
			cUp.degradeMeVMol = 0.0;
			iterateInHashTable(&a_electrodePtr->massDefectHashTable[a_massDefectBy], &cUp, iterateClenupMassDefect);
			if(initMaxMinOfClenup(&cUp)){
				//I recommend the degrading process of SCAT_BIG_MASS_DEFECT_NOW with using the aproximatic calculation for faster calculation.
				//fprintf(stderr, "DEBUG:%s:doProcess:%d cUp.cntOfBigMolNode:%d cUp.maxChangeMeVOfBigMol:%lg cUp.minChangeMeVOfBigMol:%lg\n", __FUNCTION__, doProcess, cUp.cntOfBigMolNode, cUp.maxChangeMeVOfBigMol, cUp.minChangeMeVOfBigMol);
				//fprintf(stderr, "DEBUG:%s:doProcess:%d cUp.cntOfSmallMolNode:%d\n", __FUNCTION__, doProcess, cUp.cntOfSmallMolNode);
				if(cUp.cntOfBigMolNode > DEGRADE_MAG * DEGRADE_DIV || cUp.cntOfSmallMolNode > DEGRADE_MAG * DEGRADE_DIV * 8){
					iterateInHashTable(&a_electrodePtr->massDefectHashTable[a_massDefectBy], &cUp, iterateDegreadeMassDefect);
					if(abs((int)(10000 * (cUp.degradeMeVMol - cUp.totalMeVMol) / cUp.totalMeVMol)) > 1){
						fprintf(stderr, "ERROR:%s(%d):a_massDefectBy:%d cUp.degradeMeVMol:%lg != cUp.totalMeVMol:%lg\n", __FUNCTION__, __LINE__, a_massDefectBy, cUp.degradeMeVMol, cUp.totalMeVMol);
						exit(1);
					}
					moveFromClenup(&cUp);
				}
			}
		}
		//fprintf(stderr, "DEBUG:%s:massDefectHashTable[a_massDefectBy].usedCnt:%u j:%d\n", __FUNCTION__, a_electrodePtr->massDefectHashTable[a_massDefectBy].usedCnt, j);

		if(j == 0){
			break;
		}
	}
	//fprintf(stderr, "DEBUG:%s:tableName:%s usedCnt:%d\n", __FUNCTION__, a_electrodePtr->atomHashTable.tableName, a_electrodePtr->atomHashTable.usedCnt);
	//fprintf(stderr, "DEBUG:%s(%d):}END calcPastTime-E:%lg[sec]\n", __FUNCTION__, __LINE__, calcPastTime(&tv2, &tv1));	
	//fprintf(stderr, "DEBUG:%s(%d):}END\n", __FUNCTION__, __LINE__);	
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
	//struct atomNodeConst * positiveElectronPtr;
	//struct atomNodeConst * negativeElectronPtr;
	//positiveElectronPtr = findElectronInElectrode(&e_positiveElectrode);	
	//negativeElectronPtr = findElectronInElectrode(&e_negativeElectrode);	
	
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
		registOutput(&e_negativeElectrode, SCAT_NOT_COLLIDE, MASS_DEFECT_BY_GANMMA, e_stepByStepLostingEnergyMeV, a_pastSecond * e_genedHydrogenInSpaceMol * 0.5);
		registOutput(&e_positiveElectrode, SCAT_NOT_COLLIDE, MASS_DEFECT_BY_GANMMA, e_stepByStepLostingEnergyMeV, a_pastSecond * e_genedHydrogenInSpaceMol * 0.5);
	}

	collideBulletToElectrode(&e_positiveElectrode, COLLIDE_ELECTRON, e_appliedVoltageMeV, 
		a_pastSecond * e_arrivedElectronMol, e_collideElectronRateOnElectrode, NULL, 1);

	collideBulletToElectrode(&e_negativeElectrode, COLLIDE_PROTON, e_appliedVoltageMeV, 
		a_pastSecond * e_arrivedProtonMol * protonRate, e_collideProtonRateOnElectrode, NULL, 1);
	collideBulletToElectrode(&e_negativeElectrode, COLLIDE_DEUTERIUM, e_appliedVoltageMeV,
		a_pastSecond * e_arrivedProtonMol * deuteriumRate, e_collideProtonRateOnElectrode, NULL, 1);
	collideBulletToElectrode(&e_negativeElectrode, COLLIDE_TRITIUM, e_appliedVoltageMeV,
		a_pastSecond * e_arrivedProtonMol * tritiumRate, e_collideProtonRateOnElectrode, NULL, 1);

	//fprintf(stderr, "DEBUG:%s:}END\n\n", __FUNCTION__);
}


//---------------------------------------------------------------------

#define OLD_VERSION_OF_serializeRunningCondition 1
#define VERSION_OF_serializeRunningCondition 2
#define SERIALIZE_RUNNING_CONDITION(RW) \
extern int RW ## RunningCondition(FILE * a_fp) \
{\
	SERIALIZE_VERSION_CHECK2(RW, VERSION_OF_serializeRunningCondition, OLD_VERSION_OF_serializeRunningCondition, a_fp) \
	SERIALIZE_VALUE(RW, e_rc.startTime, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.termTime, a_fp)\
	SERIALIZE_VALUE(RW, e_rc.intervalTime, a_fp)\
	if((void *)RW == (void *)fread){\
		if(version == OLD_VERSION_OF_serializeRunningCondition){\
			e_rc.numberOfBreakingIntervals = 0;\
		}else{\
			SERIALIZE_VALUE(RW, e_rc.numberOfBreakingIntervals, a_fp)\
		}\
	}else{\
		SERIALIZE_VALUE(RW, e_rc.numberOfBreakingIntervals, a_fp)\
	}\
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
	a_rc->numberOfBreakingIntervals = 0;
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
		fprintf(a_fp, "%sR%d numberOfBreakingIntervals %d\n", a_title, i++, e_rc.numberOfBreakingIntervals);
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
	struct sumMeVMol massDefectNowN[COUNT_OF_MASS_DEFECT_HASH_TABLE];
	struct sumMeVMol massDefectNowP[COUNT_OF_MASS_DEFECT_HASH_TABLE];
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
	/* sumOfMassUMol_NE = */ printElectrode(a_fp, &e_negativeElectrode, massDefectNowN, &sumOfMassUMolIni_NE, &sumOfMassUMolAdd_NE, &sumOfMassUMolSub_NE, &heatCapacity_NE);
	fprintf(a_fp, "%s -----\n\n", timeMess);
	/*sumOfMassUMol_PE = */ printElectrode(a_fp, &e_positiveElectrode, massDefectNowP, &sumOfMassUMolIni_PE, &sumOfMassUMolAdd_PE, &sumOfMassUMolSub_PE, &heatCapacity_PE);
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
		struct sumMeVMol massDefectNow[COUNT_OF_MASS_DEFECT_HASH_TABLE];
		struct sumMeVMol massDefectAll;
		struct sumMeVMol byNeutrinoAll;
		struct sumMeVMol scattered[SIZE_OF_SCATTERD];
		int i;
		for(i = 0; i < COUNT_OF_MASS_DEFECT_HASH_TABLE; ++i){
			mergeSumMeVMol(&massDefectNow[i], &massDefectNowN[i], &massDefectNowP[i]);
		}
		mergeSumMeVMol(&massDefectAll, &e_negativeElectrode.massDefectAll, &e_positiveElectrode.massDefectAll);
		mergeSumMeVMol(&byNeutrinoAll, &e_negativeElectrode.byNeutrinoAll, &e_positiveElectrode.byNeutrinoAll);
		mergeScattered(scattered, e_negativeElectrode.scattered, e_positiveElectrode.scattered);
		outputHeat = printMassDefectHeat(a_fp, "TOTAL ENERGY IN THE SYSTEM", massDefectNow, &massDefectAll, &byNeutrinoAll, scattered);
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

	if(a_BeginOrEnd == END_PRINT || e_rc.logNecleusReactions){
		fprintf(a_fp, "%s -----\n\n", timeMess);
		//fprintf(stderr, "DEBUG:%s:e_nuclearReaction.usedCnt=%d\n", __FUNCTION__, e_nuclearReaction.usedCnt);
		printAllNuclearReaction(a_fp, timeMess);
		printUnregistedAtomNumbers(a_fp);
		printUnregistedIsotopeNumbers(a_fp);
	}
	fprintf(a_fp, "%s >>>>>\n\n", timeMess);
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
		}else if(memcmp(argv[i], "-b=", 3) == 0){
			//fprintf(stderr, "DEBUG:%s:numberOfBreakingIntervals -b\n", __FUNCTION__);
			rc.numberOfBreakingIntervals = getNumber(argv[i]);
			fprintf(stderr, "DEBUG:%s:numberOfBreakingIntervals:%d\n", __FUNCTION__, rc.numberOfBreakingIntervals);
			if(rc.numberOfBreakingIntervals < 0){
				fprintf(stderr, "FATAL ERROR:Wrong numberOfBreakingIntervals:%d\n", rc.numberOfBreakingIntervals);
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
		}else if(strcmp(argv[i], "-EC") == 0){
			printAllElectronCapture(stdout);
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
-b=numberOfBreakingIntervals :\n\
   This is the number of breaking intervals with stopping the electric current.\n\
   The 'numberOfBreakingIntervals' is zero (0) or a positive integer number.\n\
   The zero (0) means that the electric current is not stopped.\n\
   The default of numberOfBreakingIntervals is 0.\n\
   If electricPulseLength is positive, \n\
   it means that the the electric current is stopped \n\
   while the following number of intervals.\n\
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
		double e_collideMidiMeV = %lg;\n\
		double e_collideMiniMeV = %lg;\n\
		int e_usebulletCrossSection = %d;\n\
		int e_useProtonScattering = %d;\n\
		int e_useNeutonScattering = %d;\n\
		int e_useAlphaScattering = %d;\n\
		int e_useElectronScattering = %d;\n\
		int e_useComptonEffect = %d;\n\
		double e_rateForAlphaParticle = %lg;\n\
		double e_rateForProtonAtBetaPlus = %lg;\n\
		double e_rateFor2ProtonAtBetaPlus = %lg;\n\
		double e_rateFor3ProtonAtBetaPlus = %lg;\n\
		double e_rateForAlphaParticleAtBetaPlus = %lg;\n\
		double e_rateForProtonAtEC = %lg;\n\
		double e_rateFor2ProtonAtEC = %lg;\n\
		double e_rateFor3ProtonAtEC = %lg;\n\
		double e_rateForAlphaParticleAtEC = %lg;\n\
		double e_rateForNeytonAtBetaMinus = %lg;\n\
		double e_rateFor2NeytonAtBetaMinus = %lg;\n\
		double e_rateFor3NeytonAtBetaMinus = %lg;\n\
		double e_rateFor4NeytonAtBetaMinus = %lg;\n\
		double e_rateForAlphaParticleAtBetaMinus = %lg;\n\
-H- : It specifies to reduce the amount of hydrogen and helium in the positive electrode\n\
    into one-1,000,000th before running the simulation.\n\
-EC : Print the property of Electron Capture in decay modes. Not simulate.\n\
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
			e_collideMidiMeV,
			e_collideMiniMeV,
			e_usebulletCrossSection,
			e_useProtonScattering,
			e_useNeutonScattering,
			e_useAlphaScattering,
			e_useElectronScattering,
			e_useComptonEffect,
			e_rateForAlphaParticle,
			e_rateForProtonAtBetaPlus,
			e_rateFor2ProtonAtBetaPlus,
			e_rateFor3ProtonAtBetaPlus,
			e_rateForAlphaParticleAtBetaPlus,
			e_rateForProtonAtEC,
			e_rateFor2ProtonAtEC,
			e_rateFor3ProtonAtEC,
			e_rateForAlphaParticleAtEC,
			e_rateForNeytonAtBetaMinus,
			e_rateFor2NeytonAtBetaMinus,
			e_rateFor3NeytonAtBetaMinus,
			e_rateFor4NeytonAtBetaMinus,
			e_rateForAlphaParticleAtBetaMinus);
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
	//exit(0);//debug
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
		if((e_rc.numberOfBreakingIntervals == 0)
		|| (loop % (e_rc.numberOfBreakingIntervals + 1) == 0)){
			pulseCurrent(e_rc.intervalTime);
			e_rc.inputTotalEnergy += (e_emittedPulseEnergyMeV * e_rc.intervalTime);
			//printSumup(e_logFp, "after pulseCurrent", loop, tail, MID_PRINT);
		}

		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_NEUTRON);//It must be here before absorbeOrDecayNeutronInElectrode()
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_NEUTRON);//It must be here before absorbeOrDecayNeutronInElectrode()

		absorbeOrDecayNeutronInElectrode(&e_negativeElectrode, e_rc.intervalTime);
		absorbeOrDecayNeutronInElectrode(&e_positiveElectrode, e_rc.intervalTime);
		//printSumup(e_logFp, "after absorbeOrDecayNeutronInElectrode", loop, tail, MID_PRINT);
		
		decayNeuclayWithoutNeutron(&e_negativeElectrode, e_rc.intervalTime);
		decayNeuclayWithoutNeutron(&e_positiveElectrode, e_rc.intervalTime);
		//printSumup(e_logFp, "after decayNeuclayWithoutNeutron", loop, tail, MID_PRINT);
		
		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_PROTON);
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_PROTON);
		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_DEUTERIUM);
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_DEUTERIUM);
		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_TRITIUM);
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_TRITIUM);
		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_ALPHA);
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_ALPHA);
		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_BETA);
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_BETA);
		scatterIn(&e_negativeElectrode, MASS_DEFECT_BY_GANMMA);
		scatterIn(&e_positiveElectrode, MASS_DEFECT_BY_GANMMA);
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