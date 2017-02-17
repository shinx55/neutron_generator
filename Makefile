all:simulate analyzelog

simulate: simulate.c hashtable.c
	gcc -Wall -Wformat -g -O0 -o simulate simulate.c hashtable.c timeformat.c

#	gcc -g -O0 -v -da -Q -o simulate simulate.c hashtable.c

analyzelog: analyzelog.c hashtable.c
	gcc -Wall -Wformat -g -O0 -o analyzelog analyzelog.c hashtable.c timeformat.c

clean:
	rm simulate.exe simulate.exe.stackdump
	rm simulate.c.1* simulate.c.2* hashtable.c.1* hashtable.c.2*
	rm analyzelog.exe analyzelog.exe.stackdump

rmlog:
	rm -f data.dat s00*.log m*.log h*.log

RUN_SIMULATE = ../../neutron_generator/simulate
RUN_ANALYZELOG = ../../neutron_generator/analyzelog

help:
	./simulate -h

ElectronCaptureList:
	mkdir -p ElectronCaptureList
	cd ElectronCaptureList ; \
	rm -f data.dat s00*.log m*.log h*.log ; \
	$(RUN_SIMULATE) -EC >ElectronCaptureList.txt

testMAX:
	mkdir -p testMAX
	cd testMAX ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.05, D=0.05, Ni=0.8, Li=0.2\";" > testMAX.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.05, D=0.05, Ni=0.8, Li=0.2\";" >> testMAX.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testMAX.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testMAX.txt ; \
	echo "double e_appliedVoltageScale = 10.0;" >>testMAX.txt ; \
	$(RUN_SIMULATE) -P=testMAX.txt ; \
	$(RUN_SIMULATE) -t=15m -P=testMAX.txt

# testNi is :
# (1) run 'sumulate' while 1 minute + many hour with default contitions.
testNiLessProton:
	mkdir -p testNiLessProton
	cd testNiLessProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testNiLessProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testNiLessProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testNiLessProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testNiLessProton.txt ; \
	$(RUN_SIMULATE) -P=testNiLessProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiLessProton.txt

testAlLessProton:
	mkdir -p testAlLessProton
	cd testAlLessProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlLessProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" >> testAlLessProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testAlLessProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testAlLessProton.txt ; \
	$(RUN_SIMULATE) -P=testAlLessProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlLessProton.txt

testNiPoorProton:
	mkdir -p testNiPoorProton
	cd testNiPoorProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" > testNiPoorProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testNiPoorProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testNiPoorProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testNiPoorProton.txt ; \
	$(RUN_SIMULATE) -P=testNiPoorProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiPoorProton.txt

testAlPoorProton:
	mkdir -p testAlPoorProton
	cd testAlPoorProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" > testAlPoorProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" >> testAlPoorProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testAlPoorProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testAlPoorProton.txt ; \
	$(RUN_SIMULATE) -P=testAlPoorProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlPoorProton.txt

testFePoorProton:
	mkdir -p testFePoorProton
	cd testFePoorProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" > testFePoorProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" >> testFePoorProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testFePoorProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testFePoorProton.txt ; \
	$(RUN_SIMULATE) -P=testFePoorProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFePoorProton.txt

testMnPoorProton:
	mkdir -p testMnPoorProton
	cd testMnPoorProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" > testMnPoorProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" >> testMnPoorProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testMnPoorProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testMnPoorProton.txt ; \
	$(RUN_SIMULATE) -P=testMnPoorProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnPoorProton.txt

testLiPoorProton:
	mkdir -p testLiPoorProton
	cd testLiPoorProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" > testLiPoorProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiPoorProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiPoorProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiPoorProton.txt ; \
	$(RUN_SIMULATE) -P=testLiPoorProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiPoorProton.txt

testNiMuchProton:
	mkdir -p testNiMuchProton
	cd testNiMuchProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" > testNiMuchProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testNiMuchProton.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testNiMuchProton.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testNiMuchProton.txt ; \
	$(RUN_SIMULATE) -P=testNiMuchProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiMuchProton.txt

testAlMuchProton:
	mkdir -p testAlMuchProton
	cd testAlMuchProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" > testAlMuchProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" >> testAlMuchProton.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testAlMuchProton.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testAlMuchProton.txt ; \
	$(RUN_SIMULATE) -P=testAlMuchProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlMuchProton.txt

testFeMuchProton:
	mkdir -p testFeMuchProton
	cd testFeMuchProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" > testFeMuchProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" >> testFeMuchProton.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testFeMuchProton.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testFeMuchProton.txt ; \
	$(RUN_SIMULATE) -P=testFeMuchProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFeMuchProton.txt

testMnMuchProton:
	mkdir -p testMnMuchProton
	cd testMnMuchProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" > testMnMuchProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" >> testMnMuchProton.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testMnMuchProton.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testMnMuchProton.txt ; \
	$(RUN_SIMULATE) -P=testMnMuchProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnMuchProton.txt

testLiMuchProton:
	mkdir -p testLiMuchProton
	cd testLiMuchProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" > testLiMuchProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiMuchProton.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testLiMuchProton.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testLiMuchProton.txt ; \
	$(RUN_SIMULATE) -P=testLiMuchProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiMuchProton.txt

testNiMuchD:
	mkdir -p testNiMuchD
	cd testNiMuchD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Ni=1.0\";" > testNiMuchD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Ni=1.0\";" >> testNiMuchD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testNiMuchD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testNiMuchD.txt ; \
	$(RUN_SIMULATE) -P=testNiMuchD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiMuchD.txt

testAlMuchD:
	mkdir -p testAlMuchD
	cd testAlMuchD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" > testAlMuchD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" >> testAlMuchD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testAlMuchD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testAlMuchD.txt ; \
	$(RUN_SIMULATE) -P=testAlMuchD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlMuchD.txt

testFeMuchD:
	mkdir -p testFeMuchD
	cd testFeMuchD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Fe=1.0\";" > testFeMuchD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Fe=1.0\";" >> testFeMuchD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testFeMuchD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testFeMuchD.txt ; \
	$(RUN_SIMULATE) -P=testFeMuchD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFeMuchD.txt

testMnMuchD:
	mkdir -p testMnMuchD
	cd testMnMuchD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Mn=1.0\";" > testMnMuchD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Mn=1.0\";" >> testMnMuchD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testMnMuchD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testMnMuchD.txt ; \
	$(RUN_SIMULATE) -P=testMnMuchD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnMuchD.txt

testLiMuchD:
	mkdir -p testLiMuchD
	cd testLiMuchD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Li=1.0\";" > testLiMuchD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Li=1.0\";" >> testLiMuchD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testLiMuchD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testLiMuchD.txt ; \
	$(RUN_SIMULATE) -P=testLiMuchD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiMuchD.txt

testNiMaxD:
	mkdir -p testNiMaxD
	cd testNiMaxD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Ni=1.0\";" > testNiMaxD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Ni=1.0\";" >> testNiMaxD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testNiMaxD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testNiMaxD.txt ; \
	echo "double e_appliedVoltageScale = 10.0;" >>testNiMaxD.txt ; \
	$(RUN_SIMULATE) -P=testNiMaxD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiMaxD.txt

testAlMaxD:
	mkdir -p testAlMaxD
	cd testAlMaxD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" > testAlMaxD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" >> testAlMaxD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testAlMaxD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testAlMaxD.txt ; \
	echo "double e_appliedVoltageScale = 4.5;" >>testAlMaxD.txt ; \
	$(RUN_SIMULATE) -P=testAlMaxD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlMaxD.txt

testFeMaxD:
	mkdir -p testFeMaxD
	cd testFeMaxD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Fe=1.0\";" > testFeMaxD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Fe=1.0\";" >> testFeMaxD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testFeMaxD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testFeMaxD.txt ; \
	echo "double e_appliedVoltageScale = 10.0;" >>testFeMaxD.txt ; \
	$(RUN_SIMULATE) -P=testFeMaxD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFeMaxD.txt

testMnMaxD:
	mkdir -p testMnMaxD
	cd testMnMaxD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Mn=1.0\";" > testMnMaxD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Mn=1.0\";" >> testMnMaxD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testMnMaxD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testMnMaxD.txt ; \
	echo "double e_appliedVoltageScale = 10.0;" >>testMnMaxD.txt ; \
	$(RUN_SIMULATE) -P=testMnMaxD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnMaxD.txt

testLiMaxD:
	mkdir -p testLiMaxD
	cd testLiMaxD ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Li=1.0\";" > testLiMaxD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Li=1.0\";" >> testLiMaxD.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testLiMaxD.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testLiMaxD.txt ; \
	echo "double e_appliedVoltageScale = 4.5;" >>testLiMaxD.txt ; \
	$(RUN_SIMULATE) -P=testLiMaxD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiMaxD.txt

testMMM: testNiLessProton testAlLessProton testNiPoorProton testAlPoorProton testFePoorProton testMnPoorProton testLiPoorProton testNiMuchProton testAlMuchProton testFeMuchProton testMnMuchProton testLiMuchProton testNiMuchD testAlMuchD testFeMuchD testMnMuchD testLiMuchD testNiMaxD testAlMaxD testFeMaxD testMnMaxD testLiMaxD

# ---------------------------------------------
testNiLessProtonHV:
	mkdir -p testNiLessProtonHV
	cd testNiLessProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testNiLessProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testNiLessProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testNiLessProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testNiLessProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testNiLessProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testNiLessProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiLessProtonHV.txt

testAlLessProtonHV:
	mkdir -p testAlLessProtonHV
	cd testAlLessProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlLessProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" >> testAlLessProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testAlLessProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testAlLessProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testAlLessProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testAlLessProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlLessProtonHV.txt

testNiPoorProtonHV:
	mkdir -p testNiPoorProtonHV
	cd testNiPoorProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" > testNiPoorProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testNiPoorProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testNiPoorProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testNiPoorProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testNiPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testNiPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiPoorProtonHV.txt

testAlPoorProtonHV:
	mkdir -p testAlPoorProtonHV
	cd testAlPoorProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" > testAlPoorProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" >> testAlPoorProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testAlPoorProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testAlPoorProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testAlPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testAlPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlPoorProtonHV.txt

testFePoorProtonHV:
	mkdir -p testFePoorProtonHV
	cd testFePoorProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" > testFePoorProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" >> testFePoorProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testFePoorProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testFePoorProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testFePoorProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testFePoorProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFePoorProtonHV.txt

testMnPoorProtonHV:
	mkdir -p testMnPoorProtonHV
	cd testMnPoorProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" > testMnPoorProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" >> testMnPoorProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testMnPoorProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testMnPoorProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testMnPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testMnPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnPoorProtonHV.txt

testLiPoorProtonHV:
	mkdir -p testLiPoorProtonHV
	cd testLiPoorProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" > testLiPoorProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiPoorProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiPoorProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiPoorProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testLiPoorProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiPoorProtonHV.txt

testNiMuchProtonHV:
	mkdir -p testNiMuchProtonHV
	cd testNiMuchProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" > testNiMuchProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testNiMuchProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testNiMuchProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testNiMuchProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testNiMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testNiMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiMuchProtonHV.txt

testAlMuchProtonHV:
	mkdir -p testAlMuchProtonHV
	cd testAlMuchProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" > testAlMuchProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Al=1.0\";" >> testAlMuchProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testAlMuchProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testAlMuchProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testAlMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testAlMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlMuchProtonHV.txt

testFeMuchProtonHV:
	mkdir -p testFeMuchProtonHV
	cd testFeMuchProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" > testFeMuchProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Fe=1.0\";" >> testFeMuchProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testFeMuchProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testFeMuchProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testFeMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testFeMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFeMuchProtonHV.txt

testMnMuchProtonHV:
	mkdir -p testMnMuchProtonHV
	cd testMnMuchProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" > testMnMuchProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Mn=1.0\";" >> testMnMuchProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testMnMuchProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testMnMuchProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testMnMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testMnMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnMuchProtonHV.txt

testLiMuchProtonHV:
	mkdir -p testLiMuchProtonHV
	cd testLiMuchProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" > testLiMuchProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiMuchProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testLiMuchProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testLiMuchProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testLiMuchProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiMuchProtonHV.txt

testNiMuchDHV:
	mkdir -p testNiMuchDHV
	cd testNiMuchDHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Ni=1.0\";" > testNiMuchDHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Ni=1.0\";" >> testNiMuchDHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testNiMuchDHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testNiMuchDHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testNiMuchDHV.txt ; \
	$(RUN_SIMULATE) -P=testNiMuchDHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNiMuchDHV.txt

testAlMuchDHV:
	mkdir -p testAlMuchDHV
	cd testAlMuchDHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" > testAlMuchDHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" >> testAlMuchDHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testAlMuchDHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testAlMuchDHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testAlMuchDHV.txt ; \
	$(RUN_SIMULATE) -P=testAlMuchDHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testAlMuchDHV.txt

testFeMuchDHV:
	mkdir -p testFeMuchDHV
	cd testFeMuchDHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Fe=1.0\";" > testFeMuchDHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Fe=1.0\";" >> testFeMuchDHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testFeMuchDHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testFeMuchDHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testFeMuchDHV.txt ; \
	$(RUN_SIMULATE) -P=testFeMuchDHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testFeMuchDHV.txt

testMnMuchDHV:
	mkdir -p testMnMuchDHV
	cd testMnMuchDHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Mn=1.0\";" > testMnMuchDHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Mn=1.0\";" >> testMnMuchDHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testMnMuchDHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testMnMuchDHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testMnMuchDHV.txt ; \
	$(RUN_SIMULATE) -P=testMnMuchDHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testMnMuchDHV.txt

testLiMuchDHV:
	mkdir -p testLiMuchDHV
	cd testLiMuchDHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Li=1.0\";" > testLiMuchDHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Li=1.0\";" >> testLiMuchDHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testLiMuchDHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testLiMuchDHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiMuchDHV.txt ; \
	$(RUN_SIMULATE) -P=testLiMuchDHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiMuchDHV.txt


testHV: testNiLessProtonHV testAlLessProtonHV testNiPoorProtonHV testAlPoorProtonHV testFePoorProtonHV testMnPoorProtonHV testLiPoorProtonHV testNiMuchProtonHV testAlMuchProtonHV testFeMuchProtonHV testMnMuchProtonHV testLiMuchProtonHV testNiMuchDHV testAlMuchDHV testFeMuchDHV testMnMuchDHV testLiMuchDHV     

testLiLessProton:
	mkdir -p testLiLessProton
	cd testLiLessProton ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Li=1.0\";" > testLiLessProton.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiLessProton.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiLessProton.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiLessProton.txt ; \
	$(RUN_SIMULATE) -P=testLiLessProton.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiLessProton.txt

testLiLessProtonHV:
	mkdir -p testLiLessProtonHV
	cd testLiLessProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Li=1.0\";" > testLiLessProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiLessProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiLessProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiLessProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiLessProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testLiLessProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiLessProtonHV.txt

testLiLessMProtonHV:
	mkdir -p testLiLessMProtonHV
	cd testLiLessMProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Li=1.0\";" > testLiLessMProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiLessMProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-11;" >>testLiLessMProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-12;" >>testLiLessMProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiLessMProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testLiLessMProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiLessMProtonHV.txt

testLiNiLessProtonHV:
	mkdir -p testLiNiLessProtonHV
	cd testLiNiLessProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Li=1.0\";" > testLiNiLessProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testLiNiLessProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiNiLessProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiNiLessProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiNiLessProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testLiNiLessProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiNiLessProtonHV.txt

testLiNiLessProtonB10:
	mkdir -p testLiNiLessProtonB
	cd testLiNiLessProtonB ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Li=1.0\";" > testLiNiLessProtonB.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testLiNiLessProtonB.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiNiLessProtonB.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiNiLessProtonB.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiNiLessProtonB.txt ; \
	$(RUN_SIMULATE) -b=10 -P=testLiNiLessProtonB.txt ; \
	$(RUN_SIMULATE) -b=10 -t=1h -P=testLiNiLessProtonB.txt

testLiNiLessProtonB100:
	mkdir -p testLiNiLessProtonB100
	cd testLiNiLessProtonB100 ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Li=1.0\";" > testLiNiLessProtonB100.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Ni=1.0\";" >> testLiNiLessProtonB100.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiNiLessProtonB100.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiNiLessProtonB100.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiNiLessProtonB100.txt ; \
	$(RUN_SIMULATE) -b=100 -P=testLiNiLessProtonB100.txt ; \
	$(RUN_SIMULATE) -b=100 -t=24h -P=testLiNiLessProtonB100.txt

testLiNiXLessProtonHV:
	mkdir -p testLiNiXLessProtonHV
	cd testLiNiXLessProtonHV ; \
	rm -f *.dat *.log test*.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testLiNiXLessProtonHV.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.1, Li=1.0\";" >> testLiNiXLessProtonHV.txt ; \
	echo "double e_emittedElectronMol = 4e-12;" >>testLiNiXLessProtonHV.txt ; \
	echo "double e_emittedProtonMol = 4e-11;" >>testLiNiXLessProtonHV.txt ; \
	echo "double e_appliedVoltageScale = 2.4;" >>testLiNiXLessProtonHV.txt ; \
	$(RUN_SIMULATE) -P=testLiNiXLessProtonHV.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiNiXLessProtonHV.txt

# ---------------- old targets ----------------
# testNistopcur is :
# (1) run 'sumulate' with stopping the high voltage charge and with reducing hydrogen in the new directry 'testNistopcur' after helitaging 'testNi'.
testNistopcur:
	mkdir -p testNistopcur
	cd testNistopcur ; \
	cp -f ../testNi/data.dat . ; \
	cp -f ../testNi/*.log . ; \
	echo "double e_emittedElectronMol = 0.0;" >testNistopcur.txt ; \
	echo "double e_emittedProtonMol = 0.0;" >>testNistopcur.txt ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcur.txt ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcur.txt ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcur.txt ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcur.txt 

# testNistopcurH is :
# (1) run 'sumulate' with stopping the high voltage charge and with reducing hydrogen in the new directry 'testNistopcurH' after helitaging 'testNi'.
testNistopcurH:
	mkdir -p testNistopcurH
	cd testNistopcurH ; \
	cp -f ../testNi/data.dat . ; \
	cp -f ../testNi/*.log . ; \
	echo "double e_emittedElectronMol = 0.0;" >testNistopcurH.txt ; \
	echo "double e_emittedProtonMol = 0.0;" >>testNistopcurH.txt ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcurH.txt -H- ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcurH.txt -H- ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcurH.txt ; \
	$(RUN_SIMULATE) -t=24h -P=testNistopcurH.txt 

# testNi_10days is :
# (1) run 'sumulate' during 31 days in the new directry 'testNi_10days' after helitaging 'testNi'.
testNi_10days:
	mkdir -p testNi_10days
	cd testNi_10days ; \
	cp -f ../testNi/data.dat . ; \
	cp -f ../testNi/*.log . ; \
	$(RUN_SIMULATE) -t=2d ; \
	$(RUN_SIMULATE) -t=2d ; \
	$(RUN_SIMULATE) -t=2d ; \
	$(RUN_SIMULATE) -t=2d 

# testNi_woCE is :
# (1) run 'sumulate' while 1 minute + some hours without ComptonEffect.
testNi_woCE:
	mkdir -p testNi_woCE
	cd testNi_woCE ; \
	rm -f data.dat s00*.log m*.log h*.log param2.dat ; \
	echo "int e_useComptonEffect = 0;" >param2.dat ; \
	$(RUN_SIMULATE) -P=param2.dat ; \
	$(RUN_SIMULATE) -t=1h -P=param2.dat ; \
	$(RUN_SIMULATE) -t=1d -P=param2.dat

# test_woH is :
# (1) run 'sumulate' while (1 minute + some hours) without hydrogen.
test_woH:
	mkdir -p test_woH
	cd test_woH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log param3.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > param3.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" >> param3.dat ; \
	echo "double e_emittedProtonMol = 0.0;" >> param3.dat ; \
	$(RUN_SIMULATE) -P=param3.dat ; \
	$(RUN_SIMULATE) -t=1h -P=param3.dat ; \
	$(RUN_SIMULATE) -t=1d -P=param3.dat 


# testNi_cur100 is :
# (1) run 'sumulate' with rich electric current in the new directry 'testNi_cur100' .
testNi_cur100:
	mkdir -p testNi_cur100
	cd testNi_cur100 ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >param4.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>param4.dat ; \
	$(RUN_SIMULATE) -P=param4.dat ; \
	$(RUN_SIMULATE) -t=1h -P=param4.dat ; \
	$(RUN_SIMULATE) -t=1d -P=param4.dat 


# testNi_rate4 is :
# (1) run 'sumulate' with high rate (* 4) in the new directry 'testNi_rate4' .
testNi_rate4:
	mkdir -p testNi_rate4
	cd testNi_rate4 ; \
	echo "double e_collideElectronRateOnElectrode = 0.8;" >param5.dat ; \
	echo "double e_collideElectronRateForMidiMeV = 0.4;" >>param5.dat ; \
	echo "double e_collideElectronRateForMiniMeV = 0.04;" >>param5.dat ; \
	echo "double e_neutronGenInSpaceProtonRate = 0.04;" >>param5.dat ; \
	echo "double e_hydrogenGenInSpaceProtonRate = 0.04;" >>param5.dat ; \
	echo "double e_collideProtonRateOnElectrode = 0.8;" >>param5.dat ; \
	$(RUN_SIMULATE) -P=param5.dat ; \
	$(RUN_SIMULATE) -t=1h -P=param5.dat ; \
	$(RUN_SIMULATE) -t=1d -P=param5.dat 


# testNi_rate025 is :
# (1) run 'sumulate' with high rate (* 0.25) in the new directry 'testNi_rate025' .
testNi_rate025:
	mkdir -p testNi_rate025
	cd testNi_rate025 ; \
	echo "double e_collideElectronRateOnElectrode = 0.05;" >param8.dat ; \
	echo "double e_collideElectronRateForMidiMeV = 0.025;" >>param8.dat ; \
	echo "double e_collideElectronRateForMiniMeV = 0.0025;" >>param8.dat ; \
	echo "double e_neutronGenInSpaceProtonRate = 0.0025;" >>param8.dat ; \
	echo "double e_hydrogenGenInSpaceProtonRate = 0.0025;" >>param8.dat ; \
	echo "double e_collideProtonRateOnElectrode = 0.05;" >>param8.dat ; \
	$(RUN_SIMULATE) -P=param8.dat ; \
	$(RUN_SIMULATE) -t=1h -P=param8.dat ; \
	$(RUN_SIMULATE) -t=1d -P=param8.dat

testNi_rateNeutron:
	mkdir -p testNi_rateNeutron
	cd testNi_rateNeutron ; \
	echo "double e_collideElectronRateOnElectrode = 0.2;" >param_rateNeutron.dat ; \
	echo "double e_collideElectronRateForMidiMeV = 0.1;" >>param_rateNeutron.dat ; \
	echo "double e_collideElectronRateForMiniMeV = 0.01;" >>param_rateNeutron.dat ; \
	echo "double e_neutronGenInSpaceProtonRate = 0.1;" >>param_rateNeutron.dat ; \
	echo "double e_hydrogenGenInSpaceProtonRate = 0.1;" >>param_rateNeutron.dat ; \
	echo "double e_collideProtonRateOnElectrode = 0.2;" >>param_rateNeutron.dat ; \
	$(RUN_SIMULATE) -P=param_rateNeutron.dat ; \
	$(RUN_SIMULATE) -t=1h -P=param_rateNeutron.dat ; \
	$(RUN_SIMULATE) -t=1d -P=param_rateNeutron.dat 

# Deuterium(D, heavy hydrogen) #1
testD:
	mkdir -p testD
	cd testD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramD.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramD.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.6, Ni=1.0\";" >> paramD.dat ; \
	$(RUN_SIMULATE) -P=paramD.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramD.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramD.dat 


testDH:
	mkdir -p testDH
	cd testDH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testDH.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH.dat ; \
	$(RUN_SIMULATE) -P=testDH.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testDH.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testDH.dat 

testDH3M:
	mkdir -p testDH3M
	cd testDH3M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH3M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testDH3M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH3M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testDH3M.dat ; \
	$(RUN_SIMULATE) -P=testDH3M.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testDH3M.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testDH3M.dat 

testD2:
	mkdir -p testD2
	cd testD2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramD2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.3, Ni=1.0\";" > paramD2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.6, Ni=1.0\";" >> paramD2.dat ; \
	$(RUN_SIMULATE) -P=paramD2.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramD2.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramD2.dat 

testDH2:
	mkdir -p testDH2
	cd testDH2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.15, D=0.15, Ni=1.0\";" > testDH2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH2.dat ; \
	$(RUN_SIMULATE) -P=testDH2.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testDH2.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testDH2.dat 

testDH23M:
	mkdir -p testDH23M
	cd testDH23M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH23M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.15, D=0.15, Ni=1.0\";" > testDH23M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH23M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testDH23M.dat ; \
	$(RUN_SIMULATE) -P=testDH23M.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testDH23M.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testDH23M.dat 

# Deuterium(D, heavy hydrogen) 1%
testD1_Ni:
	mkdir -p testD1_Ni
	cd testD1_Ni ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testD1_Ni.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testD1_Ni.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.6, D=0.006, Ni=1.0\";" >> testD1_Ni.txt ; \
	$(RUN_SIMULATE) -P=testD1_Ni.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testD1_Ni.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testD1_Ni.txt 

testD1_Fe:
	mkdir -p testD1_Fe
	cd testD1_Fe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testD1_Fe.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testD1_Fe.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.4, D=0.004, Fe=1.0\";" >> testD1_Fe.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testD1_Fe.txt ; \
	$(RUN_SIMULATE) -P=testD1_Fe.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testD1_Fe.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testD1_Fe.txt 

testD1_Al:
	mkdir -p testD1_Al
	cd testD1_Al ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testD1_Al.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testD1_Al.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.2, D=0.002, Al=1.0\";" >> testD1_Al.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testD1_Al.dat ; \
	$(RUN_SIMULATE) -P=testD1_Al.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testD1_Al.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testD1_Al.dat 

# Li 1%
testLi1_Ni:
	mkdir -p testLi1_Ni
	cd testLi1_Ni ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi1_Ni.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testLi1_Ni.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, Li=0.006, Ni=1.0\";" >> testLi1_Ni.txt ; \
	$(RUN_SIMULATE) -P=testLi1_Ni.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLi1_Ni.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLi1_Ni.txt 

testLi1_Fe:
	mkdir -p testLi1_Fe
	cd testLi1_Fe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi1_Fe.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testLi1_Fe.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Li=0.004, Fe=1.0\";" >> testLi1_Fe.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testLi1_Fe.txt ; \
	$(RUN_SIMULATE) -P=testLi1_Fe.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLi1_Fe.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLi1_Fe.txt 

testLi1_Al:
	mkdir -p testLi1_Al
	cd testLi1_Al ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi1_Al.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLi1_Al.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=0.002, Al=1.0\";" >> testLi1_Al.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLi1_Al.dat ; \
	$(RUN_SIMULATE) -P=testLi1_Al.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testLi1_Al.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testLi1_Al.dat 

# Na 1%
testNa1_Ni:
	mkdir -p testNa1_Ni
	cd testNa1_Ni ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testNa1_Ni.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testNa1_Ni.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, Na=0.006, Ni=1.0\";" >> testNa1_Ni.txt ; \
	$(RUN_SIMULATE) -P=testNa1_Ni.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNa1_Ni.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testNa1_Ni.txt 

testNa1_Fe:
	mkdir -p testNa1_Fe
	cd testNa1_Fe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testNa1_Fe.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testNa1_Fe.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Na=0.004, Fe=1.0\";" >> testNa1_Fe.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testNa1_Fe.txt ; \
	$(RUN_SIMULATE) -P=testNa1_Fe.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testNa1_Fe.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testNa1_Fe.txt 

testNa1_Al:
	mkdir -p testNa1_Al
	cd testNa1_Al ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testNa1_Al.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testNa1_Al.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Na=0.002, Al=1.0\";" >> testNa1_Al.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testNa1_Al.dat ; \
	$(RUN_SIMULATE) -P=testNa1_Al.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testNa1_Al.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testNa1_Al.dat 

# Mercury(Hg) #80
# testHg is :(1) run 'sumulate' while 1 minute + some hours with mercury(Hg).
testHg196:
	mkdir -p testHg196
	cd testHg196 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHg196.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testHg196.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, 196Hg=0.015, 198Hg=0.997\";" >> testHg196.txt ; \
	$(RUN_SIMULATE) -P=testHg196.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testHg196.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testHg196.txt 

testHgDAlcur100:
	mkdir -p testHgDAlcur100
	cd testHgDAlcur100 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHgDAlcur100.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testHgDAlcur100.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 196Hg=0.015, 198Hg=0.997\";" >> testHgDAlcur100.dat ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testHgDAlcur100.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>testHgDAlcur100.dat ; \
	$(RUN_SIMULATE) -P=testHgDAlcur100.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testHgDAlcur100.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testHgDAlcur100.dat 

testHgDAl3Mcur100:
	mkdir -p testHgDAl3Mcur100
	cd testHgDAl3Mcur100 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHgDAl3Mcur100.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testHgDAl3Mcur100.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 196Hg=0.015, 198Hg=0.997\";" >> testHgDAl3Mcur100.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testHgDAl3Mcur100.dat ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testHgDAl3Mcur100.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>testHgDAl3Mcur100.dat ; \
	$(RUN_SIMULATE) -P=testHgDAl3Mcur100.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testHgDAl3Mcur100.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testHgDAl3Mcur100.dat 

testHgDAl23Mcur100:
	mkdir -p testHgDAl23Mcur100
	cd testHgDAl23Mcur100 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHgDAl23Mcur100.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.05, D=0.05, Al=1.0\";" > testHgDAl23Mcur100.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 196Hg=0.015, 198Hg=0.997\";" >> testHgDAl23Mcur100.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testHgDAl23Mcur100.dat ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testHgDAl23Mcur100.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>testHgDAl23Mcur100.dat ; \
	$(RUN_SIMULATE) -P=testHgDAl23Mcur100.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testHgDAl23Mcur100.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testHgDAl23Mcur100.dat 


# Iron(Fe) #26 
testFe:
	mkdir -p testFe
	cd testFe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramFe.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > paramFe.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Fe=1.0\";" >> paramFe.dat ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> paramFe.dat ; \
	$(RUN_SIMULATE) -P=paramFe.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramFe.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramFe.dat 

# Aluminum(Al) #13 has a largest COP = 3.51 after 1 day.
testAl:
	mkdir -p testAl
	cd testAl ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAl.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAl.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Al=1.0\";" >> testAl.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAl.dat ; \
	$(RUN_SIMULATE) -P=testAl.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAl.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAl.dat 

testAl_cur100_10days:
	mkdir -p testAl_cur100_10days
	cd testAl_cur100_10days ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAl_cur100_10days.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAl_cur100_10days.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Al=1.0\";" >> testAl_cur100_10days.txt ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testAl_cur100_10days.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-8;" >> testAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -P=testAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testAl_cur100_10days.txt  ; \
	$(RUN_SIMULATE) -t=2d -P=testAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=2d -P=testAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=2d -P=testAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=2d -P=testAl_cur100_10days.txt

testAlD:
	mkdir -p testAlD
	cd testAlD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlD.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlD.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.2, Al=1.0\";" >> testAlD.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlD.dat ; \
	$(RUN_SIMULATE) -P=testAlD.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAlD.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAlD.dat 

testAlD2:
	mkdir -p testAlD2
	cd testAlD2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlD2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" > testAlD2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.2, Al=1.0\";" >> testAlD2.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlD2.dat ; \
	$(RUN_SIMULATE) -P=testAlD2.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAlD2.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAlD2.dat 

testAlHD:
	mkdir -p testAlHD
	cd testAlHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlHD.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD.dat ; \
	$(RUN_SIMULATE) -P=testAlHD.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAlHD.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAlHD.dat 

testAlHD2:
	mkdir -p testAlHD2
	cd testAlHD2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.05, D=0.05, Al=1.0\";" > testAlHD2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD2.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD2.dat ; \
	$(RUN_SIMULATE) -P=testAlHD2.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAlHD2.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAlHD2.dat 

testAlHD3M:
	mkdir -p testAlHD3M
	cd testAlHD3M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD3M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlHD3M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD3M.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD3M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testAlHD3M.dat ; \
	$(RUN_SIMULATE) -P=testAlHD3M.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAlHD3M.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAlHD3M.dat 

testAlHD23M:
	mkdir -p testAlHD23M
	cd testAlHD23M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD23M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.05, D=0.05, Al=1.0\";" > testAlHD23M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD23M.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD23M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> paramAlD.dat ; \
	$(RUN_SIMULATE) -P=testAlHD23M.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testAlHD23M.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testAlHD23M.dat 

#----------------------------------------------------------
testLi:
	mkdir -p testLi
	cd testLi ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramLi.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramLi.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> paramLi.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramLi.dat ; \
	$(RUN_SIMULATE) -P=paramLi.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramLi.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramLi.dat 

testLiAl:
	mkdir -p testLiAl
	cd testLiAl ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> testLiAl.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLiAl.txt ; \
	$(RUN_SIMULATE) -P=testLiAl.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiAl.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLiAl.txt 

testLiAl6:
	mkdir -p testLiAl6
	cd testLiAl6 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl6.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl6.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, 6Li=1.0\";" >> testLiAl6.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLiAl6.txt ; \
	$(RUN_SIMULATE) -P=testLiAl6.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiAl6.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLiAl6.txt 

testLiAl7:
	mkdir -p testLiAl7
	cd testLiAl7 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl7.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl7.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, 7Li=1.0\";" >> testLiAl7.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLiAl7.txt ; \
	$(RUN_SIMULATE) -P=testLiAl7.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testLiAl7.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLiAl7.txt 

testLi_Short:
	mkdir -p testLi_Short
	cd testLi_Short ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi_Short.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testLi_Short.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> testLi_Short.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLi_Short.txt ; \
	$(RUN_SIMULATE) -t=120s -i=1s -s=1 -n1 -P=testLi_Short.txt; \
	$(RUN_ANALYZELOG) -Psec -D=1 -M=1 s00_m02s00_1_1_60.log >ana.log


testLiAl_cur100_10days:
	mkdir -p testLiAl_cur100_10days
	cd testLiAl_cur100_10days ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl_cur100_10days.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl_cur100_10days.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> testLiAl_cur100_10days.txt ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testLiAl_cur100_10days.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-8;" >> testLiAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -P=testLiAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLiAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testLiAl_cur100_10days.txt  ; \
	$(RUN_SIMULATE) -t=2d -P=testLiAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=2d -P=testLiAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=2d -P=testLiAl_cur100_10days.txt ; \
	$(RUN_SIMULATE) -t=2d -P=testLiAl_cur100_10days.txt

testCa:
	mkdir -p testCa
	cd testCa ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramCa.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramCa.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ca=1.0\";" >> paramCa.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramCa.dat ; \
	$(RUN_SIMULATE) -P=paramCa.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramCa.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramCa.dat 

testNa:
	mkdir -p testNa
	cd testNa ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramNa.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramNa.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Na=1.0\";" >> paramNa.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramNa.dat ; \
	$(RUN_SIMULATE) -P=paramNa.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramNa.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramNa.dat 

testK:
	mkdir -p testK
	cd testK ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramK.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramK.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, K=1.0\";" >> paramK.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramK.dat ; \
	$(RUN_SIMULATE) -P=paramK.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramK.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramK.dat 

testMg:
	mkdir -p testMg
	cd testMg ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramMg.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramMg.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Mg=1.0\";" >> paramMg.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramMg.dat ; \
	$(RUN_SIMULATE) -P=paramMg.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramMg.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramMg.dat 

testTi:
	mkdir -p testTi
	cd testTi ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramTi.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ti=1.0\";" > paramTi.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ti=1.0\";" >> paramTi.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramTi.dat ; \
	$(RUN_SIMULATE) -P=paramTi.dat ; \
	$(RUN_SIMULATE) -t=1h -P=paramTi.dat ; \
	$(RUN_SIMULATE) -t=1d -P=paramTi.dat 

#----------------------------------------------------------
testBe:
	mkdir -p testBe
	cd testBe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testBe.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Be=1.0\";" > testBe.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Be=1.0\";" >> testBe.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testBe.dat ; \
	$(RUN_SIMULATE) -P=testBe.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testBe.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testBe.dat 

testSc:
	mkdir -p testSc
	cd testSc ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSc.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Sc=1.0\";" > testSc.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Sc=1.0\";" >> testSc.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testSc.dat ; \
	$(RUN_SIMULATE) -P=testSc.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testSc.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testSc.dat 

testV:
	mkdir -p testV
	cd testV ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testV.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"V=1.0\";" > testV.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, V=1.0\";" >> testV.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testV.dat ; \
	$(RUN_SIMULATE) -P=testV.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testV.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testV.dat 

testCr:
	mkdir -p testCr
	cd testCr ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testCr.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Cr=1.0\";" > testCr.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Cr=1.0\";" >> testCr.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testCr.dat ; \
	$(RUN_SIMULATE) -P=testCr.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testCr.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testCr.dat 

testMn:
	mkdir -p testMn
	cd testMn ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testMn.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Mn=1.0\";" > testMn.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Mn=1.0\";" >> testMn.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testMn.dat ; \
	$(RUN_SIMULATE) -P=testMn.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testMn.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testMn.dat 

testCo:
	mkdir -p testCo
	cd testCo ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testCo.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Co=1.0\";" > testCo.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Co=1.0\";" >> testCo.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testCo.dat ; \
	$(RUN_SIMULATE) -P=testCo.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testCo.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testCo.dat 

testCu:
	mkdir -p testCu
	cd testCu ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testCu.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Cu=1.0\";" > testCu.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Cu=1.0\";" >> testCu.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testCu.dat ; \
	$(RUN_SIMULATE) -P=testCu.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testCu.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testCu.dat 

testZn:
	mkdir -p testZn
	cd testZn ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testZn.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testZn.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Zn=1.0\";" >> testZn.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testZn.dat ; \
	$(RUN_SIMULATE) -P=testZn.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testZn.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testZn.dat 

testGa:
	mkdir -p testGa
	cd testGa ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testGa.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testGa.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ga=1.0\";" >> testGa.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testGa.dat ; \
	$(RUN_SIMULATE) -P=testGa.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testGa.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testGa.dat 

testGe:
	mkdir -p testGe
	cd testGe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testGe.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ge=1.0\";" > testGe.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ge=1.0\";" >> testGe.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testGe.dat ; \
	$(RUN_SIMULATE) -P=testGe.dat ; \
	$(RUN_SIMULATE) -t=1h -P=testGe.dat ; \
	$(RUN_SIMULATE) -t=1d -P=testGe.dat 

#----------------------------------------------------------
#super high voltage for AL, Fe, Ni
testSHV_AlH:
	mkdir -p testSHV_AlH
	cd testSHV_AlH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_AlH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testSHV_AlH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Al=1.0\";" >> testSHV_AlH.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testSHV_AlH.txt ; \
	echo "double e_appliedVoltageScale = 5.97;" >> testSHV_AlH.txt ; \
	$(RUN_SIMULATE) -P=testSHV_AlH.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_AlH.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_AlH.txt

testSHV_AlHD:
	mkdir -p testSHV_AlHD
	cd testSHV_AlHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_AlHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testSHV_AlHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testSHV_AlHD.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testSHV_AlHD.txt ; \
	echo "double e_appliedVoltageScale = 5.97;" >> testSHV_AlHD.txt ; \
	$(RUN_SIMULATE) -P=testSHV_AlHD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_AlHD.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_AlHD.txt

testSHV_FeH:
	mkdir -p testSHV_FeH
	cd testSHV_FeH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_FeH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testSHV_FeH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Fe=1.0\";" >> testSHV_FeH.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_FeH.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_FeH.txt ; \
	$(RUN_SIMULATE) -P=testSHV_FeH.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_FeH.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_FeH.txt

testSHV_FeHD:
	mkdir -p testSHV_FeHD
	cd testSHV_FeHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_FeHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testSHV_FeHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.2, D=0.2, Fe=1.0\";" >> testSHV_FeHD.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_FeHD.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_FeHD.txt ; \
	$(RUN_SIMULATE) -P=testSHV_FeHD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_FeHD.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_FeHD.txt

testSHV_56FeH:
	mkdir -p testSHV_56FeH
	cd testSHV_56FeH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_56FeH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"56Fe=1.0\";" > testSHV_56FeH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, 56Fe=1.0\";" >> testSHV_56FeH.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_56FeH.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_56FeH.txt ; \
	$(RUN_SIMULATE) -P=testSHV_56FeH.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_56FeH.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_56FeH.txt

testSHV_56FeHD:
	mkdir -p testSHV_56FeHD
	cd testSHV_56FeHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_56FeHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"56Fe=1.0\";" > testSHV_56FeHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.2, D=0.2, 56Fe=1.0\";" >> testSHV_56FeHD.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_56FeHD.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_56FeHD.txt ; \
	$(RUN_SIMULATE) -P=testSHV_56FeHD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_56FeHD.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_56FeHD.txt

testSHV_NiH:
	mkdir -p testSHV_NiH
	cd testSHV_NiH ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_NiH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testSHV_NiH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, Ni=1.0\";" >> testSHV_NiH.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_NiH.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_NiH.txt ; \
	$(RUN_SIMULATE) -P=testSHV_NiH.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_NiH.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_NiH.txt

testSHV_NiHD:
	mkdir -p testSHV_NiHD
	cd testSHV_NiHD ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_NiHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testSHV_NiHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testSHV_NiHD.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_NiHD.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_NiHD.txt ; \
	$(RUN_SIMULATE) -P=testSHV_NiHD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_NiHD.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_NiHD.txt

testSHV_62NiH:
	mkdir -p testSHV_62NiH
	cd testSHV_62NiH ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_62NiH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"62Ni=1.0\";" > testSHV_62NiH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, 62Ni=1.0\";" >> testSHV_62NiH.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_62NiH.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_62NiH.txt ; \
	$(RUN_SIMULATE) -P=testSHV_62NiH.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_62NiH.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_62NiH.txt

testSHV_62NiHD:
	mkdir -p testSHV_62NiHD
	cd testSHV_62NiHD ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_62NiHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"62Ni=1.0\";" > testSHV_62NiHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 62Ni=1.0\";" >> testSHV_62NiHD.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_62NiHD.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_62NiHD.txt ; \
	$(RUN_SIMULATE) -P=testSHV_62NiHD.txt ; \
	$(RUN_SIMULATE) -t=1h -P=testSHV_62NiHD.txt ; \
	$(RUN_SIMULATE) -t=1d -P=testSHV_62NiHD.txt

#----------------------------------------------------------
alltestA:testNi testNi_woCE test_woH testNi_cur100 testNi_10days

alltestStop: testNistopcur testNistopcurH

alltestC:testNi_rate4 testNi_rate025 testNi_rateNeutron

alltestD:testD testD2 testDH testDH2 testDH3M testDH23M

alltestX:testAl testAlD testAlD2 testAlHD testAlHD2 testAlHD3M testAlHD23M

alltestX2:testFe testLi testLiAl testCa testNa testK testMg testTi
 
alltestX3:testHg testHgDAlcur100 testHgDAl3Mcur100 testHgDAl23Mcur100

alltestX4: testBe testSc testV testCr testMn testCo testCu testZn testGa testGe

allSHV1: testSHV_AlH  testSHV_56FeH  testSHV_62NiH  testSHV_FeH  testSHV_NiH

allSHV2: testSHV_AlHD testSHV_56FeHD testSHV_62NiHD testSHV_FeHD testSHV_NiHD

all10: testAl_cur100_10days testLiAl_cur100_10days testNi_10days 