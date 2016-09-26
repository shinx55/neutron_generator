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

# testNi is :
# (1) run 'sumulate' while 1 minute + many hour with default contitions.
testNi:
	mkdir -p testNi
	cd testNi ; \
	rm -f data.dat s00*.log m*.log h*.log param1.dat ; \
	../simulate ; \
	../simulate -t=1h ; \
	../simulate -t=1d ; \
	../simulate -t=1d

# testNistopcur is :
# (1) run 'sumulate' with stopping the high voltage charge and with reducing hydrogen in the new directry 'testNistopcur' after helitaging 'testNi'.
testNistopcur:
	mkdir -p testNistopcur
	cd testNistopcur ; \
	cp -f ../testNi/data.dat . ; \
	cp -f ../testNi/*.log . ; \
	echo "double e_emittedElectronMol = 0.0;" >testNistopcur.txt ; \
	echo "double e_emittedProtonMol = 0.0;" >>testNistopcur.txt ; \
	../simulate -t=24h -P=testNistopcur.txt ; \
	../simulate -t=24h -P=testNistopcur.txt ; \
	../simulate -t=24h -P=testNistopcur.txt ; \
	../simulate -t=24h -P=testNistopcur.txt 

# testNistopcurH is :
# (1) run 'sumulate' with stopping the high voltage charge and with reducing hydrogen in the new directry 'testNistopcurH' after helitaging 'testNi'.
testNistopcurH:
	mkdir -p testNistopcurH
	cd testNistopcurH ; \
	cp -f ../testNi/data.dat . ; \
	cp -f ../testNi/*.log . ; \
	echo "double e_emittedElectronMol = 0.0;" >testNistopcurH.txt ; \
	echo "double e_emittedProtonMol = 0.0;" >>testNistopcurH.txt ; \
	../simulate -t=24h -P=testNistopcurH.txt -H- ; \
	../simulate -t=24h -P=testNistopcurH.txt -H- ; \
	../simulate -t=24h -P=testNistopcurH.txt ; \
	../simulate -t=24h -P=testNistopcurH.txt 

# testNi_10days is :
# (1) run 'sumulate' during 31 days in the new directry 'testNi_10days' after helitaging 'testNi'.
testNi_10days:
	mkdir -p testNi_10days
	cd testNi_10days ; \
	cp -f ../testNi/data.dat . ; \
	cp -f ../testNi/*.log . ; \
	../simulate -t=2d ; \
	../simulate -t=2d ; \
	../simulate -t=2d ; \
	../simulate -t=2d 

# testNi_woCE is :
# (1) run 'sumulate' while 1 minute + some hours without ComptonEffect.
testNi_woCE:
	mkdir -p testNi_woCE
	cd testNi_woCE ; \
	rm -f data.dat s00*.log m*.log h*.log param2.dat ; \
	echo "int e_useComptonEffect = 0;" >param2.dat ; \
	../simulate -P=param2.dat ; \
	../simulate -t=1h -P=param2.dat ; \
	../simulate -t=1d -P=param2.dat

# test_woH is :
# (1) run 'sumulate' while (1 minute + some hours) without hydrogen.
test_woH:
	mkdir -p test_woH
	cd test_woH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log param3.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > param3.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" >> param3.dat ; \
	echo "double e_emittedProtonMol = 0.0;" >> param3.dat ; \
	../simulate -P=param3.dat ; \
	../simulate -t=1h -P=param3.dat ; \
	../simulate -t=1d -P=param3.dat 


# testNi_cur100 is :
# (1) run 'sumulate' with rich electric current in the new directry 'testNi_cur100' .
testNi_cur100:
	mkdir -p testNi_cur100
	cd testNi_cur100 ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >param4.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>param4.dat ; \
	../simulate -P=param4.dat ; \
	../simulate -t=1h -P=param4.dat ; \
	../simulate -t=1d -P=param4.dat 


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
	../simulate -P=param5.dat ; \
	../simulate -t=1h -P=param5.dat ; \
	../simulate -t=1d -P=param5.dat 


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
	../simulate -P=param8.dat ; \
	../simulate -t=1h -P=param8.dat ; \
	../simulate -t=1d -P=param8.dat

testNi_rateNeutron:
	mkdir -p testNi_rateNeutron
	cd testNi_rateNeutron ; \
	echo "double e_collideElectronRateOnElectrode = 0.2;" >param_rateNeutron.dat ; \
	echo "double e_collideElectronRateForMidiMeV = 0.1;" >>param_rateNeutron.dat ; \
	echo "double e_collideElectronRateForMiniMeV = 0.01;" >>param_rateNeutron.dat ; \
	echo "double e_neutronGenInSpaceProtonRate = 0.1;" >>param_rateNeutron.dat ; \
	echo "double e_hydrogenGenInSpaceProtonRate = 0.1;" >>param_rateNeutron.dat ; \
	echo "double e_collideProtonRateOnElectrode = 0.2;" >>param_rateNeutron.dat ; \
	../simulate -P=param_rateNeutron.dat ; \
	../simulate -t=1h -P=param_rateNeutron.dat ; \
	../simulate -t=1d -P=param_rateNeutron.dat 

# Deuterium(D, heavy hydrogen) #1
testD:
	mkdir -p testD
	cd testD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramD.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramD.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.6, Ni=1.0\";" >> paramD.dat ; \
	../simulate -P=paramD.dat ; \
	../simulate -t=1h -P=paramD.dat ; \
	../simulate -t=1d -P=paramD.dat 


testDH:
	mkdir -p testDH
	cd testDH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testDH.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH.dat ; \
	../simulate -P=testDH.dat ; \
	../simulate -t=1h -P=testDH.dat ; \
	../simulate -t=1d -P=testDH.dat 

testDH3M:
	mkdir -p testDH3M
	cd testDH3M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH3M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testDH3M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH3M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testDH3M.dat ; \
	../simulate -P=testDH3M.dat ; \
	../simulate -t=1h -P=testDH3M.dat ; \
	../simulate -t=1d -P=testDH3M.dat 

testD2:
	mkdir -p testD2
	cd testD2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramD2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.3, Ni=1.0\";" > paramD2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.6, Ni=1.0\";" >> paramD2.dat ; \
	../simulate -P=paramD2.dat ; \
	../simulate -t=1h -P=paramD2.dat ; \
	../simulate -t=1d -P=paramD2.dat 

testDH2:
	mkdir -p testDH2
	cd testDH2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.15, D=0.15, Ni=1.0\";" > testDH2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH2.dat ; \
	../simulate -P=testDH2.dat ; \
	../simulate -t=1h -P=testDH2.dat ; \
	../simulate -t=1d -P=testDH2.dat 

testDH23M:
	mkdir -p testDH23M
	cd testDH23M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testDH23M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.15, D=0.15, Ni=1.0\";" > testDH23M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testDH23M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testDH23M.dat ; \
	../simulate -P=testDH23M.dat ; \
	../simulate -t=1h -P=testDH23M.dat ; \
	../simulate -t=1d -P=testDH23M.dat 

# Deuterium(D, heavy hydrogen) 1%
testD1_Ni:
	mkdir -p testD1_Ni
	cd testD1_Ni ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testD1_Ni.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testD1_Ni.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.6, D=0.006, Ni=1.0\";" >> testD1_Ni.txt ; \
	../simulate -P=testD1_Ni.txt ; \
	../simulate -t=1h -P=testD1_Ni.txt ; \
	../simulate -t=1d -P=testD1_Ni.txt 

testD1_Fe:
	mkdir -p testD1_Fe
	cd testD1_Fe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testD1_Fe.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testD1_Fe.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.4, D=0.004, Fe=1.0\";" >> testD1_Fe.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testD1_Fe.txt ; \
	../simulate -P=testD1_Fe.txt ; \
	../simulate -t=1h -P=testD1_Fe.txt ; \
	../simulate -t=1d -P=testD1_Fe.txt 

testD1_Al:
	mkdir -p testD1_Al
	cd testD1_Al ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testD1_Al.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testD1_Al.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.2, D=0.002, Al=1.0\";" >> testD1_Al.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testD1_Al.dat ; \
	../simulate -P=testD1_Al.dat ; \
	../simulate -t=1h -P=testD1_Al.dat ; \
	../simulate -t=1d -P=testD1_Al.dat 

# Li 1%
testLi1_Ni:
	mkdir -p testLi1_Ni
	cd testLi1_Ni ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi1_Ni.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testLi1_Ni.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, Li=0.006, Ni=1.0\";" >> testLi1_Ni.txt ; \
	../simulate -P=testLi1_Ni.txt ; \
	../simulate -t=1h -P=testLi1_Ni.txt ; \
	../simulate -t=1d -P=testLi1_Ni.txt 

testLi1_Fe:
	mkdir -p testLi1_Fe
	cd testLi1_Fe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi1_Fe.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testLi1_Fe.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Li=0.004, Fe=1.0\";" >> testLi1_Fe.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testLi1_Fe.txt ; \
	../simulate -P=testLi1_Fe.txt ; \
	../simulate -t=1h -P=testLi1_Fe.txt ; \
	../simulate -t=1d -P=testLi1_Fe.txt 

testLi1_Al:
	mkdir -p testLi1_Al
	cd testLi1_Al ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLi1_Al.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLi1_Al.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=0.002, Al=1.0\";" >> testLi1_Al.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLi1_Al.dat ; \
	../simulate -P=testLi1_Al.dat ; \
	../simulate -t=1h -P=testLi1_Al.dat ; \
	../simulate -t=1d -P=testLi1_Al.dat 

# Na 1%
testNa1_Ni:
	mkdir -p testNa1_Ni
	cd testNa1_Ni ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testNa1_Ni.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testNa1_Ni.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, Na=0.006, Ni=1.0\";" >> testNa1_Ni.txt ; \
	../simulate -P=testNa1_Ni.txt ; \
	../simulate -t=1h -P=testNa1_Ni.txt ; \
	../simulate -t=1d -P=testNa1_Ni.txt 

testNa1_Fe:
	mkdir -p testNa1_Fe
	cd testNa1_Fe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testNa1_Fe.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testNa1_Fe.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Na=0.004, Fe=1.0\";" >> testNa1_Fe.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testNa1_Fe.txt ; \
	../simulate -P=testNa1_Fe.txt ; \
	../simulate -t=1h -P=testNa1_Fe.txt ; \
	../simulate -t=1d -P=testNa1_Fe.txt 

testNa1_Al:
	mkdir -p testNa1_Al
	cd testNa1_Al ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testNa1_Al.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testNa1_Al.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Na=0.002, Al=1.0\";" >> testNa1_Al.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testNa1_Al.dat ; \
	../simulate -P=testNa1_Al.dat ; \
	../simulate -t=1h -P=testNa1_Al.dat ; \
	../simulate -t=1d -P=testNa1_Al.dat 

# Mercury(Hg) #80
# testHg is :(1) run 'sumulate' while 1 minute + some hours with mercury(Hg).
testHg196:
	mkdir -p testHg196
	cd testHg196 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHg196.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testHg196.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, 196Hg=0.015, 198Hg=0.997\";" >> testHg196.txt ; \
	../simulate -P=testHg196.txt ; \
	../simulate -t=1h -P=testHg196.txt ; \
	../simulate -t=1d -P=testHg196.txt 

testHgDAlcur100:
	mkdir -p testHgDAlcur100
	cd testHgDAlcur100 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHgDAlcur100.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testHgDAlcur100.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 196Hg=0.015, 198Hg=0.997\";" >> testHgDAlcur100.dat ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testHgDAlcur100.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>testHgDAlcur100.dat ; \
	../simulate -P=testHgDAlcur100.dat ; \
	../simulate -t=1h -P=testHgDAlcur100.dat ; \
	../simulate -t=1d -P=testHgDAlcur100.dat 

testHgDAl3Mcur100:
	mkdir -p testHgDAl3Mcur100
	cd testHgDAl3Mcur100 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHgDAl3Mcur100.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testHgDAl3Mcur100.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 196Hg=0.015, 198Hg=0.997\";" >> testHgDAl3Mcur100.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testHgDAl3Mcur100.dat ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testHgDAl3Mcur100.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>testHgDAl3Mcur100.dat ; \
	../simulate -P=testHgDAl3Mcur100.dat ; \
	../simulate -t=1h -P=testHgDAl3Mcur100.dat ; \
	../simulate -t=1d -P=testHgDAl3Mcur100.dat 

testHgDAl23Mcur100:
	mkdir -p testHgDAl23Mcur100
	cd testHgDAl23Mcur100 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testHgDAl23Mcur100.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.05, D=0.05, Al=1.0\";" > testHgDAl23Mcur100.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 196Hg=0.015, 198Hg=0.997\";" >> testHgDAl23Mcur100.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testHgDAl23Mcur100.dat ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testHgDAl23Mcur100.dat ; \
	echo "double e_emittedProtonMol = 0.4E-8;" >>testHgDAl23Mcur100.dat ; \
	../simulate -P=testHgDAl23Mcur100.dat ; \
	../simulate -t=1h -P=testHgDAl23Mcur100.dat ; \
	../simulate -t=1d -P=testHgDAl23Mcur100.dat 


# Iron(Fe) #26 
testFe:
	mkdir -p testFe
	cd testFe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramFe.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > paramFe.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Fe=1.0\";" >> paramFe.dat ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> paramFe.dat ; \
	../simulate -P=paramFe.dat ; \
	../simulate -t=1h -P=paramFe.dat ; \
	../simulate -t=1d -P=paramFe.dat 

# Aluminum(Al) #13 has a largest COP = 3.51 after 1 day.
testAl:
	mkdir -p testAl
	cd testAl ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAl.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAl.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Al=1.0\";" >> testAl.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAl.dat ; \
	../simulate -P=testAl.dat ; \
	../simulate -t=1h -P=testAl.dat ; \
	../simulate -t=1d -P=testAl.dat 

testAl_cur100_10days:
	mkdir -p testAl_cur100_10days
	cd testAl_cur100_10days ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAl_cur100_10days.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAl_cur100_10days.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Al=1.0\";" >> testAl_cur100_10days.txt ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testAl_cur100_10days.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-8;" >> testAl_cur100_10days.txt ; \
	../simulate -P=testAl_cur100_10days.txt ; \
	../simulate -t=1d -P=testAl_cur100_10days.txt ; \
	../simulate -t=1d -P=testAl_cur100_10days.txt  ; \
	../simulate -t=2d -P=testAl_cur100_10days.txt ; \
	../simulate -t=2d -P=testAl_cur100_10days.txt ; \
	../simulate -t=2d -P=testAl_cur100_10days.txt ; \
	../simulate -t=2d -P=testAl_cur100_10days.txt

testAlD:
	mkdir -p testAlD
	cd testAlD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlD.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlD.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.2, Al=1.0\";" >> testAlD.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlD.dat ; \
	../simulate -P=testAlD.dat ; \
	../simulate -t=1h -P=testAlD.dat ; \
	../simulate -t=1d -P=testAlD.dat 

testAlD2:
	mkdir -p testAlD2
	cd testAlD2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlD2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.1, Al=1.0\";" > testAlD2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"D=0.2, Al=1.0\";" >> testAlD2.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlD2.dat ; \
	../simulate -P=testAlD2.dat ; \
	../simulate -t=1h -P=testAlD2.dat ; \
	../simulate -t=1d -P=testAlD2.dat 

testAlHD:
	mkdir -p testAlHD
	cd testAlHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlHD.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD.dat ; \
	../simulate -P=testAlHD.dat ; \
	../simulate -t=1h -P=testAlHD.dat ; \
	../simulate -t=1d -P=testAlHD.dat 

testAlHD2:
	mkdir -p testAlHD2
	cd testAlHD2 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD2.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.05, D=0.05, Al=1.0\";" > testAlHD2.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD2.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD2.dat ; \
	../simulate -P=testAlHD2.dat ; \
	../simulate -t=1h -P=testAlHD2.dat ; \
	../simulate -t=1d -P=testAlHD2.dat 

testAlHD3M:
	mkdir -p testAlHD3M
	cd testAlHD3M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD3M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testAlHD3M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD3M.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD3M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> testAlHD3M.dat ; \
	../simulate -P=testAlHD3M.dat ; \
	../simulate -t=1h -P=testAlHD3M.dat ; \
	../simulate -t=1d -P=testAlHD3M.dat 

testAlHD23M:
	mkdir -p testAlHD23M
	cd testAlHD23M ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testAlHD23M.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.05, D=0.05, Al=1.0\";" > testAlHD23M.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testAlHD23M.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testAlHD23M.dat ; \
	echo "double e_appliedVoltageScale = 3.9;" >> paramAlD.dat ; \
	../simulate -P=testAlHD23M.dat ; \
	../simulate -t=1h -P=testAlHD23M.dat ; \
	../simulate -t=1d -P=testAlHD23M.dat 

#----------------------------------------------------------
testLi:
	mkdir -p testLi
	cd testLi ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramLi.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramLi.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> paramLi.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramLi.dat ; \
	../simulate -P=paramLi.dat ; \
	../simulate -t=1h -P=paramLi.dat ; \
	../simulate -t=1d -P=paramLi.dat 

testLiAl:
	mkdir -p testLiAl
	cd testLiAl ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> testLiAl.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLiAl.txt ; \
	../simulate -P=testLiAl.txt ; \
	../simulate -t=1h -P=testLiAl.txt ; \
	../simulate -t=1d -P=testLiAl.txt 

testLiAl6:
	mkdir -p testLiAl6
	cd testLiAl6 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl6.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl6.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, 6Li=1.0\";" >> testLiAl6.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLiAl6.txt ; \
	../simulate -P=testLiAl6.txt ; \
	../simulate -t=1h -P=testLiAl6.txt ; \
	../simulate -t=1d -P=testLiAl6.txt 

testLiAl7:
	mkdir -p testLiAl7
	cd testLiAl7 ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl7.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl7.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, 7Li=1.0\";" >> testLiAl7.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testLiAl7.txt ; \
	../simulate -P=testLiAl7.txt ; \
	../simulate -t=1h -P=testLiAl7.txt ; \
	../simulate -t=1d -P=testLiAl7.txt 

testLiAl_cur100_10days:
	mkdir -p testLiAl_cur100_10days
	cd testLiAl_cur100_10days ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testLiAl_cur100_10days.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testLiAl_cur100_10days.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Li=1.0\";" >> testLiAl_cur100_10days.txt ; \
	echo "double e_emittedElectronMol = 1.0E-8;" >>testLiAl_cur100_10days.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-8;" >> testLiAl_cur100_10days.txt ; \
	../simulate -P=testLiAl_cur100_10days.txt ; \
	../simulate -t=1d -P=testLiAl_cur100_10days.txt ; \
	../simulate -t=1d -P=testLiAl_cur100_10days.txt  ; \
	../simulate -t=2d -P=testLiAl_cur100_10days.txt ; \
	../simulate -t=2d -P=testLiAl_cur100_10days.txt ; \
	../simulate -t=2d -P=testLiAl_cur100_10days.txt ; \
	../simulate -t=2d -P=testLiAl_cur100_10days.txt

testCa:
	mkdir -p testCa
	cd testCa ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramCa.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramCa.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ca=1.0\";" >> paramCa.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramCa.dat ; \
	../simulate -P=paramCa.dat ; \
	../simulate -t=1h -P=paramCa.dat ; \
	../simulate -t=1d -P=paramCa.dat 

testNa:
	mkdir -p testNa
	cd testNa ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramNa.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramNa.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Na=1.0\";" >> paramNa.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramNa.dat ; \
	../simulate -P=paramNa.dat ; \
	../simulate -t=1h -P=paramNa.dat ; \
	../simulate -t=1d -P=paramNa.dat 

testK:
	mkdir -p testK
	cd testK ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramK.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramK.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, K=1.0\";" >> paramK.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramK.dat ; \
	../simulate -P=paramK.dat ; \
	../simulate -t=1h -P=paramK.dat ; \
	../simulate -t=1d -P=paramK.dat 

testMg:
	mkdir -p testMg
	cd testMg ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramMg.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > paramMg.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Mg=1.0\";" >> paramMg.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramMg.dat ; \
	../simulate -P=paramMg.dat ; \
	../simulate -t=1h -P=paramMg.dat ; \
	../simulate -t=1d -P=paramMg.dat 

testTi:
	mkdir -p testTi
	cd testTi ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log paramTi.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ti=1.0\";" > paramTi.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ti=1.0\";" >> paramTi.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> paramTi.dat ; \
	../simulate -P=paramTi.dat ; \
	../simulate -t=1h -P=paramTi.dat ; \
	../simulate -t=1d -P=paramTi.dat 

#----------------------------------------------------------
testBe:
	mkdir -p testBe
	cd testBe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testBe.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Be=1.0\";" > testBe.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Be=1.0\";" >> testBe.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testBe.dat ; \
	../simulate -P=testBe.dat ; \
	../simulate -t=1h -P=testBe.dat ; \
	../simulate -t=1d -P=testBe.dat 

testSc:
	mkdir -p testSc
	cd testSc ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSc.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Sc=1.0\";" > testSc.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Sc=1.0\";" >> testSc.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testSc.dat ; \
	../simulate -P=testSc.dat ; \
	../simulate -t=1h -P=testSc.dat ; \
	../simulate -t=1d -P=testSc.dat 

testV:
	mkdir -p testV
	cd testV ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testV.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"V=1.0\";" > testV.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, V=1.0\";" >> testV.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testV.dat ; \
	../simulate -P=testV.dat ; \
	../simulate -t=1h -P=testV.dat ; \
	../simulate -t=1d -P=testV.dat 

testCr:
	mkdir -p testCr
	cd testCr ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testCr.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Cr=1.0\";" > testCr.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Cr=1.0\";" >> testCr.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testCr.dat ; \
	../simulate -P=testCr.dat ; \
	../simulate -t=1h -P=testCr.dat ; \
	../simulate -t=1d -P=testCr.dat 

testMn:
	mkdir -p testMn
	cd testMn ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testMn.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Mn=1.0\";" > testMn.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Mn=1.0\";" >> testMn.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testMn.dat ; \
	../simulate -P=testMn.dat ; \
	../simulate -t=1h -P=testMn.dat ; \
	../simulate -t=1d -P=testMn.dat 

testCo:
	mkdir -p testCo
	cd testCo ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testCo.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Co=1.0\";" > testCo.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Co=1.0\";" >> testCo.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testCo.dat ; \
	../simulate -P=testCo.dat ; \
	../simulate -t=1h -P=testCo.dat ; \
	../simulate -t=1d -P=testCo.dat 

testCu:
	mkdir -p testCu
	cd testCu ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testCu.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Cu=1.0\";" > testCu.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Cu=1.0\";" >> testCu.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testCu.dat ; \
	../simulate -P=testCu.dat ; \
	../simulate -t=1h -P=testCu.dat ; \
	../simulate -t=1d -P=testCu.dat 

testZn:
	mkdir -p testZn
	cd testZn ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testZn.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testZn.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Zn=1.0\";" >> testZn.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testZn.dat ; \
	../simulate -P=testZn.dat ; \
	../simulate -t=1h -P=testZn.dat ; \
	../simulate -t=1d -P=testZn.dat 

testGa:
	mkdir -p testGa
	cd testGa ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testGa.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testGa.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ga=1.0\";" >> testGa.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testGa.dat ; \
	../simulate -P=testGa.dat ; \
	../simulate -t=1h -P=testGa.dat ; \
	../simulate -t=1d -P=testGa.dat 

testGe:
	mkdir -p testGe
	cd testGe ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testGe.dat ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ge=1.0\";" > testGe.dat ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.2, Ge=1.0\";" >> testGe.dat ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testGe.dat ; \
	../simulate -P=testGe.dat ; \
	../simulate -t=1h -P=testGe.dat ; \
	../simulate -t=1d -P=testGe.dat 

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
	../simulate -P=testSHV_AlH.txt ; \
	../simulate -t=1h -P=testSHV_AlH.txt ; \
	../simulate -t=1d -P=testSHV_AlH.txt

testSHV_AlHD:
	mkdir -p testSHV_AlHD
	cd testSHV_AlHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_AlHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Al=1.0\";" > testSHV_AlHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.1, D=0.1, Al=1.0\";" >> testSHV_AlHD.txt ; \
	echo "double e_emittedProtonMol = 0.13333E-10;" >> testSHV_AlHD.txt ; \
	echo "double e_appliedVoltageScale = 5.97;" >> testSHV_AlHD.txt ; \
	../simulate -P=testSHV_AlHD.txt ; \
	../simulate -t=1h -P=testSHV_AlHD.txt ; \
	../simulate -t=1d -P=testSHV_AlHD.txt

testSHV_FeH:
	mkdir -p testSHV_FeH
	cd testSHV_FeH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_FeH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testSHV_FeH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, Fe=1.0\";" >> testSHV_FeH.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_FeH.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_FeH.txt ; \
	../simulate -P=testSHV_FeH.txt ; \
	../simulate -t=1h -P=testSHV_FeH.txt ; \
	../simulate -t=1d -P=testSHV_FeH.txt

testSHV_FeHD:
	mkdir -p testSHV_FeHD
	cd testSHV_FeHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_FeHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Fe=1.0\";" > testSHV_FeHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.2, D=0.2, Fe=1.0\";" >> testSHV_FeHD.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_FeHD.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_FeHD.txt ; \
	../simulate -P=testSHV_FeHD.txt ; \
	../simulate -t=1h -P=testSHV_FeHD.txt ; \
	../simulate -t=1d -P=testSHV_FeHD.txt

testSHV_56FeH:
	mkdir -p testSHV_56FeH
	cd testSHV_56FeH ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_56FeH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"56Fe=1.0\";" > testSHV_56FeH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.4, 56Fe=1.0\";" >> testSHV_56FeH.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_56FeH.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_56FeH.txt ; \
	../simulate -P=testSHV_56FeH.txt ; \
	../simulate -t=1h -P=testSHV_56FeH.txt ; \
	../simulate -t=1d -P=testSHV_56FeH.txt

testSHV_56FeHD:
	mkdir -p testSHV_56FeHD
	cd testSHV_56FeHD ; \
	rm -f data.dat s00*.log m*.log h*.log D*.log testSHV_56FeHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"56Fe=1.0\";" > testSHV_56FeHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.2, D=0.2, 56Fe=1.0\";" >> testSHV_56FeHD.txt ; \
	echo "double e_emittedProtonMol = 0.26666E-10;" >> testSHV_56FeHD.txt ; \
	echo "double e_appliedVoltageScale = 10.01;" >> testSHV_56FeHD.txt ; \
	../simulate -P=testSHV_56FeHD.txt ; \
	../simulate -t=1h -P=testSHV_56FeHD.txt ; \
	../simulate -t=1d -P=testSHV_56FeHD.txt

testSHV_NiH:
	mkdir -p testSHV_NiH
	cd testSHV_NiH ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_NiH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testSHV_NiH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, Ni=1.0\";" >> testSHV_NiH.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_NiH.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_NiH.txt ; \
	../simulate -P=testSHV_NiH.txt ; \
	../simulate -t=1h -P=testSHV_NiH.txt ; \
	../simulate -t=1d -P=testSHV_NiH.txt

testSHV_NiHD:
	mkdir -p testSHV_NiHD
	cd testSHV_NiHD ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_NiHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"Ni=1.0\";" > testSHV_NiHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, Ni=1.0\";" >> testSHV_NiHD.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_NiHD.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_NiHD.txt ; \
	../simulate -P=testSHV_NiHD.txt ; \
	../simulate -t=1h -P=testSHV_NiHD.txt ; \
	../simulate -t=1d -P=testSHV_NiHD.txt

testSHV_62NiH:
	mkdir -p testSHV_62NiH
	cd testSHV_62NiH ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_62NiH.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"62Ni=1.0\";" > testSHV_62NiH.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"H=0.6, 62Ni=1.0\";" >> testSHV_62NiH.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_62NiH.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_62NiH.txt ; \
	../simulate -P=testSHV_62NiH.txt ; \
	../simulate -t=1h -P=testSHV_62NiH.txt ; \
	../simulate -t=1d -P=testSHV_62NiH.txt

testSHV_62NiHD:
	mkdir -p testSHV_62NiHD
	cd testSHV_62NiHD ; \
	rm -f data.dat s00*.log m*.log h*.log testSHV_62NiHD.txt ; \
	echo "char e_negativeElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"62Ni=1.0\";" > testSHV_62NiHD.txt ; \
	echo "char e_positiveElectrodeAtomicMols[TEXT_LEN_OF_MOLS] = \"1H=0.3, D=0.3, 62Ni=1.0\";" >> testSHV_62NiHD.txt ; \
	echo "double e_emittedProtonMol = 0.4e-10;" >> testSHV_62NiHD.txt ; \
	echo "double e_appliedVoltageScale = 10.58;" >> testSHV_62NiHD.txt ; \
	../simulate -P=testSHV_62NiHD.txt ; \
	../simulate -t=1h -P=testSHV_62NiHD.txt ; \
	../simulate -t=1d -P=testSHV_62NiHD.txt

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