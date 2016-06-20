evo : cmaes.c ThermoExample.c All_Inputs.h ThermoModel.o
	g++ -g cmaes.c ThermoExample.c ThermoModel.o -o evo 

ThermoModel.o : ThermoModel.cpp ThermoModel.h 
	g++ -c -g ThermoModel.cpp -o ThermoModel.o

read: readAllInputs.cpp readAllInputs.h
	g++ -o read readAllInputs.cpp

check : check.cpp ThermoModel.o
	g++ check.cpp ThermoModel.o -o check 

predict : getResults_and_Predict.cpp All_Inputs.h ThermoModel.o 
	g++ getResults_and_Predict.cpp ThermoModel.o -o predict 

clean:
	rm *.o
	rm read
	rm All_Inputs.h
