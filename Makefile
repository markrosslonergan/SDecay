d1:	inflight.cxx sterile_flux.h sterile_flux.cxx fourmomentum.h fourmomentum.cxx
	g++ -g -std=c++11 -c inflight.cxx -o inflight.o -I.
	g++ -g -std=c++11 -c fourmomentum.cxx -o fourmomentum.o -I.
	g++ -g -std=c++11 -c sterile_flux.cxx -o sterile_flux.o -I.
	g++ -g -std=c++11 -c resonantZprime_decay.cxx -o decay_function1.o -I.
	g++ -g -std=c++11 -c threebody_decay.cxx -o decay_function2.o -I.
	g++ -g -std=c++11 -o inflight inflight.o fourmomentum.o sterile_flux.o decay_function1.o decay_function2.o -lgomp -lnlopt -lgsl -lgslcblas
	rm *.o




#d1:	inflight.cxx sterile_flux.h sterile_flux.cxx fourmomentum.h fourmomentum.cxx
#	g++ -g -std=c++11 -c inflight.cxx -o inflight.o -I.
#	g++ -g -std=c++11 -c fourmomentum.cxx -o fourmomentum.o -I.
#	g++ -g -std=c++11 -c sterile_flux.cxx -o sterile_flux.o -I.
#	g++ -g -std=c++11 -c resonantZprime_decay.cxx -o decay_function.o -I.
#	g++ -g -std=c++11 -o inflight inflight.o fourmomentum.o sterile_flux.o decay_function.o -lgomp -lnlopt -lgsl -lgslcblas
#	rm *.o
