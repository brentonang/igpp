all:
	g++ -g -std=c++11 CompareTopo_AlongPath.cpp cget_latlon.cpp Cpath_topoSRTM1.cpp csetminmax.cpp -o igppout

cget:
	g++ -g -std=c++11 cget_latlon.cpp csetminmax.cpp -o cgetout

cpath:
	g++ -g -std=c++11 Cpath_topoSRTM1.cpp cget_latlon.cpp csetminmax.cpp -o cpathout

tester:
	g++ -g -std=c++11 test.cpp -o testoutput