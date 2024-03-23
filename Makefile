target : test_enum_u_t.o ExactSynthesis.o enum_u_t.o grid_solver.o U2_ZOmega.o quaternion.o rings.o
	g++  test_enum_u_t.o ExactSynthesis.o enum_u_t.o grid_solver.o U2_ZOmega.o quaternion.o rings.o -ltbb -o exfile

test_enum_u_t.o : tests/enum_u_t.cpp src/enum_u_t.cpp
	g++ -c tests/enum_u_t.cpp -o test_enum_u_t.o

ExactSynthesis.o : src/ExactSynthesis.cpp
	g++ -c src/ExactSynthesis.cpp

enum_u_t.o : src/enum_u_t.cpp
	g++ -c src/enum_u_t.cpp

grid_solver.o : src/grid_solver.cpp
	g++ -c src/grid_solver.cpp

U2_ZOmega.o : src/U2_ZOmega.cpp
	g++ -c src/U2_ZOmega.cpp

rings.o : src/rings.cpp
	g++ -c src/rings.cpp

quaternion.o : src/quaternion.cpp
	g++ -c src/quaternion.cpp
