CXXFLAGS = -O3 -fopenmp


bench_Prob_Unitary_Synthesis : bench_Prob_Unitary_Synthesis.o eps_net_verification.o Prob_Synthesis.o ExactSynthesis.o enum_u_t.o grid_solver.o U2_ZOmega.o quaternion.o rings.o
	g++ $^ $(CXXFLAGS) -ltbb -o benchmark/Prob_Unitary_Synthesis

bench_Prob_Unitary_Synthesis.o : benchmark/Prob_Unitary_Synthesis.cpp src/Prob_Synthesis.cpp
	g++ -c $< $(CXXFLAGS) -o $@

test_Prob_Synthesis : test_Prob_Synthesis.o eps_net_verification.o Prob_Synthesis.o ExactSynthesis.o enum_u_t.o grid_solver.o U2_ZOmega.o quaternion.o rings.o
	g++ $^ $(CXXFLAGS) -ltbb -o tests/Prob_Synthesis

test_Prob_Synthesis.o : tests/Prob_Synthesis.cpp src/Prob_Synthesis.cpp
	g++ -c $< $(CXXFLAGS) -o $@

# target : test_enum_u_t.o ExactSynthesis.o enum_u_t.o grid_solver.o U2_ZOmega.o quaternion.o rings.o
# 	g++ $^ $(CXXFLAGS) -ltbb -o tests/enum_u_t

# test_enum_u_t.o : tests/enum_u_t.cpp src/enum_u_t.cpp
# 	g++ -c $< $(CXXFLAGS) -o $@

Prob_Synthesis.o : src/Prob_Synthesis.cpp
	g++ -c $< $(CXXFLAGS)

eps_net_verification.o : src/eps_net_verification.cpp
	g++ -c $< $(CXXFLAGS)

ExactSynthesis.o : src/ExactSynthesis.cpp
	g++ -c $< $(CXXFLAGS)

enum_u_t.o : src/enum_u_t.cpp
	g++ -c $< $(CXXFLAGS)

grid_solver.o : src/grid_solver.cpp
	g++ -c $< $(CXXFLAGS)

U2_ZOmega.o : src/U2_ZOmega.cpp
	g++ -c $< $(CXXFLAGS)

rings.o : src/rings.cpp
	g++ -c $< $(CXXFLAGS)

quaternion.o : src/quaternion.cpp
	g++ -c $< $(CXXFLAGS)


clean : 
	rm quaternion.o rings.o U2_ZOmega.o grid_solver.o \
	   enum_u_t.o ExactSynthesis.o eps_net_verification.o \
	   Prob_Synthesis.o