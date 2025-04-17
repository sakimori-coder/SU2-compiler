CXXFLAGS = -O3 -fopenmp -I /home/morisaki/SU-2-compiler -lquadmath


# SU2_compiler_MPIのテストのbuild
test_SU2_compiler_MPI : test_SU2_compiler_MPI.o rings.o U2_ZOmega.o quaternion.o grid_solver.o ExactSynthesis.o
	/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $^ $(CXXFLAGS) -o tests/SU2_compiler_MPI

test_SU2_compiler_MPI.o : tests/SU2_compiler_MPI.cpp src/SU2_compiler_MPI.cpp
	/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx -c $< $(CXXFLAGS) -o $@


# ユニモジュラー行列のテストのbuild
test_Unimodular : test_Unimodular.o quaternion.o
	g++ $^ $(CXXFLAGS) -o tests/Unimodular

test_Unimodular.o : tests/Unimodular_Z.cpp src/type.hpp src/quaternion.hpp src/Unimodular_Z.hpp
	g++ -c $< $(CXXFLAGS) -o $@

probabilistic_bench : benchmark/probabilistic.cpp eps_net_verification.o Prob_Synthesis.o ExactSynthesis.o grid_solver.o U2_ZOmega.o quaternion.o rings.o
	g++ $^ $(CXXFLAGS) -o benchmark/probabilistic


# benchのテストのbuild
deterministic_bench : benchmark/deterministic.cpp rings.o U2_ZOmega.o quaternion.o grid_solver.o ExactSynthesis.o
	g++ $^ $(CXXFLAGS) -o benchmark/deterministic



test_Prob_Synthesis : test_Prob_Synthesis.o eps_net_verification.o Prob_Synthesis.o ExactSynthesis.o grid_solver.o U2_ZOmega.o quaternion.o rings.o
	g++ $^ $(CXXFLAGS) -o tests/Prob_Synthesis

test_Prob_Synthesis.o : tests/Prob_Synthesis.cpp src/Prob_Synthesis.cpp src/SU2_compiler.cpp
	g++ -c $< $(CXXFLAGS) -o $@

Prob_Synthesis.o : src/Prob_Synthesis.cpp src/Prob_Synthesis.hpp src/type.hpp src/rings.hpp src/quaternion.hpp src/U2_ZOmega.hpp src/ExactSynthesis.hpp src/eps_net_verification.hpp src/sdpa_dd.hpp src/SU2_compiler.hpp
	g++ $(CXXFLAGS) -c $< -o $@

eps_net_verification.o : src/eps_net_verification.cpp src/eps_net_verification.hpp src/type.hpp src/quaternion.hpp src/linalg.hpp
	g++ $(CXXFLAGS) -c $< -o $@






# SU2_compilerのテストのbuild
test_SU2_compiler : test_SU2_compiler.o rings.o U2_ZOmega.o quaternion.o ExactSynthesis.o SU2_compiler.o ExactSynthesis.o
	g++ $^ $(CXXFLAGS) -o tests/SU2_compiler

test_SU2_compiler.o : tests/SU2_compiler.cpp src/SU2_compiler.hpp src/type.hpp src/quaternion.hpp src/ExactSynthesis.hpp
	g++ -c $< $(CXXFLAGS) -o $@

SU2_compiler.o : src/SU2_compiler.cpp src/type.hpp src/quaternion.hpp src/rings.hpp src/U2_ZOmega.hpp src/Unimodular_Z.hpp src/ExactSynthesis.hpp
	g++ -c $< $(CXXFLAGS) -o $@

ExactSynthesis.o : src/ExactSynthesis.cpp src/ExactSynthesis.hpp src/U2_ZOmega.hpp
	g++ $(CXXFLAGS) -c $< -o $@

grid_solver.o : src/grid_solver.cpp src/grid_solver.hpp src/type.hpp src/rings.hpp
	g++ $(CXXFLAGS) -c $< -o $@

U2_ZOmega.o : src/U2_ZOmega.cpp src/U2_ZOmega.hpp src/type.hpp src/rings.hpp src/quaternion.hpp
	g++ $(CXXFLAGS) -c $< -o $@

rings.o : src/rings.cpp src/rings.hpp src/type.hpp
	g++ $(CXXFLAGS) -c $< -o $@

quaternion.o : src/quaternion.cpp src/quaternion.hpp src/type.hpp
	g++ $(CXXFLAGS) -c $< -o $@


clean : 
	rm quaternion.o rings.o U2_ZOmega.o grid_solver.o \
	   enum_u_t.o ExactSynthesis.o eps_net_verification.o \
	   Prob_Synthesis.o \
	   SU2_compiler.o