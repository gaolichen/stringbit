# Add executable called "helloDemo" that is built from the source files

SET(STRBIT_SRCS 
	BitUtility.cpp 
	Coefficient.cpp
	Hamiltonian.cpp
	HamOperator.cpp
	main.cpp
	MixedTrace.cpp
	MixState.cpp
	NormCalculator.cpp
	Polynomial.cpp
	SingleTrace.cpp
	StateCollection.cpp
	StateGenerator.cpp
	StateId.cpp
	TraceState.cpp
	ScriptGenerator.cpp)

IF(WIN32)
	SET(STRBIT_SRCS 
		${STRBIT_SRCS}
		BitUtility.h 
		Coefficient.h 
		Hamiltonian.h 
		HamOperator.h 
		MixedTrace.h 
		MixState.h 
		NormCalculator.h 
		Polynomial.h 
		SingleTrace.h 
		StateCollection.h 
		StateGenerator.h 
		StateId.h 
		TraceState.h 
		ScriptGenerator.h
	)
	
	source_group("Header Files" FILES *.h)
ENDIF(WIN32)



add_executable (strbit ${STRBIT_SRCS})


# Link the executable to the Hello library. Since the Hello library has

# public include directories we will use those link directories when building

# helloDemo

# target_link_libraries (helloDemo LINK_PUBLIC Hello)

