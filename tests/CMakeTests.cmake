add_executable (OpenSMOKE_GasKinetics tests/OpenSMOKE_GasKinetics.cpp)
target_link_libraries (OpenSMOKE_GasKinetics PRIVATE opensmoke)

add_executable (OpenSMOKE_ODESolver tests/OpenSMOKE_ODESolver.cpp)
target_link_libraries (OpenSMOKE_ODESolver PRIVATE opensmoke)

add_test (NAME OpenSMOKE_GasKinetics
          COMMAND bash -c "
            ./OpenSMOKE_GasKinetics 2> OpenSMOKE_GasKinetics.log
            diff OpenSMOKE_GasKinetics.log ../tests/OpenSMOKE_GasKinetics.ref")

add_test (NAME OpenSMOKE_ODESolver
          COMMAND bash -c "
            ./OpenSMOKE_ODESolver 2> OpenSMOKE_ODESolver.log
            diff OpenSMOKE_ODESolver.log ../tests/OpenSMOKE_ODESolver.ref")

add_test (NAME OpenSMOKE_PythonInterface
          COMMAND zsh -c "python3 ../tests/OpenSMOKE_PythonInterface.py")

