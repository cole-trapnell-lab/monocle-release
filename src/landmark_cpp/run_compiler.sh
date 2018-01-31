# /opt/local/bin/g++-mp-4.8 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -o landmark_selector.o landmark_selector.cpp
# /opt/local/bin/g++-mp-4.8 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -o random.o random.cpp
# /opt/local/bin/g++-mp-4.8 -std=c++11 -fopenmp  -o landmark_selector  landmark_selector.o random.o

/opt/local/bin/mpicc-mpich-gcc7 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -o landmark_selector.o landmark_selector.cpp
/opt/local/bin/mpicc-mpich-gcc7 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -o random.o random.cpp
/opt/local/bin/mpicc-mpich-gcc7 -std=c++11 -fopenmp  -o landmark_selector  landmark_selector.o random.o

# /opt/local/bin/g++-mp-4.8 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"landmark_selector.d" -MT"landmark_selector.d" -o "landmark_selector.o" "landmark_selector.cpp"

# /opt/local/bin/g++-mp-4.8 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"random.d" -MT"random.d" -o "random.o" "random.cpp"

# /opt/local/bin/g++-mp-4.8 -std=c++11 -fopenmp  -o "landmark_selector"  ./landmark_selector.o ./random.o  

# matlab -nodesktop
# mex libsvmread.c
# exit

