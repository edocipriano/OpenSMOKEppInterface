
g++ -w -fPIC -std=c++11 -c *.cpp \
    -I$EIGEN_LIBRARY_PATH \
    -I/usr/local/Cellar/boost/1.80.0/include \
    -I$OPENSMOKE_LIBRARY_PATH \
    -I/usr/local/Cellar/libconfig/1.7.3/include \
    -L/usr/local/Cellar/boost/1.80.0/lib \
    -L/usr/local/Cellar/libconfig/1.7.3/lib \
    -lconfig++ \
    -lboost_date_time \
    -lboost_filesystem \
    -lboost_program_options \
    -lboost_system \
    -lboost_regex \
    -lboost_timer \
    -lboost_chrono \
    -lboost_program_options \
    -lboost_container \
    -lboost_exception


g++ -std=c++11 -shared *.o -o libopensmoke.so \
    -I$EIGEN_LIBRARY_PATH \
    -I/usr/local/Cellar/boost/1.80.0/include \
    -I$OPENSMOKE_LIBRARY_PATH \
    -I/usr/local/Cellar/libconfig/1.7.3/include \
    -L/usr/local/Cellar/boost/1.80.0/lib \
    -L/usr/local/Cellar/libconfig/1.7.3/lib \
    -lconfig++ \
    -lboost_date_time \
    -lboost_filesystem \
    -lboost_program_options \
    -lboost_system \
    -lboost_regex \
    -lboost_timer \
    -lboost_chrono \
    -lboost_program_options \
    -lboost_container \
    -lboost_exception

rm *.o
