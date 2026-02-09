if [ ! -d "./bin" ]
then
  mkdir ./bin
fi

if [ ! -f "./bin/MBE_ALL" ]
then
  cd ./code || exit
  cd MBE || exit
  mkdir build
  cd build || exit
  cmake -G "MinGW Makefiles" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
  mingw32-make VERBOSE=1
  mv MBE_ALL ../../../bin/
  cd ../../../
fi
if [ ! -f "./bin/mbbp" ]
then
  cd ./code || exit
  cd cohesive_subgraph_bipartite || exit
  mkdir build
  cd build || exit
  cmake -G "MinGW Makefiles" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
  mingw32-make VERBOSE=1
  mv mbbp ../../../bin/
  cd ../../../
fi

