mkdir build
cd build
cmake .. && make
./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt ../output/sample-laser-radar-measurement-data-1.txt
./UnscentedKF ../data/sample-laser-radar-measurement-data-2.txt ../output/sample-laser-radar-measurement-data-2.txt
./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../output/obj_pose-laser-radar-synthetic-input.txt
