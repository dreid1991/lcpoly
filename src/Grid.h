#include <vector>
#include <Eigen/Dense>
#include <iostream>
template <class T>
class Grid {
public:
        std::vector<T> vals;
        Vector3i ns;
        Vector3d os;
        Vector3d ds;

        Grid(){}
        Grid(Vector3i ns_) : ns(ns_) {
            allocate();
        }
        Grid(Vector3d size, double gridSize) {
            for (int i=0; i<3; i++) {
                ns[i] = size[i] / gridSize;
                ds[i] = size[i] / ns[i];
            }

            printf("sizes %f %f %f %f \n", size[0], size[1], size[2], gridSize);
            os = {0, 0, 0};
            allocate();

        }
        Vector3i coord(Vector3d pos) {
            Vector3i x;
            auto normed = pos - os;

            for (int i=0; i<3; i++) {
                x[i] = normed[i] / ds[i];
            }
            return x;
        }
        void allocate() {
            vals = std::vector<T>(ns[0]*ns[1]*ns[2], T());
        }
        void reset_grid(T val = T()){
            fill(vals.begin(),vals.end(),val);
        }
        void reset_gridlist(){
            vals.clear();
        }
        T &operator()(Vector3i idxs) {
            return vals[idxs[0]*ns[1]*ns[2] + idxs[1]*ns[2] + idxs[2]];
        }
        T &operator()(Vector3d pos) {
            return (*this)(coord(pos));
        }
        double volCell() {
            return ds[0]*ds[1]*ds[2];
        }
};
