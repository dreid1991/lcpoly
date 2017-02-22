#include <vector>
#include <Eigen/Dense>
template <class T>
class Grid {
public:
        std::vector<T> vals;
        int ns[3];
        double os[3];
        //doub
        Grid(){}
        Grid(int nx, int ny, int nz) {
                vals = std::vector<T>(nx*ny*nz, 0);
                ns[0] = nx;
                ns[1] = ny;
                ns[2] = nz;
        }

        Grid(int nx, int ny, int nz, int vlen) {
          vals = std::vector<T>(nx*ny*nz, std::vector<int>(vlen));
        }
        void reset_grid(){
          fill(vals.begin(),vals.end(),T());
        }
        void reset_gridlist(){
          vals.clear();
        }
        T &operator()(int x, int y, int z) {
          return vals[x + ns[1]*y + ns[0]*ns[1]*z];
        }
        T &operator()(int x, int y, int z, int el) {
          return vals[x + ns[1]*y + ns[0]*ns[1]*z].push_back(el);
        }
};
