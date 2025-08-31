#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

static constexpr int    N   = 300;
static constexpr double rho = 0.25;
static constexpr double L   = std::sqrt(N / rho);
static constexpr double L2  = L / 2.0;

static constexpr std::size_t T    = 1'000'000; 
static constexpr double CAGE      = 2.0;       

static constexpr std::size_t DT_STEP        = 10;   
static constexpr std::size_t ORIGINS_TOTAL  = 100; 


// static const std::vector<std::string> TEMPS = {
//     "0.47","0.45","0.43","0.41","0.39","0.37","0.35","0.33",
//     "0.31","0.29","0.27","0.25","0.23","0.21","0.19","0.17"
// };
static const std::vector<std::string> TEMPS = {
    "0.47", "0.43", "0.41", "0.39", "0.37", "0.35", "0.33", "0.31", "0.29", "0.27", "0.25", "0.23", "0.21", "0.19", "0.17"
};

struct MMap {
    std::size_t bytes = 0;
    double* ptr = nullptr;
    int fd = -1;

    void open_map(const std::string& path) {
        fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) { perror("open"); throw std::runtime_error("open failed: " + path); }
        struct stat st{};
        if (fstat(fd, &st) != 0) { perror("fstat"); throw std::runtime_error("fstat failed: " + path); }
        bytes = static_cast<std::size_t>(st.st_size);
        void* p = ::mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
        if (p == MAP_FAILED) { perror("mmap"); throw std::runtime_error("mmap failed: " + path); }
        ptr = static_cast<double*>(p);
    }
    void close_map() {
        if (ptr) ::munmap(ptr, bytes);
        if (fd >= 0) ::close(fd);
        ptr = nullptr; fd = -1; bytes = 0;
    }
    ~MMap() { close_map(); }
};

inline double pbc(double dr) {
    if (dr >  L2) dr -= L;
    if (dr < -L2) dr += L;
    return dr;
}

inline const double& at(const double* A, std::size_t t, int i) {
    return A[t * static_cast<std::size_t>(N) + i];
}

void unwrap_positions(const double* mod, float* out) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        double acc = mod[0 * N + i];
        out[0 * N + i] = static_cast<float>(acc);
        for (std::size_t t = 1; t < T; ++t) {
            double prev = mod[(t-1) * N + i];
            double curr = mod[t * N + i];
            double step = curr - prev;
            if (step >  L2) step -= L;
            if (step < -L2) step += L;
            acc += step;
            out[t * N + i] = static_cast<float>(acc);
        }
    }
}

static inline void build_neighbors_at_t0(
    const double* XX, const double* YY, std::size_t t0,
    std::vector<int>& row_ptr, std::vector<int>& cols, std::vector<int>& deg)
{
    row_ptr.assign(N + 1, 0);
    cols.clear();
    deg.assign(N, 0);

    for (int i = 0; i < N; ++i) {
        int cnt = 0;
        const double xi = at(XX, t0, i);
        const double yi = at(YY, t0, i);
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            const double dx = pbc(xi - at(XX, t0, j));
            const double dy = pbc(yi - at(YY, t0, j));
            const double d2 = dx*dx + dy*dy;
            if (d2 <= CAGE*CAGE) ++cnt;
        }
        deg[i] = cnt;
        row_ptr[i+1] = row_ptr[i] + cnt;
    }
    cols.resize(static_cast<std::size_t>(row_ptr.back()));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        int write = row_ptr[i];
        const double xi = at(XX, t0, i);
        const double yi = at(YY, t0, i);
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            const double dx = pbc(xi - at(XX, t0, j));
            const double dy = pbc(yi - at(YY, t0, j));
            const double d2 = dx*dx + dy*dy;
            if (d2 <= CAGE*CAGE) cols[write++] = j;
        }
    }
}

static inline void make_evenly_spaced_origins(std::size_t n_possible, std::vector<std::size_t>& origins) {
    origins.clear();
    if (n_possible == 0) return;
    const std::size_t m = std::min<std::size_t>(ORIGINS_TOTAL, n_possible);
    origins.reserve(m);
    if (m == 1) {
        origins.push_back(0);
        return;
    }
    const std::size_t last = n_possible - 1;
    for (std::size_t k = 0; k < m; ++k) {
        std::size_t t0 = (last * k) / (m - 1);
        origins.push_back(t0);
    }
}

std::vector<double> msd_cage_relative(
    const float* X_unw, const float* Y_unw,
    const double* XX_mod, const double* YY_mod)
{
    const std::size_t NDTS = (T - 1) / DT_STEP;
    std::vector<double> msd(NDTS, std::numeric_limits<double>::quiet_NaN());

    for (std::size_t k = 0; k < NDTS; ++k) {
        const std::size_t dt = (k + 1) * DT_STEP;
        const std::size_t n_possible = T - dt;
        if (n_possible == 0) { msd[k] = std::numeric_limits<double>::quiet_NaN(); continue; }

        std::vector<std::size_t> origins;
        make_evenly_spaced_origins(n_possible, origins);
        const std::size_t total_origins = origins.size();

        double accum = 0.0;
        long long count = 0;
        std::size_t origins_done = 0;

        #pragma omp parallel
        {
            std::vector<int> row_ptr_local, cols_local, deg_local;
            std::vector<double> nx(N), ny(N);

            double accum_local = 0.0;
            long long count_local = 0;

            #pragma omp for schedule(dynamic)
            for (std::size_t oi = 0; oi < total_origins; ++oi) {
                const std::size_t t0 = origins[oi];

                // neighbors at this origin
                build_neighbors_at_t0(XX_mod, YY_mod, t0,
                                      row_ptr_local, cols_local, deg_local);

                // neighbor-mean displacement
                for (int i = 0; i < N; ++i) {
                    const int d = deg_local[i];
                    if (d == 0) { nx[i] = 0.0; ny[i] = 0.0; continue; }
                    double sx = 0.0, sy = 0.0;
                    const int s = row_ptr_local[i];
                    const int e = row_ptr_local[i+1];
                    for (int kk = s; kk < e; ++kk) {
                        const int j = cols_local[kk];
                        const double dx = static_cast<double>(X_unw[(t0 + dt) * N + j])
                                        - static_cast<double>(X_unw[t0 * N + j]);
                        const double dy = static_cast<double>(Y_unw[(t0 + dt) * N + j])
                                        - static_cast<double>(Y_unw[t0 * N + j]);
                        sx += dx; sy += dy;
                    }
                    const double invd = 1.0 / static_cast<double>(d);
                    nx[i] = sx * invd;
                    ny[i] = sy * invd;
                }
                for (int i = 0; i < N; ++i) {
                    if (deg_local[i] == 0) continue;
                    const double dxi = static_cast<double>(X_unw[(t0 + dt) * N + i])
                                     - static_cast<double>(X_unw[t0 * N + i]);
                    const double dyi = static_cast<double>(Y_unw[(t0 + dt) * N + i])
                                     - static_cast<double>(Y_unw[t0 * N + i]);
                    const double dxcr = dxi - nx[i];
                    const double dycr = dyi - ny[i];
                    accum_local += dxcr*dxcr + dycr*dycr;
                    ++count_local;
                }

                std::size_t done;
                #pragma omp atomic capture
                { ++origins_done; done = origins_done; }
            }

            #pragma omp atomic
            accum += accum_local;
            #pragma omp atomic
            count += count_local;
        } 
        
    } 

    return msd;
}

int main() {
    try {
        const std::size_t expected_bytes = T * static_cast<std::size_t>(N) * sizeof(double);
        const std::size_t NDTS = (T - 1) / DT_STEP;

        for (const auto& temp : TEMPS) {
            std::cerr << "[*] Processing T=" << temp << "\n";

            const std::string xbin = "/Users/leo/C++/MolecularDynamics/300bin/300_"+temp+"_X.bin";
            const std::string ybin = "/Users/leo/C++/MolecularDynamics/300bin/300_"+temp+"_Y.bin";

            MMap Xmod, Ymod;
            Xmod.open_map(xbin);
            Ymod.open_map(ybin);
            if (Xmod.bytes != expected_bytes || Ymod.bytes != expected_bytes) {
                throw std::runtime_error("Input size mismatch for " + temp +
                    " (expected " + std::to_string(expected_bytes) + " bytes per file).");
            }

            std::vector<float> X_unw(T * static_cast<std::size_t>(N));
            std::vector<float> Y_unw(T * static_cast<std::size_t>(N));
            std::cerr << "    Unwrapping X...\n";
            unwrap_positions(Xmod.ptr, X_unw.data());
            std::cerr << "    Unwrapping Y...\n";
            unwrap_positions(Ymod.ptr, Y_unw.data());

            std::cerr << "    Computing MSD_cr (DT_STEP=" << DT_STEP
                      << ", ORIGINS_TOTAL=" << ORIGINS_TOTAL
                      << ", NDTS=" << NDTS << ")...\n";
            std::vector<double> msd = msd_cage_relative(X_unw.data(), Y_unw.data(),
                                                         Xmod.ptr, Ymod.ptr);

            const std::string outbin = "/Users/leo/C++/MolecularDynamics/MSD_CR/MSDfull_"+temp+".bin";
            std::ofstream fout(outbin, std::ios::binary);
            if (!fout) throw std::runtime_error("Cannot open output: " + outbin);
            fout.write(reinterpret_cast<const char*>(msd.data()),
                       static_cast<std::streamsize>(msd.size() * sizeof(double)));
            fout.close();
            std::cerr << "    Saved: " << outbin << " (" << msd.size() << " doubles)\n";
        }

        std::cerr << "[*] All temperatures done.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
