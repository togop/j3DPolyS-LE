#include "mc_cuda.h"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>

#define EVMAX 65536
#define DIRPAD 14

#define CUDA_TRY(call)                                                         \
    do {                                                                       \
        cudaError_t _e = (call);                                               \
        if (_e != cudaSuccess) {                                               \
            std::fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, \
                         cudaGetErrorString(_e));                              \
            return -1;                                                         \
        }                                                                      \
    } while (0)

struct McShared {
    int *opp;
    float *voisxyz;
    float *costhet;
    int *voisnn;
    int *connec;
    int *nbr;
    float *boundary;
    float *loadfac;
    int *istate;
    int N;
    int naddr;
    int iku;
    int ikm;
    int ikb;
    float kb;
    float ku;
    float km;
    float Ea;
    float Ei;
    int z_loop;
    int unidirectional;
};

struct McCudaState {
    McCudaParams p;
    int Npad;
    int Apad;
    McShared sh;
    int *d_config;
    int *d_contact;
    int *d_occ;
    int *d_nint;
    int *d_nfree;
    float *d_dr;
    unsigned long long *d_rng;
    int *d_nsteps_t;
    int use_nsteps_t;
    int *d_ev_n;
    int *d_ev_id;
    int *d_ev_cnt;
};

/* ---- host table copies: Fortran column-major -> 1-based C pads ---- */

static void copy_opp(int *dst, const int *src)
{
    dst[0] = 0;
    for (int i = 1; i <= 13; i++)
        dst[i] = src[i - 1];
}

static void copy_voisxyz(float *dst, const float *src)
{
    std::memset(dst, 0, (size_t)DIRPAD * 4 * sizeof(float));
    for (int v = 1; v <= 13; v++)
        for (int c = 1; c <= 3; c++)
            dst[v * 4 + c] = src[(v - 1) * 3 + (c - 1)];
}

static void copy_costhet(float *dst, const float *src)
{
    std::memset(dst, 0, (size_t)DIRPAD * DIRPAD * sizeof(float));
    for (int j = 1; j <= 13; j++)
        for (int i = 1; i <= 13; i++)
            dst[i * DIRPAD + j] = src[(j - 1) * 13 + (i - 1)];
}

static void copy_table3(int *dst, const int *src)
{
    std::memset(dst, 0, (size_t)DIRPAD * DIRPAD * DIRPAD * sizeof(int));
    for (int c = 1; c <= 13; c++)
        for (int b = 1; b <= 13; b++)
            for (int a = 1; a <= 13; a++)
                dst[(c * DIRPAD + b) * DIRPAD + a] =
                    src[((c - 1) * 13 + (b - 1)) * 13 + (a - 1)];
}

static void copy_nbr(int *dst, const int *bittable, int naddr, int Apad)
{
    std::memset(dst, 0, (size_t)Apad * DIRPAD * sizeof(int));
    for (int a = 1; a <= naddr; a++)
        for (int iv = 2; iv <= 13; iv++)
            dst[a * DIRPAD + iv] = bittable[(a - 1) * 14 + (iv - 1)];
}

static uint64_t splitmix64_host(uint64_t *x)
{
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* ---- device RNG ---- */

__device__ __forceinline__ uint64_t rotl64(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

__device__ __forceinline__ double rnd(unsigned long long *rng)
{
    uint64_t s0 = (uint64_t)rng[0];
    uint64_t s1 = (uint64_t)rng[1];
    const uint64_t result = s0 + s1;
    s1 ^= s0;
    rng[0] = (unsigned long long)(rotl64(s0, 24) ^ s1 ^ (s1 << 16));
    rng[1] = (unsigned long long)rotl64(s1, 37);
    return (double)(result >> 11) * (1.0 / 9007199254740992.0);
}

__device__ __forceinline__ double hash_u01(unsigned long long x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return (double)(x >> 11) * (1.0 / 9007199254740992.0);
}

/* 1-based accessors (Npad / Apad already applied on the pointers) */

#define CFG1(n) (config[(n) * 2])
#define CFG2(n) (config[(n) * 2 + 1])
#define CON1(n) (contact[(n) * 3])
#define CON2(n) (contact[(n) * 3 + 1])
#define CON3(n) (contact[(n) * 3 + 2])
#define DRX(n) (dr[(n) * 3])
#define DRY(n) (dr[(n) * 3 + 1])
#define DRZ(n) (dr[(n) * 3 + 2])
#define OPP(i) (S.opp[(i)])
#define CTH(i, j) (S.costhet[(i) * DIRPAD + (j)])
#define VNN(a, b, c) (S.voisnn[((c) * DIRPAD + (b)) * DIRPAD + (a)])
#define CNC(a, b, c) (S.connec[((c) * DIRPAD + (b)) * DIRPAD + (a)])
#define VXYZ(c, v) (S.voisxyz[(v) * 4 + (c)])
#define NBR(a, iv) (S.nbr[(a) * DIRPAD + (iv)])

__device__ __forceinline__ int node_dir(const McShared &S, int en2, int iv)
{
    return (iv == 1) ? en2 : NBR(en2, iv);
}

__device__ void nint_move(const McShared &S, int *nintf, int en, int v)
{
    nintf[en] -= 1;
    nintf[v] += 1;
    for (int j = 2; j <= 13; j++) {
        nintf[NBR(en, j)] -= 1;
        nintf[NBR(v, j)] += 1;
    }
}

__device__ void record_unbind(int *ev_n, int *ev_id, int *ev_cnt, int evmax, int n,
                              int id)
{
    int c = ev_cnt[0];
    if (c < evmax) {
        ev_n[c] = n;
        ev_id[c] = id;
        ev_cnt[0] = c + 1;
    }
}

__device__ void trialmovetad(const McShared &S, int *config, int *contact, float *dr,
                             int *occ, int *nintf, unsigned long long *rng)
{
    const int N = S.N;
    const int n = 1 + (int)(N * rnd(rng));
    const int en = CFG1(n);

    if (n == 1) {
        const int en2 = CFG1(2);
        int cn2 = OPP(CFG2(1));
        int cm2;
        if (cn2 < CFG2(2)) {
            cm2 = CFG2(2);
        } else {
            cm2 = cn2;
            cn2 = CFG2(2);
        }
        int iv = 1 + (int)(11.0 * rnd(rng));
        if (iv >= cn2)
            iv += 1;
        if (iv >= cm2)
            iv += 1;

        const int v = node_dir(S, en2, iv);
        const int b = occ[v];
        if ((b == 0) || ((b == 1) && (en2 == v))) {
            const int id = CON1(n);
            cn2 = CFG2(1);
            const int cn3 = CON2(n);
            int cc = 0;
            if (cn3 != 0)
                cc = CNC(iv, OPP(cn2), cn3);
            if ((id != 0) && (cc == 0))
                return;
            float dE = CTH(OPP(iv), CFG2(2)) - CTH(cn2, CFG2(2));
            if (S.istate[n] > 0)
                dE += S.Ei * (float)(nintf[v] - nintf[en]);
            if ((CON3(n) == -1) && (id < N)) {
                if (CNC(OPP(CON2(n)), CFG2(id), 1) != 0)
                    dE -= S.Ea;
                if (CNC(OPP(cc), CFG2(id), 1) != 0)
                    dE += S.Ea;
            }
            if (CON3(n + 1) == -1) {
                if (CNC(OPP(cn2), CON2(n + 1), 1) != 0)
                    dE -= S.Ea;
                if (CNC(iv, CON2(n + 1), 1) != 0)
                    dE += S.Ea;
            }
            if (rnd(rng) < expf(-dE)) {
                occ[en] -= 1;
                occ[v] += 1;
                nint_move(S, nintf, en, v);
                CFG1(1) = v;
                CFG2(1) = OPP(iv);
                if (id != 0) {
                    CON2(n) = cc;
                    CON2(id) = OPP(cc);
                }
                DRX(n) += VXYZ(1, iv) + VXYZ(1, cn2);
                DRY(n) += VXYZ(2, iv) + VXYZ(2, cn2);
                DRZ(n) += VXYZ(3, iv) + VXYZ(3, cn2);
            }
        }
    } else if (n == N) {
        const int en2 = CFG1(N - 1);
        int cn2 = CFG2(N - 1);
        int cm2;
        if (cn2 < OPP(CFG2(N - 2))) {
            cm2 = OPP(CFG2(N - 2));
        } else {
            cm2 = cn2;
            cn2 = OPP(CFG2(N - 2));
        }
        int iv = 1 + (int)(11.0 * rnd(rng));
        if (iv >= cn2)
            iv += 1;
        if (iv >= cm2)
            iv += 1;

        const int v = node_dir(S, en2, iv);
        const int b = occ[v];
        if ((b == 0) || ((b == 1) && (en2 == v))) {
            const int id = CON1(n);
            cn2 = CFG2(N - 1);
            const int cn3 = CON2(n);
            int cc = 0;
            if (cn3 != 0)
                cc = CNC(iv, cn2, cn3);
            if ((id != 0) && (cc == 0))
                return;
            float dE = CTH(CFG2(N - 2), iv) - CTH(CFG2(N - 2), cn2);
            if (S.istate[n] > 0)
                dE += S.Ei * (float)(nintf[v] - nintf[en]);
            if ((CON3(n) == 1) && (id > 1)) {
                if (CNC(CON2(n), OPP(CFG2(id - 1)), 1) != 0)
                    dE -= S.Ea;
                if (CNC(cc, OPP(CFG2(id - 1)), 1) != 0)
                    dE += S.Ea;
            }
            if (CON3(n - 1) == 1) {
                if (CNC(cn2, CON2(n - 1), 1) != 0)
                    dE -= S.Ea;
                if (CNC(iv, CON2(n - 1), 1) != 0)
                    dE += S.Ea;
            }
            if (rnd(rng) < expf(-dE)) {
                occ[en] -= 1;
                occ[v] += 1;
                nint_move(S, nintf, en, v);
                CFG1(N) = v;
                CFG2(N - 1) = iv;
                if (id != 0) {
                    CON2(n) = cc;
                    CON2(id) = OPP(cc);
                }
                DRX(n) += VXYZ(1, iv) - VXYZ(1, cn2);
                DRY(n) += VXYZ(2, iv) - VXYZ(2, cn2);
                DRZ(n) += VXYZ(3, iv) - VXYZ(3, cn2);
            }
        }
    } else {
        const int cn2 = CFG2(n);
        const int cm2 = CFG2(n - 1);
        const int en2 = CFG1(n - 1);
        const int nm2 = n - 2;
        const int np1 = n + 1;
        if (VNN(1, cm2, cn2) > 1) {
            int iv = 1 + (int)((VNN(1, cm2, cn2) - 1) * rnd(rng));
            if (VNN(2 * iv, cm2, cn2) >= cm2)
                iv += 1;
            const int nv1 = VNN(2 * iv, cm2, cn2);
            const int nv2 = VNN(2 * iv + 1, cm2, cn2);
            const int v = node_dir(S, en2, nv1);
            const int b = occ[v];
            if ((b == 0) || ((b == 1) && ((v == en2) || (v == CFG1(np1))))) {
                const int id = CON1(n);
                const int cn3 = CON2(n);
                int cc = 0;
                if (cn3 != 0)
                    cc = CNC(nv1, cm2, cn3);
                if ((id != 0) && (cc == 0))
                    return;
                float dE;
                if (n == 2) {
                    dE = CTH(nv1, nv2) + CTH(nv2, CFG2(np1)) - CTH(cm2, cn2) -
                         CTH(cn2, CFG2(np1));
                } else if (n == N - 1) {
                    dE = CTH(CFG2(nm2), nv1) + CTH(nv1, nv2) - CTH(CFG2(nm2), cm2) -
                         CTH(cm2, cn2);
                } else {
                    dE = CTH(CFG2(nm2), nv1) + CTH(nv1, nv2) + CTH(nv2, CFG2(np1)) -
                         CTH(CFG2(nm2), cm2) - CTH(cm2, cn2) - CTH(cn2, CFG2(np1));
                }
                if (S.istate[n] > 0)
                    dE += S.Ei * (float)(nintf[v] - nintf[en]);
                if ((CON3(n) == -1) && (id < N)) {
                    if (CNC(OPP(CON2(n)), CFG2(id), 1) != 0)
                        dE -= S.Ea;
                    if (CNC(OPP(cc), CFG2(id), 1) != 0)
                        dE += S.Ea;
                } else if ((CON3(n) == 1) && (id > 1)) {
                    if (CNC(CON2(n), OPP(CFG2(id - 1)), 1) != 0)
                        dE -= S.Ea;
                    if (CNC(cc, OPP(CFG2(id - 1)), 1) != 0)
                        dE += S.Ea;
                }
                if (CON3(n + 1) == -1) {
                    if (CNC(OPP(cn2), CON2(n + 1), 1) != 0)
                        dE -= S.Ea;
                    if (CNC(OPP(nv2), CON2(n + 1), 1) != 0)
                        dE += S.Ea;
                }
                if (CON3(n - 1) == 1) {
                    if (CNC(cm2, CON2(n - 1), 1) != 0)
                        dE -= S.Ea;
                    if (CNC(nv1, CON2(n - 1), 1) != 0)
                        dE += S.Ea;
                }
                if (rnd(rng) < expf(-dE)) {
                    occ[en] -= 1;
                    occ[v] += 1;
                    nint_move(S, nintf, en, v);
                    CFG1(n) = v;
                    CFG2(n - 1) = nv1;
                    CFG2(n) = nv2;
                    if (id != 0) {
                        CON2(n) = cc;
                        CON2(id) = OPP(cc);
                    }
                    DRX(n) += VXYZ(1, nv1) - VXYZ(1, cm2);
                    DRY(n) += VXYZ(2, nv1) - VXYZ(2, cm2);
                    DRZ(n) += VXYZ(3, nv1) - VXYZ(3, cm2);
                }
            }
        }
    }
}

__device__ void trialmoveex(const McShared &S, int *config, int *contact,
                            unsigned long long *rng)
{
    const int N = S.N;
    const int n = 1 + (int)(N * rnd(rng));
    const int s = CON3(n);
    if ((s == 0) || (n == 1) || (n == N))
        return;

    const int strand = (s + 1) / 2 + 1;
    const float impermeability = fabsf(S.boundary[n * 2 + (strand - 1)]);
    const float fc = powf(1.f - impermeability, 1.f / (float)S.ikm);
    for (int i = 1; i <= S.ikm; i++) {
        if (rnd(rng) >= (double)(S.km * fc))
            return;
    }

    if (s == -1) {
        if (CNC(1, CFG2(n - 1), CON2(n)) == 0)
            return;
        if (CON1(n - 1) != 0) {
            if (!S.z_loop)
                return;
            if ((CON3(n - 1) == -1) || (n == 2))
                return;
            if (CNC(1, OPP(CFG2(n - 1)), CON2(n - 1)) == 0)
                return;
            const int con1 = CON1(n - 1);
            const int con2 = CON2(n - 1);
            const int con3 = CON3(n - 1);

            int iv = CNC(1, CFG2(n - 1), CON2(n));
            int id = CON1(n);
            CON1(n - 1) = id;
            CON2(n - 1) = iv;
            CON3(n - 1) = -1;
            CON1(id) = n - 1;
            CON2(id) = OPP(iv);

            iv = CNC(1, OPP(CFG2(n - 1)), con2);
            id = con1;
            CON1(n) = id;
            CON2(n) = iv;
            CON3(n) = S.unidirectional ? con3 : 1;
            CON1(id) = n;
            CON2(id) = OPP(iv);
        } else {
            const int iv = CNC(1, CFG2(n - 1), CON2(n));
            const int id = CON1(n);
            CON1(n - 1) = id;
            CON2(n - 1) = iv;
            CON3(n - 1) = -1;
            CON1(id) = n - 1;
            CON2(id) = OPP(iv);
            CON1(n) = 0;
            CON2(n) = 0;
            CON3(n) = 0;
        }
    } else if (s == 1) {
        if (CNC(1, OPP(CFG2(n)), CON2(n)) == 0)
            return;
        if (CON1(n + 1) != 0) {
            if (!S.z_loop)
                return;
            if ((CON3(n + 1) == 1) || (n == (N - 1)))
                return;
            /* Z-loop (+) swap: use bond n->n+1 = config(2,n), not config(2,n+1) */
            if (CNC(1, CFG2(n), CON2(n + 1)) == 0)
                return;
            const int con1 = CON1(n + 1);
            const int con2 = CON2(n + 1);
            const int con3 = CON3(n + 1);

            int iv = CNC(1, OPP(CFG2(n)), CON2(n));
            int id = CON1(n);
            CON1(n + 1) = id;
            CON2(n + 1) = iv;
            CON3(n + 1) = 1;
            CON1(id) = n + 1;
            CON2(id) = OPP(iv);

            iv = CNC(1, CFG2(n), con2);
            id = con1;
            CON1(n) = id;
            CON2(n) = iv;
            CON3(n) = S.unidirectional ? con3 : -1;
            CON1(id) = n;
            CON2(id) = OPP(iv);
        } else {
            const int iv = CNC(1, OPP(CFG2(n)), CON2(n));
            const int id = CON1(n);
            CON1(n + 1) = id;
            CON2(n + 1) = iv;
            CON3(n + 1) = 1;
            CON1(id) = n + 1;
            CON2(id) = OPP(iv);
            CON1(n) = 0;
            CON2(n) = 0;
            CON3(n) = 0;
        }
    }
}

__device__ void trialbound(const McShared &S, int *config, int *contact, int *nfree,
                           unsigned long long *rng)
{
    const int N = S.N;
    const int n = 1 + (int)(N * rnd(rng));
    const int id = CON1(n);
    if (id != 0)
        return;
    const float kbp = S.kb * powf(S.loadfac[n], 1.f / (float)S.ikb);
    for (int j = 1; j <= S.ikb; j++) {
        if (rnd(rng) >= (double)kbp)
            return;
    }
    const int d = 1 + (int)(2.0 * rnd(rng));
    if (d == 1) {
        if (n > 1) {
            if (CON1(n - 1) == 0) {
                CON1(n) = n - 1;
                CON2(n) = OPP(CFG2(n - 1));
                CON1(n - 1) = n;
                CON2(n - 1) = CFG2(n - 1);
                if (S.unidirectional) {
                    if (rnd(rng) < 0.5) {
                        CON3(n) = 2;
                        CON3(n - 1) = -1;
                    } else {
                        CON3(n) = 1;
                        CON3(n - 1) = -2;
                    }
                } else {
                    CON3(n) = 1;
                    CON3(n - 1) = -1;
                }
                nfree[0] -= 1;
            }
        }
    } else if (n < N) {
        if (CON1(n + 1) == 0) {
            CON1(n) = n + 1;
            CON2(n) = CFG2(n);
            CON1(n + 1) = n;
            CON2(n + 1) = OPP(CFG2(n));
            if (S.unidirectional) {
                if (rnd(rng) < 0.5) {
                    CON3(n) = -2;
                    CON3(n + 1) = 1;
                } else {
                    CON3(n) = -1;
                    CON3(n + 1) = 2;
                }
            } else {
                CON3(n) = -1;
                CON3(n + 1) = 1;
            }
            nfree[0] -= 1;
        }
    }
}

__device__ void trialunbound(const McShared &S, int *contact, int *nfree,
                             unsigned long long *rng, int *ev_n, int *ev_id,
                             int *ev_cnt, int evmax)
{
    const int n = 1 + (int)(S.N * rnd(rng));
    const int id = CON1(n);
    if (id > 0) {
        for (int j = 1; j <= S.iku; j++) {
            if (rnd(rng) >= (double)S.ku)
                return;
        }
        CON1(n) = 0;
        CON2(n) = 0;
        CON3(n) = 0;
        CON1(id) = 0;
        CON2(id) = 0;
        CON3(id) = 0;
        record_unbind(ev_n, ev_id, ev_cnt, evmax, n, id);
        nfree[0] += 1;
    }
}

__device__ void unbound_all(const McShared &S, int *contact, int *nfree, int *ev_n,
                            int *ev_id, int *ev_cnt, int evmax)
{
    for (int n = 1; n <= S.N; n++) {
        const int id = CON1(n);
        if (id > 0) {
            CON1(n) = 0;
            CON2(n) = 0;
            CON3(n) = 0;
            CON1(id) = 0;
            CON2(id) = 0;
            CON3(id) = 0;
            record_unbind(ev_n, ev_id, ev_cnt, evmax, n, id);
            nfree[0] += 1;
        }
    }
}

__device__ void run_steps(const McShared &S, int *config, int *contact, float *dr,
                          int *occ, int *nintf, int *nfree, unsigned long long *rng,
                          int nloc, int mode, int step0, int coherent, int *ev_n,
                          int *ev_id, int *ev_cnt, int evmax)
{
    const int N = S.N;
    if (mode == MC_CUDA_MODE_UNBOUND_ALL) {
        unbound_all(S, contact, nfree, ev_n, ev_id, ev_cnt, evmax);
        return;
    }
    if (mode == MC_CUDA_MODE_MONOMER) {
        for (int j = 0; j < nloc; j++) {
            for (int v = 1; v <= N; v++)
                trialmovetad(S, config, contact, dr, occ, nintf, rng);
        }
        return;
    }
    for (int k = 0; k < nloc; k++) {
        const int ntrial = 3 * N + nfree[0];
        const double pt = (double)N / (double)ntrial;
        const unsigned long long base =
            (unsigned long long)(step0 + k) * (unsigned long long)(4 * N + 1024);
        for (int v = 0; v < ntrial; v++) {
            const double r = coherent ? hash_u01(base + (unsigned long long)v) : rnd(rng);
            if (r < pt)
                trialmovetad(S, config, contact, dr, occ, nintf, rng);
            else if (r < 2.0 * pt)
                trialmoveex(S, config, contact, rng);
            else if (r < 3.0 * pt)
                trialunbound(S, contact, nfree, rng, ev_n, ev_id, ev_cnt, evmax);
            else
                trialbound(S, config, contact, nfree, rng);
        }
    }
}

__global__ void k_run(McShared S, int *config, int *contact, float *dr, int *occ,
                      int *nintf, int *nfree, unsigned long long *rng, int *nsteps_t,
                      int use_nsteps_t, int *ev_n, int *ev_id, int *ev_cnt, int evmax,
                      int T, int Npad, int Apad, int mode, int nsteps, int step0,
                      int coherent)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= T)
        return;
    int nloc = nsteps;
    if (use_nsteps_t)
        nloc = nsteps_t[t];
    run_steps(S, config + t * Npad * 2, contact + t * Npad * 3, dr + t * Npad * 3,
              occ + t * Apad, nintf + t * Apad, nfree + t, rng + t * 2, nloc, mode,
              step0, coherent, ev_n + t * evmax, ev_id + t * evmax, ev_cnt + t, evmax);
}

/* ---- host API ---- */

int mc_cuda_compiled(void) { return 1; }

int mc_cuda_device_count(void)
{
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess)
        return 0;
    return n < 0 ? 0 : n;
}

int mc_cuda_available(void) { return mc_cuda_device_count() > 0 ? 1 : 0; }

int mc_cuda_create(McCudaState **out, const McCudaParams *p, const int *opp,
                   const float *voisxyz, const float *costhet, const int *voisnn,
                   const int *connec, const int *bittable, const float *boundary,
                   const float *loadfac, const int *istate)
{
    if (!out || !p || p->T <= 0 || p->N <= 0 || p->naddr <= 0)
        return -1;
    *out = NULL;

    McCudaState *s = (McCudaState *)std::calloc(1, sizeof(McCudaState));
    if (!s)
        return -1;
    s->p = *p;
    s->Npad = p->N + 1;
    s->Apad = p->naddr + 1;
    s->use_nsteps_t = 0;
    s->sh.N = p->N;
    s->sh.naddr = p->naddr;
    s->sh.iku = p->iku;
    s->sh.ikm = p->ikm;
    s->sh.ikb = p->ikb;
    s->sh.kb = p->kb;
    s->sh.ku = p->ku;
    s->sh.km = p->km;
    s->sh.Ea = p->Ea;
    s->sh.Ei = p->Ei;
    s->sh.z_loop = p->z_loop;
    s->sh.unidirectional = p->unidirectional;

    int h_opp[DIRPAD];
    float h_vxyz[DIRPAD * 4];
    float h_cth[DIRPAD * DIRPAD];
    int *h_vnn = (int *)std::calloc((size_t)DIRPAD * DIRPAD * DIRPAD, sizeof(int));
    int *h_cnc = (int *)std::calloc((size_t)DIRPAD * DIRPAD * DIRPAD, sizeof(int));
    int *h_nbr = (int *)std::calloc((size_t)s->Apad * DIRPAD, sizeof(int));
    float *h_bnd = (float *)std::calloc((size_t)s->Npad * 2, sizeof(float));
    float *h_lf = (float *)std::calloc((size_t)s->Npad, sizeof(float));
    int *h_is = (int *)std::calloc((size_t)s->Npad, sizeof(int));
    if (!h_vnn || !h_cnc || !h_nbr || !h_bnd || !h_lf || !h_is) {
        std::free(h_vnn);
        std::free(h_cnc);
        std::free(h_nbr);
        std::free(h_bnd);
        std::free(h_lf);
        std::free(h_is);
        std::free(s);
        return -1;
    }
    copy_opp(h_opp, opp);
    copy_voisxyz(h_vxyz, voisxyz);
    copy_costhet(h_cth, costhet);
    copy_table3(h_vnn, voisnn);
    copy_table3(h_cnc, connec);
    copy_nbr(h_nbr, bittable, p->naddr, s->Apad);
    for (int n = 1; n <= p->N; n++) {
        h_bnd[n * 2 + 0] = boundary[(n - 1) * 2 + 0];
        h_bnd[n * 2 + 1] = boundary[(n - 1) * 2 + 1];
        h_lf[n] = loadfac[n - 1];
        h_is[n] = istate[n - 1];
    }

    const int T = p->T;
    int rc = 0;
    auto fail = [&]() {
        std::free(h_vnn);
        std::free(h_cnc);
        std::free(h_nbr);
        std::free(h_bnd);
        std::free(h_lf);
        std::free(h_is);
        mc_cuda_destroy(s);
        return -1;
    };

    if (cudaMalloc(&s->sh.opp, sizeof(int) * DIRPAD) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.voisxyz, sizeof(float) * DIRPAD * 4) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.costhet, sizeof(float) * DIRPAD * DIRPAD) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.voisnn, sizeof(int) * DIRPAD * DIRPAD * DIRPAD) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.connec, sizeof(int) * DIRPAD * DIRPAD * DIRPAD) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.nbr, sizeof(int) * s->Apad * DIRPAD) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.boundary, sizeof(float) * s->Npad * 2) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.loadfac, sizeof(float) * s->Npad) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->sh.istate, sizeof(int) * s->Npad) != cudaSuccess)
        return fail();

    if (cudaMemcpy(s->sh.opp, h_opp, sizeof(int) * DIRPAD, cudaMemcpyHostToDevice) !=
        cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.voisxyz, h_vxyz, sizeof(float) * DIRPAD * 4,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.costhet, h_cth, sizeof(float) * DIRPAD * DIRPAD,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.voisnn, h_vnn, sizeof(int) * DIRPAD * DIRPAD * DIRPAD,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.connec, h_cnc, sizeof(int) * DIRPAD * DIRPAD * DIRPAD,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.nbr, h_nbr, sizeof(int) * s->Apad * DIRPAD,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.boundary, h_bnd, sizeof(float) * s->Npad * 2,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.loadfac, h_lf, sizeof(float) * s->Npad,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();
    if (cudaMemcpy(s->sh.istate, h_is, sizeof(int) * s->Npad,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return fail();

    if (cudaMalloc(&s->d_config, sizeof(int) * T * s->Npad * 2) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_contact, sizeof(int) * T * s->Npad * 3) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_dr, sizeof(float) * T * s->Npad * 3) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_occ, sizeof(int) * T * s->Apad) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_nint, sizeof(int) * T * s->Apad) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_nfree, sizeof(int) * T) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_rng, sizeof(unsigned long long) * T * 2) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_nsteps_t, sizeof(int) * T) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_ev_n, sizeof(int) * T * EVMAX) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_ev_id, sizeof(int) * T * EVMAX) != cudaSuccess)
        return fail();
    if (cudaMalloc(&s->d_ev_cnt, sizeof(int) * T) != cudaSuccess)
        return fail();
    if (cudaMemset(s->d_ev_cnt, 0, sizeof(int) * T) != cudaSuccess)
        return fail();
    if (cudaMemset(s->d_config, 0, sizeof(int) * T * s->Npad * 2) != cudaSuccess)
        return fail();
    if (cudaMemset(s->d_contact, 0, sizeof(int) * T * s->Npad * 3) != cudaSuccess)
        return fail();
    if (cudaMemset(s->d_dr, 0, sizeof(float) * T * s->Npad * 3) != cudaSuccess)
        return fail();
    if (cudaMemset(s->d_occ, 0, sizeof(int) * T * s->Apad) != cudaSuccess)
        return fail();
    if (cudaMemset(s->d_nint, 0, sizeof(int) * T * s->Apad) != cudaSuccess)
        return fail();

    std::free(h_vnn);
    std::free(h_cnc);
    std::free(h_nbr);
    std::free(h_bnd);
    std::free(h_lf);
    std::free(h_is);
    (void)rc;
    *out = s;
    return 0;
}

int mc_cuda_upload_traj(McCudaState *s, int t, const int *config, const int *contact,
                        const int *bittable, const float *dr, int nfree,
                        int traj_index)
{
    if (!s || t < 0 || t >= s->p.T)
        return -1;
    const int N = s->p.N;
    const int naddr = s->p.naddr;
    int *h_cfg = (int *)std::calloc((size_t)s->Npad * 2, sizeof(int));
    int *h_con = (int *)std::calloc((size_t)s->Npad * 3, sizeof(int));
    float *h_dr = (float *)std::calloc((size_t)s->Npad * 3, sizeof(float));
    int *h_occ = (int *)std::calloc((size_t)s->Apad, sizeof(int));
    int *h_ni = (int *)std::calloc((size_t)s->Apad, sizeof(int));
    if (!h_cfg || !h_con || !h_dr || !h_occ || !h_ni) {
        std::free(h_cfg);
        std::free(h_con);
        std::free(h_dr);
        std::free(h_occ);
        std::free(h_ni);
        return -1;
    }
    for (int n = 1; n <= N; n++) {
        h_cfg[n * 2 + 0] = config[(n - 1) * 2 + 0];
        h_cfg[n * 2 + 1] = config[(n - 1) * 2 + 1];
        h_con[n * 3 + 0] = contact[(n - 1) * 3 + 0];
        h_con[n * 3 + 1] = contact[(n - 1) * 3 + 1];
        h_con[n * 3 + 2] = contact[(n - 1) * 3 + 2];
        h_dr[n * 3 + 0] = dr[(n - 1) * 3 + 0];
        h_dr[n * 3 + 1] = dr[(n - 1) * 3 + 1];
        h_dr[n * 3 + 2] = dr[(n - 1) * 3 + 2];
    }
    for (int a = 1; a <= naddr; a++) {
        h_occ[a] = bittable[(a - 1) * 14 + 0];
        h_ni[a] = bittable[(a - 1) * 14 + 13];
    }

    uint64_t sm = (uint64_t)(uint32_t)s->p.seed;
    sm += 0x9E3779B97F4A7C15ULL * (uint64_t)(uint32_t)(traj_index + 1);
    unsigned long long rng[2];
    rng[0] = (unsigned long long)splitmix64_host(&sm);
    rng[1] = (unsigned long long)splitmix64_host(&sm);
    if (rng[0] == 0 && rng[1] == 0)
        rng[1] = 1ULL;

    int rc = 0;
    if (cudaMemcpy(s->d_config + t * s->Npad * 2, h_cfg, sizeof(int) * s->Npad * 2,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(s->d_contact + t * s->Npad * 3, h_con, sizeof(int) * s->Npad * 3,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(s->d_dr + t * s->Npad * 3, h_dr, sizeof(float) * s->Npad * 3,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(s->d_occ + t * s->Apad, h_occ, sizeof(int) * s->Apad,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(s->d_nint + t * s->Apad, h_ni, sizeof(int) * s->Apad,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(s->d_nfree + t, &nfree, sizeof(int), cudaMemcpyHostToDevice) !=
        cudaSuccess)
        rc = -1;
    if (cudaMemcpy(s->d_rng + t * 2, rng, sizeof(rng), cudaMemcpyHostToDevice) !=
        cudaSuccess)
        rc = -1;

    std::free(h_cfg);
    std::free(h_con);
    std::free(h_dr);
    std::free(h_occ);
    std::free(h_ni);
    return rc;
}

int mc_cuda_set_nsteps_t(McCudaState *s, const int *nsteps_t)
{
    if (!s)
        return -1;
    if (!nsteps_t) {
        s->use_nsteps_t = 0;
        return 0;
    }
    if (cudaMemcpy(s->d_nsteps_t, nsteps_t, sizeof(int) * s->p.T,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return -1;
    s->use_nsteps_t = 1;
    return 0;
}

int mc_cuda_run(McCudaState *s, int mode, int nsteps, int step0, int coherent)
{
    if (!s)
        return -1;
    const int T = s->p.T;
    const int threads = 64;
    const int blocks = (T + threads - 1) / threads;
    k_run<<<blocks, threads>>>(s->sh, s->d_config, s->d_contact, s->d_dr, s->d_occ,
                               s->d_nint, s->d_nfree, s->d_rng, s->d_nsteps_t,
                               s->use_nsteps_t, s->d_ev_n, s->d_ev_id, s->d_ev_cnt,
                               EVMAX, T, s->Npad, s->Apad, mode, nsteps, step0,
                               coherent);
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess) {
        std::fprintf(stderr, "CUDA kernel launch: %s\n", cudaGetErrorString(e));
        return -1;
    }
    e = cudaDeviceSynchronize();
    if (e != cudaSuccess) {
        std::fprintf(stderr, "CUDA synchronize: %s\n", cudaGetErrorString(e));
        return -1;
    }
    return 0;
}

int mc_cuda_download_traj(McCudaState *s, int t, int *config, int *contact,
                          int *bittable, float *dr, int *nfree)
{
    if (!s || t < 0 || t >= s->p.T)
        return -1;
    const int N = s->p.N;
    const int naddr = s->p.naddr;
    int *h_cfg = (int *)std::calloc((size_t)s->Npad * 2, sizeof(int));
    int *h_con = (int *)std::calloc((size_t)s->Npad * 3, sizeof(int));
    float *h_dr = (float *)std::calloc((size_t)s->Npad * 3, sizeof(float));
    if (!h_cfg || !h_con || !h_dr) {
        std::free(h_cfg);
        std::free(h_con);
        std::free(h_dr);
        return -1;
    }
    int rc = 0;
    if (cudaMemcpy(h_cfg, s->d_config + t * s->Npad * 2, sizeof(int) * s->Npad * 2,
                   cudaMemcpyDeviceToHost) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(h_con, s->d_contact + t * s->Npad * 3, sizeof(int) * s->Npad * 3,
                   cudaMemcpyDeviceToHost) != cudaSuccess)
        rc = -1;
    if (cudaMemcpy(h_dr, s->d_dr + t * s->Npad * 3, sizeof(float) * s->Npad * 3,
                   cudaMemcpyDeviceToHost) != cudaSuccess)
        rc = -1;
    if (nfree &&
        cudaMemcpy(nfree, s->d_nfree + t, sizeof(int), cudaMemcpyDeviceToHost) !=
            cudaSuccess)
        rc = -1;
    if (rc == 0) {
        for (int n = 1; n <= N; n++) {
            if (config) {
                config[(n - 1) * 2 + 0] = h_cfg[n * 2 + 0];
                config[(n - 1) * 2 + 1] = h_cfg[n * 2 + 1];
            }
            if (contact) {
                contact[(n - 1) * 3 + 0] = h_con[n * 3 + 0];
                contact[(n - 1) * 3 + 1] = h_con[n * 3 + 1];
                contact[(n - 1) * 3 + 2] = h_con[n * 3 + 2];
            }
            if (dr) {
                dr[(n - 1) * 3 + 0] = h_dr[n * 3 + 0];
                dr[(n - 1) * 3 + 1] = h_dr[n * 3 + 1];
                dr[(n - 1) * 3 + 2] = h_dr[n * 3 + 2];
            }
        }
    }
    if (rc == 0 && bittable) {
        int *h_occ = (int *)std::calloc((size_t)s->Apad, sizeof(int));
        int *h_ni = (int *)std::calloc((size_t)s->Apad, sizeof(int));
        if (!h_occ || !h_ni)
            rc = -1;
        else {
            if (cudaMemcpy(h_occ, s->d_occ + t * s->Apad, sizeof(int) * s->Apad,
                           cudaMemcpyDeviceToHost) != cudaSuccess)
                rc = -1;
            if (cudaMemcpy(h_ni, s->d_nint + t * s->Apad, sizeof(int) * s->Apad,
                           cudaMemcpyDeviceToHost) != cudaSuccess)
                rc = -1;
            if (rc == 0) {
                for (int a = 1; a <= naddr; a++) {
                    bittable[(a - 1) * 14 + 0] = h_occ[a];
                    bittable[(a - 1) * 14 + 13] = h_ni[a];
                }
            }
        }
        std::free(h_occ);
        std::free(h_ni);
    }
    std::free(h_cfg);
    std::free(h_con);
    std::free(h_dr);
    return rc;
}

int mc_cuda_download_unbind(McCudaState *s, int t, int max_events, int *n_out,
                            int *id_out, int *count)
{
    if (!s || t < 0 || t >= s->p.T || !count)
        return -1;
    int cnt = 0;
    if (cudaMemcpy(&cnt, s->d_ev_cnt + t, sizeof(int), cudaMemcpyDeviceToHost) !=
        cudaSuccess)
        return -1;
    int ncopy = cnt;
    if (ncopy > max_events)
        ncopy = max_events;
    if (ncopy > 0 && n_out && id_out) {
        if (cudaMemcpy(n_out, s->d_ev_n + t * EVMAX, sizeof(int) * ncopy,
                       cudaMemcpyDeviceToHost) != cudaSuccess)
            return -1;
        if (cudaMemcpy(id_out, s->d_ev_id + t * EVMAX, sizeof(int) * ncopy,
                       cudaMemcpyDeviceToHost) != cudaSuccess)
            return -1;
    }
    *count = ncopy;
    cnt = 0;
    if (cudaMemcpy(s->d_ev_cnt + t, &cnt, sizeof(int), cudaMemcpyHostToDevice) !=
        cudaSuccess)
        return -1;
    return 0;
}

void mc_cuda_destroy(McCudaState *s)
{
    if (!s)
        return;
    cudaFree(s->sh.opp);
    cudaFree(s->sh.voisxyz);
    cudaFree(s->sh.costhet);
    cudaFree(s->sh.voisnn);
    cudaFree(s->sh.connec);
    cudaFree(s->sh.nbr);
    cudaFree(s->sh.boundary);
    cudaFree(s->sh.loadfac);
    cudaFree(s->sh.istate);
    cudaFree(s->d_config);
    cudaFree(s->d_contact);
    cudaFree(s->d_occ);
    cudaFree(s->d_nint);
    cudaFree(s->d_nfree);
    cudaFree(s->d_dr);
    cudaFree(s->d_rng);
    cudaFree(s->d_nsteps_t);
    cudaFree(s->d_ev_n);
    cudaFree(s->d_ev_id);
    cudaFree(s->d_ev_cnt);
    std::free(s);
}
