#ifndef MC_CUDA_H
#define MC_CUDA_H

/* C ABI for the optional CUDA Monte-Carlo backend.
 * Trajectory index t is 0-based. Fortran arrays are column-major. */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct McCudaState McCudaState;

typedef struct McCudaParams {
    int T;
    int N;
    int naddr;
    int nlef0;
    int iku;
    int ikm;
    int ikb;
    int z_loop;
    int unidirectional;
    int seed;
    float kb;
    float ku;
    float km;
    float Ea;
    float Ei;
} McCudaParams;

enum {
    MC_CUDA_MODE_MIXED = 0,
    MC_CUDA_MODE_MONOMER = 1,
    MC_CUDA_MODE_UNBOUND_ALL = 2
};

/* 1 if this binary was linked with the real CUDA implementation. */
int mc_cuda_compiled(void);

/* Number of CUDA devices, or 0 if none / driver missing. */
int mc_cuda_device_count(void);

/* 1 if at least one device is visible. */
int mc_cuda_available(void);

int mc_cuda_create(McCudaState **out,
                   const McCudaParams *p,
                   const int *opp,
                   const float *voisxyz,
                   const float *costhet,
                   const int *voisnn,
                   const int *connec,
                   const int *bittable,
                   const float *boundary,
                   const float *loadfac,
                   const int *istate);

int mc_cuda_upload_traj(McCudaState *s, int t,
                        const int *config,
                        const int *contact,
                        const int *bittable,
                        const float *dr,
                        int nfree,
                        int traj_index);

/* If nsteps_t is non-NULL, each thread uses nsteps_t[t] (length T).
 * Pass NULL to clear and use the nsteps argument of mc_cuda_run. */
int mc_cuda_set_nsteps_t(McCudaState *s, const int *nsteps_t);

int mc_cuda_run(McCudaState *s, int mode, int nsteps, int step0, int coherent);

int mc_cuda_download_traj(McCudaState *s, int t,
                          int *config,
                          int *contact,
                          int *bittable,
                          float *dr,
                          int *nfree);

/* Copy and clear the per-trajectory unbind event buffer (process.out). */
int mc_cuda_download_unbind(McCudaState *s, int t,
                            int max_events, int *n_out, int *id_out, int *count);

void mc_cuda_destroy(McCudaState *s);

#ifdef __cplusplus
}
#endif

#endif /* MC_CUDA_H */
