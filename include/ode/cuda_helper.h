#ifndef _CUDA_HELPER_H_
#define _CUDA_HELPER_H_

#include <ode/common.h>

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__cplusplus)
#define RM_STRUCT(X) \
typedef struct X X;

RM_STRUCT(dxBody)
RM_STRUCT(dxWorld)
#endif
ODE_API dReal* cuda_realMalloc(unsigned int real_n);
ODE_API void cuda_testMemcpy();

ODE_API dReal *cuda_copyToDevice(dReal *a, int n);
ODE_API dReal *cuda_copyPaddedToDevice(dReal *a, int dim);
ODE_API dReal *cuda_copyFromDevice(dReal *dev_a, dReal *a, int n);
ODE_API void cuda_freeFromDevice(dReal *dev_a);
ODE_API dReal *cuda_makeOnDevice(int n);

ODE_API dxBody *cuda_copyBodiesToDevice(dxBody *cuda_body, dxBody **body, int NUM);
ODE_API dxBody *cuda_copyWorldBodiesToDevice(dxBody *cuda_body, dxWorld *world);
ODE_API dxJoint *cuda_copyWorldJointsToDevice(dxJoint *cuda_joint, dxWorld *world);

ODE_API dxJoint **cuda_copyWorldJointsFromDevice(dxWorld *world, dxJoint *cuda_joint, dxJoint *j_buff);
ODE_API dxBody **cuda_copyBodiesFromDevice(dxBody **body, dxBody *cuda_body, int NUM, dxBody *b_buff);
ODE_API dxBody **cuda_copyWorldBodiesFromDevice(dxWorld *world, dxBody *cuda_body, int NUM, dxBody *b_buff);

ODE_API dxBody *cuda_initBodiesOnDevice(int NUM);
ODE_API dxJoint *cuda_initJointsOnDevice(int NUM);
ODE_API void cuda_free(dxBody *ptr);
ODE_API void cuda_jfree(dxJoint *ptr);
ODE_API void cuda_free2(dReal *ptr);
ODE_API void *cuda_malloc(void **dev, int NUM);

#ifdef __cplusplus
}
#endif
#endif

