// cuda_step.cu

//#include "objects.h"
//#include "joints/joint.h"
//#include <ode/odeconfig.h>
//#include "config.h"
//#include <ode/odemath.h>
//#include <ode/rotation.h>
//#include <ode/timer.h>
//#include <ode/error.h>
//#include <ode/matrix.h>
//#include "lcp.h"
#include "util.h"

#include <ode/cuda_step.h>
#include <cuda.h>
#include <ode/cuda_helper.h>
// #include <ode/cuda_demo_helper.h>
#include <ode/cuda_matrix.h>
// #include <ode/cuPrintf.cuh>

#define BLOCKSIZE 16

__device__ void printMatrixBase(const char *name, const char *fmt,
                                dReal const *a, const int h, const int w,
                                bool pad) {
  printf("%s (%s):\n", name, (pad ? "padded" : "no-pad"));
  for (int row = 0; row < h; row++) {
    for (int col = 0; col < w; col++)
      printf(fmt, a[(pad ? (row * (dPAD(w)) + col) : (row * w + col))]);
    printf("\n");
  }
  printf("\n");
}

#define show_mat(M, X, Y) printMatrixBase(#M, "%f\t", (M),(X),(Y),false)
#define show_Pmat(M, X, Y) printMatrixBase(#M, "%f\t", (M),(X),(Y),true)

__device__ void dQMultiply0 (dQuaternion qa, const dQuaternion qb, const dQuaternion qc) {
  qa[0] = qb[0]*qc[0] - qb[1]*qc[1] - qb[2]*qc[2] - qb[3]*qc[3];
  qa[1] = qb[0]*qc[1] + qb[1]*qc[0] + qb[2]*qc[3] - qb[3]*qc[2];
  qa[2] = qb[0]*qc[2] + qb[2]*qc[0] + qb[3]*qc[1] - qb[1]*qc[3];
  qa[3] = qb[0]*qc[3] + qb[3]*qc[0] + qb[1]*qc[2] - qb[2]*qc[1];
}

__device__ void dDQfromW (dReal dq[4], const dVector3 w, const dQuaternion q)
{
  dq[0] = REAL(0.5)*(- w[0]*q[1] - w[1]*q[2] - w[2]*q[3]);
  dq[1] = REAL(0.5)*(  w[0]*q[0] + w[1]*q[3] - w[2]*q[2]);
  dq[2] = REAL(0.5)*(- w[0]*q[3] + w[1]*q[0] + w[2]*q[1]);
  dq[3] = REAL(0.5)*(  w[0]*q[2] - w[1]*q[1] + w[2]*q[0]);
}

__device__ void dWtoDQ(const dVector3 w, const dQuaternion q, dReal dq[4]) {
	return dDQfromW(dq,w,q);
}

__device__ void dRfromQ (dMatrix3 R, const dQuaternion q)
{
  // q = (s,vx,vy,vz)
  dReal qq1 = 2*q[1]*q[1];
  dReal qq2 = 2*q[2]*q[2];
  dReal qq3 = 2*q[3]*q[3];
  R[(0)*4+(0)] = 1 - qq2 - qq3;
  R[(0)*4+(1)] = 2*(q[1]*q[2] - q[0]*q[3]);
  R[(0)*4+(2)] = 2*(q[1]*q[3] + q[0]*q[2]);
  R[(0)*4+(3)] = (0.0);
  R[(1)*4+(0)] = 2*(q[1]*q[2] + q[0]*q[3]);
  R[(1)*4+(1)] = 1 - qq1 - qq3;
  R[(1)*4+(2)] = 2*(q[2]*q[3] - q[0]*q[1]);
  R[(1)*4+(3)] = (0.0);
  R[(2)*4+(0)] = 2*(q[1]*q[3] - q[0]*q[2]);
  R[(2)*4+(1)] = 2*(q[2]*q[3] + q[0]*q[1]);
  R[(2)*4+(2)] = 1 - qq1 - qq2;
  R[(2)*4+(3)] = (0.0);
}

__device__ void dQtoR(const dQuaternion q, dMatrix3 R) {
	return dRfromQ(R, q);
}

// __device__ dReal dDOTpq(dReal *a, dReal *b, int p, int q) {
// 	return ((a)[0]*(b)[0] + (a)[p]*(b)[q] + (a)[2*(p)]*(b)[2*(q)]);
// }

#define dDOTpq(a,b,p,q)  (a[0]*b[0] + a[p]*b[q] + a[2*(p)]*b[2*(q)])

// __device__ dReal dDOT(dReal *a, dReal *b) {
// 	return dDOTpq(a,b,1,1);
// }


#define dDOT11(a,b) dDOTpq(a,b,1,1)
#define dDOT(a,b) dDOT11(a,b)
#define dDOT13(a,b) dDOTpq(a,b,1,3)
#define dDOT14(a,b) dDOTpq(a,b,1,4)
#define dDOT41(a,b) dDOTpq(a,b,4,1)

__device__ int dNormalize4(dVector4 a) {
  dReal l = dDOT(a,a)+a[3]*a[3];
  if (l > 0) {
    //l = dRecipSqrt(l);
	l = ((1.0f/sqrtf(l)));
    a[0] *= l;
    a[1] *= l;
    a[2] *= l;
    a[3] *= l;
	return 1;
  }
  else {
    a[0] = 1;
    a[1] = 0;
    a[2] = 0;
    a[3] = 0;
    return 0;
  }
}

/*ODE_API void cuda_dxProcessIslands (dxWorld *world, dReal stepsize, dstepper_fn_t cuda_stepper)
{
	int cuda_bodies_count = 0;
	dxBody *cuda_bodies;
	cudaMalloc((void**) &cuda_bodies, sizeof(dxBody)*world->nb);
	dxBody *bb;
	for (bb=world->firstbody;bb;bb=(dxBody*)bb->next)
		cudaMemcpy(cuda_bodies+sizeof(dxBody)*cuda_bodies_count++, bb, sizeof(dxBody), cudaMemcpyHostToDevice);
    cuda_stepper (world,cuda_bodies,cuda_bodies_count,NULL,NULL,stepsize);*/

/* special-case matrix multiplication functions */

// A = B*C  A, B, C all 3x3
// B: pad	A,C: nopad
__device__ void cuda_dMultiply0_333(dReal *A, dReal *B, dReal *C) {

	A[0] = dDOT13((B),(C)); 
	A[1] = dDOT13((B),(C+1)); 
	A[2] = dDOT13((B),(C+2)); 

	A[3] = dDOT13((B+4),(C)); 
	A[4] = dDOT13((B+4),(C+1)); 
	A[5] = dDOT13((B+4),(C+2)); 

	A[6] = dDOT13((B+8),(C)); 
	A[7] = dDOT13((B+8),(C+1)); 
	A[8] = dDOT13((B+8),(C+2)); 
}

// A = B*C^T  A, B, C all 3x3
// A: nopad,	B,C: pad
__device__ void cuda_dMultiply2_333(dReal *A, dReal *B, dReal *C) {
	A[0] = dDOT11((B),(C)); 
	A[1] = dDOT11((B),(C+4)); 
	A[2] = dDOT11((B),(C+8)); 

	A[3] = dDOT11((B+4),(C)); 
	A[4] = dDOT11((B+4),(C+4)); 
	A[5] = dDOT11((B+4),(C+8));

	A[6] = dDOT11((B+8),(C)); 
	A[7] = dDOT11((B+8),(C+4)); 
	A[8] = dDOT11((B+8),(C+8)); 
}

#define CU_ARRAY_DBG(A)                                 \
  printf("---\narray " #A "\n");                        \
  for (i = 0; i < 3; ++i) {                             \
    for (j = 0; j < 3; ++j) {                           \
      printf(#A "[%d,%d]{%f}\t", i, j, (A)[3 * i + j]); \
    }                                                   \
    printf("\n");                                       \
  }

#define MUL_PROC(B, C, OFFSET)                                                 \
  for (j = 0; j < 3; ++j) {                                                    \
    printf("B[%d]*C[%d] = %f * %f = %f\n", j + OFFSET, j, B[j + OFFSET], C[j], \
           B[j + OFFSET] * C[j]);                                              \
  }

// A = B*C  A 3x1, B 3x3, C 3x1
// A,B,C: nopad
__device__ void cuda_dMultiply0_331(dReal *A, dReal  const*B, dReal const*C) {
	A[0] = dDOT11(B, C);
	A[3] = dDOT11((B+3), C);
	A[6] = dDOT11((B+6), C);

	// int i,j;
	// CU_ARRAY_DBG(B);
	// CU_ARRAY_DBG(C);

    // MUL_PROC(B,C,0);
	// MUL_PROC(B,C,3);
	// MUL_PROC(B,C,6);
}

// A = B*C  A 1x3, B 3x3, C 3x1
// A,B,C  nopad
__device__ void cuda_dMultiplyAdd0_331(dReal *A, dReal  const*B, dReal const*C) {
	A[0] += dDOT11(B, C);
	A[1] += dDOT11((B+3), C);
	A[2] += dDOT11((B+6), C);
}

// A = B*C  A 1x3, B 1x3, C 3x3
__device__ void cuda_dMultiply0_133(dReal *A, dReal *B, dReal *C) {
	A[0] = dDOT13((B),(C));
	A[1] = dDOT13((B),(C+1));
	A[2] = dDOT13((B),(C+2));
}


// a -= b cross c
__device__ void cuda_dCross(dReal *a, dReal *b, dReal *c) {
	a[0] -= ((b)[1]*(c)[2] - (b)[2]*(c)[1]);
	a[1] -= ((b)[2]*(c)[0] - (b)[0]*(c)[2]);
	a[2] -= ((b)[0]*(c)[1] - (b)[1]*(c)[0]);
}

// A = B*C  A pxr, B pxq, C qxr
__device__ void naiveMatMultiply(dReal *A, dReal *B, dReal *C, int p, int q, int r) {
  int i, j, k;
  for (j = 0; j < p; ++j)
    for (k = 0; k < r; ++k) {
      for (i = 0; i < q; ++i) {
        A[k + j * r] += B[i] * C[k + r * i];
      }
    }

  // for (i = 0; i < p; i++) {
  // 	for (j = 0; j < r; j++) {
  // 		for (k = 0; k < q; k++) {
  // 			A[i*r + j] += (B[i*q + k])*(C[k*r + j]);
  // 		}
  // 	}
  // }
}

__device__ dReal cuda_sinc(dReal x)
{
	// if |x| < 1e-4 then use a taylor series expansion. this two term expansion
	// is actually accurate to one LS bit within this range if double precision
	// is being used - so don't worry!
	if (fabs(x) < 1.0e-4) return (1.0) - x*x*(0.166666666666666666667);
	else return sinf(x)/x;
}


__device__ void disp_bodyd(dxBody *body) {
  int i;
  printf("flags: %d\n", body->flags);
  printf("mass: %f\n", body->mass);
  printf("InvMass: %f\n", body->invMass);
  printf("posr:\n");
  printf("\tpos: (%f,%f,%f,%f)\n", body->posr.pos[0], body->posr.pos[1],
         body->posr.pos[2], body->posr.pos[3]);
  for (i = 0; i < 3; ++i)
    printf("\tR[%d]: (%f,%f,%f,%f)\n", i, (body->posr.R + i * 4)[0],
           (body->posr.R + i * 4)[1], (body->posr.R + i * 4)[2],
           (body->posr.R + i * 4)[3]);
}

// for debug only
__global__ void cuda_step_none(dxBody *body, int nb, dReal stepsize, dReal g1,
                               dReal g2, dReal g3) {}


// #define _CUDA_DBG
#if defined(_CUDA_DBG)
#define _CUDA_DBG_DO(DO) DO
#else
#define _CUDA_DBG_DO(DO)
#endif


//****************************************************************************
// the slow, but sure way
// note that this does not do any joint feedback!

// given lists of bodies and joints that form an island, perform a first
// order timestep.
//
// `body' is the body array, `nb' is the size of the array.
// `_joint' is the body array, `nj' is the size of the array.
 __global__ void cuda_step(dxBody *body, int nb, dxJoint *joint, int nj, dReal stepsize, dReal g1, dReal g2, dReal g3)
{
	dVector3 gravity; 
	gravity[0] = g1;
	gravity[1] = g2;
	gravity[2] = g3;
	int i,j,k;

	dReal I[3*3], invI[3*3];

	int bid = threadIdx.x + blockDim.x * blockIdx.x;
	if (bid >= nb) { return; }

	// for all bodies, compute the inertia tensor and its inverse in the global
	// frame, and compute the rotational force and add it to the torque
	// accumulator.
	// @@@ check computation of rotational force.

	//dSetZero (I,3*nb*4);
	//dSetZero (invI,3*nb*4);

	_CUDA_DBG_DO(printf("[IN CUDA]================Before: Body[%d]: \n", bid));
	_CUDA_DBG_DO(disp_bodyd(body+bid));

	dReal tmp[9];
#if defined(_CUDA_DBG)
	show_mat(tmp,3,3);
	printf("before compute inertia tensor \n");
	show_Pmat(body[bid].mass.I,3,3);
	show_Pmat(body[bid].posr.R,3,3);
	printf("=== compute inertia tensor \n");
#endif
    // compute inertia tensor in global frame
    cuda_dMultiply2_333(tmp, body[bid].mass.I, body[bid].posr.R);
    cuda_dMultiply0_333(I, body[bid].posr.R, tmp);

	_CUDA_DBG_DO(printf("I after compute inertia tensor \n")); 
	_CUDA_DBG_DO(show_mat(I,3,3);)
	_CUDA_DBG_DO(printf("tmp after compute inertia tensor \n");)
	_CUDA_DBG_DO( show_mat(tmp,3,3);)

if (body[bid].flags & dxBodyGyroscopic) {
    // compute inverse inertia tensor in global frame
    cuda_dMultiply2_333(tmp, body[bid].invI, body[bid].posr.R);
    cuda_dMultiply0_333(invI, body[bid].posr.R, tmp);
#if defined(_CUDA_DBG)
	_CUDA_DBG_DO(printf("tmp after compute inverse inertia tensor \n"); )
	_CUDA_DBG_DO(show_mat(tmp,3,3);)
#endif
	_CUDA_DBG_DO(show_mat(I,3,3));
	_CUDA_DBG_DO(show_mat(body[bid].avel,3,1));
	for(i=0;i<9;++i)tmp[i]=0;

    // compute rotational force
    cuda_dMultiply0_331(tmp, I, body[bid].avel);
	_CUDA_DBG_DO(printf("tmp after compute rotational force \n"));
	_CUDA_DBG_DO(show_mat(tmp,3,3));
    cuda_dCross(body[bid].tacc, body[bid].avel, tmp);
}

	// add the gravity force to all bodies

    if ((body[bid].flags & dxBodyNoGravity)==0) {
		body[bid].facc[0] += body[bid].mass.mass * gravity[0];
		body[bid].facc[1] += body[bid].mass.mass * gravity[1];
		body[bid].facc[2] += body[bid].mass.mass * gravity[2];
    }

	/// Joint relavant implement
	// here
	/// Joint relavant implement

	// create (6*nb,6*nb) inverse mass matrix `invM', and fill it with mass
	// parameters  
	dReal invM[6*6];

	for(i = 0; i < 6*6; i++) invM[i] = 0;



    invM[0] = body[bid].invMass;
    invM[6+1] = body[bid].invMass;
    invM[2*6+2] = body[bid].invMass;
    for (j = 3; j < 6; j++) for (k = 3; k < 6; k++) {
			invM[j*6+k] = invI[j*6+k];
		}
  

	// assemble some body vectors: fe = external forces, v = velocities
	dReal fe[6];
	dReal v[6];

	//dSetZero (fe,n6);
	//dSetZero (v,n6);

	// calculate fe(external force)
    for (j = 0; j < 3; j++) fe[j] = body[bid].facc[j];
    for (j = 0; j < 3; j++) {
		fe[3+j] = body[bid].tacc[j];
		_CUDA_DBG_DO( printf("update body[%d].tacc[%d] = %f\n",bid,j,body[bid].tacc[j]);)
	}

	// this will be set to the velocity update. vnew means cforce
	dReal vnew[6];

	// here should calculate vnew(cforce) according to joints restriction
	// temperarily just set it to zero
	for(i = 0; i < 6; i++) vnew[i] = 0;
	
	// add fe to vnew(cforce)
	for(i = 0; i < 6; i++) vnew[i] += fe[i];
	// // multiply cforce by stepsize
	for(i = 0; i < 6; i++) vnew[i] *= stepsize;



	_CUDA_DBG_DO(show_mat(invM, 6,6);)
	_CUDA_DBG_DO(show_mat(fe, 6,1);)
	_CUDA_DBG_DO(show_mat(vnew, 1, 6));
	// no constraints
	// naiveMatMultiply(vnew, invM, fe, 6, 6, 1);


	// for (i = 0; i < 6; i++) {
	// 	vnew[i] = v[i] + stepsize * vnew[i];
	// 	_CUDA_DBG_DO(printf("update vnew[%d] = %f + %f * %f\n", i, v[i], stepsize, vnew[i]));
	// }

	// apply the velocity update to the bodies

	for (j = 0; j < 3; j++) {
		body[bid].lvel[j] += vnew[j] * (body[bid].invMass);
		_CUDA_DBG_DO( printf("update body[%d].lvel[%d] to %f\n", bid, j, vnew[j]);)
	}
	// add invM * cforce to the body velocity
	cuda_dMultiplyAdd0_331(body[bid].avel, invI, (vnew+3));
    // for (j = 0; j < 3; j++) {
	// 	body[bid].avel[j] += vnew[3 + j];
	// 	_CUDA_DBG_DO(printf("update body[%d].avel[%d] to %f\n", bid, j, vnew[3+j]));
	// }

	// ! impl of void dxStepBody (dxBody *b, dReal h) in CUDA
	// update the position and orientation from the new linear/angular velocity
	// (over the given timestep)
	//dxBody *b = &(body[bid]);

	// cap the angular velocity
	if (body[bid].flags & dxBodyMaxAngularSpeed) {
        const dReal max_ang_speed = body[bid].max_angular_speed;
        const dReal aspeed = dDOT( body[bid].avel, body[bid].avel );
        if (aspeed > max_ang_speed*max_ang_speed) {
			const dReal coef = max_ang_speed/sqrtf(aspeed);
			//dOPEC(b.avel, *=, coef); // multiply vector by scalar coef
			body[bid].avel[0] *= coef;
			body[bid].avel[1] *= coef;
			body[bid].avel[2] *= coef;
        }
	}
	// end of angular velocity cap

	dReal h = stepsize;

	// handle linear velocity
	for (j=0; j<3; j++) body[bid].posr.pos[j] += h * body[bid].lvel[j];

	if (body[bid].flags & dxBodyFlagFiniteRotation) {
		dVector3 irv;	// infitesimal rotation vector
		dQuaternion q;	// quaternion for finite rotation

		if (body[bid].flags & dxBodyFlagFiniteRotationAxis) {
			// split the angular velocity vector into a component along the finite
			// rotation axis, and a component orthogonal to it.
			dVector3 frv;		// finite rotation vector
			dReal k = dDOT (body[bid].finite_rot_axis,body[bid].avel);
			frv[0] = body[bid].finite_rot_axis[0] * k;
			frv[1] = body[bid].finite_rot_axis[1] * k;
			frv[2] = body[bid].finite_rot_axis[2] * k;
			irv[0] = body[bid].avel[0] - frv[0];
			irv[1] = body[bid].avel[1] - frv[1];
			irv[2] = body[bid].avel[2] - frv[2];

			// make a rotation quaternion q that corresponds to frv * h.
			// compare this with the full-finite-rotation case below.
			h *= REAL(0.5);
			dReal theta = k * h;
			q[0] = cosf(theta);
			dReal s = cuda_sinc(theta) * h;
			q[1] = frv[0] * s;
			q[2] = frv[1] * s;
			q[3] = frv[2] * s;
		}
		else {
			// make a rotation quaternion q that corresponds to w * h
			dReal wlen = sqrtf (body[bid].avel[0]*body[bid].avel[0] + body[bid].avel[1]*body[bid].avel[1] +
								body[bid].avel[2]*body[bid].avel[2]);
			h *= REAL(0.5);
			dReal theta = wlen * h;
			q[0] = cosf(theta);
			dReal s = cuda_sinc(theta) * h;
			q[1] = body[bid].avel[0] * s;
			q[2] = body[bid].avel[1] * s;
			q[3] = body[bid].avel[2] * s;
		}

		// do the finite rotation
		dQuaternion q2;
		dQMultiply0 (q2,q,body[bid].q);
		for (j=0; j<4; j++) body[bid].q[j] = q2[j];

		// do the infitesimal rotation if required
		if (body[bid].flags & dxBodyFlagFiniteRotationAxis) {
			dReal dq[4];
			dWtoDQ (irv,body[bid].q,dq);
			for (j=0; j<4; j++) body[bid].q[j] += h * dq[j];
		}
	}
	else {
		// the normal way - do an infitesimal rotation
		dReal dq[4];
		dWtoDQ (body[bid].avel,body[bid].q,dq);
		for (j=0; j<4; j++) body[bid].q[j] += h * dq[j];
	}

	// normalize the quaternion and convert it to a rotation matrix
	dNormalize4 (body[bid].q);
	dQtoR (body[bid].q,body[bid].posr.R);

	// damping
	if (body[bid].flags & dxBodyLinearDamping) {
		const dReal lin_threshold = body[bid].dampingp.linear_threshold;
        const dReal lin_speed = dDOT( body[bid].lvel, body[bid].lvel );
        if ( lin_speed > lin_threshold) {
			const dReal k = 1 - body[bid].dampingp.linear_scale;
			//dOPEC(b.lvel, *=, k);
			body[bid].lvel[0] *= k;
			body[bid].lvel[1] *= k;
			body[bid].lvel[2] *= k;
        }
	}
	if (body[bid].flags & dxBodyAngularDamping) {
		const dReal ang_threshold = body[bid].dampingp.angular_threshold;
        const dReal ang_speed = dDOT( body[bid].avel, body[bid].avel );
        if ( ang_speed > ang_threshold) {
			const dReal k = 1 - body[bid].dampingp.angular_scale;
			//dOPEC(b.avel, *=, k);
			body[bid].avel[0] *= k;
			body[bid].avel[1] *= k;
			body[bid].avel[2] *= k;
        }
	}

	// zero all force accumulators
    body[bid].facc[0] = 0;
    body[bid].facc[1] = 0;
    body[bid].facc[2] = 0;
    body[bid].facc[3] = 0;
    body[bid].tacc[0] = 0;
    body[bid].tacc[1] = 0;
    body[bid].tacc[2] = 0;
    body[bid].tacc[3] = 0;

	_CUDA_DBG_DO( printf("================After: Body[%d]: \n", bid);)
	_CUDA_DBG_DO(disp_bodyd(body+bid);)

}

ODE_API void cuda_dInternalStepIsland_x1 (dxWorld *world, dxBody *cuda_body, int nb, dxJoint *_joint, int nj, dReal stepsize)
{

	cuda_step<<<nb, 1>>>(cuda_body, world->nb, _joint, world->nj, stepsize, world->gravity[0], world->gravity[1], world->gravity[2]);
	// cuda_step<<<1, 1>>>(cuda_body, world->nb, stepsize, world->gravity[0], world->gravity[1], world->gravity[2]);

	//cuda_step<<<BLOCKSIZE/nb, 256>>>(cuda_body, world->nb, stepsize, world->gravity[0], world->gravity[1], world->gravity[2]);
}

 ODE_API void cuda_dxProcessIslands(dxWorld *world, dxBody *cuda_body, dxJoint *cuda_joint, dReal stepsize, dstepper_fn_t stepper)
{
	const int block_size = BLOCKSIZE;
	dim3 dimBlock(block_size, block_size);
	dim3 dimGrid(block_size, block_size);

	cuda_dInternalStepIsland_x1(world, cuda_body, world->nb, cuda_joint, world->nj, stepsize);
}

