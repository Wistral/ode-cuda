// sample6.cpp  by Kosei Demura 2005-2008
// My web site is http://demura.net
// This program uses the Open Dynamics Engine (ODE) by Russell Smith.
// The ODE web site is http://ode.org/

#include <drawstuff/drawstuff.h>
#include <ode/ode.h>
#include <iostream>

#include "../src/objects.h"
#include "joints/joint.h"

#ifdef dDOUBLE
#define dsDrawSphere dsDrawSphereD
#endif

// #define USECUDA

static dWorldID world;
static dSpaceID space;
static dGeomID ground;
static dJointGroupID contactgroup;
static int flag = 0;
dsFunctions fn;

const dReal radius = 0.2;
const dReal mass = 1.0;

typedef struct {
  dBodyID body;
  dBodyID bbuf;

  dJointID joint;
  dJointID jbuf;
  dGeomID geom;
} MyObject;
MyObject ball, ball_cuda;

static void nearCallback(void *data, dGeomID o1, dGeomID o2) {
  const int N = 10;
  dContact contact[N];

  int isGround = ((ground == o1) || (ground == o2));

  int n = dCollide(o1, o2, N, &contact[0].geom, sizeof(dContact));

  if (isGround) {
    if (n >= 1)
      flag = 1;
    else
      flag = 0;
    for (int i = 0; i < n; i++) {
      contact[i].surface.mode = dContactBounce;
      contact[i].surface.mu = dInfinity;
      contact[i].surface.bounce = 0.9;  // (0.0~1.0) restitution parameter
      contact[i].surface.bounce_vel =
          0.0;  // minimum incoming velocity for bounce
      dJointID c = dJointCreateContact(world, contactgroup, &contact[i]);
      dJointAttach(c, dGeomGetBody(contact[i].geom.g1),
                   dGeomGetBody(contact[i].geom.g2));
    }
  }
}

static void simLoop(int pause) {
  const dReal *pos, *R;
  static unsigned int sc;
  // printf("Simloop: %d\n", sc++);
  std::cerr << "Simloop: " << sc++ << " \n";

  flag = 0;
  dSpaceCollide(space, 0, &nearCallback);

#if defined(USECUDA)
  cuda_copyWorldBodiesToDevice(ball_cuda.body, world);
  // cuda_copyBodiesToDevice(ball_cuda.body, &ball.body, 1);
  cuda_dxProcessIslands(world, ball_cuda.body, ball_cuda.joint, 0.1, NULL);
  cuda_copyWorldBodiesFromDevice(world, ball_cuda.body, world->nb,
                                 ball_cuda.bbuf);
  // cuda_copyBodiesFromDevice2(world, ball_cuda.body, 1, ball.body);
#else
  dWorldStep(world, 0.1);
#endif

  dJointGroupEmpty(contactgroup);

  if (flag == 0)
    dsSetColor(1.0, 0.0, 0.0);
  else
    dsSetColor(0.0, 0.0, 1.0);
  pos = dBodyGetPosition(ball.body);
  R = dBodyGetRotation(ball.body);

  std::cerr << "ball pos: " << pos[0] << ", " << pos[1] << ", " << pos[2]
            << "\n";
  dsDrawSphere(pos, R, radius);
}

void start() {
  static float xyz[3] = {0.0, -3.0, 1.0};
  static float hpr[3] = {90.0, 0.0, 0.0};
  dsSetViewpoint(xyz, hpr);

#if defined(USECUDA)
  ball_cuda.body = cuda_initBodiesOnDevice(1);
  ball_cuda.joint = cuda_initJointsOnDevice(1);

  ball_cuda.bbuf = (dBodyID)malloc(sizeof(dxBody) * 1);
  ball_cuda.jbuf = (dJointID)malloc(sizeof(dxJoint) * 1);
#endif
}

void prepDrawStuff() {
  fn.version = DS_VERSION;
  fn.start = &start;
  fn.step = &simLoop;
  fn.command = NULL;
  fn.stop = NULL;
  fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;
}

int main(int argc, char *argv[]) {
  dReal x0 = 0.0, y0 = 0.0, z0 = 2.0;
  dMass m1;

  prepDrawStuff();

  dInitODE();
  world = dWorldCreate();
  space = dHashSpaceCreate(0);
  contactgroup = dJointGroupCreate(0);

  dWorldSetGravity(world, 0, 0, -0.5);

  // Create a ground
  ground = dCreatePlane(space, 0, 0, 1, 0);

  // Create a ball
  ball.body = dBodyCreate(world);
  dMassSetZero(&m1);
  dMassSetSphereTotal(&m1, mass, radius);
  dBodySetMass(ball.body, &m1);
  dBodySetPosition(ball.body, x0, y0, z0);

  ball.geom = dCreateSphere(space, radius);
  dGeomSetBody(ball.geom, ball.body);

  dsSimulationLoop(argc, argv, 352, 288, &fn);

  dWorldDestroy(world);
  dCloseODE();

  return 0;
}
