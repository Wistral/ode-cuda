/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

#include "ode/ode.h"
#include "objects.h"
#include "joints/joint.h"
#include "util.h"

#define ALLOCA dALLOCA16

// Auto disabling
void dInternalHandleAutoDisabling (dxWorld *world, dReal stepsize)
{
	dxBody *bb;
	for ( bb=world->firstbody; bb; bb=(dxBody*)bb->next )
	{
		// don't freeze objects mid-air (patch 1586738)
		if ( bb->firstjoint == NULL ) continue;

		// nothing to do unless this body is currently enabled and has
		// the auto-disable flag set
		if ( (bb->flags & (dxBodyAutoDisable|dxBodyDisabled)) != dxBodyAutoDisable ) continue;

		// if sampling / threshold testing is disabled, we can never sleep.
		if ( bb->adis.average_samples == 0 ) continue;

		//
		// see if the body is idle
		//
		
#ifndef dNODEBUG
		// sanity check
		if ( bb->average_counter >= bb->adis.average_samples )
		{
			dUASSERT( bb->average_counter < bb->adis.average_samples, "buffer overflow" );

			// something is going wrong, reset the average-calculations
			bb->average_ready = 0; // not ready for average calculation
			bb->average_counter = 0; // reset the buffer index
		}
#endif // dNODEBUG

		// sample the linear and angular velocity
		bb->average_lvel_buffer[bb->average_counter][0] = bb->lvel[0];
		bb->average_lvel_buffer[bb->average_counter][1] = bb->lvel[1];
		bb->average_lvel_buffer[bb->average_counter][2] = bb->lvel[2];
		bb->average_avel_buffer[bb->average_counter][0] = bb->avel[0];
		bb->average_avel_buffer[bb->average_counter][1] = bb->avel[1];
		bb->average_avel_buffer[bb->average_counter][2] = bb->avel[2];
		bb->average_counter++;
		
		// buffer ready test
		if ( bb->average_counter >= bb->adis.average_samples )
		{
			bb->average_counter = 0; // fill the buffer from the beginning
			bb->average_ready = 1; // this body is ready now for average calculation
		}

		int idle = 0; // Assume it's in motion unless we have samples to disprove it.

		// enough samples?
		if ( bb->average_ready )
		{
			idle = 1; // Initial assumption: IDLE

			// the sample buffers are filled and ready for calculation
			dVector3 average_lvel, average_avel;

			// Store first velocity samples
			average_lvel[0] = bb->average_lvel_buffer[0][0];
			average_avel[0] = bb->average_avel_buffer[0][0];
			average_lvel[1] = bb->average_lvel_buffer[0][1];
			average_avel[1] = bb->average_avel_buffer[0][1];
			average_lvel[2] = bb->average_lvel_buffer[0][2];
			average_avel[2] = bb->average_avel_buffer[0][2];
			
			// If we're not in "instantaneous mode"
			if ( bb->adis.average_samples > 1 )
			{
				// add remaining velocities together
				for ( unsigned int i = 1; i < bb->adis.average_samples; ++i )
				{
					average_lvel[0] += bb->average_lvel_buffer[i][0];
					average_avel[0] += bb->average_avel_buffer[i][0];
					average_lvel[1] += bb->average_lvel_buffer[i][1];
					average_avel[1] += bb->average_avel_buffer[i][1];
					average_lvel[2] += bb->average_lvel_buffer[i][2];
					average_avel[2] += bb->average_avel_buffer[i][2];
				}

				// make average
				dReal r1 = dReal( 1.0 ) / dReal( bb->adis.average_samples );

				average_lvel[0] *= r1;
				average_avel[0] *= r1;
				average_lvel[1] *= r1;
				average_avel[1] *= r1;
				average_lvel[2] *= r1;
				average_avel[2] *= r1;
			}

			// threshold test
			dReal av_lspeed, av_aspeed;
			av_lspeed = dDOT( average_lvel, average_lvel );
			if ( av_lspeed > bb->adis.linear_average_threshold )
			{
				idle = 0; // average linear velocity is too high for idle
			}
			else
			{
				av_aspeed = dDOT( average_avel, average_avel );
				if ( av_aspeed > bb->adis.angular_average_threshold )
				{
					idle = 0; // average angular velocity is too high for idle
				}
			}
		}

		// if it's idle, accumulate steps and time.
		// these counters won't overflow because this code doesn't run for disabled bodies.
		if (idle) {
			bb->adis_stepsleft--;
			bb->adis_timeleft -= stepsize;
		}
		else {
			// Reset countdowns
			bb->adis_stepsleft = bb->adis.idle_steps;
			bb->adis_timeleft = bb->adis.idle_time;
		}

		// disable the body if it's idle for a long enough time
		if ( bb->adis_stepsleft <= 0 && bb->adis_timeleft <= 0 )
		{
			bb->flags |= dxBodyDisabled; // set the disable flag

			// disabling bodies should also include resetting the velocity
			// should prevent jittering in big "islands"
			bb->lvel[0] = 0;
			bb->lvel[1] = 0;
			bb->lvel[2] = 0;
			bb->avel[0] = 0;
			bb->avel[1] = 0;
			bb->avel[2] = 0;
		}
	}
}


//****************************************************************************
// body rotation

// return sin(x)/x. this has a singularity at 0 so special handling is needed
// for small arguments.
static inline dReal sinc (dReal x)
{
  // if |x| < 1e-4 then use a taylor series expansion. this two term expansion
  // is actually accurate to one LS bit within this range if double precision
  // is being used - so don't worry!
  if (dFabs(x) < 1.0e-4) return REAL(1.0) - x*x*REAL(0.166666666666666666667);
  else return dSin(x)/x;
}


// given a body b, apply its linear and angular rotation over the time
// interval h, thereby adjusting its position and orientation.
void dxStepBody (dxBody *b, dReal h)
{
  // cap the angular velocity
  if (b->flags & dxBodyMaxAngularSpeed) {
        const dReal max_ang_speed = b->max_angular_speed;
        const dReal aspeed = dDOT( b->avel, b->avel );
        if (aspeed > max_ang_speed*max_ang_speed) {
                const dReal coef = max_ang_speed/dSqrt(aspeed);
                dOPEC(b->avel, *=, coef);
        }
  }
  // end of angular velocity cap


  int j;

  // handle linear velocity
  for (j=0; j<3; j++) b->posr.pos[j] += h * b->lvel[j];

  if (b->flags & dxBodyFlagFiniteRotation) {
    dVector3 irv;	// infitesimal rotation vector
    dQuaternion q;	// quaternion for finite rotation

    if (b->flags & dxBodyFlagFiniteRotationAxis) {
      // split the angular velocity vector into a component along the finite
      // rotation axis, and a component orthogonal to it.
      dVector3 frv;		// finite rotation vector
      dReal k = dDOT (b->finite_rot_axis,b->avel);
      frv[0] = b->finite_rot_axis[0] * k;
      frv[1] = b->finite_rot_axis[1] * k;
      frv[2] = b->finite_rot_axis[2] * k;
      irv[0] = b->avel[0] - frv[0];
      irv[1] = b->avel[1] - frv[1];
      irv[2] = b->avel[2] - frv[2];

      // make a rotation quaternion q that corresponds to frv * h.
      // compare this with the full-finite-rotation case below.
      h *= REAL(0.5);
      dReal theta = k * h;
      q[0] = dCos(theta);
      dReal s = sinc(theta) * h;
      q[1] = frv[0] * s;
      q[2] = frv[1] * s;
      q[3] = frv[2] * s;
    }
    else {
      // make a rotation quaternion q that corresponds to w * h
      dReal wlen = dSqrt (b->avel[0]*b->avel[0] + b->avel[1]*b->avel[1] +
			  b->avel[2]*b->avel[2]);
      h *= REAL(0.5);
      dReal theta = wlen * h;
      q[0] = dCos(theta);
      dReal s = sinc(theta) * h;
      q[1] = b->avel[0] * s;
      q[2] = b->avel[1] * s;
      q[3] = b->avel[2] * s;
    }

    // do the finite rotation
    dQuaternion q2;
    dQMultiply0 (q2,q,b->q);
    for (j=0; j<4; j++) b->q[j] = q2[j];

    // do the infitesimal rotation if required
    if (b->flags & dxBodyFlagFiniteRotationAxis) {
      dReal dq[4];
      dWtoDQ (irv,b->q,dq);
      for (j=0; j<4; j++) b->q[j] += h * dq[j];
    }
  }
  else {
    // the normal way - do an infitesimal rotation
    dReal dq[4];
    dWtoDQ (b->avel,b->q,dq);
    for (j=0; j<4; j++) b->q[j] += h * dq[j];
  }

  // normalize the quaternion and convert it to a rotation matrix
  dNormalize4 (b->q);
  dQtoR (b->q,b->posr.R);

  // notify all attached geoms that this body has moved
  for (dxGeom *geom = b->geom; geom; geom = dGeomGetBodyNext (geom))
    dGeomMoved (geom);

  // notify the user
  if (b->moved_callback)
    b->moved_callback(b);


  // damping
  if (b->flags & dxBodyLinearDamping) {
        const dReal lin_threshold = b->dampingp.linear_threshold;
        const dReal lin_speed = dDOT( b->lvel, b->lvel );
        if ( lin_speed > lin_threshold) {
                const dReal k = 1 - b->dampingp.linear_scale;
                dOPEC(b->lvel, *=, k);
        }
  }
  if (b->flags & dxBodyAngularDamping) {
        const dReal ang_threshold = b->dampingp.angular_threshold;
        const dReal ang_speed = dDOT( b->avel, b->avel );
        if ( ang_speed > ang_threshold) {
                const dReal k = 1 - b->dampingp.angular_scale;
                dOPEC(b->avel, *=, k);
        }
  }

}

//****************************************************************************
// island processing

// this groups all joints and bodies in a world into islands. all objects
// in an island are reachable by going through connected bodies and joints.
// each island can be simulated separately.
// note that joints that are not attached to anything will not be included
// in any island, an so they do not affect the simulation.
//
// this function starts new island from unvisited bodies. however, it will
// never start a new islands from a disabled body. thus islands of disabled
// bodies will not be included in the simulation. disabled bodies are
// re-enabled if they are found to be part of an active island.
void dxProcessIslands (dxWorld *world, dReal stepsize, dstepper_fn_t stepper)
{
  dxBody *b,*bb,**body;
  dxJoint *j,**joint;

  // nothing to do if no bodies
  if (world->nb <= 0) return;

  // handle auto-disabling of bodies
  dInternalHandleAutoDisabling (world,stepsize);

  // make arrays for body and joint lists (for a single island) to go into
  body = (dxBody**) ALLOCA (world->nb * sizeof(dxBody*));
  joint = (dxJoint**) ALLOCA (world->nj * sizeof(dxJoint*));
  int bcount = 0;	// number of bodies in `body'
  int jcount = 0;	// number of joints in `joint'

  // set all body/joint tags to 0
  for (b=world->firstbody; b; b=(dxBody*)b->next) b->tag = 0;
  for (j=world->firstjoint; j; j=(dxJoint*)j->next) j->tag = 0;

  // allocate a stack of unvisited bodies in the island. the maximum size of
  // the stack can be the lesser of the number of bodies or joints, because
  // new bodies are only ever added to the stack by going through untagged
  // joints. all the bodies in the stack must be tagged!
  int stackalloc = (world->nj < world->nb) ? world->nj : world->nb;
  dxBody **stack = (dxBody**) ALLOCA (stackalloc * sizeof(dxBody*));

  for (bb=world->firstbody; bb; bb=(dxBody*)bb->next) {
    // get bb = the next enabled, untagged body, and tag it
    if (bb->tag || (bb->flags & dxBodyDisabled)) continue;
    bb->tag = 1;

    // tag all bodies and joints starting from bb.
    int stacksize = 0;
    b = bb;
    body[0] = bb;
    bcount = 1;
    jcount = 0;

    // traverse and tag all body's joints, add untagged connected bodies to
    // stack
    auto traverse_joints = [&]() {
      for (dxJointNode *n = b->firstjoint; n; n = n->next) {
        if (!n->joint->tag && n->joint->isEnabled()) {
          n->joint->tag = 1;
          joint[jcount++] = n->joint;
          if (n->body && !n->body->tag) {
            n->body->tag = 1;
            stack[stacksize++] = n->body;
          }
        }
      }
    };

    traverse_joints();
    while (stacksize > 0) {
      b = stack[--stacksize];	// pop body off stack
      body[bcount++] = b;	// put body on body list

      traverse_joints();

      dIASSERT(stacksize <= world->nb);
      dIASSERT(stacksize <= world->nj);
    }

    // now do something with body and joint lists
    stepper (world,body,bcount,joint,jcount,stepsize);

    // what we've just done may have altered the body/joint tag values.
    // we must make sure that these tags are nonzero.
    // also make sure all bodies are in the enabled state.
    int i;
    for (i=0; i<bcount; i++) {
      body[i]->tag = 1;
      body[i]->flags &= ~dxBodyDisabled;
    }
    for (i=0; i<jcount; i++) joint[i]->tag = 1;
  }

  // if debugging, check that all objects (except for disabled bodies,
  // unconnected joints, and joints that are connected to disabled bodies)
  // were tagged.
# ifndef dNODEBUG
  for (b=world->firstbody; b; b=(dxBody*)b->next) {
    if (b->flags & dxBodyDisabled) {
      if (b->tag) dDebug (0,"disabled body tagged");
    }
    else {
      if (!b->tag) dDebug (0,"enabled body not tagged");
    }
  }
  for (j=world->firstjoint; j; j=(dxJoint*)j->next) {
    if ( (( j->node[0].body && (j->node[0].body->flags & dxBodyDisabled)==0 ) ||
          (j->node[1].body && (j->node[1].body->flags & dxBodyDisabled)==0) )
         && 
         j->isEnabled() ) {
      if (!j->tag) dDebug (0,"attached enabled joint not tagged");
    }
    else {
      if (j->tag) dDebug (0,"unattached or disabled joint tagged");
    }
  }
# endif
}



