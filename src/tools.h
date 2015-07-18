/**
 * @file 	tools.h
 * @brief 	Tools for creating distributions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef TOOLS_H
#define TOOLS_H
#include "particle.h"
struct reb_simulation;



/**
 * This function sets up a Plummer sphere.
 * @param _N Number of particles in the plummer sphere.
 * @param M Total mass of the cluster.
 * @param R Characteristic radius of the cluster.
 */

void reb_tools_init_plummer(struct reb_simulation* r, int _N, double M, double R);

/**
 * Initialize a particle on an orbit in the xy plane.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param omega Pericenter of the particle.
 * @param f true anomaly of the particle.
 */
struct reb_particle reb_tools_init_orbit2d(double G, double M, double m, double a, double e, double omega, double f);

/**
 * Initialize a particle on a 3D orbit.  See Fig. 2.13 of Murray & Dermott Solar System Dynamics for diagram.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param i inclination of the particle to the reference plane.
 * @param Omega Longitude of the ascending node of the particle.
 * @param omega argument of pericenter of the particle.
 * @param f true anomaly of the particle.
 */

struct reb_particle reb_tools_init_orbit3d(double G, double M, double m, double a, double e, double i, double Omega, double omega, double f);

/* 
 * Init the MEGNO particles
 **/
void reb_tools_megno_init(struct reb_simulation* const r, double delta);

/*
 * Returns the current value of <Y>
 **/
double reb_tools_megno(struct reb_simulation* r);

/*
 * Returns the largest Lyapunov characteristic number (LCN), or maximal Lyapunov exponent
 **/
double reb_tools_lyapunov(struct reb_simulation* r);

/*
 * Returns deltad/delta (Note, there is a typo in Gozdziewski et al 2001).
 **/

double reb_tools_megno_deltad_delta(struct reb_simulation* const r);

/*
 * Update MEGNO after a successful timestep by adding dY (=ddelta/delta*dt)
 **/
void reb_tools_megno_update(struct reb_simulation* r, double dY);


/**
 * Init random number generator based on time and process id.
 */
void reb_tools_init_srand();

#endif 	// TOOLS_H
