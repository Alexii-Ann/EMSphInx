/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                     *
 * Copyright (c) 2019-2019, De Graef Group, Carnegie Mellon University *
 * All rights reserved.                                                *
 *                                                                     *
 * Author: William C. Lenthe                                           *
 *                                                                     *
 * This package is free software; you can redistribute it and/or       *
 * modify it under the terms of the GNU General Public License as      *
 * published by the Free Software Foundation; either version 2 of the  *
 * License, or (at your option) any later version.                     *
 *                                                                     *
 * This program is distributed in the hope that it will be useful,     *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
 * GNU General Public License for more details.                        *
 *                                                                     *
 * You should have received a copy of the GNU General Public License   *
 * along with this program; if not, check the Free Software Foundation *
 * website: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.html>   *
 *                                                                     *
 *                                                                     *
 * Interested in a commercial license? Contact:                        *
 *                                                                     *
 * Center for Technology Transfer and Enterprise Creation              *
 * 4615 Forbes Avenue, Suite 302                                       *
 * Pittsburgh, PA 15213                                                *
 *                                                                     *
 * phone. : 412.268.7393                                               *
 * email  : innovation@cmu.edu                                         *
 * website: https://www.cmu.edu/cttec/                                 *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _LAUE_DETECTOR_H_
#define _LAUE_DETECTOR_H_

#include <algorithm>
#include <memory>
#include <array>

#include "H5Cpp.h"

#include "util/fft.hpp"
#include "util/image.hpp"
#include "idx/base.hpp"
#include "nml.hpp"

namespace emsphinx {

	namespace laue {

		//helper for back projecting from an Laue detector -> sphere
		template <typename Real>
		struct BackProjector : public emsphinx::BackProjector<Real> {
			GeomParam             geo    ;//geometry
			const size_t          dim    ;//side length of spherical grid
			image::Rescaler<Real> sclr   ;//pattern rescaler
			std::vector<Real>     rPat   ;//space for pattern as real
			fft::vector<Real>     sWrk   ;//work space for pattern rescaler
			std::vector<Real>     omeg   ;//solid angle of each pixel (or 0 if there is no contribution from the detector)
			Real                  omegSm ;//sum of omeg
			size_t                bPas[2];//bandpass cutoffs

			//@brief  : construct a back projector
			//@param g: detector geometry to back project from
			//@param d: side length of square legendre grid to back project onto
			//@param f: additional scale factor for detector resizing
			//@param b: bandpass parameters
			BackProjector(GeomParam g, const size_t d, const Real f, uint16_t const * b);

			//@brief    : unproject a pattern from the detector to a square grid
			//@param pat: pattern to back project to unit sphere
			//@param sph: location to write back projected pattern
			//@param iq : location to write pattern image quality (NULL to skip computation)
			//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
			Real unproject(Real * const pat, Real * const sph, Real * const iq = NULL);

			//@brief : get a copy of the stored back projector
			//@return: unique pointer to copy of current projector
			std::unique_ptr<emsphinx::BackProjector<Real> > clone() const {return std::unique_ptr<BackProjector>(new BackProjector(*this));}

			private:

				//@brief: do actual backproject of an appropriately sized pattern
				//@param pat: rescaled pattern to back project to unit sphere
				//@param sph: location to write back projected pattern
				void doUnproject(Real const * const pat, Real * const sph);
		};

	}//laue

}//emsphinx

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "constants.hpp"
#include "sht/square_sht.hpp"
#include "sht/sht_xcorr.hpp"
#include "xtal/rotations.hpp"
#include "xtal/quaternion.hpp"

namespace emsphinx {

	namespace laue {

		////////////////////////////////////////////////////////////////////////
		//                           BackProjector                            //
		////////////////////////////////////////////////////////////////////////

		//@brief  : construct a back projector
		//@param g: detector geometry to back project from
		//@param d: side length of square legendre grid to back project onto
		//@param f: additional scale factor for detector resizing
		//@param b: bandpass parameters 
		template <typename Real>
		BackProjector<Real>::BackProjector(GeomParam g, const size_t d, const Real f, uint16_t const * b) :
			geo(g.rescale(f)),
			dim(d),
			sclr(g.Ny, g.Nz, f, fft::flag::Plan::Patient),
			rPat(std::max<size_t>(g.Ny * g.Nz, sclr.wOut * sclr.hOut)),
			sWrk(sclr.allocateWork()),
			omeg(d * d) {
				//copy bandpass cutoffs
				bPas[0] = b[0];
				bPas[1] = b[1];
				
				//build mask of which points back project onto the detector
				//this is equivalent to back projecting an image of 1
				std::vector<Real> temp(sclr.wOut * sclr.hOut, 1);
				doUnproject(temp.data(), omeg.data());//unproject onto solid angle

				//now scale omega points by solid angle of spherical pixel size
				std::vector<Real> omegaRing = square::solidAngles<Real>(dim, square::Layout::Legendre);
				for(size_t i = 0; i < omeg.size(); i++) omeg[i] *= omegaRing[square::ringNum(dim, i)];
				omegSm = std::accumulate(omeg.begin(), omeg.end(), Real(0));
			}

		//@brief    : unproject a pattern from the detector to a square grid
		//@param pat: pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
		//@param iq : location to write pattern image quality (NULL to skip computation)
		template <typename Real>
		Real BackProjector<Real>::unproject(Real * const pat, Real * const sph, Real * const iq) {
			//rescale pattern
			Real vIq = sclr.scale(pat, rPat.data(), sWrk.data(), true, bPas[0], NULL != iq);//rescale (TODO this is where bandpass could be added in)
			if(NULL != iq) *iq = vIq;//save iq value if needed

			//do the back projection
			doUnproject(rPat.data(), sph);

			//make weighted mean 0
			const Real vMean = std::inner_product(omeg.begin(), omeg.end(), sph, Real(0)) / omegSm;
			Real stdev(0);
			std::transform(omeg.begin(), omeg.end(), sph, sph,
				[vMean,&stdev](const Real o, const Real s){
					if(o > Real(0)) {
						const Real vs = s - vMean;
						stdev += vs * vs * o;
						return vs;
					} else {
						return Real(0);
					}
				}
			);
			stdev = std::sqrt(stdev / omegSm);

			//now make standard deviation 1 and compute sqrt(sumsq(intensity))
			Real var(0);
			std::transform(omeg.begin(), omeg.end(), sph, sph,
				[stdev](const Real o, const Real s){
					if(o > Real(0)) {
						const Real vs = s / stdev;
						return vs;
					} else {
						return Real(0);
					}
				}
			);
			return 1;
		}

		//@brief: do actual backproject of an appropriately sized pattern
		//@param pat: rescaled pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		template <typename Real>
		void BackProjector<Real>:: doUnproject(Real const * const pat, Real * const sph) {
			//zero out sphere
			std::fill(sph, sph + dim * dim, Real(0));

			const double L = geo.sampletodetector;
			const double kk0 = geo.VoltageL * 5067730.76;// k = E / (hbar * c), 1 keV / (hbar * c) = 5067730.76 mm^-1
			const double kk1 = geo.VoltageH * 5067730.76;// k = E / (hbar * c), 1 keV / (hbar * c) = 5067730.76 mm^-1
			const double kl0 = kk0 + L;
			const double kl1 = kk1 + L;
			xtal::Quat<double> yquat(0.7071067811865475, 0, -0.7071067811865475, 0);

			//now loop over pattern back projecting
			image::LineWeight<double> lineWeight;
			for(size_t j = 0; j < geo.Nz; j++) {//loop over rows
				const double y = double(j) - double(geo.Nz)/2 * geo.ps;
				for(size_t i = 0; i < geo.Ny; i++) {//loop over cols
					const double x = double(i) - double(geo.Ny)/2 * geo.ps;
					// get the azimuthal angle phi from x and y
					const double phi = std::atan2(y, x) - 1.570796326794897;
					const double cPhi =  std::cos(phi / 2);
					const double sPhi = -std::sin(phi / 2);
					const double r2 = x * x + y * y;
					const double r = std::sqrt(r2);
					const double p0 = std::sqrt( kl0 * kl0 + r2 );
					const double p1 = std::sqrt( kl1 * kl1 + r2 );
					const double q0 = 1.0 / std::sqrt( ( p0 - kl0 ) * p0 * 2 );
					const double q1 = 1.0 / std::sqrt( ( p1 - kl1 ) * p0 * 2 );

					//compute sphere direction for kk0 and kk1
					double n0[3] = {(kl0 - p0) * q0, 0, r * q0};
					double n1[3] = {(kl1 - p1) * q1, 0, r * q1};

					// these are the normals in the azimuthal plane (x,z); next we need to apply the rotation by phi 
					// around x to bring the vector into the correct location
					// also rotate these unit vectors by 90° around the y-axis so that they fall in along the equator
					xtal::Quat<double> quat(cPhi, sPhi, 0, 0);
					quat.rotateVector(n0, n0);
					quat.rotateVector(n1, n1);
					yquat.rotateVector(n0, n0);
					yquat.rotateVector(n1, n1);
					
					// and project both points onto the Lambert square
					double X0, Y0, X1, Y1;
					square::lambert::sphereToSquare(n0[0], n0[1], n0[2], X0, Y0);
					square::lambert::sphereToSquare(n1[0], n1[1], n1[2], X1, Y1);

					// finally distribute intensity across square image
					const Real vPat = pat[j * geo.Ny + i];
					lineWeight.bilinearCoeff(X0, Y0, X1, Y1, dim, dim);
					for(const std::pair< size_t, Real>& p : lineWeight.pairs) sph[p.first] += p.second * vPat;
				}
			}
		}

	}//laue

}//emsphinx

#endif//_LAUE_DETECTOR_H_
