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
		/*
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
		*/

			//make (unweighted) standard deviation 1
			Real stdev = std::sqrt( std::inner_product(sph, sph + dim * dim, sph, Real(0)) ) / dim;
			std::for_each(sph, sph + dim * dim,
				[stdev](Real& v){
					v /= stdev;
				}
			);

			//copy to southern hemisphere with 
			for(size_t j = 0; j < dim; j++) {
				for(size_t i = 0; i < dim; i++) {
					sph[dim * dim + (dim - 1 - j) * dim + (dim - 1 - i)] = sph[dim * j + i];
				}
			}

			// static bool once = true;
			// std::ofstream os("back.raw", std::ios::out | std::ios::binary);
			// os.write((char*)sph, 2 * dim * dim * sizeof(Real));
			return 1;
		}

		//@brief: do actual backproject of an appropriately sized pattern
		//@param pat: rescaled pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		template <typename Real>
		void BackProjector<Real>:: doUnproject(Real const * const pat, Real * const sph) {
			//zero out sphere
			std::fill(sph, sph + dim * dim, Real(0));

			//define some constants for consistency with fortran code
			const double kk1 = geo.VoltageL / 1.23984193;// this is mislabeled in EMsoft
			const double kk2 = geo.VoltageH / 1.23984193;// this is mislabeled in EMsoft
			const double delta = geo.ps;
			const double L = geo.sampletodetector;
			const double kl1 = kk1 + L;
			const double kl2 = kk2 + L;
			xtal::Quat<double> yquat(0.7071067811865475, 0, -0.7071067811865475, 0);

			//now loop over pattern back projecting
			image::LineWeight<double> lineWeight;
			for(size_t iy = 0; iy < geo.Nz; iy++) {//loop over rows
				const double py = (double(iy) - double(geo.Nz)/2) * delta;
				for(size_t ix = 0; ix < geo.Ny; ix++) {//loop over cols
					const double px = (double(ix) - double(geo.Ny)/2) * delta;

					const Real vPat = pat[iy * geo.Ny + ix];
					if(vPat < 1e-2) continue; //don't bother back projecting empty pixels

					// get the azimuthal angle phi from x and y
					const double phi = std::atan2(px, py) - 1.570796326794897;
					const double cPhi =  std::cos(phi / 2);
					const double sPhi = -std::sin(phi / 2);
					xtal::Quat<double> quat(cPhi, sPhi, 0, 0);
					const double r2 = px * px + py * py;
					const double r = std::sqrt(r2);
					const double p1 = std::sqrt( kl1 * kl1 + r2 );
					const double p2 = std::sqrt( kl2 * kl2 + r2 );
					const double q1 = 1.0 / std::sqrt( ( p1 - kl1 ) * p1 * 2 );
					const double q2 = 1.0 / std::sqrt( ( p2 - kl2 ) * p2 * 2 );

					//compute sphere direction for kk0 and kk1
					double n1[3] = {(kl1 - p1) * q1, 0, r * q1};
					double n2[3] = {(kl2 - p2) * q2, 0, r * q2};

					// these are the normals in the azimuthal plane (x,z); next we need to apply the rotation by phi 
					// around x to bring the vector into the correct location
					// also rotate these unit vectors by 90° around the y-axis so that they fall in along the equator
					double rn1[3], rn2[3];
					quat.rotateVector(n1, rn1);
					quat.rotateVector(n2, rn2);
					yquat.rotateVector(rn1, rn1);
					yquat.rotateVector(rn2, rn2);

					if(std::signbit(rn1[2])) {
						for(size_t j = 0; j < 3; j++) rn1[j] = -rn1[j];
					}
					if(std::signbit(rn2[2])) {
						for(size_t j = 0; j < 3; j++) rn2[j] = -rn2[j];
					}

					// and project both points onto the Lambert square
					double X1, Y1, X2, Y2;
					square::lambert::sphereToSquare(rn1[0], rn1[1], rn1[2], X1, Y1);
					square::lambert::sphereToSquare(rn2[0], rn2[1], rn2[2], X2, Y2);

					// finally distribute intensity across square image
					lineWeight.bilinearCoeff(X1, Y1, X2, Y2, dim, dim);
					for(const std::pair< size_t, Real>& p : lineWeight.pairs) sph[p.first] += p.second * vPat;
					
					//use this to do bounds back projection like the Fortran version
					// const size_t idx1 = std::round(Y1 * (dim-1)) * dim + std::round(X1 * (dim-1));
					// const size_t idx2 = std::round(Y2 * (dim-1)) * dim + std::round(X2 * (dim-1));
					// sph[idx1] += vPat;
					// sph[idx2] += vPat;
				}
			}
		}

	}//laue

}//emsphinx

#endif//_LAUE_DETECTOR_H_
