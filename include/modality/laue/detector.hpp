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
			struct Constants;
			const std::shared_ptr<const Constants> pLut;//read only values (can be shared across threads)

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

		};

		//helper struct to hold read only constants needed by BackProjector (for sharing across threads)
		template <typename Real>
		struct BackProjector<Real>::Constants {
			const bool                        flip;//do input patterns need to be vertically flipped
			GeomParam                         geo ;
			size_t                            dim ;//side length of spherical grid
			std::vector<Real>                 zLat;//ring lattitudes
			std::vector<Real>                 sNrm;//legendre normals

			//@brief    : construct a back projector
			//@param geo: detector geometry to back project from
			//@param dim: side length of square legendre grid to back project onto
			//@param fct: additional scale factor for detector resizing
			//@param bp : bandpass parameters
			//@param qu : quaternion to rotate detector location by
			Constants(GeomParam geo, const size_t dim, const Real fct, uint16_t const * bp, Real const * const qu = NULL);
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

		//@brief    : construct a back projector
		//@param geo: detector geometry to back project from
		//@param dim: side length of square legendre grid to back project onto
		//@param fct: additional scale factor for detector resizing
		//@param bp : bandpass parameters
		//@param qu : quaternion to rotate detector location by
		template <typename Real>
		BackProjector<Real>::Constants::Constants(GeomParam geo, const size_t dim, const Real fct, uint16_t const * bp, Real const * const qu) :
			geo(geo),
			flip(false),
			dim(dim)
		{
			//build zLat and sNrm
			zLat = square::cosLats<Real>(dim, square::Layout::Legendre);
			sNrm = square::normals<Real>(dim, square::Layout::Legendre);
		}

		//@brief  : construct a back projector
		//@param g: detector geometry to back project from
		//@param d: side length of square legendre grid to back project onto
		//@param f: additional scale factor for detector resizing
		//@param b: bandpass parameters 
		template <typename Real>
		BackProjector<Real>::BackProjector(GeomParam g, const size_t d, const Real f, uint16_t const * b) :
			pLut(std::make_shared<const Constants>(g, d, f, b, (Real const * const)NULL)) {}

		//@brief    : unproject a pattern from the detector to a square grid
		//@param pat: pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
		//@param iq : location to write pattern image quality (NULL to skip computation)
		template <typename Real>
		Real BackProjector<Real>::unproject(Real * const pat, Real * const sph, Real * const iq) {
			const Real uTru = pLut->geo.Dy + pLut->geo.ps * pLut->geo.Ny / 2;
			const Real vTru = pLut->geo.Dz + pLut->geo.ps * pLut->geo.Nz / 2;

			std::fill(sph, sph + pLut->dim * pLut->dim, Real(0));//zero out north hemisphere

			//loop over interior pixels
			for(size_t j = 1; j < pLut->geo.Nz - 1; j++) {
				Real const * row[3] = {
					pat + (j-1) * pLut->geo.Ny,//previous row
					pat +  j    * pLut->geo.Ny,//current row
					pat + (j+1) * pLut->geo.Ny,//next row
				};

				for(size_t i = 1; i < pLut->geo.Ny - 1; i++) {
					//get pixel value
					const Real v = row[1][i];
					if(v <= Real(0)) continue;//we're only interested in positive intensity

					//note: currently for flat top peaks every pixel will be found
					//correcting this requires a second pass and a working space
					if(v > row[1][i-1] &&//greater than left
					   v > row[1][i+1] &&//greater than right
					   v > row[0][i  ] &&//greater than above
					   v > row[2][i  ] ) {//greater than above
						//pixel is a local maxima, fit a peak
						//f(x,y) = A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) +c(y-y0)^2 ) )
						//a = cos^2(theta) / (2 sigmax^2) + sin^2(theta) / (2 sigmay^2)
						//b =-sin(2*theta) / (4 sigmax^2) + sin(2*theta) / (4 sigmay^2)
						//c = sin^2(theta) / (2 sigmax^2) + cos^2(theta) / (2 sigmay^2)
						Real params[6] = {
							Real(i),//x0
							Real(j),//y0
							1,//sigmax
							1,//sigmay
							0,//theta
							v,//amplitude
						};

						//TODO: refine gaussian peak estimate

						//TODO: skip peaks with bad fits (threshold on amplitude, and/or standard deviations)

						//back project peak center to spherical direction
						Real xyz[3] = {
							  pLut->geo.sampletodetector,
							-(pLut->geo.ps * params[0] - uTru),
							  pLut->geo.ps * params[1] - vTru ,
						};

						Real mag = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
						for(size_t k = 0; k < 3; k++) xyz[k] /= mag; // Ku
						xyz[0] -= 1;//subtract beam direction
						mag = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
						for(size_t k = 0; k < 3; k++) xyz[k] /= mag; // renormalize

						// rotate 90 @ y to put ring around equator
						Real qq[3] = {xyz[2], xyz[1], -xyz[0]};//to back project around equator
						// Real qq[3] = {xyz[0], xyz[1], xyz[2]};//to back project around beam direction

						//find spherical grid point closest to pack projected direction
						size_t inds[4];
						square::legendre::boundingInds(pLut->dim, pLut->zLat.data(), qq, inds);

						//compute dot product of bounding indices
						Real dp[4] = {0};
						for(size_t k = 0; k < 4; k++) dp[k] = std::inner_product(qq, qq + 3, pLut->sNrm.begin() +3*inds[k], Real(0));
						for(size_t k = 0; k < 4; k++) sph[inds[k]] += dp[k] * params[5];

						//distribute intensity around 3x3 pixel grid using von mises fisher distribution

					}
				}
			}

			//fill in southern hemisphere with inversion symmetry
			Real * const sh = sph + pLut->dim * pLut->dim;
			for(size_t j = 0; j < pLut->dim; j++) {
				for(size_t i = 0; i < pLut->dim; i++) {
					sh[(pLut->dim - 1 - j) * pLut->dim + (pLut->dim - 1 - i)] = sph[pLut->dim * j + i];
				}
			}
			return 0;
		}

	}//laue

}//emsphinx

#endif//_LAUE_DETECTOR_H_
