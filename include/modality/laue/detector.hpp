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
			std::vector<Real>                      rPat;//space for pattern as real
			fft::vector<Real>                      sWrk;//work space for pattern rescaler
			std::vector<Real>                      iVal;//space for intermediate interpolation result

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
			std::vector< image::BiPix<Real> > iPts;//interpolation coefficients
			std::vector<Real>                 omeg   ;//solid angles of iPts (relative to average spherical pixel size)
			Real                              omgW   ;//sum of omeg (in window)
			Real                              omgS   ;//sum of omeg (over sphere)
			image::Rescaler<Real>             sclr   ;//pattern rescaler
			const bool                        flip   ;//do input patterns need to be vertically flipped
			GeomParam                         geo    ;
			size_t                            bPas[2];//bandpass cutoffs

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
			sclr(geo.Ny, geo.Nz, geo.scaleFactor(dim) * fct, fft::flag::Plan::Patient),//we're going to be doing lots of 
			geo(geo),
			flip(false)
		{
			//save band pass values
			bPas[0] = bp[0]; bPas[1] = bp[1];

			//compute normals and solid angles of square legendre grid
			std::vector<Real> xyz(dim * dim * 3);
			square::legendre::normals(dim, xyz.data());
			std::vector<Real> omegaRing = square::solidAngles<Real>(dim, square::Layout::Legendre);

			//rescale detector
			const Real sclFct = std::sqrt( Real(geo.Ny * geo.Nz) / (sclr.wOut * sclr.hOut));
			GeomParam g = geo.rescale(sclFct);

			//determine weights north hemisphere
			image::BiPix<Real> p;
			for(size_t i = 0; i < xyz.size() / 3; i++) {//loop over square legendre grid points
				//get direction and apply rotatation
				Real n[3];
				if(NULL != qu) {
					xtal::quat::rotateVector(qu, xyz.data() + 3 * i, n);
				} else {
					std::copy(xyz.data() + 3 * i, xyz.data() + 3 * i + 3, n);
				}

				//save interpolation weights if it hits the detector
				if(g.interpolatePixel(n, &p, flip)) {
					p.idx = i;
					iPts.push_back(p);
					omeg.push_back(omegaRing[square::ringNum(dim, i)]);
				}
			}

			//determine weights south hemisphere
			for(size_t i = 0; i < xyz.size() / 3; i++) {//loop over square legendre grid points
				//split into x and y indices
				const size_t y = i / dim;
				if(y == 0 || y+1 == dim) continue;//dont' double count equator
				const size_t x = i - dim * y;
				if(x == 0 || x+1 == dim) continue;//dont' double count equator
				xyz[3*i+2] = -xyz[3*i+2];//move normal to southern hemisphere

				//get direction and apply rotatation
				Real n[3];
				if(NULL != qu) {
					xtal::quat::rotateVector(qu, xyz.data() + 3 * i, n);
				} else {
					std::copy(xyz.data() + 3 * i, xyz.data() + 3 * i + 3, n);
				}

				//save interpolation weights if it hits the detector
				if(g.interpolatePixel(n, &p, flip)) {
					p.idx = i;
					iPts.push_back(p);
					omeg.push_back(omegaRing[square::ringNum(dim, i)]);
				}
			}

			//accumulate solid angle sums
			omgW = std::accumulate(omeg.cbegin(), omeg.cend(), Real(0));
			omgS = 0;
			for(size_t j = 0; j < dim; j++) {
				for(size_t i = 0; i < dim; i++) {
					const Real& o = omegaRing[square::ringNum(dim, j * dim + i)];
					const bool eq = j == 0 || i == 0 || j == dim-1 || i == dim-1;//are we on the equator
					omgS += o;
					if(!eq) omgS += o;//don't double count equator
				}
			}
		}

		//@brief  : construct a back projector
		//@param g: detector geometry to back project from
		//@param d: side length of square legendre grid to back project onto
		//@param f: additional scale factor for detector resizing
		//@param b: bandpass parameters 
		template <typename Real>
		BackProjector<Real>::BackProjector(GeomParam g, const size_t d, const Real f, uint16_t const * b) :
			pLut(std::make_shared<const Constants>(g, d, f, b, (Real const * const)NULL)),
			rPat(std::max<size_t>(g.Ny * g.Nz, pLut->sclr.wOut * pLut->sclr.hOut)),
			sWrk(pLut->sclr.allocateWork()),
			iVal(pLut->iPts.size()) {}

		//@brief    : unproject a pattern from the detector to a square grid
		//@param pat: pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
		//@param iq : location to write pattern image quality (NULL to skip computation)
		template <typename Real>
		Real BackProjector<Real>::unproject(Real * const pat, Real * const sph, Real * const iq) {
			//rescale pattern
			Real vIq = pLut->sclr.scale(pat, rPat.data(), sWrk.data(), true, 0, NULL != iq);//rescale TODO: this is where bandpass can be implemented (replace 0 with bPas[0], update scaler to allow high pass)
			if(NULL != iq) *iq = vIq;//save iq value if needed

			//interpolate pattern intensities and compute weighted mean
			image::BiPix<Real> const * const pIpt = pLut->iPts.data();
			for(size_t i = 0; i < pLut->iPts.size(); i++) iVal[i] = pIpt[i].interpolate(rPat.data());//determine interpolated values
			const Real mean = std::inner_product(iVal.cbegin(), iVal.cend(), pLut->omeg.cbegin(), Real(0)) / pLut->omgW;//compute weighted average

			//make mean zero and compute weighted standard deviation
			Real stdev = 0;
			for(size_t i = 0; i < iVal.size(); i++) {
				iVal[i] -= mean;
				stdev += iVal[i] * iVal[i] * pLut->omeg[i];
			}
			stdev = std::sqrt(stdev / pLut->omgW);

			if(Real(0) == stdev) {//the back projected image had no contrast
				//copy to output grid making value uniformly 1
				for(size_t i = 0; i < pLut->iPts.size(); i++) sph[pIpt[i].idx] = Real(1);
				return 0;
			} else {
				//make standard deviation 1
				for(size_t i = 0; i < iVal.size(); i++) iVal[i] /= stdev;

				//copy to output grid making mean 0
				// for(size_t i = 0; i < pLut->iPts.size(); i++) sph[pIpt[i].idx] = iVal[i];
				//NOTE: this is to prevent a band of negative values around the back projected spots
				//TODO: eliminate rescaling entirely by preprocessing laue pattern
				for(size_t i = 0; i < pLut->iPts.size(); i++) sph[pIpt[i].idx] = std::max(iVal[i], Real(0));

				//compute sum if needed
				static const Real var = std::sqrt(pLut->omgW / pLut->omgS * emsphinx::Constants<Real>::pi * Real(4));//standard deviation within window (since we normalized to 1)
				return var;
			}

		}

	}//laue

}//emsphinx

#endif//_LAUE_DETECTOR_H_
