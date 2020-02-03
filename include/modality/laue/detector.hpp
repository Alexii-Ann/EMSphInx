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

		//@brief: helper for describing the geometry of a detector
		template <typename Real>
		struct Geometry {
			GeomParam geo;

			//@brief    : bilinearlly interpolate a pixel value for the given direction
			//@param n  : unit direction to interpolate value for in the sample reference frame (pattern center is in +z direction)
			//@param pix: [optional] location to store interpolation pixel
			//@return   : true / false if the direction lies inside / outside the detector
			bool interpolatePixel(Real const * const n, image::BiPix<Real> * const pix = NULL) const;

			//@brief        : approximate the solid angle captured by a detector
			//@param gridRes: side length of square lambert projection to grid over
			//@return       : fractional solid angle i.e [0,1] not [0,4*pi]
			Real solidAngle(size_t gridRes) const;

			//@brief      : create a new detector geometry by rescaling the pixels size while maintaining solid angle
			//@param scale: rescaling factor, must lie in [0, min(w, h)] with values < 1 corresponding to increasing pixel density
			//@note       : this is analogous to continuous camera binning with e.g. scale = 4.0 corresponding to 4x4 camera binning
			Geometry<Real> rescale(const Real scale) const;

			//@brief    : compute the rescale factor required to make average detector pixel size the same as an average pixel on a square spherical grid
			//@param dim: side length of square spherical grid to match pixel size to
			//@return   : geometry scale factor to make average pixel sizes the same
			Real scaleFactor(const size_t dim) const;

		};

		/*
		//helper for back projecting from an Laue detector -> sphere
		template <typename Real>
		struct BackProjector : public emsphinx::BackProjector<Real> {
			struct Constants;//helper struct to hold read only constants
			const std::shared_ptr<const Constants> pLut;//read only values (can be shared across threads)
			std::vector<Real>                      rPat;//space for pattern as real
			fft::vector<Real>                      sWrk;//work space for pattern rescaler
			std::vector<Real>                      iVal;//space for intermediate interpolation result

			//@brief    : construct a back projector
			//@param geo: detector geometry to back project from
			//@param dim: side length of square legendre grid to back project onto
			//@param fct: additional scale factor for detector resizing
			//@param qu : quaternion to rotate detector location by
			BackProjector(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu = NULL);

			//@brief    : unproject a pattern from the detector to a square grid
			//@param pat: pattern to back project to unit sphere
			//@param sph: location to write back projected pattern
			//@param iq : location to write pattern image quality (NULL to skip computation)
			//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
			Real unproject(Real * const pat, Real * const sph, Real * const iq = NULL);

			//@brief : get a copy of the stored back projector
			//@return: unique pointer to copy of current projector
			std::unique_ptr<emsphinx::BackProjector<Real> > clone() const {return std::unique_ptr<BackProjector>(new BackProjector(*this));}

			//@brief    : build a mask of spherical pixels covered by the detector
			//@param sph: location to write mask
			void mask(Real * const sph);
		};


		//helper struct to hold read only constants needed by BackProjector (for sharing across threads)
		template <typename Real>
		struct BackProjector<Real>::Constants {
			std::vector< image::BiPix<Real> > iPts;//interpolation coefficients
			std::vector<Real>                 omeg;//solid angles of iPts (relative to average spherical pixel size)
			Real                              omgW;//sum of omeg (in window)
			Real                              omgS;//sum of omeg (over sphere)
			image::Rescaler<Real>             sclr;//pattern rescaler

			//@brief    : construct a back projector
			//@param geo: detector geometry to back project from
			//@param dim: side length of square legendre grid to back project onto
			//@param fct: additional scale factor for detector resizing
			//@param qu : quaternion to rotate detector location by
			Constants(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu = NULL);
		};
		*/

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
		//                              Geometry                              //
		////////////////////////////////////////////////////////////////////////

		//@brief    : bilinearlly interpolate a pixel value for the given direction
		//@param n  : unit direction to interpolate value for in the sample reference frame (pattern center is in +z direction)
		//@param pix: [optional] location to store interpolation pixel
		//@return   : true / false if the direction lies inside / outside the detector
		template <typename Real>
		bool Geometry<Real>::interpolatePixel(Real const * const n, image::BiPix<Real> * const pix) const {
			//sanity check
			if(std::signbit(n[2])) return false;//can't project through sample

			/*
			//compute angle between detector and sample normal
			static const Real degToRad = emsphinx::Constants<Real>::pi / 180;
			const Real alpha = (Real(90) - sTlt + dTlt) * degToRad;//angle between sample and detector is 90 with no tilts, increasing sample tilt decreases the angle, increasing camera tilt increases the angle
			if(std::fabs(alpha) > emsphinx::Constants<Real>::pi_2) throw std::runtime_error("pattern center not on detector");//beyond +/-90 the detector is behind the sample

			//compute sin / cos of angle and distance to detector
			const Real sA = std::sin(alpha);
			const Real cA = std::sqrt(Real(1) - sA * sA);//cos(alpha)
			const Real d =  sDst / (n[0] * sA + n[2] * cA);
			if(std::signbit(d)) return false;//negative distance

			//compute offset from pattern center in microns
			const Real x =                   n[1]  * d;
			const Real y = (sA * n[2] - cA * n[0]) * d;

			//convert from microns to fractional position (0->100%) and check if we're inside the detector
			const Real X = (cX + x / pX) / w + Real(0.5);
			      Real Y = (cY + y / pY) / h + Real(0.5);

			// check if we're inside the detector
			if(std::signbit(X) || std::signbit(Y) || X > 1 || Y > 1) return false;//outside of frame

			//check against circular mask if needed
			if(circ) {
				const Real dX = (X - Real(0.5)) * w;//horizontal distance from center in pixels
				const Real dY = (Y - Real(0.5)) * h;//vertical distance from center in pixels
				const Real r = Real(std::min(w, h)) / 2;//radius of circular mask
				if(r * r < dX * dX + dY * dY) return false;
			}

			//if we made it this far we're in frame, compute relative contribution from neighboring pixels
			if(pix != NULL) pix->bilinearCoeff(X, Y, w, h);
			return true;
			*/
			return false;
		}

		//@brief        : approximate the solid angle captured by a detector
		//@param gridRes: side length of square lambert projection to grid over
		//@return       : fractional solid angle i.e [0,1] not [0,4*pi]
		template <typename Real>
		Real Geometry<Real>::solidAngle(size_t gridRes) const {
			//loop over northern hemisphere
			Real XY[2], xyz[3];
			size_t goodPoints = 0;
			for(size_t j = 0; j <= gridRes; j++) {//loop over rows of square projection
				XY[1] = Real(j) / gridRes;
				for(size_t i = 0; i <= gridRes; i++) {//loop over columns of square projection
					XY[0] = Real(i) / gridRes;
					square::lambert::squareToSphere(XY[0], XY[1], xyz[0], xyz[1], xyz[2]);
					if(interpolatePixel(xyz)) ++goodPoints;//count this point if it falls on the detector
				}
			}
			const size_t totalPoints = (gridRes * gridRes + (gridRes - 2) * (gridRes - 2));//total number of points in both hemispheres (don't double count equators)
			return Real(goodPoints) / totalPoints;//count fraction of grid points on detector
		}

		//@brief      : create a new detector geometry by rescaling the pixels size while maintaining solid angle
		//@param scale: rescaling factor, must lie in [0, min(w, h)] with values < 1 corresponding to increasing pixel density
		//@note       : this is analogous to continuous camera binning with e.g. scale = 4.0 corresponding to 4x4 camera binning
		template <typename Real>
		Geometry<Real> Geometry<Real>::rescale(const Real scale) const {
			const size_t wNew = (size_t) std::round(Real(geo.Ny) / scale);
			const size_t hNew = (size_t) std::round(Real(geo.Nz) / scale);
			if(wNew == 0 || hNew == 0) throw std::runtime_error("cannot rescale detector to less than 1 pixel");
			
			//copy everything and update size
			Geometry<Real> geom(*this);
			geom.geo.Ny = wNew;
			geom.geo.Nz = hNew;
			geom.ps *= scale;
			return geom;
		}

		//@brief    : compute the rescale factor required to make average detector pixel size the same as an average pixel on a square spherical grid
		//@param dim: side length of square spherical grid to match pixel size to
		//@return   : geometry scale factor to make average pixel sizes the same
		template <typename Real>
		Real Geometry<Real>::scaleFactor(const size_t dim) const {
			const size_t sqrPix = dim * dim * 2 - (dim - 1) * 4;//number of pixels in square grid
			const size_t detPix = geo.Ny * geo.Nz;//number of pixels on detector
			return std::sqrt(solidAngle(501) * sqrPix / detPix);
		}

		////////////////////////////////////////////////////////////////////////
		//                           BackProjector                            //
		////////////////////////////////////////////////////////////////////////
/*
		template <typename Real>
		BackProjector<Real>::Constants::Constants(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu) :
			sclr(geo.w, geo.h, geo.scaleFactor(dim) * fct, fft::flag::Plan::Patient)//we're going to be doing lots of 
		{
			//compute normals and solid angles of square legendre grid
			std::vector<Real> xyz(dim * dim * 3);
			square::legendre::normals(dim, xyz.data());
			std::vector<Real> omegaRing = square::solidAngles<Real>(dim, square::Layout::Legendre);

			//rescale detector
			Geometry<Real> g = geo.rescale(sclr.wOut, sclr.hOut);

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
				if(g.interpolatePixel(n, &p)) {
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
				if(g.interpolatePixel(n, &p)) {
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

		//@brief    : construct a back projector
		//@param geo: detector geometry to back project from
		//@param dim: side length of square legendre grid to back project onto
		//@param fct: additional scale factor for detector resizing
		//@param qu : quaternion to rotate detector location by
		template <typename Real>
		BackProjector<Real>::BackProjector(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu) :
			pLut(std::make_shared<const Constants>(geo, dim, fct, qu)),
			rPat(std::max(geo.w * geo.h, pLut->sclr.wOut * pLut->sclr.hOut)),
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
			Real vIq = pLut->sclr.scale(pat, rPat.data(), sWrk.data(), true, 0, NULL != iq);//rescale
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
				for(size_t i = 0; i < pLut->iPts.size(); i++) sph[pIpt[i].idx] = iVal[i];

				//compute sum if needed
				static const Real var = std::sqrt(pLut->omgW / pLut->omgS * emsphinx::Constants<Real>::pi * Real(4));//standard deviation within window (since we normalized to 1)
				return var;
			}

		}

		//@brief    : build a mask of spherical pixels covered by the detector
		//@param sph: location to write mask
		template <typename Real>
		void BackProjector<Real>::mask(Real * const sph) {
			for(const image::BiPix<Real>& p: pLut->iPts) sph[p.idx] = 1;
		}
	*/

	}//laue

}//emsphinx

#endif//_LAUE_DETECTOR_H_
