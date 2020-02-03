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

		//helper for back projecting from an Laue detector -> sphere
		template <typename Real>
		struct BackProjector : public emsphinx::BackProjector<Real> {
			GeomParam             geo ;//geometry
			const size_t          dim ;//side length of spherical grid
			image::Rescaler<Real> sclr;//pattern rescaler

			//TODO: add in bandpass filtering parameters
			BackProjector(GeomParam g, const size_t d) : geo(g), dim(d), sclr(geo.Ny, geo.Nz, 1.0, fft::flag::Plan::Patient) {}

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
			// if(std::signbit(n[2])) return false;//can't project through sample
/*
L sample = detector distance
kk = wave number

[X,Y] --> [x,0,z]
        phi = datan2(X, Y) - pi/2
        r = sqrt(X*X+Y*Y)
        r2 = r*r
        p = sqrt((kk+L)**2+r2)
        q = 1.0/sqrt(2.0*p**2-2.0*(kk+L)*p)
        x = (kk + L - p) * q
        y = 0.0
        z = r * q

[x,0,z] --> [X,Y]
	
		q = z / r
		p = (kk + L)
*/







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

		//@brief    : unproject a pattern from the detector to a square grid
		//@param pat: pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
		//@param iq : location to write pattern image quality (NULL to skip computation)
		template <typename Real>
		Real BackProjector<Real>::unproject(Real * const pat, Real * const sph, Real * const iq) {
			//rescale pattern
			// Real vIq = sclr.scale(pat, rPat.data(), sWrk.data(), true, 0, NULL != iq);//rescale (TODO this is where bandpass could be added in)
			// if(NULL != iq) *iq = vIq;//save iq value if needed

			//zero out sphere
			std::fill(sph, sph + dim * dim * 2, Real(0));

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

					//finally distribute intensity across square image
					const Real vPat = pat[j * geo.Ny + i];
					lineWeight.bilinearCoeff(X0, Y0, X1, Y1, dim, dim);
					for(const std::pair< size_t, Real>& p : lineWeight.pairs) sph[p.first] += p.second * vPat;
				}
			}

			/*
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
			*/

		}

	/*
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
