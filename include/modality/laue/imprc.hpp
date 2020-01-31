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

#ifndef _LAUE_IM_PRC_
#define _LAUE_IM_PRC_

#include "util/gaussian.hpp"
#include "util/ahe.hpp"

namespace emsphinx {

	namespace laue {

		//@brief: class to hold pattern processing details
		template <typename Real>
		class PatternProcessor : public ImageProcessor<Real> {
			size_t nPix ;//Nnumber of pixels

			public:
				//@brief: set processing parameters
				//@param w: image width in pixels
				//@param h: image height in pixels
				void setSize(const size_t w, const size_t h);

				//@brief    : process an image out of place (background subtract and/or AHE)
				//@param im : image to process
				//@param buf: location to write floating point processed image
				void process(uint8_t  const * const im, Real * const buf);
				void process(uint16_t const * const im, Real * const buf);
				void process(float    const * const im, Real * const buf);

				//@brief : get size of target image to process
				//@return: size of input image in pixels
				size_t numPix() const {return nPix;}

				//@brief : get a copy of the stored image processor
				//@return: unique pointer to copy of current processor
				std::unique_ptr<ImageProcessor<Real> > clone() const {return std::unique_ptr<PatternProcessor>(new PatternProcessor(*this));}
		};

	}//laue

}//emsphinx

namespace emsphinx {

	namespace laue {
    
        //@brief: set processing parameters
				//@param w: image width in pixels
				//@param h: image height in pixels
				template <typename Real> 
        void setSize(const size_t w, const size_t h){
          nPix = w * h;
        }
    
        //@brief    : process an image out of place (background subtract and/or AHE)
				//@param im : image to process
				//@param buf: location to write floating point processed image
        template <typename Real>
				void process(uint8_t  const * const im, Real * const buf) {
          for(size_t i = 0; i < nPix; i++){
            buf[i] = (Real) im[i];
          }
        }
        template <typename Real>
				void process(uint16_t  const * const im, Real * const buf) {
          for(size_t i = 0; i < nPix; i++){
            buf[i] = (Real) im[i];
          }
        }
        template <typename Real>
				void process(float  const * const im, Real * const buf) {
          for(size_t i = 0; i < nPix; i++){
            buf[i] = (Real) im[i];
          }
        }
        
  }//laue
}//emsphinx



#endif//_LAUE_IM_PRC_
