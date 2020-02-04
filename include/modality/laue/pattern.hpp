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

#ifndef _LAUE_PATTERN_H_
#define _LAUE_PATTERN_H_

#include <vector>
#include <string>
#include <mutex>

#include "idx/base.hpp"

namespace emsphinx {

	namespace laue {

		//@brief: abstract base class to hold patterns for indexing
		//@note: there aren't currently large datasets to index so I'll just read all the patterns into memory up front
		class PatternFile : public ImageSource {
			public:
		
				//@brief : get patterns held in file
				//@return: number of patterns
				size_t numPat() const {return num;}

				//@brief    : extract the next batch of patterns into a buffer
				//@param out: location to write
				//@param cnt: maximum number of patterns to extract
				//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
				//@note     : implementations should be thread safe (only for other extract calls)
				std::vector<size_t> extract(char * const out, const size_t cnt) const;

				//@brief     : read experimental patterns from a file
				//@param name: h5 file to read from
				//@param dset: path to the dataset within the file
				void read(const std::string name, const std::string dset);

			protected:
				size_t               num;//number of patterns
				size_t               byt;//bytes per pattern (w * h * bytes/pix)
				std::vector<char>    buf;//raw pattern buffer
				mutable std::mutex   mut;//mutex for thread safe access to extraction
				mutable char const * ptr;//pointer to next pattern
				mutable size_t       idx;//index of next pattern
		};

	}//laue

}//emsphinx

////////////////////////////////////////////////////////////////////////
//                        Non-abstract Classes                        //
////////////////////////////////////////////////////////////////////////


#include "H5Cpp.h"

#include <algorithm>
#include <numeric>
#include <cctype>
#include <sstream>
#include <stack>
#include <cmath>

namespace emsphinx {
	
	namespace laue {

		//@brief     : read experimental patterns from a file
		//@param name: h5 file to read from
		//@param dset: path to the dataset within the file
		//@return    : shared pointer to a pattern file
		void PatternFile::read(const std::string name, const std::string dset) {
			try {
				//open the dataset and get shape
				H5::H5File file = H5::H5File(name, H5F_ACC_RDONLY);//open the file
				H5::DataSet dSet = file.openDataSet(dset);//get the dataset we're after
				hsize_t dims[3];
				if(3 != dSet.getSpace().getSimpleExtentNdims()) throw std::runtime_error("hdf pattern dataset must be 3D");
				dSet.getSpace().getSimpleExtentDims(dims);//read extent in each dimension

				//determine pixel type
				Bits fBits = ImageSource::Bits::UNK;
				H5::DataType type = dSet.getDataType();
				if     (type == H5::DataType(H5::PredType::NATIVE_UINT8 )) fBits = ImageSource::Bits::U8 ;
				else if(type == H5::DataType(H5::PredType::NATIVE_UCHAR )) fBits = ImageSource::Bits::U8 ;
				else if(type == H5::DataType(H5::PredType::NATIVE_UINT16)) fBits = ImageSource::Bits::U16;
				else if(type == H5::DataType(H5::PredType::NATIVE_FLOAT )) fBits = ImageSource::Bits::F32;
				if(ImageSource::Bits::UNK == fBits) throw std::runtime_error("only uint8, uint16, and float hdf patterns are supported");

				//save meta data
				bits = fBits  ;
				w    = dims[2];
				h    = dims[1];
				num  = dims[0];

				//allocate memory and read
				buf = std::vector<char>(imBytes() * num);
				ptr = buf.data();
				idx = 0;//next pattern is first pattern

				if     (type == H5::DataType(H5::PredType::NATIVE_UINT8 )) dSet.read((void*)buf.data(), H5::PredType::NATIVE_UINT8 );
				else if(type == H5::DataType(H5::PredType::NATIVE_UCHAR )) dSet.read((void*)buf.data(), H5::PredType::NATIVE_UCHAR );
				else if(type == H5::DataType(H5::PredType::NATIVE_UINT16)) dSet.read((void*)buf.data(), H5::PredType::NATIVE_UINT16);
				else if(type == H5::DataType(H5::PredType::NATIVE_FLOAT )) dSet.read((void*)buf.data(), H5::PredType::NATIVE_FLOAT);
				else throw std::logic_error("unsupported hdf5 data type");
			} catch (H5::Exception& e) {//convert hdf5 exceptions into std::exception
				std::ostringstream ss;
				ss << "H5 error attempting to read Laue patterns:\n";
				ss << "\tfile - " << name << '\n';
				ss << "\tdset - " << dset  << '\n';
				ss << "\tfunc - " << e.getCFuncName() << '\n';
				ss << "detailed message:\n" << e.getCDetailMsg();
				throw std::runtime_error(ss.str());
			}
		}

		//@brief    : extract the next batch of patterns into a buffer
		//@param ptr: location to write
		//@param cnt: maximum number of patterns to extract
		//@return   : index of first and last pattern extract (exlusive so second - first patterns were extracted)
		//@note     : implementations should be thread safe (only for other extract calls)
		std::vector<size_t> PatternFile::extract(char * const out, const size_t cnt) const {
			std::lock_guard<std::mutex> lock(mut);//only 1 thread can get patterns at once
			const size_t numExt = std::min(numPat() - idx, cnt);//determine how many patterns we can extract
			char const * const ptrNew = ptr + numExt * imBytes();//get end of range to copy
			std::copy(ptr, ptrNew, out);//copy patterns to output
			std::vector<size_t> res(numExt);//build vector to hold extraction index list
			std::iota(res.begin(), res.end(), idx);//fill index list
			idx += numExt;//update index of next patterns
			ptr = ptrNew;//update pointer to next pattern
			if(idx == numPat()) ptr = 0;//we've run out of patterns
			return res;
		}

	}//laue

}//emsphinx

#endif//_LAUE_PATTERN_H_
