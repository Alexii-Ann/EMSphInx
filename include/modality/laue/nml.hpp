
#ifndef _LAUE_NML_
#define _LAUE_NML_

#include <string>

#include "H5Cpp.h"
#include "idx/roi.h"
#include "util/nml.hpp"

namespace emsphinx {

	namespace laue {
		//@brief: class to encapsulate geometry parameters
		struct GeomParam {
			/* Coordinate system : 
			 * O(0, 0, 0) is the origin of the lab frame located at source point
			 * X is the X-ray direction, normal to the detector
			 * Y is horizontal on the detector
			 * Z is vertical on the detector
			 */
			double       Dy, Dz     ;//position of through beam on detector in mm relative to detector center
			double       Lw, Lh     ;//width/height of X-Ray slits in mm
			double   Lx, Ly, Lz     ;//position of slit opening center in mm
			double   Sx             ;//position of sample
			double   VoltageL       ;//lower  X-Ray energy in kV
			double   VoltageH       ;//higher X-Ray energy in kV
			double   samplethickness;//thickness of sample in mm
			double   ps             ;//pixel size in mm
			uint16_t    Ny, Nz      ;//detector size in pixels

			//@brief: initialize / reset values with defaults
			void clear();

			//@brief: sanity check some basic values
			//@note : throws an exception if check fails
			void sanityCheck() const;

			//@brief    : read geometry data from an h5 file
			//@param grp: hdf5 group to read data from
			void read(H5::Group grp);
		};

		//@brief: encapsulation of all the parameters needed for Laue indexing
		struct Namelist {

			std::string              ipath          ;//input file path
			std::string              dataFile       ;//pattern file (H5 with geometry and Laue patterns)
			std::string              dataPath       ;//path to experimental data folder within dataFile
			std::string              masterFile     ;//master pattern files

			GeomParam                geoParam       ;//experimental geometry

			uint16_t                 bandPass[2]    ;//bandpass filter bounds
			RoiSelection             mask_ib        ;//mask of the incident beam and other pixels to ignore

			int32_t                  bw             ;//spherical harmonic bandwidth to index with
			bool                     normed         ;//should normalized cross correlation be used instead of unnormalized
			bool                     refine         ;//should newton's method refinement be used instead of subpixel interpolation
			int32_t                  nThread        ;//number of threads to index with
			int32_t                  batchSize      ;//number of patterns to dispatch to a thread at once

			std::string              opath          ;//output file path
			std::string              outputFile     ;//output HDF file

			std::string              namelistFile   ;//namelist file used for building (actual file contents)

			//@brief: initialize / reset values with defaults
			void clear();

			//@brief: sanity check some basic values
			//@note : throws an exception if check fails
			void sanityCheck() const;

			//@brief    : parse indexing values from a namelist file
			//@param nml: file to parse (actual contents)
			//@return   : warning string (list of unused namelist values)
			std::string from_string(std::string nml);

			//@brief    : parse indexing values from a namelist file
			//@param nml: file to parse (namelist object)
			//@return   : warning string (list of unused namelist values)
			std::string parse_nml(nml::NameList& nml);

			//@brief    : convert a namelist to a string
			//@param nml: namelist name (for fortran)
			//@return   : namelist file string
			std::string to_string(std::string nml = "EMSphInx") const;
		};

	}//laue

}//emsphinx

#include <sstream>
#include <cmath>

namespace emsphinx {

	namespace laue {

		//@brief: initialize / reset values with defaults
		void GeomParam::clear() {
			     Dy = Dz    = NAN;
			     Lw = Lh    = NAN;
			Sx              = NAN;
			VoltageL        = NAN;
			VoltageH        = NAN;
			samplethickness = NAN;
			ps              = NAN;
			Ny = Nz = 0;
		}

		//@brief: sanity check some basic values
		//@note : throws an exception if check fails
		void GeomParam::sanityCheck() const {
			const double detW = ps * Ny;
			const double detH = ps * Nz;
			if(std::fabs(Dy) > detW / 2 || std::fabs(Dz) > detH / 2) throw std::runtime_error("through beam not on detector");
			if(Ny < 2 || Ny > 16384 || 
			   Nz < 2 || Nz > 16384) throw std::runtime_error("unreasonable pattern dimension (outside [2, 16384] pix)");
			if(Lw < 0.1 || Lw > 100.0) throw std::runtime_error("unreasonable slit width (outside [0.1, 100] mm)");
			if(Lh < 0.1 || Lh > 100.0) throw std::runtime_error("unreasonable slit height (outside [0.1, 100] mm)");
			if(Lx < 5.0 || Lx > 1000.0) throw std::runtime_error("unreasonable slit X (outside [5, 1000] mm)");
			if(Ly < -100.0 || Ly > 100.0) throw std::runtime_error("unreasonable slit Y (outside [-100, 100] mm)");
			if(Lz < -100.0 || Lz > 100.0) throw std::runtime_error("unreasonable slit Z (outside [-100, 100] mm)");
			if(Sx < 5.0 || Sx > 1000.0) throw std::runtime_error("unreasonable sample X (outside [5, 1000] mm)");
			if(VoltageL > VoltageH) throw std::runtime_error("voltage upper bound must be >= lower bound");
			if(VoltageL < 10.0 || VoltageL > 100.0) throw std::runtime_error("unreasonable voltage (outside [10, 100] kV)");
			if(VoltageH < 10.0 || VoltageH > 100.0) throw std::runtime_error("unreasonable voltage (outside [10, 100] kV)");
			if(samplethickness < 0.1 || samplethickness > 100.0) throw std::runtime_error("unreasonable sample thickness (outside [0.1, 100] mm)");
			if(ps > 10.0 || ps < 0.0001) throw std::runtime_error("unreasonable pixel size (outside of [0.1, 10000] um");
		}

		//@brief    : read geometry data from an h5 file
		//@param grp: hdf5 group to read data from
		void GeomParam::read(H5::Group grp) {
			clear();
			try {
				H5::DataSpace ds(H5S_SCALAR);//single element data space
				grp.openDataSet("Dy"             ).read(&Dy             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Dz"             ).read(&Dz             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Lw"             ).read(&Lw             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Lh"             ).read(&Lh             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Lx"             ).read(&Lx             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Ly"             ).read(&Ly             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Lz"             ).read(&Lz             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Sx"             ).read(&Sx             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("VoltageL"       ).read(&VoltageL       , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("VoltageH"       ).read(&VoltageH       , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("samplethickness").read(&samplethickness, H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("ps"             ).read(&ps             , H5::PredType::NATIVE_DOUBLE, ds);
				grp.openDataSet("Ny"             ).read(&Ny             , H5::PredType::NATIVE_UINT16, ds);
				grp.openDataSet("Nz"             ).read(&Nz             , H5::PredType::NATIVE_UINT16, ds);
			} catch (H5::Exception& e) {//convert hdf5 exceptions into std::exception
				std::ostringstream ss;
				ss << "H5 error attempting to read Laue geometry:\n";
				ss << "\tfunc - " << e.getCFuncName() << '\n';
				ss << "detailed message:\n" << e.getCDetailMsg();
				throw std::runtime_error(ss.str());
			}
		}

		//@brief: initialize / reset values with defaults
		void Namelist::clear() {
			//clear inputs
			ipath = dataFile = dataPath = masterFile = "";

			//clear geometry
			geoParam.clear();

			//clear image processing
			bandPass[0] = 0x0000;
			bandPass[1] = 0xFFFF;
			mask_ib.clear();

			//clear indexing parameters
			bw = -1; normed = false; refine = false; nThread = 0; batchSize = 0;

			//clear outputs
			opath = outputFile = "";
		}

		//@brief: sanity check some basic values
		//@note : throws an exception if check fails
		void Namelist::sanityCheck() const {
			if(dataFile.empty()) throw std::runtime_error("no experimental input file");
			if(dataPath.empty()) throw std::runtime_error("no data path");
			if(masterFile.empty()) throw std::runtime_error("no master pattern");

			if(bandPass[0] >= bandPass[1]) throw std::runtime_error("first bandpass cutoff must be less than second");
			geoParam.sanityCheck();

			if(bw < 16 || bw > 512) throw std::runtime_error("unreasonable bandwidth (should be [16, 512])");
			if(nThread < 0) throw std::runtime_error("negative thread count");
			if(batchSize < 0) throw std::runtime_error("negative batch size");
			if(outputFile.empty()) throw std::runtime_error("missing output data file");
		}

		//@brief    : parse indexing values from a namelist file
		//@param nml: file to parse (actual contents)
		//@return   : warning string (list of unused namelist values)
		std::string Namelist::from_string(std::string nml) {
			//start by wrapping input as stream and parsing with namelist reader
			clear();
			namelistFile = nml;//save a copy of the file
			std::istringstream iss(nml);
			nml::NameList nameList;
			nameList.read(iss); 
			return parse_nml(nameList);
		}

		//@brief    : parse indexing values from a namelist file
		//@param nml: file to parse (namelist object)
		//@return   : warning string (list of unused namelist values)
		std::string Namelist::parse_nml(nml::NameList& nml) {
			//parse inputs first
			try {ipath      =         nml.getString("ipath"     );} catch (...) {ipath   .clear();}//if ipath isn't found we'll just use cwd
			     dataFile   = ipath + nml.getString("dataFile"  );
			     dataPath   =         nml.getString("dataPath"  );
			     masterFile = ipath + nml.getString("masterfile");

			//read geometry
			H5::Group expDir = H5::H5File(dataFile, H5F_ACC_RDONLY).openGroup(dataPath.c_str());
			geoParam.read(expDir.openGroup("Header"));

			//parse image processing
			std::vector<int> dims = nml.getInts("bandPass");
			if(2 != dims.size()) throw std::runtime_error("bandPass should have 2 values");
			bandPass[0] = dims[0]; bandPass[1] = dims[1];

			//parse indexing parameters
			bw        = (size_t) nml.getInt ("bw"       );//what bandwidth should be used, if 2*bw-1 is product of small primes it is a good candidate for speed (fft is significant fraction of time): 32,38,41,53,63,68,74,88,95,113,123,158
			normed    = (size_t) nml.getBool("normed"   );//should normalized or unnormalized cross correlation be used
			refine    = (size_t) nml.getBool("refine"   );//should refinement be used
			nThread   = (size_t) nml.getInt ("nthread"  );
			batchSize = (size_t) nml.getInt ("batchsize");//number of patterns per work item (should be large enough that the task is significant but small enough that there are enough jobs for load balancing)

			//parse outputs
			try {opath      = nml.getString("opath"     );} catch (...) {opath     .clear();}//if ipath isn't found we'll just use cwd
			     dataFile   = nml.getString("datafile"  );

			//check for unused inputs
			sanityCheck();
			if(!nml.fullyParsed()) return nml.unusedTokens();
			return "";
		}

		//@brief    : convert a namelist to a string
		//@param nml: namelist name (for fortran)
		//@return   : namelist file string
		std::string Namelist::to_string(std::string nml) const {
			std::ostringstream ss;
			ss << " &" << nml << "\n";//this if for the fortran version
			ss << "!#################################################################\n";
			ss << "! Input Files\n";
			ss << "!#################################################################\n";
			ss << "\n";
			if(!ipath.empty()) {
				ss << "! input path, empty for current working directory\n";
				ss << " ipath      = '" << ipath << "',\n";//ignored for fortran version
				ss << "\n";
			}
			ss << "! experimental data file (hdf5) with patterns to index and geometry (relative to ipath)\n";
			ss << " dataFile    = '" << dataFile << "',\n";
			ss << "\n";
			ss << "! path to experimental data folder within dataFile\n";
			ss << " dataPath    = '" << dataPath << "',\n";
			ss << "\n";
			ss << "! master pattern file (*.sht) with phase to index (relative to ipath)\n";
			ss << " masterfile = " << masterFile << ",\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Pattern Processing\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! bandpass bounds (upper, lower) for frequency filtering\n";
			ss << " bandPass    = " << bandPass[0] << ", " << bandPass[1] << ",\n";
			ss << "\n";
			ss << "! mask of incident beam and any other pixels to ignore when back projecting\n";
			ss << "! 0 (or omitted) to back project the entire detector\n";
			ss << "! x0, y0, dx, dy for a (dx, dy) rectangle starting at pixel (x0, y0)\n";
			ss << "! prepend 'e' for ellipse instead of rectangle\n";
			ss << "! prepend 'i' to specify included instead of excluded region\n";
			ss << "! multiple points for polygon region\n";
			ss << " mask_ib   = " << mask_ib.to_string() << ",\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Indexing Parameters\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! spherical harmonic bandwidth to be used (2*bw-1 should be a product of small primes for speed)\n";
			ss << "! some reasonable values are: 53, 63, 68, 74, 88, 95, 113, 122, 123, 158, 172, 188, 203, 221, 263, 284, 313\n";
			ss << "! a nice range for parameter studies is 53, 68, 88, 113, 158, 203, 263, 338 (~a factor of 1.3 between each)\n";
			ss << "! any value is now pretty fast since the transform is zero padded to the nearest fast FFT size\n";
			ss << " bw         = " << bw << ",\n";
			ss << "\n";
			ss << "! should normalized / unnormalized spherical cross correlation be used?\n";
			ss << "! normalization is more robust for (esp. for lower symmetries) but is slower\n";
			ss << " normed     = ." << (normed ? "TRUE" : "FALSE") << ".,\n";
			ss << "\n";
			ss << "! should newton's method orientation refinement be used?\n";
			ss << "! normalization is more robust for (esp. for lower symmetries) but is slower\n";
			ss << " refine     = ." << (refine ? "TRUE" : "FALSE") << ".,\n";
			ss << "\n";
			ss << "! number of work threads\n";
			ss << "! 0 (or omitted) to multithread with an automatic number of threads\n";
			ss << "! 1 for serial threading\n";
			ss << "! N to multithread with N threads\n";
			ss << " nthread    = " << nThread<< ",\n";
			ss << "\n";
			ss << "! number of patterns to index per work itme (ignored for single threading)\n";
			ss << "! should be large enough to make the task significant compared to thread overhead\n";
			ss << "! should be small enough to enable enough work items for load balancing\n";
			ss << "! should be small enough so nthread * batchsize patterns can be held in memory\n";
			ss << "! 0 (or omitted) to estimate a reasonable value based on speed\n";
			ss << " batchsize  = " << batchSize << ",\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Output Files\n";
			ss << "!#################################################################\n";
			ss << "\n";
			if(!opath.empty()) {
				ss << "! output path, empty for current working directory\n";
				ss << " opath      = '" << opath << "',\n";//ignored for fortran version
				ss << "\n";
			}
			ss << "! output orientation map name relative to opath [must be hdf5 type]\n";
			ss << " outputFile   = '" << outputFile << "',\n"; //should be hdf
			ss << "\n";
			return ss.str();
		}

	}//laue

}//emsphinx


#endif//_LAUE_NML_
