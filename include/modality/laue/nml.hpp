
#ifndef _LAUE_NML_
#define _LAUE_NML_

#include <vector>
#include <string>

namespace emsphinx {

	namespace laue {

		//@brief: encapsulation of all the parameters needed for Laue indexing
		struct Namelist {

			std::string              ipath          ;//input file path
			std::string              patFile        ;//pattern file [txt]
			std::string              masterFile     ;//master pattern files

			int32_t                  patDims[2]     ;//image width and height in pixels
      int32_t                  bandpass[2]    ;//bandpass filter bounds
      std::vector<uint8_t>     mask_ib        ;//mask of the incident beam and other pixels to ignore
  
  /* Coordinate system : 
   * O(0, 0, 0) is the origin of the lab frame located at source point
   * X is the X-ray direction, normal to the detector
   * Y is horizontal on the detector
   * Z is vertical on the detector
  */
  
      double                   kv[2]          ;//source bandwidth in kV
      double                   Lx, Ly, Lz     ;//slits center in the lab frame, mm
      double                   Lw, Lh         ;//slits opening aera, mm
      double                   Sx, t, vsize   ;//sample position, thickness and voxel size, in mm
      uint16_t                 Ny, Nz         ;//detector size in pixels
      double                   psize          ;//pixel size, in mm
      double                   Px, Py, Pz     ;//Px is the perpendicular distance to the detector, Py, Pz the vertors in the detector plane from the det center to the optical axis intersection, in mm

			int32_t                  bw             ;//spherical harmonic bandwidth to index with
			bool                     normed         ;//should normalized cross correlation be used instead of unnormalized
			bool                     refine         ;//should newton's method refinement be used instead of subpixel interpolation
			int32_t                  nThread        ;
			int32_t                  batchSize      ;

			std::string              opath          ;//output file path
			std::string              dataFile       ;//output HDF file
		};

	}//laue

}//emsphinx

#endif//_LAUE_NML_
