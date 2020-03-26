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
 
#include <iostream>

#include "idx/master.hpp"
#include "idx/indexer.hpp"

#include "modality/laue/nml.hpp"
#include "modality/laue/pattern.hpp"
#include "modality/laue/imprc.hpp"
#include "modality/laue/detector.hpp"


#include "util/timer.hpp"
#include "idx/master.hpp"
#include "sht/sht_xcorr.hpp"
#include "xtal/rotations.hpp"

//@brief     : get the extension of a file name
//@param name: file name to get extension of
//@return    : extension (all lower case)
std::string getFileExt(const std::string name) {
	size_t pos = name.find_last_of(".");//find the last '.' in the name
	if(std::string::npos == pos) return "";//handle files with no extension
	std::string ext = name.substr(pos+1);//extract the file extension
	std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});//convert to lowercase
	return ext;
}

//@brief   : read a list of xyz normals from an input stream
//@param is: input stream to read comma separated normals from
//@return  : list of normals
template <typename Real>
std::vector<Real> readNormals(std::istream& is) {
	Real n[3];
	char delim;//space to read common
	std::vector<Real> normals;
	while(is >> n[0] >> delim >> n[1] >> delim >> n[2]) {
		normals.push_back(n[0]);
		normals.push_back(n[1]);
		normals.push_back(n[2]);
	}
	return normals;
}

//@brief         : read a list of xyz normals from a csv
//@param fileName: csv to read normals from
//@return        : list of normals
template <typename Real>
std::vector<Real> readNormalFile(std::string fileName) {
	std::ifstream is(fileName);
	return readNormals<Real>(is);
}

//@brief       : convert from a list of vectors to a spherical image by placing a delta function at each point
//@param hkls  : location to read unit vectors from as x,y,z,x,y,z,...,
//@param numDir: number of normals to read
//@param dim   : side length of square lambert grid to write image to
//@param nh    : location to write north hemisphere of spherical image
//@param sh    : location to write south hemisphere of spherical image
template <typename Real>
void fillLambert(Real const * const hkls, const size_t numDir, const size_t dim, Real * const nh, Real * const sh) {
	for(size_t i = 0; i < numDir; i++) {//loop over peak locations
		//project peak from sphere to square lambert
		Real x, y;
		Real const * const n = hkls + 3 * i;
		// emsphinx::square::lambert::sphereToSquare(n[0], n[1], n[2], x, y);
		emsphinx::square::lambert::sphereToSquare(-n[0], -n[1], n[2], x, y);

		//spread peak over 4 pixels (reverse bilinear interpolate)
		image::BiPix<Real> pix;
		pix.bilinearCoeff(x, y, dim, dim);
		Real * const ptr = std::signbit(n[2]) ? nh : sh;
		for(size_t j = 0; j < 4; j++) ptr[pix.inds[j]] += pix.wgts[j];
	}
}

int main(int argc, char *argv[]) {

	const bool debug = false;//should extra debug info be saved

try {

	////////////////////////////////////////////////////////////////////////
	//                          Parse Arguments                           //
	////////////////////////////////////////////////////////////////////////

	//check argument count
	if(2 != argc) {
		std::cout << "useage: ";
		std::cout << "\tindex using a nml file : " << argv[0] << " input.nml\n";
		std::cout << "\tgenerate a template nml: " << argv[0] << " -t\n";
		return EXIT_FAILURE;
	}

	//check for template request
	const std::string nmlName(argv[1]);
	if(0 == nmlName.compare("-t")) {
		std::ofstream os(std::string("IndexLaue") + ".nml");//create programname.nml (in the future argv[0] could be used with std::filesystem to remove full path)
		emsphinx::laue::Namelist nml;
		os << nml.to_string();
		return EXIT_SUCCESS;
	}

	//read nml and parse
	emsphinx::laue::Namelist nml;
	{
		std::ifstream is(nmlName);
		std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
		std::string warning = nml.from_string(str);
		if(!warning.empty()) {
			std::cout << "\n * * * * * * warning: some namelist parameters weren't used: " << warning << " * * * * * * \n" << std::endl;
		}
	}
	const size_t dim = nml.bw % 2 == 0 ? nml.bw + 3 : nml.bw + 1;

	////////////////////////////////////////////////////////////////////////
	//                    Read Inputs / Build Indexers                    //
	////////////////////////////////////////////////////////////////////////

	//blank master pattern
	emsphinx::MasterSpectra<double> spec;//(mp, nml.bw);

	if("txt" == getFileExt(nml.masterFile)) {
		//create master pattern from directions if needed
		std::vector<double> masterDirs;//spherical directions to place delta functions at for
		masterDirs = readNormalFile<double>(nml.masterFile);
		const size_t numMas = masterDirs.size() / 3;

		//now build a spherical grid (square lambert) and build image from normals
		emsphinx::MasterPattern<double> mp(dim);
		fillLambert(masterDirs.data(), numMas, dim, mp.nh.data(), mp.sh.data());
		mp.toLegendre();

		//save out the constructed master pattern to a raw binary file
		if(debug) {
			std::ofstream os("mp.raw", std::ios::out | std::ios::binary);
			os.write((char*)mp.nh.data(), mp.nh.size() * sizeof(double));
			os.write((char*)mp.sh.data(), mp.sh.size() * sizeof(double));
		}

		//convert to square legendre harmonics
		spec = emsphinx::MasterSpectra<double>(mp, nml.bw);
	} else {
		//read master pattern
		spec = emsphinx::MasterSpectra<double>(nml.masterFile);

	}

	//read experimental patterns
	emsphinx::laue::PatternFile pats;
	pats.read(nml.dataFile, nml.dataPath + "/EMData/Laue/LauePatterns");
  
	//make sure patterns match expected shape
	if(nml.geoParam.Ny != pats.width() || nml.geoParam.Nz != pats.height()) throw std::runtime_error("patterns don't match expected shape");

	//create pattern processor
	std::unique_ptr< emsphinx::laue::PatternProcessor<double> > lauePrc(new emsphinx::laue::PatternProcessor<double>());
	lauePrc->setSize(nml.geoParam.Ny, nml.geoParam.Nz);
	std::unique_ptr< emsphinx::ImageProcessor<double> > prc(std::move(lauePrc));

	//create back projector
	std::unique_ptr< emsphinx::laue::BackProjector<double> > lauePrj(new emsphinx::laue::BackProjector<double>(nml.geoParam, dim, 1.0, nml.bandPass));
	std::unique_ptr< emsphinx::BackProjector<double> > prj(std::move(lauePrj));

	//build correlator for phase
	std::vector< std::unique_ptr< emsphinx::sphere::PhaseCorrelator<double> > > corrs;

	if(nml.normed) {
		//need to buid up spectra of function squared and of window
		//also need to normalize back projected pattern properly
		throw std::runtime_error("normalized cross correlation not yet implemented for Laue");
		/*
		corrs.push_back( std::unique_ptr< emsphinx::sphere::NormalizedCorrelator<double> >(new emsphinx::sphere::NormalizedCorrelator<double>(
			nml.bw,
			spec.data(),
			//flm2
			spec.pointGroup().zMirror(),
			spec.pointGroup().zRot(),
			//mlm
		) ) );
		*/

	} else {
		spec.resize(nml.bw);
		std::shared_ptr< std::vector< std::complex<double> > > flm = std::make_shared< std::vector< std::complex<double> > >(spec.data(), spec.data() + nml.bw * nml.bw);
		corrs.push_back( std::unique_ptr< emsphinx::sphere::UnNormalizedCorrelator<double> >(new emsphinx::sphere::UnNormalizedCorrelator<double>(
			nml.bw,
			flm,
			spec.pointGroup().zMirror(),
			spec.pointGroup().zRot()
		) ) );

	}

	//for now just make a single indexer
	emsphinx::Indexer<double> idx(nml.bw, std::move(prc), std::move(prj), corrs);

	////////////////////////////////////////////////////////////////////////
	//                             Print Info                             //
	////////////////////////////////////////////////////////////////////////

	std::cout << "Running program \"" << argv[0] << "\"\n";
	std::cout << "\tCompiled From        : " << __FILE__ << '\n';
	std::cout << "\tGit Branch           : " << emsphinx::GitBranch << '\n';
	std::cout << "\tCommit Hash          : " << emsphinx::GitHash << '\n';
	std::cout << "\tVersion String       : " << emsphinx::Version << '\n';
	std::cout << '\n';
	std::cout << "Geometry\n";
	std::cout << "\tthrough beam position: " << nml.geoParam.Dy << ' ' << nml.geoParam.Dz << '\n';
	std::cout << "\tslit dimensions      : " << nml.geoParam.Lw << ' ' << nml.geoParam.Lh << '\n';
	std::cout << "\tslit position        : " << nml.geoParam.Lx << ' ' << nml.geoParam.Ly << ' ' << nml.geoParam.Lz << '\n';
	std::cout << "\tsample position      : " << nml.geoParam.Sx << '\n';
	std::cout << "\tsample detector dist : " << nml.geoParam.sampletodetector << '\n';
	std::cout << "\tvoltage range        : " << nml.geoParam.VoltageL << ' ' << nml.geoParam.VoltageH << '\n';
	std::cout << "\tsample thickness     : " << nml.geoParam.samplethickness << '\n';
	std::cout << "\tdetector size        : " << nml.geoParam.Ny << 'x' << nml.geoParam.Nz << " (pixel size = " << nml.geoParam.ps << ")\n";
	std::cout << '\n';
	std::cout << "Indexing patterns from \"" << nml.dataFile << ":" << "/EMData/Laue/LauePatterns\n";
	std::cout << "\tTotal Patterns       : " << pats.numPat() << '\n';
	std::cout << "\tPattern bitdepth     : " << pats.pixBytes() * 8 << '\n';
	std::cout << "\tdetector ROI Mask    : " << ( !nml.mask_ib.hasShape() ? "entire detector" : nml.mask_ib.to_string() )<< '\n';
	std::cout << "\tbandpass             : " << nml.bandPass[0] << ' ' << nml.bandPass[1]  << '\n';
	std::cout << "\n";
	std::cout << "Against master pattern:\n";
	std::cout << "\tFile Name            : " << nml.masterFile << '\n';
	std::cout << "\tPoint Group          : " << spec.pointGroup().name() << '\n';
	std::cout << "\tZ Rotational Symmetry: " << (int)spec.pointGroup().zRot() << '\n';
	std::cout << "\tEquatorial Mirror    : " << (spec.pointGroup().zMirror() ? "yes" : "no") << '\n';
	std::cout << "\n";
	std::cout << "Indexing with\n";
	std::cout << "\tBandwidth            : " << nml.bw << '\n';
	std::cout << "\tSide Length          : " << fft::fastSize(uint32_t(2 * nml.bw - 1)) << '\n';
	// std::cout << "\tThread Count         : " << idxData.threadCount << '\n';
	// std::cout << "\tBatch Size           : " << nml.batchSize << '\n';
	std::cout << "\tThread Count         : " << 1 << '\n';
	std::cout << "\tBatch Size           : " << 1 << '\n';
	std::cout.flush();

	////////////////////////////////////////////////////////////////////////
	//                            Do Indexing                             //
	////////////////////////////////////////////////////////////////////////

	//allocate workspace, output data, and get pointer to data as different data types
	std::vector< emsphinx::Result<double> > results(pats.numPat());
	std::vector<char> imBuff(pats.imBytes());
	uint8_t  * const p8  = (uint8_t *) imBuff.data();
	uint16_t * const p16 = (uint16_t*) imBuff.data();
	float    * const p32 = (float   *) imBuff.data();

	//loop over images indexing
	for(size_t i = 0; i < pats.numPat(); i++) {
		pats.extract(imBuff.data(), 1);//get next image
		switch(pats.pixelType()) {
			case emsphinx::ImageSource::Bits::U8 :
				idx.indexImage(p8 , &results[i], 1, nml.refine);
			break;

			case emsphinx::ImageSource::Bits::U16:
				idx.indexImage(p16, &results[i], 1, nml.refine);
			break;

			case emsphinx::ImageSource::Bits::F32:
				idx.indexImage(p32, &results[i], 1, nml.refine);
			break;

			case emsphinx::ImageSource::Bits::UNK:
				throw std::runtime_error("unknown laue image bitdepth");
			break;
		}

		//write out the back projected pattern
		if(debug) {
			std::ofstream os("sph.raw", std::ios::out | std::ios::binary);//open a file called "out.raw" in cwd, open for binary write
			std::cout << i << '\n';
			os.write((char*)idx.sph.data(), idx.sph.size() * sizeof(double));//write sph.size() doubles to the output file from sph.data()
		}

		//write out the cross correlation grid
		if(debug) {
			std::ofstream os("out.raw", std::ios::out | std::ios::binary);//open a file called "out.raw" in cwd, open for binary write
			fft::vector<double> const& xc = idx.xc[0]->getXC();
			os.write((char*)xc.data(), xc.size() * sizeof(double));//write sph.size() doubles to the output file from sph.data()
		}
	}

	////////////////////////////////////////////////////////////////////////
	//                            Save Outputs                            //
	////////////////////////////////////////////////////////////////////////

	//correct for y rotation
	xtal::Quat<double> yquat(std::sqrt(0.5), 0.0, std::sqrt(0.5) * emsphinx::pijk, 0.0);  	
	for(size_t i = 0; i < pats.numPat(); i++) {
		xtal::quat::mul(results[i].qu, yquat.data(), results[i].qu);//undo 90 degree rotation @ y for back projection
	}

	//for now just print to the screen
	for(size_t i = 0; i < pats.numPat(); i++) {
		//std::cout << results[i].corr << ": " << results[i].qu[0] << ' ' << results[i].qu[1] << ' ' << results[i].qu[2] << ' ' << results[i].qu[3] << '\n'; 
		std::cout << "(" << results[i].qu[0] << ", " << results[i].qu[1] << ", " << results[i].qu[2] << ", " << results[i].qu[3] << "),\n"; 
	}
	std::cout << '\n';
	/*
	//also print euler angles
	for(size_t i = 0; i < pats.numPat(); i++) {
		double eu[3];
		xtal::qu2eu(results[i].qu, eu);
		for(size_t j = 0; j < 3; j++) eu[j] *= 57.2957795131;
		std::cout << eu[0] << ' ' << eu[1] << ' ' << eu[2] << '\n';
	}
	*/
	//done
	return 0;

} catch (std::exception& e) {
	std::cout << e.what() << '\n';
	return EXIT_FAILURE;
}

}
