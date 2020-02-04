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

int main(int argc, char *argv[]) {

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

	////////////////////////////////////////////////////////////////////////
	//                    Read Inputs / Build Indexers                    //
	////////////////////////////////////////////////////////////////////////

	//read master pattern
	emsphinx::MasterSpectra<double> spec(nml.masterFile);

	//read experimental patterns
	emsphinx::laue::PatternFile pats;
	pats.read(nml.dataFile, nml.dataPath + "/Data/LauePatterns");

	//make sure patterns match expected shape
	if(nml.geoParam.Ny != pats.width() || nml.geoParam.Nz != pats.height()) throw std::runtime_error("patterns don't match expected shape");

	//create pattern processor
	std::unique_ptr< emsphinx::laue::PatternProcessor<double> > lauePrc(new emsphinx::laue::PatternProcessor<double>());
	lauePrc->setSize(nml.geoParam.Ny, nml.geoParam.Nz);
	std::unique_ptr< emsphinx::ImageProcessor<double> > prc(std::move(lauePrc));

	//create back projector
	const size_t dim = nml.bw % 2 == 0 ? nml.bw + 3 : nml.bw + 1;
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
	std::cout << "Indexing patterns from \"" << nml.dataFile << ":" << nml.dataPath << "/Data/LauePatterns\n";
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
	}

	////////////////////////////////////////////////////////////////////////
	//                            Save Outputs                            //
	////////////////////////////////////////////////////////////////////////

	//for now just print to the screen
	for(size_t i = 0; i < pats.numPat(); i++) {
		std::cout << results[i].corr << ": " << results[i].qu[0] << ' ' << results[i].qu[1] << ' ' << results[i].qu[2] << ' ' << results[i].qu[3] << '\n'; 
	}
	std::cout << '\n';

	//also print euler angles
	for(size_t i = 0; i < pats.numPat(); i++) {
		double eu[3];
		xtal::qu2eu(results[i].qu, eu);
		for(size_t j = 0; j < 3; j++) eu[j] *= 57.2957795131;
		std::cout << eu[0] << ' ' << eu[1] << ' ' << eu[2] << '\n';
	}

	//done
	return 0;

} catch (std::exception& e) {
	std::cout << e.what() << '\n';
	return EXIT_FAILURE;
}

}
