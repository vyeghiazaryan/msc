/***
 * millipede: SpleenIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_SPLEENIDENTIFIER3D
#define H_MILLIPEDE_SPLEENIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class SpleenIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MIN_MEAN_GREY_VALUE = 148;
	static const int MAX_MEAN_GREY_VALUE = 175;
	static const int MIN_VOXELS_PER_SLICE = 200;
	static const int MAX_VOXELS_PER_SLICE = 7000;

	//#################### CONSTRUCTORS ####################
public:
	SpleenIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool is_candidate(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const;
};

}

#endif
