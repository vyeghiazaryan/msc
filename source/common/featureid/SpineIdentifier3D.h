/***
 * millipede: SpineIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_SPINEIDENTIFIER3D
#define H_MILLIPEDE_SPINEIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class SpineIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MIN_MEAN_GREY_VALUE = 180;
	static const int MIN_VOXELS_PER_SLICE = 800;
	static const int MAX_VOXELS_PER_SLICE = 6000;
	static const double MIN_ASPECT_RATIO_XY = 0.25;
	static const double MAX_ASPECT_RATIO_XY = 4;

	//#################### CONSTRUCTORS ####################
public:
	SpineIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool is_spine(const PFNodeID& node, const BranchProperties& properties) const;
};

}

#endif
